! used in route_run_standalone and from within dynatop2001
! 1. Reads in river
!% Toby Dunne
!% April 2016
module dta_route_processing
    use dta_riv_tree_node
    use dta_utility
    implicit none

    type route_river_info_type
        logical:: is_enabled

        ! from flow_conn and point_con
        type(riv_tree_node), allocatable :: node_list(:)
        ! from river_data
        !% riv_id, area, dist, section_dist, slope, elevation
        double precision, allocatable, dimension(:,:) :: river_data


        ! only write output for node_type=NODE_TYPE_GAUGE
        integer :: flow_output_count
        integer, allocatable, dimension(:) :: flow_output_indexes
        character(64), allocatable, dimension(:) :: flow_output_headers

        character(1024) :: output_file_prefix

    end type

    ! routing time delay - processed from:
    ! * route_river_info_type
    ! * timestep
    ! * velocity params
    type route_time_delay_hist_type
        double precision :: timestep
        integer ::v_mode
        double precision, allocatable ::v_param(:)

        ! each row represents the histogram for the corresponding node_list item
        double precision, allocatable, dimension(:,:) :: route_hist_table
        ! each row represents the histogram index of flow output (:,1) for the corresponding node_list item
        ! (:,2) is the max index of the row histogram
        integer, allocatable, dimension(:,:) :: node_hist_indexes

        ! this is where all the magic happens
        ! flow is added to this table each timestep and is shifted one step left
        double precision, allocatable, dimension(:,:) :: flow_matrix

    end type


contains
    !DEMf_morpeth2_flow_conn.txt
    !DEMf_morpeth2_flow_point.txt
    !DEMf_morpeth2_river_data.txt
    subroutine route_processing_read_info(river_data_file, flow_conn_file, riv)
        use dta_utility
        implicit none
        !from route_river.f90
        ! _river_data.txt
        character(1024) river_data_file
        ! from route_tree.f90
        ! _flow_conn.txt
        character(1024) flow_conn_file
        type(route_river_info_type) riv

        ! locals
        ! _flow_point.txt
        character(1024) :: flow_point_file

        riv%is_enabled = .true.

        flow_point_file = flow_conn_file(1:len_trim(flow_conn_file)-14)//'_flow_point.txt'

        print *, 'read flow tree:', &
            trim(flow_conn_file), &
            trim(flow_point_file)

        call riv_tree_read(riv%node_list, &
            0, 0.0d0, 0.0d0, 0.0d0, &
            flow_conn_file, flow_point_file)

        call read_numeric_list(river_data_file, 8, 1, riv%river_data)

        call route_processing_init(riv)

    end subroutine

    ! GC - same routine as above but specifically for dynatop to use file IDs rather than filenames

    subroutine route_processing_read_info_dyna(riv)

        use dta_utility

        implicit none
        !from route_river.f90
        ! _river_data.txt
        !character(1024) river_data_file
        ! from route_tree.f90
        ! _flow_conn.txt
        !character(1024) flow_conn_file
        type(route_river_info_type) riv

        riv%is_enabled = .true.

        call riv_tree_read_dyna(riv%node_list, &
            0, 0.0d0, 0.0d0, 0.0d0)

        ! File ID is the same as file ID given in dyna_project.f90 GC 13/06
        !print *, 'read river data:', trim(river_data_file)
        call read_numeric_list_fid(502, 8, 1, riv%river_data)

        call route_processing_init(riv)

    end subroutine route_processing_read_info_dyna

    subroutine route_processing_init(riv)
        use dta_utility
        implicit none
        type(route_river_info_type) riv

        ! locals
        integer :: i

        ! only write output for node_type=gauge
        ! count the gauge points
        riv%flow_output_count = 0
        do i=1,size(riv%node_list)
            if(riv%node_list(i)%node_type == NODE_TYPE_GAUGE) then
                riv%flow_output_count = riv%flow_output_count + 1
            endif
        end do

        ! allocate record of indexes
        allocate(riv%flow_output_indexes(riv%flow_output_count))
        allocate(riv%flow_output_headers(riv%flow_output_count))


        riv%flow_output_count = 0
        ! list of indexes to store flow output data
        do i=1,size(riv%node_list)
            if(riv%node_list(i)%node_type == NODE_TYPE_GAUGE) then

                riv%flow_output_count = riv%flow_output_count + 1

                riv%flow_output_headers(riv%flow_output_count) = ''
                write(riv%flow_output_headers(riv%flow_output_count), '(I0)') riv%node_list(i)%gauge_id

                riv%flow_output_indexes(riv%flow_output_count) = riv%flow_output_count
            endif
        end do

    end subroutine route_processing_init

    subroutine route_processing_build_hist(riv, tdh)
        use dta_riv_tree_node
        implicit none
        type(route_river_info_type) riv
        type(route_time_delay_hist_type) tdh

        !locals
        integer :: i, j
        ! location of the river segment in the river_data grid
        ! (:,1) start_i
        ! (:,2) end_i
        integer, dimension(size(riv%node_list),2) :: riv_data_indexes
        double precision :: headwater_max_delay
        double precision :: river_point_max_delay
        double precision :: full_delay, min_dist_outlet
        double precision :: max_dist, min_dist
        double precision :: max_elev, min_elev
        double precision :: reach_dist, reach_slope
        double precision :: cell_dist_a, cell_dist_b, cell_dist, cell_slope
        integer :: riv_start_i, riv_end_i
        integer :: hist_start_i, hist_end_i
        double precision :: hist_cursor_a, hist_cursor_b
        integer :: total_timesteps
        double precision, allocatable :: delay_hist(:)
        integer :: tdh_bin_count

        if(allocated(tdh%node_hist_indexes).eqv..false.) then
            allocate(tdh%node_hist_indexes(size(riv%node_list),2))
        endif

        min_dist_outlet = minval(riv%river_data(:,7))

        do i=1,size(riv%node_list)

            riv_start_i = 0
            riv_end_i = size(riv%river_data,1)
            do j=1,size(riv%river_data,1)
                if(riv_start_i == 0) then
                    if (nint(riv%river_data(j,1)) == riv%node_list(i)%gauge_id) then
                        ! first row in the river data matching this id
                        riv_start_i = j
                    endif
                else
                    if(nint(riv%river_data(j,1)) /= riv%node_list(i)%gauge_id) then
                        ! last index of this matching this id
                        riv_end_i = j-1
                        exit
                    endif
                endif
            end do
            if(riv_start_i == 0) then
                print *, 'river_id has no cells', riv%node_list(i)%gauge_id
                tdh%node_hist_indexes(i,1) = 0
                tdh%node_hist_indexes(i,2) = 0
                riv_data_indexes(i,1) = 0
                riv_data_indexes(i,2) = 0
                riv%node_list(i)%reach_delay = 0
                cycle
            endif

            riv_data_indexes(i,1) = riv_start_i
            riv_data_indexes(i,2) = riv_end_i

            ! column 4 is secion distance
            max_dist = maxval(riv%river_data(riv_start_i:riv_end_i,3))
            ! column 7 is secion ds_distance
            min_dist = minval(riv%river_data(riv_start_i:riv_end_i,7))
            ! column 5 is elevation
            max_elev = maxval(riv%river_data(riv_start_i:riv_end_i,5))
            ! column 8 is ds_elevation
            min_elev = minval(riv%river_data(riv_start_i:riv_end_i,8))

            reach_dist = max_dist - min_dist
            reach_slope = (max_elev - min_elev) / reach_dist

            ! reach delay is time from upstream cell to downstream cell
            ! delay on this reach only
            riv%node_list(i)%reach_delay = calculate_cell_delay(reach_dist, reach_slope, tdh)

            !reach_length_list(i) = reach_dist ! used to distribute percentage of delay
            !reach_delay_list(i) = reach_delay

            ! remove delay downstream of the min outlet
            ! don't need to track this in the histogram
            min_dist = min_dist - min_dist_outlet

            ! slope incorrect for total_downstream_delay (not used in constant velocity)
            riv%node_list(i)%total_downstream_delay = calculate_cell_delay(min_dist, reach_slope, tdh)

        end do

        !% this is the total maximum time to the furtherest upstream river point.
        !% equal to river_point_max_delay plus max(reach_delay)
        headwater_max_delay = 1
        !% this is the maximum distance to any gauge point
        !% this will determine the warm-up time (the first time that any water can
        !% reach the outlet from the furtherest input point)
        river_point_max_delay = 1

        do i=1,size(riv%node_list)

            !% full distance from furtherest river cell to the outlet
            full_delay = riv%node_list(i)%reach_delay + riv%node_list(i)%total_downstream_delay

            ! +1 to make the index one based rather than zero based
            tdh%node_hist_indexes(i,1) = floor(riv%node_list(i)%total_downstream_delay) + 1
            tdh%node_hist_indexes(i,2) = ceiling(full_delay) + 1

            if(riv%node_list(i)%total_downstream_delay > river_point_max_delay) then
                ! just recorded to see how many timesteps the flow will reach the bottom
                river_point_max_delay = riv%node_list(i)%total_downstream_delay
            endif
            if (full_delay > headwater_max_delay) then
                headwater_max_delay = full_delay
            endif

        end do

        !total_timesteps = ceiling(headwater_max_delay / tdh%timestep) + 1
        total_timesteps = ceiling(headwater_max_delay) + 1
        !warmup_timesteps = ceiling(river_point_max_delay / tdh%timestep) + 1

        !print *, 'warm-up timesteps = ', warmup_timesteps
        !print *, 'flow timesteps = ', total_timesteps

        if(allocated(tdh%route_hist_table)) then
            deallocate(tdh%route_hist_table)
            deallocate(tdh%flow_matrix)
        endif

        ! Allocate the histogram table and the flow matrix
        allocate(tdh%route_hist_table(size(riv%node_list), total_timesteps))
        allocate(tdh%flow_matrix(size(riv%node_list), total_timesteps))

        !print*, 'point_count      =', size(tdh%route_hist_table,1) ! rows
        !print*, 'route_step_count =', size(tdh%route_hist_table,2) ! columns

        tdh%route_hist_table(:,:) = 0
        tdh%flow_matrix(:,:) = 0

        do i=1,size(riv%node_list)

            !indexes into the
            riv_start_i = riv_data_indexes(i,1)
            riv_end_i = riv_data_indexes(i,2)

            if(riv_start_i == 0) then
                cycle
            endif


            !% Create the histogram for the flow input from dynamic topmodel
            !% flow is distributed along the river reach by the area accumulation

            !indexes the time delay cell to write
            hist_start_i = tdh%node_hist_indexes(i,1)
            hist_end_i = tdh%node_hist_indexes(i,2)

            tdh_bin_count = (hist_end_i - hist_start_i) + 1

            !print*, 'hist cells: ',tdh_bin_count, ' delay: ',  &
            !    riv%node_list(i)%reach_delay, hist_start_i, hist_end_i

            allocate(delay_hist(tdh_bin_count))
            delay_hist(:) = 0

            ! loop over each cell in the input (river distance)
            do j=riv_start_i,riv_end_i

                ! cell_dist = dist - ds_dist
                ! length of this cell
                cell_dist = riv%river_data(j,3) - riv%river_data(j,7)

                ! slope based calculation will be wrong - not used in constant velocity version
                ! need to use slope for the entire downstream
                ! need to determine if sum of slope delays needs to be calculated
                ! or if average slope from this point to reach outlet is sufficient
                cell_slope = riv%river_data(j,6)

                cell_dist_a = riv%river_data(j,4)
                cell_dist_b = cell_dist_a + cell_dist

                ! how long flow takes to the outlet of this reach from the bottom of this cell
                hist_cursor_a = calculate_cell_delay(cell_dist_a, cell_slope, tdh)
                ! how long flow takes to the outlet of this reach from the top of this cell
                hist_cursor_b = calculate_cell_delay(cell_dist_b, cell_slope, tdh)

                ! add the area contribution to the histogram at the cursor locations
                call hist_add_value(delay_hist, hist_cursor_a, hist_cursor_b, riv%river_data(j,2), tdh_bin_count)

            end do

            ! normalise to 1
            delay_hist = delay_hist / sum(delay_hist)
            !write the histogram to the histogram mask table
            tdh%route_hist_table(i,hist_start_i:hist_end_i) = delay_hist

            deallocate(delay_hist)
        end do

    end subroutine route_processing_build_hist

    subroutine hist_add_value(delay_hist, hist_cursor_a, hist_cursor_b, value, tdh_bin_count)
        implicit none
        double precision :: delay_hist(*)
        double precision :: hist_cursor_a
        double precision :: hist_cursor_b
        double precision :: value
        integer :: tdh_bin_count

        double precision :: fraction_value
        integer :: hist_a
        integer :: hist_b

        ! get the dest cell index (+1 to make this 1 based indexing)
        hist_a = floor(hist_cursor_a) + 1
        hist_b = floor(hist_cursor_b) + 1
        if(hist_a < 0) then
            print *, 'hist_a < 0'
            stop
        endif

        !print *, 'add ',hist_a,hist_b,tdh_bin_count

        if(hist_b > tdh_bin_count) then
            print *, 'delay error - sum of delays greater than bin_count',hist_b,tdh_bin_count
            hist_b = tdh_bin_count
            return
        endif

        if(hist_a == hist_b) then
            ! if the delay is smaller than a single bin
            ! scale the accumulation by the fraction in this cell
            fraction_value = (hist_cursor_b - hist_cursor_a)
            delay_hist(hist_a:hist_b) = delay_hist(hist_a:hist_b) + value * fraction_value
        else
            !    fill with area  accumulation for the cell

            ! a more accurate method would be to alter the first and last cell amounts by the fractional delay
            ! e.g. hist_cursor_a = 0.6
            !      hist_cursor_b = 2.3
            ! this delay will be put into 3 cells
            ! - first cell delay is 0.6 the 'fractional accumulation' is 0.4
            ! - middle cell(s) 100% of accumulation
            ! - last cell delay is 0.3 the 'fractional accumulation' is 0.3
            ! |______wwww|wwwwwwwwww|www_______|
            !
            ! since a delay of 0.6 cannot be split as we have fixed bin size
            ! (1 - 0.6) * accumulation will be routed at 0 delay

            ! middle cell(s)
            if( (hist_a+1) <= (hist_b-1) ) then
                delay_hist((hist_a+1):(hist_b-1)) = delay_hist((hist_a+1):(hist_b-1)) + value
            endif

            ! first cell (+1 to convert from one based index)
            fraction_value = (1 - (hist_cursor_a - hist_a + 1))
            delay_hist(hist_a) = delay_hist(hist_a) + (value * fraction_value)
            ! last cell (+1 to convert from one based index)
            fraction_value = (hist_cursor_b - hist_b + 1)
            delay_hist(hist_b) = delay_hist(hist_b) + (value * fraction_value)
        endif
    end subroutine


    function calculate_cell_delay(dist, slope, tdh) result(reach_delay)
        implicit none
        double precision :: dist
        double precision :: reach_delay
        double precision :: slope
        type(route_time_delay_hist_type) :: tdh

        if(tdh%v_mode == V_MODE_CONSTANT) then
            !% param is velocity
            ! v_param(1) unit is metres travelled per timestep
            ! v = d/t
            ! t = d/v
            reach_delay = dist / tdh%v_param(1)
        elseif (tdh%v_mode == V_MODE_SLOPE) then
            ! slope delay
            ! get a reach_delay based on total slope
            ! when calculating histogram, and allocate by percentage dist for each cell
            ! v_param(1) when multiplied by slope should give metres travelled per timestep
            reach_delay = dist / (tdh%v_param(1) * slope)
        else
            print *, 'unknown velocity mode'
            stop
        endif
    end function

    subroutine add_accumulation_test()
        implicit none
        double precision, allocatable :: delay_hist(:)
        integer :: tdh_bin_count

        tdh_bin_count = 6

        allocate(delay_hist(tdh_bin_count))

        delay_hist(:) = 0
        call hist_add_value(delay_hist, 2.2d0, 2.4d0, 1.0d0, tdh_bin_count)

        delay_hist(:) = 0
        call hist_add_value(delay_hist, 2.2d0, 3.1d0, 1.0d0, tdh_bin_count)

        delay_hist(:) = 0
        call hist_add_value(delay_hist, 2.2d0, 4.2d0, 1.0d0, tdh_bin_count)

        delay_hist(:) = 0
        call hist_add_value(delay_hist, 2.2d0, 5.2d0, 1.0d0, tdh_bin_count)


        deallocate(delay_hist)

    end subroutine add_accumulation_test

    subroutine build_hist_test()
        implicit none
        type(route_river_info_type) :: riv
        type(route_time_delay_hist_type) :: tdh

        character(1024) river_data_file
        character(1024) flow_conn_file

        integer :: j

        tdh%timestep = 3600.0d0
        tdh%v_mode = V_MODE_CONSTANT
        allocate(tdh%v_param(1))
        tdh%v_param(1) = 30.133d0

        river_data_file = '54020000_riv_data.dat'
        flow_conn_file = '54020000_flow_conn.txt'



        call route_processing_read_info(river_data_file, flow_conn_file, riv)

        do j=1,size(riv%river_data, 1)
           if (riv%river_data(j,7)<0.01) then
              riv%river_data(j,7) = riv%river_data(j,3)-50
           end if
        end do


        tdh%v_param(1) = 20.133d0
        call route_processing_build_hist(riv, tdh)

        !filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, write_precision
        call write_ascii_grid('hist_v20.txt', tdh%route_hist_table, &
            size(tdh%route_hist_table, 2), size(tdh%route_hist_table, 1) , &
            0.0d0, 0.0d0, 0.0d0, 0.0d0, 6)


        tdh%v_param(1) = 60.133d0
        call route_processing_build_hist(riv, tdh)

        !filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, write_precision
        call write_ascii_grid('hist_v60.txt', tdh%route_hist_table, &
            size(tdh%route_hist_table, 2), size(tdh%route_hist_table, 1) , &
            0.0d0, 0.0d0, 0.0d0, 0.0d0, 6)


        tdh%v_param(1) = 80.133d0
        call route_processing_build_hist(riv, tdh)

        !filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, write_precision
        call write_ascii_grid('hist_v80.txt', tdh%route_hist_table, &
            size(tdh%route_hist_table, 2), size(tdh%route_hist_table, 1) , &
            0.0d0, 0.0d0, 0.0d0, 0.0d0, 6)

        tdh%v_param(1) = 100.133d0
        call route_processing_build_hist(riv, tdh)

        !filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, write_precision
        call write_ascii_grid('hist_v100.txt', tdh%route_hist_table, &
            size(tdh%route_hist_table, 2), size(tdh%route_hist_table, 1) , &
            0.0d0, 0.0d0, 0.0d0, 0.0d0, 6)


    end subroutine build_hist_test


    subroutine shift_flow_test()
        implicit none
        ! 4 timesteps
        ! 3 nodes
        integer flow_matrix(3, 4)
        integer route_step
        integer point_count
        integer i, j

        point_count = size(flow_matrix,1) ! rows
        route_step = size(flow_matrix,2) ! columns

        print*, 'route_step =', route_step
        print*, 'point_count =', point_count

        print *, 'start with trace on right column'

        flow_matrix(:,:) = 1
        flow_matrix(:,route_step) = 9

        do i=1,point_count
            print *, flow_matrix(i,:)
        end do

        do j=1,4
            print *, 'shift left'
            ! shift the matrix by one timestep
            ! this in effect moves all flows downstream by one timestep
            flow_matrix(:,1) = 0
            flow_matrix = cshift(flow_matrix, 1, 2)
            do i=1,point_count
                print *, flow_matrix(i,:)
            end do

        end do

    end subroutine shift_flow_test

    subroutine route_processing_test()

        call add_accumulation_test()

        call build_hist_test()

    end subroutine route_processing_test

    subroutine route_process_run_timeseries(riv, tdh, &
        node_to_flow_mapping, &
        in_timestep_count, &
        timeseries_flow_input, &
        timeseries_flow_output)
        use dta_riv_tree_node
        implicit none
        type(route_river_info_type) riv
        type(route_time_delay_hist_type) tdh
        ! number of rows to use within timeseries_flow_input
        integer :: in_timestep_count
        ! each row maps to a row in timeseries_flow_input
        ! the value is an index to a node in the riv%node_list
        integer :: node_to_flow_mapping(:)
        double precision :: timeseries_flow_input(:,:)
        double precision :: timeseries_flow_output(:,:)

        ! locals

        integer :: tt

        !% matrix to put the total flow sum for each point
        ! same dimensions as timeseries_flow_input
        timeseries_flow_output(:,:) = 0

        !% add a trace
            !%for ii=1:size(flow_data,2)
            !%route_flow_matrix(:,size(route_flow_matrix,2)) =1;
            !%end
        !%imagesc(flow_data);

        if(size(timeseries_flow_output(1,:))/=size(riv%flow_output_indexes)) then
            print*,'route_process_run_timeseries: error size mismatch'
            print*,'timeseries_flow_output: ', size(timeseries_flow_output(1,:))
            print*,'flow_output_indexes',size(riv%flow_output_indexes)
            stop
        endif

        do tt=1,in_timestep_count

            call route_process_run_step(riv, tdh, &
                node_to_flow_mapping, &
                timeseries_flow_input(tt,:), &
                timeseries_flow_output(tt,:))

        end do

    end subroutine route_process_run_timeseries


    subroutine route_process_run_step(riv, tdh, &
        node_to_flow_mapping, &
        step_flow_input, &
        step_flow_output)

        use dta_riv_tree_node
        implicit none
        type(route_river_info_type) riv
        type(route_time_delay_hist_type) tdh
        ! number of rows to use within timeseries_flow_input
        ! each row maps to a row in timeseries_flow_input
        ! the value is an index to a node in the riv%node_list
        integer :: node_to_flow_mapping(:)
        double precision, intent(in) :: step_flow_input(:)
        double precision, intent(out) :: step_flow_output(:)

        ! locals
        integer :: i, j
        double precision :: q_sum, q_value
        integer :: node_index
        integer :: sum_index
        integer :: upstream_index
        integer :: full_route

        ! Flag to enable full routing (i.e. accumulate all the river branches.
        !  You will want this flag set to zero is you only want the routing per reach
        ! i.e. LISFLOOD simulation - GC to code into dynaTOP.

        full_route = 1

        if(size(step_flow_output)/=size(riv%flow_output_indexes)) then
            print *, size(step_flow_output)
            print *, size(riv%flow_output_indexes)
            print*,'route_process_run_step: error size mismatch'
            stop
        endif


        do i=1,size(node_to_flow_mapping)
            ! read from the timeseries_flow_input and write to the mapped node_index row

            node_index = node_to_flow_mapping(i)

            ! multiply the flow input by the route_hist_table to distribute
            ! the flow across the river reach
            tdh%flow_matrix(node_index,:) = tdh%flow_matrix(node_index,:) &
                + tdh%route_hist_table(node_index,:) &
                * step_flow_input(i)
        end do

        ! sum up the flows for the selected indexes
        ! selected output indexes are stored in flow_output_indexes
        do i=1,size(riv%flow_output_indexes)
            node_index = riv%flow_output_indexes(i)

            !print *, 'node:', riv%node_list(node_index)%gauge_id, riv%node_list(node_index)%upstream_tree_count

            q_sum = 0.0d0
            ! sum_index is the point in time when flow passes this point
            ! need to check which rivers to use from the node's upstream_tree
            ! the read index will be the same for each row
            sum_index = tdh%node_hist_indexes(node_index,1)

            ! add the flow from this reach
            q_value = tdh%flow_matrix(node_index,sum_index)
            q_sum = q_sum + q_value

            ! add flow for full river network - todo don't perform this step if reach only mode
            if (full_route.eq.1) then

                do j=1,riv%node_list(node_index)%upstream_tree_count
                    upstream_index = riv%node_list(node_index)%upstream_tree_indexes(j)

                    !row flow_matrix row is a tributary, add it to the total
                    q_value = tdh%flow_matrix(upstream_index,sum_index)
                    q_sum = q_sum + q_value
                end do

            end if

            step_flow_output(i) = q_sum
        end do

        ! debugging - print out each timestep
        !        print*, 'step -------------'
        !        do i=1,size(tdh%flow_matrix,2)
        !        print*,tdh%flow_matrix(:,i)
        !        end do

        ! set the downstream flow to 0 (cshift wraps the first colum around)
        tdh%flow_matrix(:,1) = 0
        ! shift the matrix by one timestep
        ! this in effect moves all flows downstream by one timestep
        !print*, 'flow_matrix', size(tdh%flow_matrix,1), size(tdh%flow_matrix,2)
        tdh%flow_matrix = cshift(tdh%flow_matrix, 1, 2)

    end subroutine  route_process_run_step

end module dta_route_processing
