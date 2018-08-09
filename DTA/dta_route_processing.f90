! used in route_run_standalone and from within dynatop2001
! 1. Reads in river
!% Toby Dunne
!% April 2016
module dta_route_processing
    use dta_riv_tree_node
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
        integer :: i

        riv%is_enabled = .true.

        flow_point_file = flow_conn_file(1:len_trim(flow_conn_file)-14)//'_flow_point.txt'

        print *, 'read flow tree:', &
            trim(flow_conn_file), &
            trim(flow_point_file)

        call riv_tree_read(riv%node_list, &
            0, 0.0d0, 0.0d0, 0.0d0, &
            flow_conn_file, flow_point_file)



        print *, 'read river data:', trim(river_data_file)
        call read_numeric_list(river_data_file, 6, 1, riv%river_data)
        print *, 'river cell count: ', size(riv%river_data,1)


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

        ! locals
        ! _flow_point.txt
        !character(1024) :: flow_point_file
        integer :: i

        riv%is_enabled = .true.

        print *, 'reading flow tree'

        call riv_tree_read_dyna(riv%node_list, &
            0, 0.0d0, 0.0d0, 0.0d0)

        ! File ID is the same as file ID given in dyna_project.f90 GC 13/06
        !print *, 'read river data:', trim(river_data_file)
        call read_numeric_list_fid(502, 6, 1, riv%river_data)
        print *, 'river cell count: ', size(riv%river_data,1)

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

    end subroutine route_processing_read_info_dyna

    subroutine route_processing_build_hist(riv,   tdh    )
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
        double precision :: full_delay
        double precision :: max_dist
        double precision :: reach_delay
        double precision :: max_elev, min_elev, slope
        double precision :: reach_delay_list(size(riv%node_list))
        double precision :: reach_length_list(size(riv%node_list))
        integer :: riv_start_i, riv_end_i
        integer :: hist_start_i, hist_end_i
        integer :: total_timesteps
        double precision, allocatable :: delay_hist(:)
        integer :: bin_count, target_bin
        double precision :: bin_width, dist_value
        integer :: sea_count

        if(allocated(tdh%node_hist_indexes).eqv..false.) then
            allocate(tdh%node_hist_indexes(size(riv%node_list),2))
        endif

        sea_count = 0
        do i=1,size(riv%node_list)
            !% node type 1 is the outlet - start counting from outlets
            if(riv%node_list(i)%downstream_index == 0) then
            !if(riv%node_list(i)%node_type == NODE_TYPE_SEA) then
                call riv_node_river_delay(riv%node_list, i, 0.0d0, tdh%v_mode, tdh%v_param)
                sea_count = sea_count + 1
            endif
        end do

        if(sea_count == 0) then
           print *, 'Could not find any stream outlets'
           stop
        endif

        !% this is the total maximum time to the furtherest upstream river point.
        !% equal to river_point_max_delay plus max(reach_delay)
        headwater_max_delay = 1
        !% this is the maximum distance to any gauge point
        !% this will determine the warm-up time (the first time that any water can
        !% reach the outlet from the furtherest input point)
        river_point_max_delay = 1

        print *, 'river cells: ', size(riv%river_data,2)

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
                reach_length_list(i) = 0
                reach_delay_list(i) = 0
                tdh%node_hist_indexes(i,1) = 0
                tdh%node_hist_indexes(i,2) = 0
                cycle
            endif

            riv_data_indexes(i,1) = riv_start_i
            riv_data_indexes(i,2) = riv_end_i

            ! column 4 is secion distance
            max_dist = maxval(riv%river_data(riv_start_i:riv_end_i,4))

            if(tdh%v_mode == V_MODE_CONSTANT) then
                !% param is velocity
                ! v_param(1) unit is metres travelled per timestep
                ! v = d/t
                ! t = d/v
                reach_delay = max_dist / tdh%v_param(1)
            elseif (tdh%v_mode == V_MODE_SLOPE) then
                ! slope delay
                ! get a reach_delay based on total slope
                ! when calculating histogram, and allocate by percentage dist for each cell
                 ! v_param(1) when multiplied by slope should give metres travelled per timestep

                ! column 5 is elevation
                max_elev = maxval(riv%river_data(riv_start_i:riv_end_i,5))
                min_elev = minval(riv%river_data(riv_start_i:riv_end_i,5))

                slope = max_elev - min_elev / max_dist
                reach_delay = max_dist / (tdh%v_param(1) * slope)
            else
                print *, 'unknown velocity mode'
                stop
            endif

            reach_length_list(i) = max_dist ! used to distribute percentage of delay
            reach_delay_list(i) = reach_delay

            !% full distance from furtherest river cell to the outlet
            full_delay = reach_delay + riv%node_list(i)%total_downstream_delay

            tdh%node_hist_indexes(i,1) = floor(riv%node_list(i)%total_downstream_delay/tdh%timestep) + 1
            tdh%node_hist_indexes(i,2) = ceiling(full_delay/tdh%timestep) + 1

            if(riv%node_list(i)%total_downstream_delay > river_point_max_delay) then
                ! just recorded to see how many timesteps the flow will reach the bottom
                river_point_max_delay = riv%node_list(i)%total_downstream_delay
            endif
            if(full_delay > headwater_max_delay) then
                headwater_max_delay = full_delay
            endif
        end do
        total_timesteps = ceiling(headwater_max_delay / tdh%timestep) + 1
        !warmup_timesteps = ceiling(river_point_max_delay / tdh%timestep) + 1

        !print *, 'warm-up timesteps = ', warmup_timesteps
        print *, 'flow timesteps = ', total_timesteps

        if(allocated(tdh%route_hist_table)) then
            deallocate(tdh%route_hist_table)
            deallocate(tdh%flow_matrix)
        endif

        ! Allocate the histogram table and the flow matrix
        allocate(tdh%route_hist_table(size(riv%node_list), total_timesteps))
        allocate(tdh%flow_matrix(size(riv%node_list), total_timesteps))

        print*, 'point_count      =', size(tdh%route_hist_table,1) ! rows
        print*, 'route_step_count =', size(tdh%route_hist_table,2) ! columns

        tdh%route_hist_table(:,:) = 0
        tdh%flow_matrix(:,:) = 0

        do i=1,size(riv%node_list)

            !indexes into the
            riv_start_i = riv_data_indexes(i,1)
            riv_end_i = riv_data_indexes(i,2)

            !area_values_section = river_data(start_i:end_i,2)
            !dist_values_section = river_data(start_i:end_i,4)

            max_dist = reach_length_list(i)
            reach_delay = reach_delay_list(i)

            !% Create the histogram for the flow input from dynamic topmodel
            !% flow is distributed along the river reach by the area accumulation
            !% time delay is proportional to distance

            !indexes the time delay cell to write
            hist_start_i = tdh%node_hist_indexes(i,1)
            hist_end_i = tdh%node_hist_indexes(i,2)

            bin_count = (hist_end_i - hist_start_i)+1

            !bin_count = ceiling(reach_delay / timestep)

            ! width in distance units
            bin_width = max_dist / bin_count
            if(int(bin_width) == 0) then
                bin_width = 1
            endif

            allocate(delay_hist(bin_count))
            delay_hist(:) = 0

            ! count accumulation at each distance
            ! this divides the distance into equal bins
            ! - works constant velocity and average slope options
            do j=riv_start_i,riv_end_i
                ! column 4 is river section distance
                dist_value = riv%river_data(j,4)
                target_bin = floor(dist_value / bin_width) + 1
                if(target_bin > bin_count) then ! the value equal to the max dist should be put in the max bin
                    target_bin = bin_count
                endif
                !% add the area contributing to this river distance
                delay_hist(target_bin) = delay_hist(target_bin) + riv%river_data(j,2)
            end do


            ! normalise to 1
            delay_hist = delay_hist / sum(delay_hist)
            !write the histogram to the histogram mask table
            tdh%route_hist_table(i,hist_start_i:hist_end_i) = delay_hist

            deallocate(delay_hist)
        end do

!debugging print the time delay histogram
!        print*, 'tdh -------------'
!        !do i=1,size(riv%node_list)
!        print*,( riv%node_list(i)%gauge_id, i=1,size(riv%node_list) )
!
!        print*,( riv%node_list(i)%total_downstream_delay, i=1,size(riv%node_list) )
!        print*,''
!        !end do
!
!        do i=1,size(tdh%route_hist_table,2)
!        print*,tdh%route_hist_table(:,i)
!        end do

    ! Morpeth test
    !% add an additional row for the ungauged outlet
    !template = zeros(size(flow_data,1),1);
    !flow_data = [template flow_data];


    !timestep_count = size(flow_data,1);
    !point_count = size(flow_data,2);
    !% matrix to put the total flow sum for each point
    !q_total = zeros(timestep_count, point_count);



    end subroutine route_processing_build_hist

    subroutine route_processing_test()
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

        if(size(step_flow_output)/=size(riv%flow_output_indexes)) then
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

            do j=1,riv%node_list(node_index)%upstream_tree_count
                upstream_index = riv%node_list(node_index)%upstream_tree_indexes(j)

                !row flow_matrix row is a tributary, add it to the total
                q_value = tdh%flow_matrix(upstream_index,sum_index)
                q_sum = q_sum + q_value
            end do
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
