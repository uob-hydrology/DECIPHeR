!% Toby Dunne
!% Apr 2016
program route_run_standalone
    use dta_utility
    use dta_riv_tree_node
    use dta_route_processing
    implicit none

    !from route_river.f90
    character(1024) river_data_file
    ! from route_tree.f90
    character(1024) flow_conn_file
    character(1024) flow_point_file
    character(1024) timeseries_flow_input_file

    type(time_type) :: start_time, end_time, run_start_time

    character(1024) :: tmp_char
    CHARACTER(len=1024) :: arg
    integer :: i

    type(route_river_info_type) :: riv
    type(route_time_delay_hist_type) :: tdh
    !
    double precision, allocatable, dimension(:,:) :: timeseries_flow_input
    !double precision, allocatable, dimension(:,:) :: timeseries_flow_input2

    double precision, allocatable, dimension(:,:) :: timeseries_flow_output

    ! for each column in the flow input file
    ! which river reach node does this add its contribution
    ! in this case we have an extra node which is the sea node which doesn't have a contribution
    integer, allocatable :: node_to_flow_mapping(:)

    logical :: input_is_valid, run_test

    integer :: n_points, n_steps
    ! not used
    double precision:: xllcorner, yllcorner, cellsize, double_nodata

    !    integer test(3,4)
    !    test(:,:) = 0
    !    test(1,:) = 1
    !    test(:,1) = 2
    !    do i=1,size(test,1)
    !        print *, test(i,:)
    !     end do
    !           2           1           1           1
    !           2           0           0           0
    !           2           0           0           0

    !call route_processing_test()
    !stop


    CALL timer_get(run_start_time)

    tdh%timestep = 3600.0d0
    tdh%v_mode = V_MODE_CONSTANT
    allocate(tdh%v_param(1))
    tdh%v_param(1) = 0.133d0

    river_data_file = ''
    flow_conn_file= ''
    flow_point_file= ''
    timeseries_flow_input_file = ''

    input_is_valid = .true.
    run_test = .false.


    i = 0
    print *, '--- route_run_standalone ---'
    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-river_data')) then
            CALL get_command_argument(i+1, river_data_file)
        elseif (are_equal(arg, '-tree')) then
            CALL get_command_argument(i+1, flow_conn_file)
        elseif (are_equal(arg, '-timeseries_in')) then
            CALL get_command_argument(i+1, timeseries_flow_input_file)
        elseif (are_equal(arg, '-test')) then
            run_test = .true.
        endif
        i = i + 1
    enddo

    if(check_file_arg(river_data_file,'-river_data').eqv..false.) then
        input_is_valid = .false.
    endif
    if(check_file_arg(flow_conn_file,'-tree').eqv..false.) then
        input_is_valid = .false.
    endif

    if(run_test.eqv. .true.) then
        call route_processing_test()
        stop
    endif

    if(check_file_arg(timeseries_flow_input_file,'-timeseries_in').eqv..false.) then
        input_is_valid = .false.
    endif


    if(input_is_valid .eqv. .false.) then
        print *, ' command options '
        print *, '-river_data <xxx_river_data.txt> river data (from route_river_file.e)'
        print *, '-tree <xxx_flow_conn.txt> the routing tree file'
        print *, '         note: xxx_flow_point.txt will also be read '
        print *, '-timeseries_in'
        stop
    endif

    print *, 'river data file:', trim(river_data_file)
    !print *, 'flow connection file:', trim(flow_conn_file)
    !print *, 'flow point file:', trim(flow_point_file)

    print *, 'flow timeseries input: ', trim(timeseries_flow_input_file)

    call route_processing_read_info(river_data_file, flow_conn_file, riv)

    CALL timer_get(start_time)
    print *, 'route_processing_build_hist'

    call route_processing_build_hist(riv, tdh)

    CALL timer_get(end_time)
    call timer_print('route_processing_build_hist', start_time, end_time)

    ! reading the timeseries data in ascii grid format
    ! todo replace with reading from dynamic topmode format
    call read_ascii_grid(timeseries_flow_input_file, timeseries_flow_input, &
        n_points, n_steps, xllcorner, yllcorner, cellsize, double_nodata)

    ! each input point will map to a node
    allocate(node_to_flow_mapping(n_points))

    !timeseries_flow_input_file
    if(n_points /= size(riv%node_list)) then
        ! one column too short, this is because input doesn't have data for the sea outlet
        ! add a zero row
        if((n_points + 1) == size(riv%node_list)) then

            do i=1,n_points
                node_to_flow_mapping(i) = i+1
            end do

            !!call move_alloc(timeseries_flow_input, timeseries_flow_input2)
            !!allocate(timeseries_flow_input(n_steps, size(riv%node_list)))

            !! set first column to 0
            !!timeseries_flow_input(:,1) = 0
            !! copy back the flow data
            !!timeseries_flow_input(:,2:size(riv%node_list)) = timeseries_flow_input2
            !deallocate(timeseries_flow_input2)
        else

            print *, 'timeseries_flow_input does not match the _flow_conn.txt'
            print *, ' timeseries data must have one column for each node in the _flow_conn.txt'
            print *, ' each row in the timeseries data is one timestep'
            stop

        endif
    else
        ! sizes match do 1 to 1 mapping
        print *, 'flow file points matches nodes'
        do i=1,n_points
            node_to_flow_mapping(i) = i
        end do
    endif

    ! allocate storage for flow output
    allocate(timeseries_flow_output(size(timeseries_flow_input,1), riv%flow_output_count))
    timeseries_flow_output(:,:) = 0

    CALL timer_get(start_time)
    print *, 'route_process_run_timeseries'

    call route_process_run_timeseries(riv, &
        tdh, &
        node_to_flow_mapping, &
        size(timeseries_flow_input,1), &
        timeseries_flow_input, &
        timeseries_flow_output)

    CALL timer_get(end_time)
    call timer_print('route_process_run_timeseries', start_time, end_time)

    tmp_char = timeseries_flow_input_file(1:len_trim(timeseries_flow_input_file)-4)//'_routed.txt'

    print *, 'write output: ', trim(tmp_char)
    print *, 'Velocity', tdh%v_param(1)

    !timeseries_flow_output = timeseries_flow_output * 1000000


    call write_numeric_list(tmp_char, riv%flow_output_headers, timeseries_flow_output, 6)


    CALL timer_get(end_time)
    call timer_print('route_run_standalone', run_start_time, end_time)

    stop



end program route_run_standalone
