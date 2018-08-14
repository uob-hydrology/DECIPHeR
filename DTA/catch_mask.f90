!% Catchment cutting
!% input: List of station id, easting, northing and reference area
!% river file with each river branch having a unique id (see headwater_run.m)
!%   used to find candidate starting points
!% sinkfilled dem - with no flat areas used
!%   used to walk upstream to find all in catchment cells
!%
!% output series of grids representing the in catchment area,
!% each cell contains the flow length to outlet in m (same unit of cellsize/dx)
!% - expanded search if no decent matches
!%
!% Toby Dunne
!% Apr 2016
program catch_mask
    use dta_utility
    use dta_catch_mask
    use dta_riv_tree_node
    implicit none

    type(time_type) :: end_time, run_start_time
    character(1024) in_dem_file
    character(1024) input_mask_dirs
    character(1024) flow_conn_file
    character(1024) flow_point_file
    integer, allocatable :: manual_outlet_id_list(:)
    character(1024), allocatable:: manual_outlet_id_list_char(:)
    logical, allocatable:: manual_start(:)
    integer :: i, j, ioerr, fp
    integer :: start_count
    character(1024) arg, tmp_char
    logical input_is_valid

    character(1024), allocatable :: mask_dirs(:)

    logical is_start_catchment
    integer :: ncols, nrows
    double precision :: xllcorner, yllcorner, cellsize
    double precision :: nodata
    integer, allocatable, dimension(:,:) :: mask_grid

    type(riv_tree_node), allocatable :: node_list(:)

    start_count = 0
    fp=999
    !call utility_test

    in_dem_file = ''
    input_mask_dirs = ''
    flow_conn_file = ''

    input_is_valid = .true.

    tmp_char = 'DTA_catch_mask.log'
    open(fp, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(fp,*) '--- catch_mask.f90 ---'
    write(fp,*) ''
    print *, '--- Starting catch_mask ---'

    i=0

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-base')) then
            CALL get_command_argument(i+1, in_dem_file)
        elseif (are_equal(arg, '-mask_dirs')) then
            CALL get_command_argument(i+1, input_mask_dirs)
        elseif (are_equal(arg, '-tree')) then
            CALL get_command_argument(i+1, flow_conn_file)
        elseif (are_equal(arg, '-start_outlet_ids')) then
            CALL get_command_argument(i+1, arg)
            call split_comma(arg, manual_outlet_id_list_char)
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_dem_file,'-dem').eqv..false.) then
        input_is_valid = .false.
    endif
    if (len_trim(input_mask_dirs) == 0) then
        print *, '-mask_dirs not specified'
        input_is_valid = .false.
    endif
    if(check_file_arg(flow_conn_file,'-tree').eqv..false.) then
        input_is_valid = .false.
    endif

    if(allocated(manual_outlet_id_list_char)) then
        print*, 'Manual start points:'
        allocate(manual_outlet_id_list(size(manual_outlet_id_list_char)))
        allocate(manual_start(size(manual_outlet_id_list_char)))
        manual_start(:) = .false.

        do i=1,size(manual_outlet_id_list_char)
            read (manual_outlet_id_list_char(i),*) manual_outlet_id_list(i)
            print*,manual_outlet_id_list(i)
        end do
        deallocate(manual_outlet_id_list_char)
    endif

    if (input_is_valid .eqv. .false.) then
        print *, 'command options'
        print *, '-base <file.asc> select dem ascii grid file (only header used)'
        print *, '-mask_dirs <dir1,dir2> list of dirs to find catchment masks'
        print *, '-tree <xxx_flow_conn.txt> the routing tree file'
        print *, '         note: xxx_flow_point.txt will also be read'
        print *, ' optional:'
        print *, '-start_outlet_ids <id> start building mask from outlet id(s) and all upstream'
        print *, '         the default behaviour (when start_gauge_id not given) is to build from'
        print *, '         all outlet points i.e. points that have no downstream point'
        print *, '         multiple can be given separated by commas'
        print *, ' e.g.:'
        print *, ' catch_mask.e -dem dem.asc -mask_dirs catch3,sea2 -tree dem_flow_conn.txt'
        stop
    endif

    CALL timer_get(run_start_time)

    call split_comma(input_mask_dirs, mask_dirs)

    call read_ascii_grid_header(in_dem_file, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

    allocate(mask_grid(ncols, nrows))
    mask_grid(:,:) = 0

    flow_point_file = flow_conn_file(1:len_trim(flow_conn_file)-14)//'_flow_point.txt'

    print *, 'read flow tree: ', &
        trim(flow_conn_file), &
        '', trim(flow_point_file)
    print *, 'Mask Directory : ', trim(input_mask_dirs)

    write(fp,*) 'Mask Directory : ', trim(input_mask_dirs)
    write(fp,*) 'read flow tree: ', &
        trim(flow_conn_file), &
        '', trim(flow_point_file)
    write(fp,*) ''

    call riv_tree_read(node_list, &
        nrows, xllcorner, yllcorner, cellsize, &
        flow_conn_file, flow_point_file)


    do i=1,size(node_list)

        is_start_catchment = .false.

        if(allocated(manual_outlet_id_list)) then
            do j=1,size(manual_outlet_id_list)
                if(node_list(i)%gauge_id == manual_outlet_id_list(j)) then
                    is_start_catchment = .true.
                    manual_start(j) = .true.
                endif
            end do
        elseif(node_list(i)%downstream_index == 0) then
            !% start from the sea outlet end point, any point with no downstream
            is_start_catchment = .true.
        endif

        if(is_start_catchment) then
            !write(fp,*) 'START: ',node_list(i)%gauge_id,'type: ', node_list(i)%node_type
            start_count = start_count + 1
            call add_catch_mask(nrows, ncols, mask_dirs, &
                node_list, i, &
                mask_grid, xllcorner, yllcorner, cellsize, fp)
        endif

    end do

    if(start_count ==0) then
        print*,'ERROR no start points matched'
        stop
    endif
    if(allocated(manual_outlet_id_list)) then
        do i=1,size(manual_outlet_id_list)
            if(manual_start(i).eqv..false.) then
                write(999,*) 'WARNING: start point not used: ', manual_outlet_id_list(i)
            endif
        end do
    endif

    arg = in_dem_file(1:len_trim(in_dem_file)-4)//'_mask.asc'

    print*,'write: ', trim(arg)
    write(fp,*) ''
    write(fp,*) 'write: ', trim(arg)
    call write_ascii_grid_int(arg, mask_grid, &
        ncols, nrows, &
        xllcorner, yllcorner, cellsize, 0)

    deallocate(node_list)
    deallocate(mask_grid)

    CALL timer_get(end_time)
    !call timer_print('catch_mask', run_start_time, end_time)

    write(fp,*) ''
    write(fp,*) 'Successfully finished catch_mask.f90'
    print *, '--- Finished catch_mask ---'
    close(fp)

    stop


end program catch_mask





