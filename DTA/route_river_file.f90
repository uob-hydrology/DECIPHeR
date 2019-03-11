!% Toby Dunne
!% Apr 2016
! extract the information needed to do the river routing
! basically takes the grid river data and puts it into list form ready to create histograms
program route_river_file
    use dta_utility
    use dta_routing_file
    use dta_riv_tree_node
    use dta_rivers
    implicit none

    character(1024) in_dem_file
    character(1024) in_riv_id_file

    character(1024) flow_conn_file
    character(1024) flow_point_file
    character(1024) in_area_file
    character(1024) in_riv_dist_file

    character(1024) dem_output_prefix

    type(time_type) :: start_time, end_time, run_start_time

    CHARACTER(len=1024) :: arg
    integer :: i, ioerr

    integer, allocatable, dimension(:,:) :: riv_id_grid
    double precision, allocatable, dimension(:,:) :: area_grid, dem_grid
    type(riv_tree_node), allocatable :: node_list(:)
    type(point_list_type), allocatable, dimension(:) :: river_point_lists
    integer :: total_river_cells


    double precision, allocatable, dimension(:,:) :: riv_dist_grid
    logical, allocatable, dimension(:,:) :: riv_mask_grid

    integer :: sea_outlet_count
    type(point_type), allocatable :: input_sea_outlets(:)

    integer :: ncols, nrows
    double precision :: xllcorner, yllcorner, cellsize
    double precision :: double_nodata
    integer :: int_nodata
    character(1024) :: tmp_char
    !% riv_id, area, dist, section_dist, slope, elevation
    double precision, allocatable, dimension(:,:) :: river_data
    character(64) :: col_headers(8)

    logical :: input_is_valid

    col_headers(1) = 'riv_id'
    col_headers(2) = 'area'
    col_headers(3) = 'dist'
    col_headers(4) = 'section_dist'
    col_headers(5) = 'elevation'
    col_headers(6) = 'slope'
    col_headers(7) = 'ds_dist'
    col_headers(8) = 'ds_elevation'


    CALL timer_get(run_start_time)

    in_dem_file= ''
    in_riv_id_file= ''

    flow_conn_file= ''
    flow_point_file= ''
    in_area_file= ''
    in_riv_dist_file = ''


    !tmp_char = ''

    input_is_valid = .true.

    tmp_char = 'DTA_route_river_file.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- route_river_file.f90 ---'
    write(999,*) ''
    print *, '--- Starting route_river_file ---'

    i=0

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, in_dem_file)
        elseif (are_equal(arg, '-river')) then
            CALL get_command_argument(i+1, in_riv_id_file)
        elseif (are_equal(arg, '-tree')) then ! optional
            CALL get_command_argument(i+1, flow_conn_file)
        elseif (are_equal(arg, '-area')) then ! optional
            CALL get_command_argument(i+1, in_area_file)
        elseif (are_equal(arg, '-riv_dist')) then ! optional
            CALL get_command_argument(i+1, in_riv_dist_file)
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_dem_file,'-dem').eqv..false.) then
        input_is_valid = .false.
    endif
    if(check_file_arg(in_riv_id_file,'-river').eqv..false.) then
        input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'command options '
        print *, '-dem <dem.asc>   select dem ascii grid file'
        print *, '-river <dem_riv_id_check.asc> labelled river ascii grid file'
        print *, ''
        print *, ' optional specify inputs (by default names from dem prefix)'
        print *, '-tree <dem_flow_conn.txt> override the routing tree file'
        print *, '         note: xxx_flow_point.txt will also be read '
        print *, '-area <dem_riv_mask_area.asc> from atb_wfp.e processed as with river file'
        print *, '-riv_dist <dem_riv_dist.asc> river labelled with distance'
        stop
    endif

    ! output based on dem filename
    dem_output_prefix = in_dem_file(1:len_trim(in_dem_file)-4)
    if (len_trim(flow_conn_file) == 0) then
        flow_conn_file = trim(dem_output_prefix)//'_flow_conn.txt'
    endif
    flow_point_file = flow_conn_file(1:len_trim(flow_conn_file)-14)//'_flow_point.txt'
    if (len_trim(in_area_file) == 0) then
        in_area_file = trim(dem_output_prefix)//'_riv_mask_area.asc'
    endif

    if (len_trim(in_riv_dist_file) == 0) then
        in_riv_dist_file = trim(dem_output_prefix)//'_riv_dist.asc'
    endif

    print *, 'dem file: ', trim(in_dem_file)
    print *, 'river id file: ', trim(in_riv_id_file)
    print *, 'flow connection file: ', trim(flow_conn_file)
    print *, 'flow point file: ', trim(flow_point_file)
    print *, 'area file: ', trim(in_area_file)
    print *, 'river dist file: ', trim(in_riv_dist_file)

    write(999,*) 'dem file:', trim(in_dem_file)
    write(999,*) 'river id file: ', trim(in_riv_id_file)
    write(999,*) 'flow connection file: ', trim(flow_conn_file)
    write(999,*) 'flow point file: ', trim(flow_point_file)
    write(999,*) 'area file: ', trim(in_area_file)
    write(999,*) 'river dist file: ', trim(in_riv_dist_file)
    write(999,*) ''

    if(file_exists(flow_conn_file).eqv..false.) then
        print*,'Missing: ', trim(flow_conn_file)
        stop
    endif
    if(file_exists(flow_point_file).eqv..false.) then
        print*,'Missing: ', trim(flow_point_file)
        stop
    endif
    if(file_exists(in_area_file).eqv..false.) then
        print*,'Missing: ', trim(in_area_file)
        stop
    endif
    if(file_exists(in_riv_dist_file).eqv..false.) then
        print*,'Missing: ', trim(in_riv_dist_file)
        stop
    endif

    !print *, 'read river id grid: ', trim(in_riv_id_file)
    write(999,*) 'read river id grid: ', trim(in_riv_id_file)
    CALL timer_get(start_time)
    call read_ascii_grid_int(in_riv_id_file, riv_id_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, int_nodata)
    CALL timer_get(end_time)
    !call timer_print('read river grid', start_time, end_time)

    write(999,*) 'read flow tree:', &
        trim(flow_conn_file), '', &
        trim(flow_point_file)

    call riv_tree_read(node_list, &
        nrows, xllcorner, yllcorner, cellsize, &
        flow_conn_file, flow_point_file)

    allocate(river_point_lists(size(node_list)))
    !print *, 'routing_file_find_rivers'
    call routing_file_find_rivers(nrows, ncols, node_list, riv_id_grid, &
        river_point_lists, total_river_cells)

    ! calculate river_dist - need to know the outlets
    ! river_dist should be already calculated, however quicker to recalculate than read in

    ! get the sea outlets from the node_list
    sea_outlet_count = 0
    do i=1,size(node_list)
        if(node_list(i)%node_type == NODE_TYPE_SEA) then
            sea_outlet_count = sea_outlet_count + 1
        endif
    end do
    allocate(input_sea_outlets(sea_outlet_count))

    sea_outlet_count = 0
    do i=1,size(node_list)
        if(node_list(i)%node_type == NODE_TYPE_SEA) then
            sea_outlet_count = sea_outlet_count + 1
            input_sea_outlets(sea_outlet_count)%y =node_list(i)%row
            input_sea_outlets(sea_outlet_count)%x =node_list(i)%col
        endif
    end do

    ! convert riv_id_grid to a mask then deallocate
    allocate(riv_mask_grid(nrows,ncols))
    riv_mask_grid(:,:) = .false.
    where(riv_id_grid/=0) riv_mask_grid = .true.
    deallocate(riv_id_grid)

!    allocate(riv_dist_grid(nrows,ncols))
!
!    ! dem_grid not used as the sea_outlets are pre-computed
!    allocate(dem_grid(0,0))
!    call river_find_dist(ncols, nrows, &
!        riv_mask_grid, &
!        dem_grid, &
!        input_sea_outlets, &
!        riv_dist_grid, sea_outlets)
!    deallocate(riv_mask_grid)
!
!    ! convert distance to metres
!    riv_dist_grid = riv_dist_grid * cellsize

    ! end of river distance calculation
    print *, 'read river dist grid: ', trim(in_riv_dist_file)
    CALL timer_get(start_time)
    call read_ascii_grid(in_riv_dist_file, riv_dist_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read river dist grid', start_time, end_time)

    write(999,*) 'read area grid: ', trim(in_area_file)
    CALL timer_get(start_time)
    call read_ascii_grid(in_area_file, area_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read dem grid', start_time, end_time)

    write(999,*) 'read dem: ', trim(in_dem_file)
    CALL timer_get(start_time)
    call read_ascii_grid(in_dem_file, dem_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)

    !print *, 'routing_file_grids_to_list'
    call routing_file_grids_to_list(nrows, ncols, cellsize, node_list, &
        river_point_lists, total_river_cells, &
        riv_dist_grid, area_grid, dem_grid, river_data)

    tmp_char = trim(dem_output_prefix)//'_river_data.txt'

    print *, 'write ', trim(tmp_char)
    write(999,*) 'write ', trim(tmp_char)
    call write_numeric_list(tmp_char, col_headers, river_data, 6)

    CALL timer_get(end_time)


    write(999,*) ''
    write(999,*) 'Successfully finished route_river_file.f90'
    print *, '--- Finished route_river_file ---'
    close(999)

    !call timer_print('route_river', run_start_time, end_time)

    stop


end program route_river_file





