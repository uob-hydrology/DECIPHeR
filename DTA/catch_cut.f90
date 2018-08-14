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
program catch_cut
    use dta_utility
    use dta_catch_cut
    use dta_riv_tree_node
    implicit none

    character(1024) in_dem_file

    character(1024) list_input_file

    character(1024) in_river_file

    character(1024) output_dir
    character(1024) tmp_char
    character(1024) result_file
    character(10) output_prefix

    type(time_type) :: start_time, end_time, run_start_time
    integer ioerr

    CHARACTER(len=1024) :: arg

    double precision, allocatable, dimension(:,:) :: dem_grid

    integer(1), allocatable :: flow_dir_grid (:,:)
    double precision, allocatable :: slope_grid (:,:)

    integer :: catchment_id
    integer :: node_type

    integer point_data_col_id
    integer point_data_col_x
    integer point_data_col_y
    integer point_data_col_type ! used when reading defined points

    double precision, allocatable, dimension(:,:) :: point_data

    integer i

    integer :: ncols, nrows
    double precision :: xllcorner, yllcorner, cellsize
    double precision :: double_nodata

    integer row, col
    integer starty, startx

    character(len=30) stop_reason
    integer candidate_out_of_bounds
    double precision dist

    integer x, y

    ! can convert to mask by selecting where value is not 0
    double precision, allocatable, dimension(:,:) :: mask_grid

    ! cropped copy of the mask for writing
    double precision, allocatable, dimension(:,:) :: small_mask_grid
    integer :: small_mask_nrows, small_mask_ncols
    double precision :: small_mask_xll, small_mask_yll

    integer :: miny, maxy, minx, maxx
    integer :: area_count
    double precision :: calc_area_m2
    double precision :: calc_area_km2
    double precision :: outlet_x_east, outlet_y_north

    logical :: input_is_valid
    logical :: write_catchment_mask

    CALL timer_get(run_start_time)

    in_dem_file = ''
    in_river_file = ''
    list_input_file = ''

    output_dir = ''
    tmp_char = ''

    input_is_valid = .true.

    tmp_char = 'DTA_catch_cut.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- catch_cut.f90 ---'
    write(999,*) ''
    print *, '--- Starting catch_cut ---'

    i=0

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, in_dem_file)
        elseif (are_equal(arg, '-points')) then
            CALL get_command_argument(i+1, list_input_file)
        elseif (are_equal(arg, '-out')) then
            CALL get_command_argument(i+1, output_dir)
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_dem_file,'-dem').eqv..false.) then
        input_is_valid = .false.
    endif

    if(check_file_arg(list_input_file,'-points').eqv..false.) then
        input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'command options '
        print *, '-dem <file.asc>    select dem ascii grid file'
        print *, '-points <file.txt> list of outlet points to build catchments from'
        print *, '-out <dir>         directory to write mask (flow length) files'
        print *, ''
        print *, 'Input data description'
        print *, 'dem: sink filled dem (ascii grid)'
        print *, 'points: tab delimited (first line header) containing: '
        print *, '         Node_ID   x_easting   y_northing  Type [subsequent columns ignored]'
        print *, ''
        print *, ' e.g.:'
        print *, ' mkdir masks'
        print *, ' catch_cut.e -dem dem.asc -points dem_flow_point.txt -out masks'
        stop
    endif

    if(len_trim(output_dir) > 0) then
        ! if the user didn't add the trailing / add it now
        if(output_dir(len_trim(output_dir):len_trim(output_dir)) /= '/') then
            output_dir = trim(output_dir) // '/'
        endif
        if(output_dir(1:1) == '~') then
            print *, 'sorry, output dir starting with ~ not supported'
            stop
        endif
    endif


    ! input file is in format
    ! tab separated, 4 columns
    ! input file is in format (tab separates)
    ! X and Y are eastings and northings
    ! Node_Id    X   Y   Type
    point_data_col_id = 1
    point_data_col_x = 2
    point_data_col_y = 3
    point_data_col_type = 4

    print *, 'read point list: ', trim(list_input_file)
    write(999,*) 'read point list: ', trim(list_input_file)
    CALL timer_get(start_time)
    call read_numeric_list(list_input_file, 4, 1, point_data)
    CALL timer_get(end_time)
    !call timer_print('read point list', start_time, end_time)

    print *, 'read dem grid: ', trim(in_dem_file)
    write(999,*) 'read dem grid: ', trim(in_dem_file)
    CALL timer_get(start_time)
    call read_ascii_grid(in_dem_file, dem_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read dem grid', start_time, end_time)

    allocate (flow_dir_grid(nrows,ncols))
    allocate (slope_grid(nrows,ncols))

    print *, 'calc_flow_direction'
    write(999,*) ''
    write(999,*) 'calculating flow direction'
    CALL timer_get(start_time)
    call calc_flow_direction(nrows, ncols, dem_grid, cellsize, flow_dir_grid, slope_grid)
    CALL timer_get(end_time)
    !call timer_print('calc_flow_direction', start_time, end_time)

    ! no longer need dem
    deallocate(dem_grid)


    !    tmp_char = in_dem_file(1:len_trim(in_dem_file)-4)//'_sfd_slope.asc'
    !
    !    print *, 'write:', trim(tmp_char)
    !    call write_ascii_grid(tmp_char, slope_grid, &
    !        ncols, nrows, &
    !        xllcorner, yllcorner, cellsize, -9999.0d0, 5)
    !
    deallocate(slope_grid)

    ! array to write the catchment mask (flow lengths)
    allocate (mask_grid(nrows,ncols))

    ! no reference area, use this column to print type
    ! sea/river/reach

    output_prefix = 'other_'
    result_file = trim(output_dir) // 'mask_info.txt'

    print *, 'results file: ', trim(result_file)
    write(999,*) 'results file: ', trim(result_file)

    ! result file
    ! this file contains all the catchment outlet points on the river cells
    open(100, file = result_file, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(result_file)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write (100, 97) 'catchment', tab,&
        'x_easting', tab, &
        'y_northing', tab,&
        'node_type', tab, &
        'calc_area_m2', tab, &
        'out_of_bound',tab, &
        'reason'

    ! write header to terminal
    !write (*, 97) 'catchment', tab, &
    !    'x_easting', tab, &
    !    'y_northing', tab, &
    !    'node_type', tab, &
    !    'calc_area_m2', tab, &
    !    'out_of_bound', tab, &
    !    'reason'


    ! format labels are a bit cryptic
    ! 97 is 13 strings (9 headers and 8 tabs)
    ! 98 is 1 int, 2 (tab floats), tab, int, tab, float, tab, int, tab, string
97  format ( 13A )
98  format ( I0,2(A,F0.1),A,I0,A,F0.1,A,I0,A,A)

    do i=1,size(point_data,1)
        catchment_id = nint(point_data(i, point_data_col_id))

        call NorthingEastingToRowCol( point_data(i,point_data_col_y), &
            point_data(i,point_data_col_x), &
            nrows, xllcorner, yllcorner, cellsize, row, col);

        starty = row
        startx = col

        ! only process outlet points inside the dem bounds
        if(starty > 1 .and. starty < nrows &
            .and. &
            startx > 1 .and. startx < ncols) then

            write_catchment_mask = .true.
            x = startx
            y = starty

            node_type = nint(point_data(i,point_data_col_type))

            if(node_type == NODE_TYPE_GAUGE) then
                stop_reason = 'input_point_gauge'
                output_prefix = 'gauge_'
            elseif(node_type == NODE_TYPE_SEA) then
                stop_reason = 'input_point_sea'
                output_prefix = 'sea_'
            elseif(node_type == NODE_TYPE_RIVER) then
                stop_reason = 'input_point_reach'
                output_prefix = 'reach_'
            else
                stop_reason = 'input_point_other'
                output_prefix = 'other_'
            endif

            calc_area_m2 = 0
            calc_area_km2 = 0

            if((write_catchment_mask.eqv..true.) .and. x > 0 .and. y > 0) then

                call RowColToNorthingEasting(y, x, &
                    nrows, xllcorner, yllcorner, cellsize, &
                    outlet_y_north, outlet_x_east, .true.)

                dist = sqrt(real((starty-y)*(starty-y)+(startx-x)*(startx-x)))

                !call calc_catch_cut_mfd( nrows, ncols, dem_grid, &
                !    y, x, &
                !    mask_grid, miny, maxy, minx, maxx, area_count)
                call calc_catch_cut_flow_dir( nrows, ncols, flow_dir_grid, &
                    y, x, &
                    mask_grid, miny, maxy, minx, maxx, area_count)

                calc_area_m2 = area_count * (cellsize*cellsize)  ! area m^2
                calc_area_km2 = calc_area_m2 / (1000*1000)       ! area km^2

                write(tmp_char, '(A,I0,''.asc'')') &
                    trim(output_prefix), &
                    catchment_id

                tmp_char = trim(output_dir) // tmp_char

                ! begin save the small catchment mask

                ! set flag to indicate catchment is out of bounds
                if  (miny == 1 .or. maxy == nrows &
                    .or. minx == 1 .or. maxx == ncols) then
                    candidate_out_of_bounds = 1
                else
                    candidate_out_of_bounds = 0
                endif

                !% expand the bounday to include a 10 cell border
                miny = max(1,     miny - 10)
                maxy = min(nrows, maxy + 10)
                minx = max(1,     minx - 10)
                maxx = min(ncols, maxx + 10)

                small_mask_nrows = maxy - miny + 1
                small_mask_ncols = maxx - minx + 1

                allocate(small_mask_grid(small_mask_nrows, small_mask_ncols))

                small_mask_grid = mask_grid(miny:maxy,minx:maxx)

                ! convert unit length to flow to distance in the unit of cellsize
                small_mask_grid = small_mask_grid * cellsize

                ! convert the maxy, minx from the location of the cropped
                ! mask into northing easting lower left
                call RowColToNorthingEasting(maxy, minx, &
                    nrows, xllcorner, yllcorner, cellsize, &
                    small_mask_yll, small_mask_xll, .false. )

                call write_ascii_grid(tmp_char, small_mask_grid, &
                    small_mask_ncols, small_mask_nrows, &
                    small_mask_xll, small_mask_yll, cellsize, 0.0d0, 1);

                deallocate(small_mask_grid)
            ! end save the small catchment mask
            else
                candidate_out_of_bounds = 0
            endif
            !write (*, 98) catchment_id, tab, &
            !    outlet_x_east, tab, &
            !    outlet_y_north, tab, &
            !    node_type, tab, &
            !    calc_area_m2, tab,&
            !    candidate_out_of_bounds, tab,&
            !    trim(stop_reason)

            write (100, 98) catchment_id, tab, &
                outlet_x_east, tab, &
                outlet_y_north, tab, &
                node_type, tab, &
                calc_area_m2, tab,&
                candidate_out_of_bounds, tab,&
                trim(stop_reason)

            flush (100)
        else
            print *, 'point out of dem bounds'
        endif
    enddo
    close (100)

    CALL timer_get(end_time)
    !call timer_print('catch_cut', run_start_time, end_time)
    write(999,*) ''
    write(999,*) 'Successfully finished catch_cut.f90'
    print *, '--- Finished catch_cut ---'
    close(999)

    stop


end program catch_cut





