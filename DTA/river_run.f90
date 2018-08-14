!% Toby Dunne
!% April 2016
program river_run
    use dta_utility
    use dta_rivers
    implicit none

    CHARACTER(len=1024) :: arg
    logical :: input_is_valid
    integer :: i
    character(1024) in_dem_file
    character(1024) hw_input_file

    character(1024) tmp_char

    double precision, allocatable, dimension(:,:) :: dem_grid
    integer, allocatable, dimension(:,:) :: riv_label_grid

    double precision, allocatable, dimension(:,:) :: headwater_input
    type(point_type), allocatable, dimension(:) :: headwater_points

    integer :: ncols, nrows, ioerr
    double precision :: xllcorner, yllcorner, cellsize
    double precision :: double_nodata
    type(time_type) :: start_time, end_time, run_start_time

    CALL timer_get(run_start_time)

    in_dem_file = ''
    hw_input_file = ''

    input_is_valid = .true.

    i = 0

    tmp_char = 'DTA_river_run.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- river_run.f90 ---'
    write(999,*) ''
    print *, '--- Starting river_run ---'

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, in_dem_file)
        elseif (are_equal(arg, '-headwater')) then
            CALL get_command_argument(i+1, hw_input_file)
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_dem_file,'-dem').eqv..false.) then
        input_is_valid = .false.
    endif
    if(check_file_arg(hw_input_file,'-headwater').eqv..false.) then
        input_is_valid = .false.
    endif


    if(input_is_valid .eqv. .false.) then
        print *, 'command options '
        print *, '-dem <file.asc>   select dem sinkfilled ascii grid file'
        print *, '-headwater <file.txt> list of headwaters id, easting northing'
        stop
    endif

    write(999,*) 'Command Options Read In'
    write(999,*) 'in_dem: ', trim(in_dem_file)
    write(999,*) 'in_headwater: ', trim(hw_input_file)
    write(999,*) ''

    print *, 'in_dem: ', trim(in_dem_file)
    print *, 'in_headwater: ', trim(hw_input_file)

    print *, 'read dem grid: ', trim(in_dem_file)
    write(999,*) 'read dem grid: ', trim(in_dem_file)

    CALL timer_get(start_time)
    call read_ascii_grid(in_dem_file, dem_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read dem grid', start_time, end_time)

    print *, 'read headwater list: ', trim(hw_input_file)
    write(999,*) 'read headwater list: ', trim(hw_input_file)
    call read_numeric_list(hw_input_file, 2, 1, headwater_input)

    ! Convert northing easting to row col (point)
    allocate (headwater_points(size(headwater_input,1)))
    do i=1,size(headwater_input,1)

        call NorthingEastingToRowCol(headwater_input(i,2), headwater_input(i,1), &
            nrows, xllcorner, yllcorner, cellsize, &
            headwater_points(i)%y, headwater_points(i)%x)

    end do
    deallocate(headwater_input)

    allocate(riv_label_grid(nrows, ncols))

    CALL timer_get(start_time)

    print *, 'Creating river layer from ', size(headwater_input, 1), ' headwaters'
    write(999,*) 'Creating river layer from ', size(headwater_input, 1), ' headwaters'

    call river_single_flow(nrows, ncols, dem_grid, headwater_points, riv_label_grid, cellsize)
    CALL timer_get(end_time)

        tmp_char = in_dem_file(1:len_trim(in_dem_file)-4)//'_riv.asc'

        print*,'Write: ', trim(tmp_char)
        write(999,*) 'Write: ', trim(tmp_char)
        call write_ascii_grid_int(tmp_char, riv_label_grid, &
            ncols, nrows, &
            xllcorner, yllcorner, cellsize, 0);

        CALL timer_get(end_time)
        !call timer_print('river_run', run_start_time, end_time)

!allocate(riv_mask_grid(nrows,ncols))

!    riv_mask_grid = .false.
!    where(riv_label_grid > 0) riv_mask_grid = .true.

!    deallocate(riv_label_grid)

    print *, '--- Finished river_run ---'
    write(999,*) ''
    write(999,*) 'Successfully finished river_run.f90'
    close(999)

    stop


end program river_run
