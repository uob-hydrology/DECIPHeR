!% Toby Dunne
!% Aug 2016
program river_thresholds
    use dta_utility
    implicit none

    character(1024) in_area_file
    character(1024) in_atb_file
    character(1024) out_riv_file

    type(time_type) :: start_time, end_time, run_start_time

    CHARACTER(len=1024) :: arg, tmp_char
    logical :: input_is_valid
    integer :: i, ioerr

    double precision, allocatable, dimension(:,:) :: area_grid
    double precision, allocatable, dimension(:,:) :: atb_grid

    integer, allocatable, dimension(:,:) :: riv_grid

    integer :: ncols, nrows
    double precision :: xllcorner, yllcorner, cellsize
    double precision :: double_nodata

    double precision :: area_threshold
    double precision :: atb_threshold

    CALL timer_get(run_start_time)

    in_area_file = ''
    in_atb_file = ''

    ! default value based on visualisation (for a 50m dem)
    area_threshold = 1250000
    atb_threshold = 14

    out_riv_file = 'river_thresholds.asc'

    input_is_valid = .true.

    i = 0

    tmp_char = 'DTA_river_thresholds.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- river_thresholds.f90 ---'
    write(999,*) ''
    print *, '--- Starting river_thresholds ---'

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-atb')) then
            CALL get_command_argument(i+1, in_atb_file)
        elseif (are_equal(arg, '-area')) then
            CALL get_command_argument(i+1, in_area_file)

        elseif (are_equal(arg, '-out_riv')) then
            CALL get_command_argument(i+1, out_riv_file)

        elseif (are_equal(arg, '-area_thresh')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) area_threshold

        elseif (are_equal(arg, '-atb_thresh')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) atb_threshold

        endif
        i = i + 1
    enddo

    if(len_trim(in_atb_file) == 0 .and. len_trim(in_area_file) == 0) then
        print *, 'Must specify -atb and/or -area'
        input_is_valid = .false.
    else
        if(len_trim(in_atb_file) /= 0) then
            if(check_file_arg(in_atb_file,'-atb').eqv..false.) then
                input_is_valid = .false.
            endif
        endif
        if(len_trim(in_area_file) /= 0) then
            if(check_file_arg(in_area_file,'-area').eqv..false.) then
                input_is_valid = .false.
            endif
        endif
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'command options '
        print *, 'choose either or both atb and area'
        print *, '-atb  <dem_atb.asc>  from atb_wfp.e processed as without river file'
        print *, '-area <dem_area.asc> from atb_wfp.e processed as without river file'
        print *, ''
        print *, 'optional'
        print *, '-atb_thresh 14 (default 14)'
        print *, '-area_thresh 1250000 (default 1250000)'
        stop
    endif

    write(999,*) 'Command Options Read In'
    write(999,*) 'atb  file: ', trim(in_atb_file)
    write(999,*) 'area file: ', trim(in_area_file)
    write(999,*) ''
    write(999,*) 'Using the following atb and area thresholds'
    write(999,*) 'atb  threshold: ', atb_threshold
    write(999,*) 'area threshold: ', area_threshold

    print *, 'atb  file: ', trim(in_atb_file)
    print *, 'area file: ', trim(in_area_file)
    print *, 'atb  threshold: ', atb_threshold
    print *, 'area threshold: ', area_threshold
    print *, ''

    write(999,*) ''

    if(len_trim(in_area_file) > 0)then
        print *, 'read area grid: ', trim(in_area_file)
        write(999,*) 'read area grid: ', trim(in_area_file)
        CALL timer_get(start_time)
        call read_ascii_grid(in_area_file, area_grid, &
            ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
        CALL timer_get(end_time)
        !call timer_print('read area grid', start_time, end_time)

        allocate(riv_grid(nrows,ncols))
        riv_grid(:,:) = 0

!print*, 'area',size(area_grid,1),size(area_grid,2)
!print*, 'riv_grid',size(riv_grid,1),size(riv_grid,2)
        where(area_grid > area_threshold) riv_grid = 1
        deallocate(area_grid)
    endif

    if(len_trim(in_atb_file) > 0)then
        print *, 'read atb grid: ', trim(in_atb_file)
        write(999,*) 'read atb grid: ', trim(in_atb_file)
        CALL timer_get(start_time)
        call read_ascii_grid(in_atb_file, atb_grid, &
            ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
        CALL timer_get(end_time)
        !call timer_print('read atb grid', start_time, end_time)

        if(allocated(riv_grid).eqv..false.) then
            allocate(riv_grid(nrows,ncols))
            riv_grid(:,:) = 0
        endif

        where(atb_grid > atb_threshold) riv_grid = riv_grid + 2

        deallocate(atb_grid)
    endif

    write(999,*) ''
    write(999,*) 'Write: ', trim(out_riv_file)

    print *, ''
    print *, 'Write: ', trim(out_riv_file)
    call write_ascii_grid_int(out_riv_file, riv_grid, &
        ncols, nrows, &
        xllcorner, yllcorner, cellsize, 0);


    CALL timer_get(end_time)

    !call timer_print('route_river', run_start_time, end_time)

    print *, '--- Finished river_thresholds ---'
    write(999,*) 'Successfully finished river_thresholds.f90'
    close(999)

    stop


end program river_thresholds
