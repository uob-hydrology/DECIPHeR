!% Calculate the Topographic wetness index
!% intput DEM must have no pits and no flat areas
!% i.e. all cells must flow downstream
!%
!% The topographic index of a every cell is dependent on the
!% connected upstream cells.
!% In the first stage only peak cells can be processed.
!% After which an A* based algorithm is used to search
!% for subsequent cells that can be processed
!%
!% Toby Dunne
!% Apr 2016
program atb_wfp
    use dta_utility
    use dta_atb_wfp
    use dta_qsort
    implicit none

    integer :: ncols
    integer :: nrows
    double precision, allocatable, dimension(:,:) :: dem
    integer :: riv_ncols
    integer :: riv_nrows
    double precision, allocatable, dimension(:,:) :: riv
    integer :: catch_mask_ncols
    integer :: catch_mask_nrows
    integer, allocatable, dimension(:,:) :: catch_mask

    double precision :: xllcorner
    double precision :: yllcorner
    double precision :: cellsize
    double precision :: nodata, nodata_riv
    integer :: nodata_int
    integer :: i
    integer :: j
    integer :: allocateStatus, ioerr

    type(point_type), allocatable, dimension(:) :: peak_list
    type(point_type), allocatable, dimension(:) :: bound_list

    double precision, allocatable, dimension(:,:) :: a
    double precision, allocatable, dimension(:,:) :: c
    double precision, allocatable, dimension(:,:) :: slope

    character(1024) in_dem_file

    ! mask file - single flow direction boundary masks
    character(1024) in_mask_file
    ! riv file - <=0 not river cell, > 0 river cell
    character(1024) in_riv_file

    character(1024) output_dir
    character(1024) tmp_char
    character(8) char_riv
    character(8) char_mask

    type(time_type) :: start_time, end_time, run_start_time

    CHARACTER(len=1024) :: arg
    logical :: input_is_valid

    CALL timer_get(run_start_time)

    input_is_valid = .true.

    ! these will remain 0 unless allocated
    riv_ncols = 0
    riv_nrows = 0
    catch_mask_ncols = 0
    catch_mask_nrows = 0

    in_dem_file = ''
    in_mask_file = ''
    in_riv_file = ''
    char_mask = ''
    char_riv = ''

    output_dir = ''

    i = 0

    tmp_char = 'DTA_atb_wfp.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- atb_wfp.f90 ---'
    write(999,*) ''
    print *, '--- Starting atb_wfp ---'

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, in_dem_file)
        else if (are_equal(arg, '-mask')) then
            CALL get_command_argument(i+1, in_mask_file)
        else if (are_equal(arg, '-river')) then
            CALL get_command_argument(i+1, in_riv_file)
        else if (are_equal(arg, '-out')) then
            CALL get_command_argument(i+1, output_dir)
        endif
        i = i + 1
    end do

    if(check_file_arg(in_dem_file,'-dem').eqv..false.) then
       input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'command options'
        print *, '-dem  select input file'
        print *, '     e.g. atb_wfp.e -dem dem.asc'
        print *, ' optional arguments:'
        print *, '-mask  select input mask file grid cells labelled with catchment'
        print *, '-river select input river file grid >0 is river cell'
        print *, ''
        print *, ' output will be: dem_atb.asc'
        print *, '               : dem_area.asc'
        print *, '               : dem_slope_mfd.asc'
        print *, ''
        stop
    endif

    write(999,*) 'in_dem_file: ', trim(in_dem_file)
    write(999,*) 'in_mask_file: ', trim(in_mask_file)
    write(999,*) 'in_riv_file: ', trim(in_riv_file)
    write(999,*) ''

    print *, 'in_dem_file: ', trim(in_dem_file)
    print *, 'in_mask_file: ', trim(in_mask_file)
    print *, 'in_riv_file: ', trim(in_riv_file)
    print *, ''

    print *, 'read_ascii_grid dem'
    write(999,*) 'read_ascii_grid dem'

    CALL timer_get(start_time)

    call read_ascii_grid(in_dem_file, dem, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

    CALL timer_get(end_time)
    !call timer_print('read_ascii_grid', start_time, end_time)

    if(len_trim(in_mask_file) > 0) then
        print *, 'read_ascii_grid mask'
        write(999,*) 'read_ascii_grid mask'
        char_mask = '_mask'
        CALL timer_get(start_time)

        call read_ascii_grid_int(in_mask_file, &
            catch_mask, catch_mask_ncols, catch_mask_nrows,&
            xllcorner, yllcorner, cellsize, nodata_int)

        CALL timer_get(end_time)
        !call timer_print('read_ascii_grid', start_time, end_time)

        if(catch_mask_ncols /= ncols .or. catch_mask_nrows /= nrows) then
            print *, 'Error: mask grid size different from dem'
            stop
        endif
    endif

    if(len_trim(in_riv_file) > 0) then
        print *, 'read_ascii_grid river'
        write(999,*) 'read_ascii_grid river'
        char_riv = '_riv'
        CALL timer_get(start_time)

        call read_ascii_grid(in_riv_file, riv, riv_ncols, riv_nrows, xllcorner, yllcorner, cellsize, nodata_riv)
        CALL timer_get(end_time)
        !call timer_print('read_ascii_grid', start_time, end_time)

        if(riv_ncols /= ncols .or. riv_nrows /= nrows) then
            print *, 'Error: river grid size different from dem'
            stop
        endif
    endif

    allocate(a(nrows, ncols), stat = AllocateStatus)
    if (allocatestatus /= 0) then
        stop "*** Not enough memory (a)***"
    endif
    allocate(c(nrows, ncols), stat = AllocateStatus)
    if (allocatestatus /= 0) then
        stop "*** Not enough memory (c)***"
    endif
    allocate(slope(nrows, ncols), stat = AllocateStatus)
    if (allocatestatus /= 0) then
        stop "*** Not enough memory (slope)***"
    endif

    if(nint(nodata) /= -9999) then
        CALL timer_get(start_time)
        do i=1,nrows
            do j=1,ncols
                if(dem(i,j) < -90 .or. &
                    abs(dem(i,j) - nodata) < 0.0001) then
                    dem(i,j) = -9999
                endif
            enddo
        enddo
        CALL timer_get(end_time)
        !call timer_print('update nodata', start_time, end_time)
        nodata = -9999
    endif

    CALL timer_get(start_time)
    call find_peaks(nrows, ncols, dem, peak_list)
    CALL timer_get(end_time)
    !call timer_print('find peaks', start_time, end_time)
    !print '(A, I0, '' / '', I0)', 'peaks found: ', size(peak_list), nrows*ncols

    ! add all the catchment boundary cells as peaks
    ! doesn't matter if too many start points, as long as they are sorted
    if(allocated(catch_mask)) then

        call find_mask_bounds(nrows, ncols, catch_mask, bound_list, peak_list)
        !print '(A, I0)', 'mask boundary cells: ', size(bound_list) - size(peak_list)
        deallocate(peak_list)
        call move_alloc(bound_list, peak_list)

    endif

    CALL timer_get(start_time)
    !print *, 'sort peaks', size(peak_list)
    call QsortC(peak_list, dem)
    CALL timer_get(end_time)
    !call timer_print('sort peaks', start_time, end_time)

    CALL timer_get(start_time)
    print *, 'calculating accumulated area, slope and topographic index'
    write(999,*) 'calculating accumulated area, slope and topographic index'

    call calc_atb(nrows, ncols, &
        riv_nrows, riv_ncols, &
        catch_mask_nrows, catch_mask_ncols, &
        dem, riv, catch_mask, peak_list, a, c, slope, cellsize)
    CALL timer_get(end_time)
    !call timer_print('calc atb', start_time, end_time)

    CALL timer_get(start_time)
    !print *, 'logarithm'
    do i=1,nrows
        do j=1,ncols
            if(c(i,j)>-9000) then
                if(c(i,j)>0) then
                    c(i,j) = log(c(i,j))
                    ! GC - ensures topographic index is never zero - can happen with small grid cells
                    if(c(i,j)<0) then
                    c(i,j) = 0.01
                    endif
                else
                    c(i,j) = 0
                endif
            else
                c(i,j) = -9999
                !a(i,j) = -9999
            endif
            if(dem(i,j) < -9000) then
                a(i,j) = -9999
            endif
        enddo
    enddo
    CALL timer_get(end_time)
    !call timer_print('logarithm', start_time, end_time)

    write(999,*) ''

    tmp_char = in_dem_file(1:len_trim(in_dem_file)-4) //trim(char_riv)//trim(char_mask)// '_atb.asc'
    print *, ''
    print *, 'write:', trim(tmp_char)
    write(999,*) 'writing: ', trim(tmp_char)
    CALL timer_get(start_time)
    call write_ascii_grid(tmp_char, c, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, 3)
    CALL timer_get(end_time)
    !call timer_print('write', start_time, end_time)

    tmp_char = in_dem_file(1:len_trim(in_dem_file)-4) //trim(char_riv)//trim(char_mask)// '_area.asc'
    print *, 'write:', trim(tmp_char)
    write(999,*) 'writing: ', trim(tmp_char)
    CALL timer_get(start_time)
    call write_ascii_grid(tmp_char, a, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, 3)
    CALL timer_get(end_time)
    !call timer_print('write', start_time, end_time)

    tmp_char = in_dem_file(1:len_trim(in_dem_file)-4) //trim(char_riv)//trim(char_mask)// '_mfd_slope.asc'
    print *, 'write:', trim(tmp_char)
    write(999,*) 'writing: ', trim(tmp_char)
    CALL timer_get(start_time)
    call write_ascii_grid(tmp_char, slope, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, 6)
    CALL timer_get(end_time)
    !call timer_print('write', start_time, end_time)

    CALL timer_get(end_time)
    print *, '--- Finished atb_wfp ---'
    write(999,*) 'Successfully finished atb_wfp.f90'
    close(999)
    !call timer_print('atb_wfp', run_start_time, end_time)

end program atb_wfp
