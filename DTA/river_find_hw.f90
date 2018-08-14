!% Toby Dunne
!% April 2016
program river_find_hw
    use dta_utility
    use dta_rivers
    implicit none

    CHARACTER(len=1024) :: arg
    logical :: input_is_valid
    integer :: i, x, y, cell_count, x_start, y_start
    character(1024) in_river_file
    character(1024) in_dem_file
    character(1024) in_outlets_file
    double precision :: search_dist
    double precision :: move_downstream_dist

    ! filenames for the searches on river file saved with river file prefix
    character(1024) output_riv_file_prefix
    character(1024) tmp_char

    double precision, allocatable, dimension(:,:) :: dem_grid
    double precision, allocatable, dimension(:,:) :: riv_grid
    logical, allocatable, dimension(:,:) :: riv_mask_grid
    logical, allocatable, dimension(:,:) :: sea_mask_grid
    integer, allocatable, dimension(:,:) :: riv_label_grid

    double precision, allocatable, dimension(:,:) :: riv_dist_grid

    type(point_type), allocatable, dimension(:) :: sea_outlets
    type(point_type), allocatable, dimension(:) :: headwater_pass1
    type(point_type), allocatable, dimension(:) :: headwater_pass2

    double precision, allocatable, dimension(:,:) :: outlets_input_data

    integer :: ncols, nrows, ioerr
    double precision :: xllcorner, yllcorner, cellsize
    double precision :: double_nodata
    type(time_type) :: start_time, end_time, run_start_time

    CALL timer_get(run_start_time)

    in_river_file = ''
    in_dem_file = ''
    in_outlets_file = ''

    search_dist = 500.0d0
    move_downstream_dist = 100.0d0

    input_is_valid = .true.

    i = 0

    tmp_char = 'DTA_river_find_hw.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- river_find_hw.f90 ---'
    write(999,*) ''
    print *, '--- Starting river_find_hw ---'

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, in_dem_file)
        elseif (are_equal(arg, '-river')) then
            CALL get_command_argument(i+1, in_river_file)
        elseif (are_equal(arg, '-outlets')) then
            CALL get_command_argument(i+1, in_river_file)
        elseif (are_equal(arg, '-search')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) search_dist
        elseif (are_equal(arg, '-move_downstream')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) move_downstream_dist
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_river_file, '-river').eqv..false.) then
        input_is_valid = .false.
    endif

    if (len_trim(in_outlets_file) == 0 .and. len_trim(in_dem_file) == 0) then
        print*,'must set either -dem or -outlets'
        input_is_valid = .false.
    elseif (len_trim(in_outlets_file) /= 0 .and. len_trim(in_dem_file) /= 0) then
        print*,'cannot set both -dem and -outlets'
        input_is_valid = .false.
    else
        if (len_trim(in_dem_file) > 0) then
            if(file_exists(in_dem_file).eqv..false.) then
                print *, '-dem file not found: ', trim(in_dem_file)
                input_is_valid = .false.
            endif
        endif

        if (len_trim(in_outlets_file) > 0) then
            if(file_exists(in_outlets_file).eqv..false.) then
                print *, '-outlets file not found: ', trim(in_outlets_file)
                input_is_valid = .false.
            endif
        endif
    endif


    if(input_is_valid .eqv. .false.) then
        print *, 'command options '
        print *, '-river <file.asc> select river ascii grid file'
        print *, '                  used as a mask any value > 0 is river cell'
        print *, '                  e.g. from a rasterized externally provided river'
        print *, ''
        print *, 'set one of -dem or -outlets:'
        print *, ''
        print *, '-dem <file.asc>   dem used as sea mask (sea = nodata)'
        print *, '                  all river cells will be joined to sea'
        print *, '                  outlets will be found where river meets the sea'
        print *, ''
        print *, '-outlets <file.txt> file with eastings and northings of outlets'
        print *, '                    if outlets are set river will be used as provided with no'
        print *, '                    joining to the sea'
        stop
    endif

    write(999,*) 'Command Options Read In'
    write(999,*) 'river  file: ', trim(in_river_file)
    write(999,*) ''
    write(999,*) 'Using the following search and move downstream distance'
    write(999,*) 'search distance', search_dist
    write(999,*) 'move downstream distance', move_downstream_dist

    print *, 'in_river: ', trim(in_river_file)
    print *, 'search distance', search_dist
    print *, 'move downstream distance', move_downstream_dist

    output_riv_file_prefix = in_river_file(1:(len_trim(in_river_file)-4))

    write(999,*) 'read river grid:', trim(in_river_file)
    print *, 'read river grid:', trim(in_river_file)

    CALL timer_get(start_time)
    call read_ascii_grid(in_river_file, riv_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read river grid', start_time, end_time)

    ! convert the river file to a logical mask
    allocate(riv_mask_grid(nrows, ncols))
    riv_mask_grid = riv_grid > 0.001
    deallocate(riv_grid)

    write(999,*) 'number of river cells:', count(riv_mask_grid)

    if(len_trim(in_dem_file) > 0) then
        ! DEM MODE
        write(999,*) 'read dem grid:', trim(in_dem_file)
        print *, 'read dem grid:', trim(in_dem_file)
        CALL timer_get(start_time)
        call read_ascii_grid(in_dem_file, dem_grid, &
            ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
        CALL timer_get(end_time)
        !call timer_print('read dem grid', start_time, end_time)

        !        dem_grid(1,:) = -99
        !        dem_grid(nrows,:) = -99
        !        dem_grid(:,1) = -99
        !        dem_grid(:,ncols) = -99


        allocate(sea_mask_grid(nrows, ncols))
        sea_mask_grid = dem_grid < -90
        !deallocate(dem_grid)

        print *, 'Finding Outlets'
        write(999,*) 'Finding Outlets'

        call river_find_outlets(nrows, ncols, riv_mask_grid, sea_mask_grid, sea_outlets )

        write(999,*) 'outlets found: ', size(sea_outlets)

        !write (tmp_char,'(A,A)') &
        !    trim(output_riv_file_prefix), &
        !    '_outlets.txt'

        ! outlets (write output not really required, maybe useful to visualise)
        !call write_point_list(tmp_char, sea_outlets, nrows, xllcorner, yllcorner, cellsize)

        ! finished with the sea mask
        deallocate(sea_mask_grid)
    else
        ! OUTLETS MODE

        call read_numeric_list(in_outlets_file, 2, 1, outlets_input_data)
        ! Convert northing easting to row col (point)
        allocate (sea_outlets(size(outlets_input_data,1)))
        do i=1,size(outlets_input_data,1)

            call NorthingEastingToRowCol(outlets_input_data(i,2), outlets_input_data(i,1), &
                nrows, xllcorner, yllcorner, cellsize, &
                sea_outlets(i)%y, sea_outlets(i)%x)

        end do

        write(999,*) 'outlets from file: ', size(sea_outlets)

        deallocate(outlets_input_data)
    endif

    print *, 'Headwater Pass 1'
    write(999,*) 'Headwater Pass 1'

    allocate(riv_dist_grid(nrows,ncols))
    riv_dist_grid(:,:) = -99
    CALL timer_get(start_time)
    call river_find_headwater_pass1(nrows, ncols, &
        riv_mask_grid, sea_outlets, riv_dist_grid, headwater_pass1)
    CALL timer_get(end_time)
    !call timer_print('headwater_pass1', start_time, end_time)

    if(len_trim(in_dem_file) > 0) then
        ! DEM MODE

        cell_count = count(abs(riv_dist_grid + 1) < 0.0001)
        write(999,*)'Found river cells not linked to sea: ', cell_count
        ! count non joined rivers (distance == -1)
        if(cell_count >0) then
            write(999,*) 'Joining river cells and re-processing'

            ! need to join rivers and reprocess

            ! sea_outlets will have to be found again
            deallocate(sea_outlets)
            ! headwater_pass1 will have to be found again
            deallocate(headwater_pass1)

            allocate(riv_label_grid(nrows,ncols))
            riv_label_grid(:,:) = 0

            where(dem_grid < -90) riv_label_grid = -99
            where(riv_mask_grid) riv_label_grid = 1

            do y=1,nrows
                do x=1,ncols
                    if(abs(riv_dist_grid(y, x) + 1) < 0.0001) then
                        call river_find_lowest(nrows, ncols, dem_grid, riv_dist_grid, x, y, x_start, y_start)
                        call river_single_flow_point(nrows, ncols, dem_grid, cellsize, riv_label_grid, &
                        x_start, y_start, 1, &
                        .false.)
                    endif
                enddo
            enddo

            ! apply joined sections to mask
            where ( riv_label_grid > 0) riv_mask_grid = .true.
            where(dem_grid < -90) riv_mask_grid = .false.

            deallocate(riv_label_grid)

            allocate(sea_mask_grid(nrows, ncols))
            sea_mask_grid = dem_grid < -90
            deallocate(dem_grid)

            !print *, 'sea cells', count(sea_mask_grid)

            call river_find_outlets(nrows, ncols, riv_mask_grid, sea_mask_grid, sea_outlets )

            !print *, 'outlets found: ', size(sea_outlets)

            !write (tmp_char,'(A,A)') &
            !    trim(output_riv_file_prefix), &
            !    '_outlets2.txt'

            ! outlets (write output not really required, maybe useful to visualise)
            !call write_point_list(tmp_char, sea_outlets, nrows, xllcorner, yllcorner, cellsize)

            ! finished with the sea mask
            deallocate(sea_mask_grid)

            riv_dist_grid(:,:) = -99

            CALL timer_get(start_time)
            call river_find_headwater_pass1(nrows, ncols, &
                riv_mask_grid, sea_outlets, riv_dist_grid, headwater_pass1)
            CALL timer_get(end_time)
            !call timer_print('headwater_pass1', start_time, end_time)

            cell_count = count(abs(riv_dist_grid + 1) < 0.0001)
            write(999,*) 'Found river cells not linked to sea: ', cell_count
            if(cell_count > 0) then
                write(999,*) 'These cells flow outside the catchment boundaries'
                tmp_char = trim(output_riv_file_prefix)//'_dist.asc'
                write(999,*) 'cells can be checked where dist == -1: ', trim(tmp_char)
            endif

        endif
    endif

    ! save headwater list (pass 1) (not really required)
    !write (tmp_char,'(A,A)') &
    !    trim(output_riv_file_prefix), &
    !    '_HW_pass1.txt'
    !call write_point_list(tmp_char, headwater_pass1, nrows, xllcorner, yllcorner, cellsize)

    tmp_char = trim(output_riv_file_prefix)//'_dist.asc'
    print*, 'Write: ', trim(tmp_char)
    write(999,*) ''
    write(999,*) 'Write: ', trim(tmp_char)
    call write_ascii_grid(tmp_char, riv_dist_grid, &
        ncols, nrows, &
        xllcorner, yllcorner, cellsize, -99.0d0, 3);

    CALL timer_get(start_time)
        print *, 'Headwater Pass 2'
    write(999,*) 'Headwater Pass 2'
    call river_find_headwater_pass2(nrows, ncols, &
        riv_dist_grid, headwater_pass1, cellsize, &
        search_dist, move_downstream_dist, &
        headwater_pass2)
    CALL timer_get(end_time)
    !call timer_print('headwater_pass2', start_time, end_time)

    ! save headwater list (pass 2)
    print *, 'Writing Headwater List'
    write(999,*) 'Writing Headwater List'
    write (tmp_char,'(A,A,I0,''m_'',I0,''m.txt'')') &
        trim(output_riv_file_prefix), &
        '_HW_',  &
        nint(search_dist), &
        nint(move_downstream_dist)
    call write_point_list(tmp_char, headwater_pass2, nrows, xllcorner, yllcorner, cellsize)

    CALL timer_get(end_time)
    !call timer_print('river_find_hw', run_start_time, end_time)

    print *, '--- Finished river_find_hw ---'
    write(999,*) 'Successfully finished river_find_hw.f90'
    close(999)

    stop


end program river_find_hw
