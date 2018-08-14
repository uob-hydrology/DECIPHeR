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
    implicit none

    type(time_type) :: end_time, run_start_time
    character(1024) in_mask_file
    character(1024) in_riv_id_file
    integer, allocatable, dimension(:,:) :: mask_grid
    integer, allocatable, dimension(:,:) :: riv_id_grid

    integer :: ncols, nrows, riv_ncols, riv_nrows
    double precision :: xllcorner, yllcorner, cellsize, riv_xllcorner, riv_yllcorner, riv_cellsize
    integer :: nodata, riv_nodata

    logical input_is_valid
    integer ioerr
    integer i, x, y
    character(1024) arg
    character(1024) tmp_char
    integer :: bad_count
    double precision :: x_easting, y_northing

    in_mask_file = ''
    in_riv_id_file = ''

    input_is_valid = .true.

    tmp_char = 'DTA_mask_check.log'
    open(999, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(999,*) '--- mask_check.f90 ---'
    write(999,*) ''
    print *, '--- Starting mask_check ---'

    i=0

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-mask')) then
            CALL get_command_argument(i+1, in_mask_file)
        elseif (are_equal(arg, '-riv_id')) then
            CALL get_command_argument(i+1, in_riv_id_file)
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_mask_file,'-mask').eqv..false.) then
        input_is_valid = .false.
    endif
    if(check_file_arg(in_riv_id_file,'-riv_id').eqv..false.) then
        input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'command options '
        print *, '-mask <file.asc> mask file as created by catch_mask'
        print *, '-riv_id <file.asc> river mask file'
        print *, ''
        print *, 'e.g.'
        print *, 'mask_check.e -mask dem_mask.asc -riv_id dem_riv_id.asc'
        stop
    endif

    CALL timer_get(run_start_time)

    write(999,*) 'read : ', trim(in_mask_file)
    write(999,*) 'read : ', trim(in_riv_id_file)
    write(999,*) ''

    print *, 'read : ', trim(in_mask_file)
    print *, 'read : ', trim(in_riv_id_file)

    call read_ascii_grid_int(in_mask_file, mask_grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
    call read_ascii_grid_int(in_riv_id_file, riv_id_grid, &
        riv_ncols, riv_nrows, riv_xllcorner, riv_yllcorner, riv_cellsize, riv_nodata)

    if (ncols/=riv_ncols .or. &
        nrows /= riv_nrows .or. &
        abs(xllcorner - riv_xllcorner) > 0.001 .or. &
        abs(yllcorner - riv_yllcorner) > 0.001 .or. &
        abs(cellsize - riv_cellsize) > 0.001) then

        print*,' ascii grids do sizes do not match'
        stop
    endif

61  FORMAT(f0.1,A,f0.1,4(A,I0))

    bad_count = 0
    do y=1,nrows
        do x=1,ncols
            if(riv_id_grid(y,x) /= riv_nodata) then
                if(mask_grid(y,x) == nodata) then
                    ! river is out side the mask, this is ok - this just means that the mask is a subset
                    ! and the river is outside the main catchment
                    ! set the river to zero
                    !riv_id_grid(y,x) = riv_nodata
                    !print *, 'NOTE: river outside of mask has been removed'
                else
                    if(riv_id_grid(y,x) /= mask_grid(y,x) ) then

                        call RowColToNorthingEasting(y, x, &
                            nrows, xllcorner, yllcorner, cellsize, &
                            x_easting, y_northing, .true.)

                        !print *, 'WARNING: bad river id (logged)',x_easting, y_northing, x, y, mask_grid(y,x), riv_id_grid(y,x)
                        if (bad_count.lt.1) then
                        write(999,*) 'x_easting',tab, 'y_northing',tab, 'x',tab,'y',tab,'expect_id',tab,'actual_id'
                        end if

                        write (999, 61) x_easting, tab, &
                                        y_northing, tab, &
                                        x, tab, &
                                        y, tab, &
                                                       mask_grid(y,x),tab, &
                                                       riv_id_grid(y,x)
                        riv_id_grid(y,x) = riv_nodata

                        bad_count = bad_count + 1

                    endif
                endif
            endif
        end do
    end do

    if(bad_count > 0) then
        print*,'River cells found on non matching mask', bad_count
        print*,'See mas_check.log'
    endif

    tmp_char = in_riv_id_file(1:len_trim(in_riv_id_file)-4)//'_check.asc'

    print*,'write: ', trim(tmp_char)
    write(999,*) 'write: ', trim(tmp_char)
    call write_ascii_grid_int(tmp_char, riv_id_grid, &
        ncols, nrows, &
        xllcorner, yllcorner, cellsize, riv_nodata)

    CALL timer_get(end_time)
    !call timer_print('catch_mask', run_start_time, end_time)

    write(999,*) ''
    write(999,*) 'Successfully finished mask_check.f90'
    print *, '--- Finished mask_check ---'
    close(999)

    stop


end program catch_mask





