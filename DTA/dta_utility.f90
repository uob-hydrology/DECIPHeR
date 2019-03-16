! useful utility functions for processing spatial grids
! reading, writing and unit conversions
!% Toby Dunne
!% Apr 2016
module dta_utility
    implicit none
    character(1) , parameter :: tab = char(9)

    ! Flow out of cell
    ! 32 | 64 | 127
    ! 16 |    |   1
    !  8 |  4 |   2
    !integer, parameter, dimension(8) :: direction_out = (/32,64,128,16,1,8,4,2/)
    integer(kind=1), parameter, dimension(8) :: dir_out = int((/32,64,127,16,1,8,4,2/),1)
    ! Flow in to cell
    !   2 |  4 |  8
    !   1 |    | 16
    ! 127 | 64 | 32
    !integer, parameter, dimension(8) :: direction_in = (/2,4,8,1,16,128,64,32/)
    integer(kind=1), parameter, dimension(8) :: dir_in = int((/2,4,8,1,16,127,64,32/),1)

    logical :: clock_init = .false.
    integer :: sys_clock_rate
    integer :: sys_clock_max

    type point_type
        integer :: x
        integer :: y
    end type point_type

    type point_list_type
        type(point_type), allocatable, dimension(:) :: list
    end type point_list_type

    type time_type
        real    :: timer_cpu_time
        integer :: timer_system_time
    end type time_type

contains

    subroutine timer_get(timer_value)
        type(time_type) :: timer_value

        if(clock_init.eqv..false.) then
            ! First initialize the system_clock
            call system_clock(count_rate=sys_clock_rate)
            call system_clock(count_max=sys_clock_max)
        endif

        call SYSTEM_CLOCK(timer_value%timer_system_time)
        call CPU_TIME(timer_value%timer_cpu_time)

    end subroutine timer_get

    function timer_diff_wall(timer_start, timer_end)
        type(time_type) :: timer_start, timer_end
        real :: timer_diff_wall

        timer_diff_wall = (timer_end%timer_system_time &
            - timer_start%timer_system_time) &
            / real(sys_clock_rate)


    end function timer_diff_wall

    subroutine timer_print(message, timer_start, timer_end)
        type(time_type) :: timer_start, timer_end
        character(*) :: message

        print '(A,'' CPU sec: '',F0.1,'' Took sec: '',F0.1)', trim(message),&
            timer_end%timer_cpu_time - timer_start%timer_cpu_time, &
            (timer_end%timer_system_time - timer_start%timer_system_time)/real(sys_clock_rate)

    end subroutine timer_print

    subroutine timer_estimate(processed_count, total, start_time, now_time, last_print_time)
        integer :: processed_count, total
        type(time_type) :: start_time, now_time, last_print_time
        real :: elapsedTime, estimatedTime, remainingTime

1502 format ('progress:',I0,'(',F0.1,'%) elapsed:',F0.1,'s est.:',F0.1,'s remaining:',F0.1,'s')

        if(timer_diff_wall(last_print_time, now_time) > 15) then
            last_print_time = now_time
            elapsedTime = timer_diff_wall(start_time, now_time)
            estimatedTime = elapsedTime * real(total) / real(processed_count)
            remainingTime = estimatedTime - elapsedTime

            print 1502, &
                total - processed_count, &
                real(processed_count)/real(total) * 100.0, &
                elapsedTime  , &
                estimatedTime, &
                remainingTime
        endif

    end subroutine

    subroutine NorthingEastingToRowCol(northing, easting, nrows, xllcorner, yllcorner, cellsize, row, col)
        implicit none
        double precision, intent(in):: northing, easting
        integer, intent(in):: nrows
        double precision, intent(in):: xllcorner, yllcorner, cellsize
        integer, intent(out) :: row, col

        row =  (nrows)- floor((northing - yllcorner)/cellsize)
        col = floor((easting - xllcorner)/cellsize) + 1

    end subroutine

    subroutine RowColToNorthingEasting(row, col, nrows, xllcorner, yllcorner, cellsize, northing, easting, cell_centre)
        implicit none
        integer, intent(in) :: row, col
        integer, intent(in) :: nrows
        double precision, intent(in) :: xllcorner, yllcorner, cellsize
        double precision, intent(out) :: northing, easting
        logical, intent(in):: cell_centre

        northing = yllcorner + ((nrows - row) * cellsize)
        easting = xllcorner + ((col-1) * cellsize )

        if(cell_centre) then
            northing = northing + cellsize / 2
            easting = easting + cellsize / 2
        endif

    end subroutine


    subroutine match_coordinates( cellsize, &
        a_ncols, a_nrows, a_xll, a_yll, &
        b_ncols, b_nrows, b_xll, b_yll, &
        ! outputs
        b_start_y, b_start_x, b_end_y, b_end_x, &
        a_start_y, a_start_x, a_end_y, a_end_x )
        double precision, intent(in) :: cellsize
        integer, intent(in) :: a_ncols, a_nrows
        double precision, intent(in) :: a_xll, a_yll
        integer, intent(in) :: b_ncols, b_nrows
        double precision, intent(in) :: b_xll, b_yll
        integer, intent(out) :: b_start_y, b_start_x, b_end_y, b_end_x
        integer, intent(out) :: a_start_y, a_start_x, a_end_y, a_end_x

        !%MATCH_COORDINATES locate cell coordinates of a in cell coordinates of b
        !% Used for two grids with different xll yll corners
        !% cell size must be the same for both grids

        !% convert coordinates of a to b
        call NorthingEastingToRowCol(a_yll, a_xll, b_nrows, b_xll, b_yll, cellsize, a_end_y, a_start_x)

        !%convert zero based to one based
        !%a_end_y = a_end_y + 1;
        !%a_start_x = a_start_x;

        a_start_y = a_end_y - (a_nrows-1)
        a_end_x = a_start_x + (a_ncols-1)

        b_start_y = max(a_start_y, 1)
        b_start_x = max(a_start_x, 1)

        b_end_y = min(a_end_y, b_nrows)
        b_end_x = min(a_end_x, b_ncols)


    end subroutine



    subroutine list_add_item_reserve(list, item, item_count, reserve)
        implicit none
        integer, intent(inout), allocatable, dimension(:) :: list
        integer, intent(in) :: item
        integer, intent(inout) :: item_count
        integer :: reserve

        !local
        integer, allocatable, dimension(:) :: list_copy
        integer :: reserve_value
        reserve_value = reserve

        item_count = item_count + 1
        if(reserve_value < item_count) then
            reserve_value = item_count
        endif

        if(item_count == 1) then
            allocate(list(reserve_value))
        elseif(size(list) < reserve_value) then
            ! allocate a new list with enough space and copy the existing values
            call move_alloc(list, list_copy)

            allocate(list(reserve_value))
            list(1:size(list_copy)) = list_copy
            deallocate(list_copy)
        endif

        list(item_count) = item
    end subroutine

    subroutine list_add_item(list, item, item_count)
        implicit none
        integer, intent(inout), allocatable, dimension(:) :: list
        integer, intent(in) :: item
        integer, intent(inout) :: item_count

        !just add the item and don't reserve any more space
        call list_add_item_reserve(list, item, item_count, 0 )

    end subroutine list_add_item

    subroutine list_remove_at(list, remove_index, item_count)
        implicit none
        integer, intent(inout), allocatable, dimension(:) :: list
        integer, intent(in) :: remove_index
        integer, intent(inout) :: item_count
        !locals
        integer :: i

        do i=(remove_index+1),item_count
            list(i-1) = list(i)
        end do
        item_count = item_count - 1
        if(item_count == 0) then
            deallocate (list)
        endif

    end subroutine list_remove_at


    subroutine write_point_list(filename, points, nrows, xllcorner, yllcorner, cellsize)
        implicit none
        character(1024) :: filename
        type(point_type) :: points(:)
        integer :: nrows
        double precision :: xllcorner
        double precision :: yllcorner
        double precision :: cellsize

        integer :: i
        double precision :: easting
        double precision :: northing

        open (10, file = filename, status = 'unknown')

        write (10,97) 'x_easting',tab,'y_northing',tab,'x_grid',tab,'y_grid'

        do i=1,size(points)
            call RowColToNorthingEasting(points(i)%y, points(i)%x, &
                nrows, xllcorner, yllcorner, cellsize, northing, easting, .true.)

            write (10,98) &
                easting,tab,  &
                northing,tab, &
                points(i)%x,tab, &
                points(i)%y

        enddo
        close (10)

        ! format is for header 4 header and 3 separators
97      format (7A)
        ! format is float tab float tab int tab int
98      format (F0.1,A,F0.1,A,I0,A,I0)

    end subroutine write_point_list

    function file_exists(filename) result(exists)
        character(len=*) ::filename
        logical:: exists
        integer::stat
        open (10, file = filename, status = 'old', iostat = stat)
        if(stat == 0) then
            exists = .true.
        else
            exists=.false.
        endif
        close (10)
    end function file_exists

    function check_file_arg(filename, arg) result(input_is_valid)
        character(len=*) :: filename
        character(len=*) :: arg
        logical ::input_is_valid
        input_is_valid = .true.
        if (len_trim(filename) == 0) then
            print *, trim(arg),' not specified'
            input_is_valid = .false.
        else if(file_exists(filename).eqv..false.) then
            print *, trim(arg),' file not found: ', trim(filename)
            input_is_valid = .false.
        endif

    end function check_file_arg

    ! similar to dlmread
    subroutine read_numeric_list(filename, columns, header_rows, data)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        integer :: columns
        integer :: header_rows
        double precision, allocatable, dimension(:,:) :: data

        integer :: i
        ! n: number of rows found
        integer :: n
        ! io: result of read
        integer :: ioerr
        character(4096) :: tmp

        open (99, file = filename, iostat=ioerr, status = 'old')
        if(ioerr/=0)then
            print *, 'File not found: ', trim(filename)
            stop
        endif


        ! skip header lines
        do i=1,header_rows
            read (99, '(A)') tmp
        enddo

        ! read through file, counting lines
        n = 0
        do
            read(99, *,iostat=ioerr) tmp
            if (ioerr/=0) exit
            n = n + 1
        end do

        allocate(data(n,columns))

        ! go back to the start of the file
        rewind 99

         ! skip header lines
        do i=1,header_rows
            read (99, '(A)') tmp
        enddo

        do i=1,n
            read (99, *) data(i,:)
        enddo

        close (99)


    end subroutine

    ! GC's alternative read_numeric_list to deal with dynatop using file IDs rather than the filename
    subroutine read_numeric_list_fid(fid, columns, header_rows, data)
        implicit none
        ! Argument declares
        integer :: fid
        integer :: columns
        integer :: header_rows
        double precision, allocatable, dimension(:,:) :: data

        integer :: i
        ! n: number of rows found
        integer :: n
        ! io: result of read
        integer :: io
        character(4096) :: tmp

        !open (99, file = filename, status = 'old')

        ! skip header lines
        do i=1,header_rows
            read (fid, '(A)') tmp
        enddo

        ! read through file, counting lines
        n = 0
        do
            read(fid, *,iostat=io) tmp
            if (io/=0) exit
            n = n + 1
        end do

        allocate(data(n,columns))

        ! go back to the start of the file
        rewind fid

         ! skip header lines
        do i=1,header_rows
            read (fid, '(A)') tmp
        enddo

        do i=1,n
            read (fid, *) data(i,:)
        enddo

        close (fid)


    end subroutine read_numeric_list_fid

    ! similar to dlmwrite
    subroutine write_numeric_list(filename, col_headers, data, data_precision)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        character(64) :: col_headers(:)
        integer :: data_precision
        double precision, dimension(:,:) :: data

        !locals
        character(65*size(col_headers)) :: header_line

        integer :: i, ioerr
        integer :: i_start, i_end
        character(64) :: tmp
        character(256) :: tmp_fmt
        character(16) :: valueFormat

        write (valueFormat, '(A,I0)') 'F0.', data_precision

        header_line=''
        do i=1,size(col_headers)
            if(i>1) then
                i_start = len_trim(header_line)+1
                i_end = i_start + 1
                header_line(i_start:i_end) = tab
            endif
            tmp = col_headers(i)
            i_start = len_trim(header_line)+1
            i_end = i_start + len_trim(tmp)
            header_line(i_start:i_end) = tmp(1:len_trim(tmp))
        end do

        open (99, file = filename, status = 'unknown', iostat=ioerr)

        if(ioerr/=0) then
            print *,ioerr
            print*,'error opening output file: ', trim(filename)
            print*,'ensure the directory exists and correct write permissions are set'
            stop
        endif

        write (99,'(A)') trim(header_line)

        ! build format based on number of columns and the required precision
        if(size(data,2)-1 > 0)then

            write(tmp_fmt, '(A,A,A,I0,A,A,A)') &
                '(1(', &
                trim(valueFormat), &
                ') ', &
                size(data,2)-1, &
                '(1x,', trim(valueFormat), &
                '))'
        else
            write(tmp_fmt, '(A,A,A)') &
                '(1(', &
                trim(valueFormat), &
                '))'
        endif

        !print *, tmp_fmt
        do i=1,size(data,1)
            write (99, tmp_fmt) data(i,:)
        enddo

        close (99)


    end subroutine

    subroutine read_ascii_grid_header(filename, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        integer :: ncols
        integer :: nrows
        double precision :: xllcorner
        double precision :: yllcorner
        double precision :: cellsize
        double precision :: nodata
        ! Local declares
        character(256) param_name
        integer :: ioerr

        open (99, file = filename, iostat=ioerr, status = 'old')
        if(ioerr/=0)then
            print *, 'File not found: ', trim(filename)
            stop
        endif

        read (99, *) param_name, ncols
        if(.NOT.are_equal(param_name, 'ncols')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, nrows
        if(.NOT.are_equal(param_name, 'nrows')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, xllcorner
        if(.NOT.are_equal(param_name, 'xllcorner')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, yllcorner
        if(.NOT.are_equal(param_name, 'yllcorner')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, cellsize
        if(.NOT.are_equal(param_name, 'cellsize')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, nodata

        if(are_equal(param_name, 'nodata') .or. are_equal(param_name, 'NODATA_value')) then
        else
            nodata = -9999
            ! move back in file to re-read as grid data
            BACKSPACE (99)
        endif

        close (99)
    end subroutine read_ascii_grid_header

    subroutine read_ascii_grid(filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        integer :: ncols
        integer :: nrows
        double precision :: xllcorner
        double precision :: yllcorner
        double precision :: cellsize
        double precision :: nodata
        double precision, allocatable, dimension(:,:) :: grid
        ! Local declares
        character(256) param_name
        integer :: i, ioerr

        open (99, file = filename, iostat=ioerr, status = 'old')
        if(ioerr/=0)then
            print *, 'File not found: ', trim(filename)
            stop
        endif

        read (99, *) param_name, ncols
        if(.NOT.are_equal(param_name, 'ncols')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, nrows
        if(.NOT.are_equal(param_name, 'nrows')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, xllcorner
        if(.NOT.are_equal(param_name, 'xllcorner')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, yllcorner
        if(.NOT.are_equal(param_name, 'yllcorner')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, cellsize
        if(.NOT.are_equal(param_name, 'cellsize')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, nodata

        if(are_equal(param_name, 'nodata') .or. are_equal(param_name, 'NODATA_value')) then
        else
            nodata = -9999
            ! move back in file to re-read as grid data
            BACKSPACE (99)
        endif

        !print *, 'nodata ', nodata
        !print *, 'nrows/ncols ', nrows, ncols

        allocate(grid(nrows, ncols));

        !print *, 'grid size ', size(grid)

        do i=1,nrows
            !print *, 'line: ' , i
            !write (99, gridformat) grid(i,:)
            read (99, *) grid(i,:)
        enddo

        !read (99, *) grid

        close (99)
    end subroutine read_ascii_grid

    subroutine read_ascii_grid_int(filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        integer :: ncols
        integer :: nrows
        double precision :: xllcorner
        double precision :: yllcorner
        double precision :: cellsize
        integer :: nodata
        integer, allocatable, dimension(:,:) :: grid
        ! Local declares
        character(256) param_name
        integer :: i, ioerr

        open (99, file = filename, iostat=ioerr, status = 'old')
        if(ioerr/=0)then
            print *, 'File not found: ', trim(filename)
            stop
        endif

        read (99, *) param_name, ncols
        if(.NOT.are_equal(param_name, 'ncols')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, nrows
        if(.NOT.are_equal(param_name, 'nrows')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, xllcorner
        if(.NOT.are_equal(param_name, 'xllcorner')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, yllcorner
        if(.NOT.are_equal(param_name, 'yllcorner')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, cellsize
        if(.NOT.are_equal(param_name, 'cellsize')) call exit_error('invalid ascii grid' // filename)

        read (99, *) param_name, nodata

        if(are_equal(param_name, 'nodata') .or. are_equal(param_name, 'NODATA_value')) then
        else
            nodata = -9999
            ! move back in file to re-read as grid data
            BACKSPACE (99)
        endif

        !print *, 'nodata ', nodata
        !print *, 'nrows/ncols ', nrows, ncols

        allocate(grid(nrows, ncols));

        !print *, 'grid size ', size(grid)

        do i=1,nrows
            !print *, 'line: ' , i
            !write (99, gridformat) grid(i,:)
            read (99, *) grid(i,:)
        enddo

        !read (99, *) grid

        close (99)
    end subroutine read_ascii_grid_int


    subroutine write_ascii_grid(filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, write_precision)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        integer :: ncols
        integer :: nrows
        double precision :: xllcorner
        double precision :: yllcorner
        double precision :: cellsize
        double precision :: nodata
        double precision, dimension(nrows, ncols) :: grid
        integer, intent(in) :: write_precision
        ! Local declares
        integer :: i
        character(30) :: gridformat
        character(10) :: valueformat
        character(10) :: valueformat2

        double precision :: max_value
        integer :: record_length
        integer :: record_length_min
        integer :: record_length_max



        ! must have enough space for the record length (RECL=)
        ! 1 record is one row
        ! record is made up of ncols * maximum space that a printed number takes up
        !
        ! each number will be n.m digits
        ! first make the format string from the write precision
        ! write out the largest value and measure the size of the string
        ! add 3 fortran needs this for some reason (1 for delimeter, 2 other)
        !
        ! write out with 3 decimal places
        !write_precision = 3

        ! set up the format string for the required precision
        ! eg.  f0.3
        !       0 means number of digits not limited
        !        .3 means number of decimal places
        write(valueformat, '(A,I0)') 'f0.', write_precision

        write(valueformat2, '(A,A,A)') '(', trim(valueformat), ')'

        ! find the number of digits required for the largest number in the grid
        ! this is needed to give enough space in the the record length for write
        max_value = maxval(grid)
        write(gridformat, valueformat2) max_value
        record_length_max = len_trim(gridformat)

        write(gridformat, valueformat2) -9999.0
        record_length_min = len_trim(gridformat)

        record_length = max(record_length_min, record_length_max)


        ! add '3' more space (seem fortran needs this, reason unknown)
        record_length = record_length + 3

        !print *, valueformat

        ! set up the format for writing a row
        ! first value without leading space
        ! then ncols-1 values with leading spaces
        !
        ! e.g. resulting format string with 500 cols precision 3
        ! gridformat = '(1(f0.3) 499(1x,f0.3))'

        write(gridformat, '(A,A,A,I0,A,A,A)') &
            '(1(', &
            trim(valueFormat), &
            ') ', &
            ncols-1, &
            '(1x,', trim(valueFormat), &
            '))'

        !print *, gridformat
        !print *, record_length*ncols

        open (99, file = filename, status = 'unknown', RECL=(record_length*ncols))

        write (99, 60) 'ncols ', ncols
        write (99, 60) 'nrows ', nrows
        write (99, 61) 'xllcorner ', xllcorner
        write (99, 61) 'yllcorner ', yllcorner
        write (99, 61) 'cellsize ', cellsize
        write (99, 61) 'NODATA_value ', nodata

        do i=1,nrows
            write (99, gridformat) grid(i,:)
        enddo

60      FORMAT(A, i0)
61      FORMAT(A, f0.3)
        !62      FORMAT(A, f0.3)

        close (99)
    end subroutine write_ascii_grid

    subroutine write_ascii_grid_int(filename, grid, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
        implicit none
        ! Argument declares
        character(len=*) :: filename
        integer :: ncols
        integer :: nrows
        double precision :: xllcorner
        double precision :: yllcorner
        double precision :: cellsize
        integer :: nodata
        integer, dimension(nrows, ncols) :: grid
        ! Local declares
        integer :: i
        character(30) :: gridformat
        character(10) :: valueformat
        character(10) :: valueformat2

        integer :: max_value
        integer :: record_length
        integer :: record_length_min
        integer :: record_length_max

        ! must have enough space for the record length (RECL=)
        ! 1 record is one row
        ! record is made up of ncols * maximum space that a printed number takes up
        !
        ! each number will be n digits
        ! first make the format string from the write precision
        ! write out the largest value and measure the size of the string
        ! add 3 fortran needs this for some reason (1 for delimeter, 2 other)
        valueformat = 'I0'
        valueformat2 = '(I0)'

        ! find the number of digits required for the largest number in the grid
        ! this is needed to give enough space in the the record length for write
        max_value = maxval(grid)
        write(gridformat, valueformat2) max_value
        record_length_max = len_trim(gridformat)

        write(gridformat, valueformat2) -9999
        record_length_min = len_trim(gridformat)

        record_length = max(record_length_min, record_length_max)


        ! add '3' more space (seem fortran needs this reason unknown)
        record_length = record_length + 3

        !print *, valueformat

        ! set up the format for writing a row
        ! first value without leading space
        ! then ncols-1 values with leading spaces
        !
        ! e.g. resulting format string with 500 cols precision 3
        ! gridformat = '(1(f0.3) 499(1x,f0.3))'

        write(gridformat, '(A,A,A,I0,A,A,A)') &
            '(1(', &
            trim(valueFormat), &
            ') ', &
            ncols-1, &
            '(1x,', trim(valueFormat), &
            '))'

        !print *, gridformat
        !print *, record_length*ncols

        open (99, file = filename, status = 'unknown', RECL=(record_length*ncols))

        write (99, 60) 'ncols ', ncols
        write (99, 60) 'nrows ', nrows
        write (99, 61) 'xllcorner ', xllcorner
        write (99, 61) 'yllcorner ', yllcorner
        write (99, 61) 'cellsize ', cellsize
        write (99, 62) 'NODATA_value ', nodata

        do i=1,nrows
            write (99, gridformat) grid(i,:)
        enddo

60      FORMAT(A, i0)
61      FORMAT(A, f0.3)
62      FORMAT(A, i0)

        close (99)
    end subroutine write_ascii_grid_int

    subroutine compare_masks_logical(a_ncols, a_nrows, &
        b_ncols, b_nrows, &
        a_mask, b_mask, &
        a_xll, a_yll, &
        b_xll, b_yll, &
        cellsize, &
        score, intersect_count, out_count)
        !%COMPARE_MASKS compare masks - may be different lower left, same cell size
        !%   calculates coordinate offsets and counts intersection of both masks
        !% Toby Dunne 2016/03/18
        !%test
        !%filenameA = 'C:\toby\LocalWork\WaterAndHealth\dta_test\mask_a.asc';
        !%filenameB = 'C:\toby\LocalWork\WaterAndHealth\dta_test\mask_b.asc';

        implicit none
        integer, intent(in) :: a_ncols, a_nrows, b_ncols, b_nrows
        logical, intent(in) :: a_mask(a_nrows, a_ncols)
        logical, intent(in) :: b_mask(b_nrows, b_ncols)
        double precision, intent(in) :: a_xll, a_yll, b_xll, b_yll
        double precision, intent(in) :: cellsize

        double precision, intent(out) :: score
        integer, intent(out) :: intersect_count, out_count

        ! Locals
        integer :: areaA, areaB
        !northing of top of a
        double precision :: a_y_top

        integer :: a_row_top, a_col_lef, a_col_rig, a_row_bot
        integer :: b_start_x, b_start_y, b_end_x, b_end_y
        integer :: yy, xx, y_a, x_a

        areaA = count(a_mask)
        areaB = count(b_mask)

        a_y_top = a_yll + cellsize * a_nrows

        ! convert coordinates of a to b
        ! a_row_top is the top  of a, in the coordinates of b
        ! a_col_lef is the left of a, in the coordinates of b
        call NorthingEastingToRowCol(a_y_top, a_xll, b_nrows, b_xll, b_yll, cellsize, &
            a_row_top, a_col_lef)

        a_row_top = a_row_top + 1
        a_col_lef = a_col_lef + 1

        a_col_rig = a_col_lef + a_ncols-1
        a_row_bot = a_row_top + a_nrows-1

        b_start_y = max(a_row_top, 1)
        b_start_x = max(a_col_lef, 1)

        b_end_y = min(a_row_bot, b_nrows)
        b_end_x = min(a_col_rig, b_ncols)

        intersect_count = 0

        do xx=b_start_x,b_end_x
            do yy=b_start_y,b_end_y
                y_a = yy - a_row_top+1
                x_a = xx - a_col_lef+1
                if(a_mask(y_a,x_a) .and. b_mask(yy,xx)) then
                    intersect_count = intersect_count + 1
                endif
            end do
        end do

        out_count = areaA + areaB - 2 * intersect_count
        if(intersect_count>0) then
            score = real(out_count)/intersect_count
        else
            score = areaA + areaB
        endif


    end subroutine compare_masks_logical


    function string_tolower( string ) result (new)
        implicit none
        character(len=*)           :: string

        character(len=len(string)) :: new

        integer                    :: i
        integer                    :: k
        integer                    :: length

        length = len(string)
        new    = string
        do i = 1,len(string)
            k = iachar(string(i:i))
            if ( k >= iachar('A') .and. k <= iachar('Z') ) then
                k = k + iachar('a') - iachar('A')
                new(i:i) = achar(k)
            endif
        enddo
    end function string_tolower


    function are_equal(string_a, string_b) result(equal)
        implicit none
        character(len=*) :: string_a
        character(len=*) :: string_b
        logical :: equal

        character(len=len(string_a)) :: lower_a
        character(len=len(string_b)) :: lower_b

        !print *, 'checkA ', len(string_a), '(', trim(string_a) , ')'
        !print *, 'checkB ', len(string_b), '(', trim(string_b) , ')'

        equal = .true.

        if(len_trim(string_a) == len_trim(string_b)) then
            lower_a = string_tolower(string_a)
            lower_b = string_tolower(string_b)

            if(lower_a /= lower_b) then
                !print *, 'not equal (', lower_a , ') : (', lower_b, ')'
                equal = .false.
            endif
        else
            !print *, 'not equal (', string_a , ') : (', string_b, ')'
            equal = .false.
        endif

        return

    end function

    subroutine split_comma(input, split_list)
        implicit none
        character(len=*) :: input
        integer :: pos
        integer :: field_count
        character(len(input)) tmp
        character(1024), allocatable :: split_list(:)

        field_count = 1

        pos = 1

        tmp = input

        do while (pos > 0)
            pos = index(tmp, ",")
            if(pos > 0) then
                tmp = tmp(pos+1:len_trim(tmp))
                field_count = field_count + 1
            endif
        end do

        allocate(split_list(field_count))

        pos=1
        field_count = 1
        tmp = input
        do while (pos > 0)
            pos = index(tmp, ",")
            if(pos > 0) then
                !split_list(field_count) = ''
                split_list(field_count) = tmp(1:pos-1)
                tmp = tmp(pos+1:len_trim(tmp))

                field_count = field_count + 1
            else
                split_list(field_count) = tmp
            endif
        end do
    end subroutine

    subroutine exit_error(message)
        implicit none
        character(len=*) :: message
        print *, message
        stop
    end subroutine exit_error

    subroutine utility_test
        implicit none
        logical result
        integer, allocatable :: list(:)
        integer :: item_count
        character(1024) :: input_string
        character(1024), allocatable :: split_list(:)

        item_count = 0

        call list_add_item(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)
        call list_add_item(list, 2, item_count)
        print *, item_count, 'values', list(1:item_count)
        call list_add_item(list, 3, item_count)
        print *, item_count, 'values', list(1:item_count)
        call list_add_item(list, 4, item_count)
        print *, item_count, 'values', list(1:item_count)


        call list_remove_at(list, 4, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_remove_at(list, 3, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_remove_at(list, 2, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_remove_at(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_add_item(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)
        call list_add_item(list, 2, item_count)
        print *, item_count, 'values', list(1:item_count)
        call list_add_item(list, 3, item_count)
        print *, item_count, 'values', list(1:item_count)
        call list_add_item(list, 4, item_count)
        print *, item_count, 'values', list(1:item_count)


        call list_remove_at(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_remove_at(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_remove_at(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)

        call list_remove_at(list, 1, item_count)
        print *, item_count, 'values', list(1:item_count)



        result = are_equal('a', 'a')
        print *, ' a a ', result
        result = are_equal('a', 'A')
        print *, ' a A ', result

        result = are_equal('all', 'All')
        print *, ' all All ', result

        result = are_equal('all', 'bll')
        print *, ' all bll ', result

        result = are_equal('all', 'allb')
        print *, ' all allb ', result


        input_string = 'test'
        print*,'----------'
        call split_comma(input_string, split_list)
        if(size(split_list) /= 1) then
            print*,'Error split list 1'
        else
            print*,'test=',trim(split_list(1))
        endif
        deallocate(split_list)

        print*,'----------'
        input_string = 'test,test2'
        call split_comma(input_string, split_list)
        if(size(split_list) /= 2) then
            print*,'Error split list 2'
        else
            print*,'test=',trim(split_list(1))
            print*,'test2=',trim(split_list(2))
        endif
        deallocate(split_list)
        print*,'----------'


    end subroutine utility_test


    subroutine calc_flow_direction( nrows, ncols, dem, cellsize, flow_dir_grid, slope_grid)
        implicit none
        integer, intent(in) :: nrows
        integer, intent(in) :: ncols
        double precision, intent(in) :: dem (nrows, ncols)
        double precision, intent(in) :: cellsize
        integer(1), intent(inout) :: flow_dir_grid (nrows, ncols)
        double precision, intent(inout) :: slope_grid(nrows, ncols)

        double precision :: dist_cardinal
        double precision :: dist_ordinal

        integer :: x, y
        integer :: max_slope_index
        double precision :: max_slope
        type(point_type) :: flow_point

        dist_cardinal = 1 * cellsize
        dist_ordinal = sqrt(2.0) * cellsize

        do x = 1,ncols
            do y = 1,nrows
                call cell_flow_dir(ncols, nrows, dem, dist_cardinal, dist_ordinal, x, y, max_slope_index, flow_point, max_slope)
                if(max_slope_index > 0) then
                    flow_dir_grid(y,x) = dir_out(max_slope_index)
                    slope_grid(y,x) = max_slope
                else
                    flow_dir_grid(y,x) = 0
                    slope_grid(y,x) = -9999
                endif
            end do
        end do


    end subroutine

    !
    ! max_slope_index
    ! |---|---|---|
    ! | 1 | 2 | 3 |
    ! | 4 | 0 | 5 |
    ! | 6 | 7 | 8 |
    ! |---|---|---|
    subroutine cell_flow_dir(ncols, nrows, dem, dist_cardinal, dist_ordinal, x, y, max_slope_index, flow_point, max_slope)
        implicit none
        integer, intent(in) :: ncols
        integer, intent(in) :: nrows
        double precision, intent(in) :: dem(nrows, ncols)
        double precision :: dist_cardinal
        double precision :: dist_ordinal
        integer :: x
        integer :: y
        integer, intent(out) :: max_slope_index
        type(point_type), intent(out) :: flow_point
        double precision, intent(out) :: max_slope

        double precision :: dist

        integer :: neighbour_index
        double precision :: slope
        integer :: x_n, y_n

        max_slope_index = 0
        flow_point%x = x
        flow_point%y = y

        if(dem(y,x) > -90) then
            max_slope = -999999
            neighbour_index = 0
            do y_n = y-1,y+1
                do x_n = x-1,x+1
                    if (y_n == y .and. x_n == x) then
                        cycle
                    endif
                    neighbour_index = neighbour_index + 1
                    if(y_n < 1 .or. y_n > nrows  &
                        .or. x_n <1 .or. x_n > ncols) then
                        cycle
                    endif
                    if (dem(y_n, x_n) < -90) then
                        cycle
                    endif

                    if(y_n == y .or. x_n == x) then
                        dist = dist_cardinal
                    else
                        dist = dist_ordinal
                    endif

                    slope = (dem(y,x) - dem(y_n, x_n)) / dist

                    if (slope > max_slope) then
                        max_slope = slope
                        max_slope_index = neighbour_index
                        flow_point%x = x_n
                        flow_point%y = y_n
                    endif
                end do
            end do
        endif
    end subroutine



end module dta_utility
