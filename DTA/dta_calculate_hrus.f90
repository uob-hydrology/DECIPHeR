module dta_calculate_hrus

    implicit none

contains

    subroutine reclass_grid(data_in, nodata, data_out, reclass_table)

        implicit none

        integer, intent(in), dimension(:,:) :: data_in
        integer, intent(in) :: nodata
        integer, allocatable, dimension(:,:), intent(out) :: data_out
        integer, allocatable, dimension(:,:), intent(out) :: reclass_table

        ! Locals
        logical, dimension(size(data_in, 1), size(data_in, 2)) :: mask
        integer :: y, x, num_unique
        ! Get all the unique elements from data_in
        num_unique = 0
        mask = .false.

        do y = 1, size(data_in, 1)
            do x = 1, size(data_in, 2)
                if ((mask(y,x).eqv..false.).and.(data_in(y,x).ne.nodata) &
                    .and.(data_in(y,x).ne.0)) then
                    num_unique = num_unique+1
                    where (data_in.eq.data_in(y,x))
                        mask = .true.
                    end where
                end if
            end do
        end do

        allocate(reclass_table(num_unique, 2))
        allocate(data_out(size(data_in,1), size(data_in,2)))
        data_out=0
        num_unique = 0

        do y = 1, size(data_in, 1)
            do x = 1, size(data_in, 2)
                if ((data_out(y,x).eq.0).and.(data_in(y,x).ne.nodata) &
                    .and.(data_in(y,x).ne.0)) then
                    num_unique = num_unique+1
                    reclass_table(num_unique,1) = num_unique
                    reclass_table(num_unique,2) = data_in(y,x)
                    where (data_in.eq.data_in(y,x))
                        data_out = num_unique
                    end where
                end if
            end do
        end do

    end subroutine reclass_grid

    subroutine calc_topo_class_frac(hru_class, j, nodata, fn_log)

        use dta_qsort
        use dta_preprocess

        implicit none

        type(hru_class_struct) :: hru_class(8)
        double precision, intent(in) :: nodata
        integer, intent(in) :: j, fn_log

        ! locals
        integer :: f
        double precision :: csum_frac
        double precision, allocatable, dimension(:) :: data_list

        data_list = pack(hru_class(j)%data, hru_class(j)%data.ne.nodata)

        ! Sort the data_list into ascending order
        call QsortG(data_list)

        csum_frac = 0

        ! Allocate the class array and set it to 0 for the catahcment
        allocate (hru_class(j)%class_array(size(hru_class(j)%data, 1), size(hru_class(j)%data, 2)))
        hru_class(j)%class_array = 0

        write(fn_log, *) hru_class(j)%num_frac, ' classes for ', hru_class(j)%class_type

        do f = 1, hru_class(j)%num_frac

            csum_frac = hru_class(j)%frac(f) + csum_frac
            hru_class(j)%frac_val(f) = data_list(ceiling(csum_frac * size(data_list, 1)))

            if(f.eq.1) then

                where ((hru_class(j)%data.ne.nodata).and.(hru_class(j)%data.le.hru_class(j)%frac_val(f)))
                    hru_class(j)%class_array = f
                end where

            else

                where ((hru_class(j)%data.gt.hru_class(j)%frac_val(f-1)).and.(hru_class(j)%data.le.hru_class(j)%frac_val(f)))
                    hru_class(j)%class_array = f
                end where

            endif

        end do

        write(fn_log,*) 'Fraction values ', hru_class(j)%frac_val

        if(csum_frac.ne.1) then

            print *, csum_frac
            print *, 'Fractions do not sum up to 1 for this class - '
            print *, hru_class(j)%class_type
            STOP

        endif

    end subroutine calc_topo_class_frac


    subroutine add_class_hru_array(hru_class_array, hru_class_meta, class_array, rivs)

        implicit none

        integer, allocatable, intent(in out) :: hru_class_array(:, :)
        integer, allocatable, intent(in out) :: hru_class_meta(:, :)
        integer, dimension(:,:), intent(in) :: class_array, rivs

        ! locals
        integer :: x, y, class_count, num_hru_classes
        integer :: ix(2)
        double precision, allocatable, dimension(:,:) :: class_all
        integer, allocatable, dimension(:,:) :: tmp_array, tmp_meta

        allocate(class_all(maxval(hru_class_array), maxval(class_array)))
        class_all = 0

        do y = 1, size(hru_class_array, 1)
            do x = 1, size(hru_class_array, 2)
                ix(1) = hru_class_array(y,x)
                ix(2) = class_array(y,x)
                ! Skip cells outside the catchment or river cells
                if ((sum(ix).lt.1).or.(rivs(y,x).gt.0)) then
                    cycle
                ! error check to make sure all catchment outlines match
                elseif (any(ix.eq.0)) then
                    print *, 'WARNING: catchment outlines do not match'
                    print *, 'cant calculate hrus for full catchment outline'
                    print *, ix(1), ix(2), y, x
                else
                    class_all(ix(1), ix(2)) = class_all(ix(1), ix(2)) + 1
                end if
            end do
        end do

        num_hru_classes = count(class_all.gt.0)

        allocate(tmp_array(size(hru_class_array, 1), size(hru_class_array, 2)))
        allocate(tmp_meta(num_hru_classes, size(hru_class_meta, 2)+1))
        tmp_array = 0
        tmp_meta = 0
        class_count = 0

        do y = 1, size(class_all, 1)
            do x = 1, size(class_all, 2)
                if (class_all(y, x).gt.0) then
                    class_count = class_count + 1;
                    tmp_meta(class_count, 1) = class_count
                    tmp_meta(class_count, 2:size(hru_class_meta, 2)) = hru_class_meta(y, 2:size(hru_class_meta, 2))
                    tmp_meta(class_count, size(hru_class_meta,2)+1) = x
                    where ((hru_class_array.eq.y).and.(class_array.eq.x).and.(rivs.eq.0))
                        tmp_array = class_count
                    end where
                end if
            end do
        end do

        hru_class_array = tmp_array
        deallocate(tmp_array)

        deallocate(hru_class_meta)
        allocate(hru_class_meta(size(tmp_meta, 1), size(tmp_meta, 2)))
        hru_class_meta = tmp_meta

    end subroutine add_class_hru_array

    subroutine write_hru_meta(hru_class_meta, hru_properties, hru_class, temp_filename, gauge)

        use dta_preprocess
        use dta_qsort

        implicit none

        type(hru_class_struct) :: hru_class(8)
        integer, allocatable, intent(in) :: hru_class_meta(:, :)
        double precision, allocatable, intent(in) :: hru_properties(:, :)

        integer, allocatable, dimension(:) :: atb_ix
        integer, allocatable, dimension(:,:) :: hru_class_meta_file
        character(len=64), allocatable, dimension(:) :: col_headers
        character(len=1024) :: temp_filename, line_format, gauge
        integer :: num_layers, num_hrus
        integer :: i, k, ix, ioerr


        allocate (atb_ix(maxval(hru_class_meta(:,1))))
        atb_ix = hru_class_meta(:, 1)

        ! HRU's are sorted by increasing mean atb
        call Qsortix(atb_ix, hru_properties(:, 2))

        num_layers = 2

        do i = 1, 7
            if (len_trim(hru_class(i)%fn_end).gt.0) then
                if (i.lt.4) then
                    num_layers=num_layers+1
                else
                    num_layers=num_layers+2
                end if
            end if
        end do

        num_hrus = maxval(hru_class_meta(:, 1))

        allocate(hru_class_meta_file(num_hrus, num_layers+1))
        allocate(col_headers(num_layers+3))

        col_headers(1) = 'HRU_ID'
        col_headers(2) = 'NUM_CELLS'
        col_headers(3) = 'MEAN_ATB'
        col_headers(4) = 'GAUGE_ID'
        col_headers(size(hru_class_meta, 2)+3) = 'GAUGE_NUM'

        hru_class_meta_file(:, 1:2) = hru_class_meta(:, 1:2)

        do k = 1, num_hrus
            hru_class_meta_file(k, size(hru_class_meta, 2) + 1) = hru_class(8)%meta(hru_class_meta(k,2), 2)
        end do

        num_layers = 2

        if (len_trim(hru_class(1)%fn_end).gt.0) then
            num_layers = num_layers+1
            col_headers(num_layers+2) = 'AREA_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
        end if

        if (len_trim(hru_class(2)%fn_end).gt.0) then
            num_layers = num_layers+1
            col_headers(num_layers+2) = 'SLOPE_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
        end if

        if (len_trim(hru_class(3)%fn_end).gt.0) then
            num_layers = num_layers+1
            col_headers(num_layers+2) = 'ELEV_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
        end if

        ix = 1

        if (len_trim(hru_class(4)%fn_end).gt.0) then

            num_layers = num_layers+1
            col_headers(num_layers+2) = 'RAIN_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
            ix = ix+1
            col_headers(size(hru_class_meta, 2)+ix+2) = 'RAIN_NUM'
            do k = 1, num_hrus
                hru_class_meta_file(k, size(hru_class_meta, 2) + ix) = hru_class(4)%meta(hru_class_meta(k,num_layers), 2)
            end do

        end if

        if (len_trim(hru_class(5)%fn_end).gt.0) then

            num_layers = num_layers+1
            col_headers(num_layers+2) = 'PET_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
            ix = ix+1
            col_headers(size(hru_class_meta, 2)+ix+2) = 'PET_NUM'
            do k = 1, num_hrus
                hru_class_meta_file(k, size(hru_class_meta, 2) + ix) = hru_class(5)%meta(hru_class_meta(k,num_layers), 2)
            end do

        end if

        if (len_trim(hru_class(6)%fn_end).gt.0) then

            num_layers = num_layers+1
            col_headers(num_layers+2) = 'MODELSTRUCT_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
            ix = ix+1
            col_headers(size(hru_class_meta, 2)+ix+2) = 'MODELSTRUCT_NUM'
            do k = 1, num_hrus
                hru_class_meta_file(k, size(hru_class_meta, 2) + ix) = hru_class(6)%meta(hru_class_meta(k,num_layers), 2)
            end do

        end if

        if (len_trim(hru_class(7)%fn_end).gt.0) then

            num_layers = num_layers+1
            col_headers(num_layers+2) = 'PARAM_ID'
            hru_class_meta_file(:, num_layers) = hru_class_meta(:, num_layers)
            ix = ix+1
            col_headers(size(hru_class_meta, 2)+ix+2) = 'PARAM_NUM'
            do k = 1, num_hrus
                hru_class_meta_file(k, size(hru_class_meta, 2) + ix) = hru_class(7)%meta(hru_class_meta(k,num_layers), 2)
            end do

        end if

        hru_class_meta_file(:, 2:size(hru_class_meta_file, 2)) = hru_class_meta_file(atb_ix, 2:size(hru_class_meta_file, 2))

        open (997, file = temp_filename, status = 'unknown', iostat=ioerr)

        if(ioerr/=0) then
            print *,ioerr
            print*,'error opening output file: ', trim(temp_filename)
            print*,'ensure the directory exists and correct write permissions are set'
            stop
        endif

            ! Print out header file
        write(997, '(A, A)') '! HRU Meta File for Gauge ', trim(gauge)
        write(997, '(A, A, I0)') 'HRUs', char(9), size(hru_class_meta_file, 1)
        write(997, '(A, A, I0)') 'NUM_COLS', char(9), size(col_headers)
        write(997, '(A)') '!'


        do i = 1, size(col_headers)
            if (i.lt.size(col_headers)) then
                write(997, '(A, 2X)', advance = "no") trim(col_headers(i))
            else
                write(997, '(A, 2X)') trim(col_headers(i))
            end if
        end do

        write(line_format,'(A,I0,A)') '(I0, A, I0, A F7.3, ', size(hru_class_meta_file, 2)-1, '(I9))'

        do i = 1, size(hru_class_meta_file, 1)
            write(997, '(I0, A, I0, A, F0.5)', advance = "no") hru_class_meta_file(i, 1), char(9), &
                int(hru_properties(atb_ix(i), 1)), char(9), hru_properties(atb_ix(i), 2)
            do k = 1, size(hru_class_meta_file, 2)-1
                if (k.lt.size(hru_class_meta_file, 2)-1) then
                    write(997, '(A, I0)', advance = "no") char(9), hru_class_meta_file(i, k+1)
                else
                    write(997, '(A, I0)') char(9), hru_class_meta_file(i, k+1)
                end if
            end do
        end do

        close(997)

    end subroutine write_hru_meta

    subroutine sharing_bw_hrus(hru_class_array, hru_share_percent, riv_share_percent, &
        hru_properties, atb, rivs, dem, area, nodata)

        use dta_preprocess

        implicit none

        integer, intent(in) :: hru_class_array(:, :)
        integer, intent(in) :: rivs(:, :)
        double precision, intent(in) :: nodata
        double precision, intent(in) :: dem(:, :), atb(:, :), area(:,:)
        double precision, allocatable, intent(out) :: hru_share_percent(:,:)
        double precision, allocatable, intent(out) :: hru_properties(:, :)
        double precision, allocatable, intent(out) :: riv_share_percent(:,:)

        ! Locals
        integer :: y, x, xl, yl, i
        double precision :: wfp(3,3)
        double precision ::routefrac, tanB
        double precision, allocatable :: hru_share(:,:)
        double precision, allocatable :: riv_class_share(:,:)

        allocate (hru_share(maxval(hru_class_array), maxval(hru_class_array)))
        hru_share = 0
        allocate (hru_share_percent(maxval(hru_class_array), maxval(hru_class_array)))
        hru_share_percent = 0
        allocate (riv_class_share(maxval(hru_class_array), maxval(rivs)))
        riv_class_share = 0
        allocate (riv_share_percent(maxval(hru_class_array), maxval(rivs)))
        riv_share_percent = 0
        allocate (hru_properties(maxval(hru_class_array), 3))
        hru_properties = 0

        !  Loop through each cell of the whole array

        do y = 2, size(dem, 1) -1
            do x = 2, size(dem, 2)-1

                wfp = 0

                ! Check the cell isn't outside the catchment and it isn't a river cell
                if ((dem(y,x).ne.nodata).and.(rivs(y,x).eq.0))then

                    ! Look to all the neighbours of that cell
                    do yl = -1, 1
                        do xl = -1, 1

                            if((xl==0.and.yl==0)) then
                                cycle
                            endif

                            if ((x+xl) < 1 .or. (y+yl) < 1 .or. (x+xl) > size(dem, 2).or. (y+yl) > size(dem, 1) ) then
                                cycle
                            endif

                            if ((dem(y+yl,x+xl)/=nodata).and.(dem(y+yl,x+xl) < dem(y,x))) then

                                if(xl == 0 .or. yl == 0) then
                                    routefrac = 0.5
                                else
                                    routefrac = 0.25
                                endif

                                wfp(yl+2, xl+2) = (dem(y, x) - dem(y+yl, x+xl)) * routefrac

                            endif

                        end do
                    end do

                    tanB = sum(wfp)

                    do yl = -1, 1
                        do xl = -1, 1

                            if((xl==0.and.yl==0)) then
                                cycle
                            endif

                            if((wfp(yl+2,xl+2).gt.0).and.(rivs(y+yl,x+xl).eq.0))then
                                ! Multiply the weighted flow proportions by the accumulated area for that cell and sum it for each class share
                                hru_share(hru_class_array(y,x), hru_class_array(y+yl,x+xl)) = &
                                    hru_share(hru_class_array(y,x), hru_class_array(y+yl,x+xl)) &
                                    + ((wfp(yl+2,xl+2)/tanB)*area(y,x))
                            elseif ((wfp(yl+2,xl+2).gt.0).and.(rivs(y+yl,x+xl).gt.0)) then
                                ! Keep a note of the calculation above for each river ID in each class - used below!
                                riv_class_share(hru_class_array(y,x), int(rivs(y+yl,x+xl))) = &
                                    riv_class_share(hru_class_array(y,x), int(rivs(y+yl,x+xl))) &
                                    + (wfp(yl+2,xl+2)/tanB)*area(y,x)
                            endif

                        end do
                    end do

                end if

            end do
        end do

        do y = 1, maxval(hru_class_array)

            hru_share_percent(y,:) = hru_share(y,:) / (sum(hru_share(y, :)) + sum(riv_class_share(y, :)))

        end do

        ! Calculate the percentage sharing into each of the rivers for each class - gets written out into DTA dynatop
        do y = 1, maxval(hru_class_array)
            do x = 1, int(maxval(rivs))

                riv_share_percent(y,x) = riv_class_share(y,x) / (sum(hru_share(y, :)) + sum(riv_class_share(y, :)))

            end do
        end do

        !  Calculate the mean atb for each hru class - this is used to sort the classes

        do i = 1, maxval(hru_class_array)

            ! This section gets the eight values for the first liune of each HRU in the file output.

            ! Calculate the percentage of cells in each HRU
            hru_properties(i, 1) = dble(count(hru_class_array.eq.i))
            ! Calculate mean atb
            hru_properties(i, 2) = sum(atb, hru_class_array.eq.i) / count(hru_class_array.eq.i)
            ! Number of classes this HRU is sharing to
            hru_properties(i, 3) = count((hru_share_percent(i,:)).gt.0)

        end do

    end subroutine sharing_bw_hrus


    subroutine write_dta_file_hru_share(temp_filename, hru_properties, hru_share_percent, riv_share_percent, &
        hru_class_array, hru_class_array_atb, hru_class_meta, hru_class, gauge, catch_area)

        use dta_qsort
        use dta_preprocess

        implicit none

        type(hru_class_struct) :: hru_class(8)
        double precision, intent(in) :: hru_share_percent(:,:), hru_properties(:,:), riv_share_percent(:,:)
        double precision, intent(in) :: catch_area
        integer, intent(in) :: hru_class_array(:, :)
        integer, allocatable, intent(in) :: hru_class_meta(:, :)
        integer, allocatable, intent(out), dimension(:,:) :: hru_class_array_atb

        ! Local
        integer :: i, k, l, ioerr, num_class_layers
        integer, allocatable, dimension(:) :: atb_ix
        double precision :: sum_cont
        character(len=1024) :: temp_filename, gauge
        double precision :: riv_share_frac(size(riv_share_percent, 2))

        allocate (atb_ix(maxval(hru_class_array)))

        allocate (hru_class_array_atb(size(hru_class_array, 1), size(hru_class_array, 2)))
        hru_class_array_atb = 0

        open(997, file = temp_filename, status="unknown", action="write", iostat=ioerr)

        if(ioerr/=0) then
            print*,'error opening output file: ', trim(temp_filename)
            print*,'ensure the directory exists and correct write permissions are set'
            stop
        endif

        ! Print out header file
        write(997, '(A, A)') 'HRU Flux File for Gauge ', trim(gauge)
        write(997, '(A, A, I0)') 'HRUs', char(9), maxval(hru_class_array)
        write(997, '(A, A, F0.2)') 'Catchment_Area', char(9), catch_area


        num_class_layers = 1
        do i = 1, 7
            if (len_trim(hru_class(i)%fn_end).gt.0) then
                num_class_layers = num_class_layers + 1
            end if
        end do
        write(997, '(A, A, I0)') 'Classifying_Layers', char(9), num_class_layers
        write(997, '(A, A, I0)') 'Gauges', char(9), maxval(hru_class_meta(:, 2))

        num_class_layers = 2
        if (len_trim(hru_class(1)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'Area_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if
        if (len_trim(hru_class(2)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'Slope_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if
        if (len_trim(hru_class(3)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'Elevation_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if
        if (len_trim(hru_class(4)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'Rainfall_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if
        if (len_trim(hru_class(5)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'PET_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if
        if (len_trim(hru_class(6)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'Model_Structure_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if
        if (len_trim(hru_class(7)%fn_end).gt.0) then
            num_class_layers = num_class_layers + 1
            write(997, '(A, A, I0)') 'Parameter_Classes', char(9), maxval(hru_class_meta(:, num_class_layers))
        end if

        ! HRU's are sorted by increasing mean atb
        do i = 1, maxval(hru_class_array)
            atb_ix(i) = i
        end do

        call Qsortix(atb_ix, hru_properties(:, 2))

        ! HRU to HRU FLUXES first

        write(997, *) ''
        write(997, *) 'HRU FLUXES'
        write(997, *) ''

        do i = 1, maxval(hru_class_array)
            sum_cont = 0
            write (997, 73) 'HRU ', i, ' Sharing to ', int(hru_properties(atb_ix(i), 3)), ' HRUs'

            do k = 1, maxval(hru_class_array)

                if (hru_share_percent(atb_ix(i), atb_ix(k)).gt.0) then

                    sum_cont = sum_cont + hru_share_percent(atb_ix(i), atb_ix(k))
                    write (997, 75) k, char(9), hru_share_percent(atb_ix(i), atb_ix(k)), char(9), sum_cont

                end if

            end do

            where (hru_class_array.eq.atb_ix(i))
                hru_class_array_atb = i
            end where

            sum_cont = sum_cont + sum(riv_share_percent(atb_ix(i), :))
            write(997, 77) 'RIVS', char(9), sum(riv_share_percent(atb_ix(i), :)), char(9), sum_cont

        end do

        write(997, *) ''
        write(997, *) 'HRU RIV FLUXES'
        write(997, *) ''

        do i = 1, maxval(hru_class_array)
            sum_cont = 0
            riv_share_frac = 0

            do l = 1, size(riv_share_percent, 2)
                if (riv_share_percent(atb_ix(i), l).gt.0) then
                    riv_share_frac(l) = riv_share_percent(atb_ix(i), l) / sum(riv_share_percent(atb_ix(i), :))
                else
                    riv_share_frac(l) = 0
                end if
            end do

            write (997, 73) 'HRU ', i, ' Sharing to ', count(riv_share_percent(atb_ix(i), :).gt.0), ' RIVS'

            do l = 1, size(riv_share_frac, 1)
                if (riv_share_frac(l).gt.0) then
                    sum_cont = sum_cont+riv_share_frac(l)
                    write(997, 76) l, char(9), riv_share_frac(l), char(9), sum_cont
                end if
            end do

            where (hru_class_array.eq.atb_ix(i))
                hru_class_array_atb = i
            end where

        end do

        close(997)

73      FORMAT(A, i0, A, i0, A)
75      FORMAT(i0, A, f0.6, A, f0.6)
76      FORMAT(i0, A, f0.6, A, f0.6)
77      FORMAT(A, A, f0.6, A, f0.6)

    end subroutine write_dta_file_hru_share


end module dta_calculate_hrus
