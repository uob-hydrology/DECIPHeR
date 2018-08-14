program calculate_hrus

    ! GC 8th April 2016
    ! GC 30th October 2017

    ! This code :
    ! 1.  Subsets your principal basin from the original DTA output (which may be at national scale)
    ! 2.  Creates your

    ! STEP 1.  Get all the data for a specific catchment

    ! Use subroutines in modules
    use dta_utility
    use dta_qsort
    use dta_calculate_hrus
    use dta_preprocess

    implicit none

    ! Declare variables

    type(hru_class_struct) :: hru_class(8)
    double precision, allocatable, dimension(:,:) :: dem, atb, area, catch_data
    double precision, allocatable, dimension(:,:) :: hru_share_percent, riv_share_percent, hru_properties
    integer, allocatable, dimension(:,:) :: hru_class_array, hru_class_array_atb, rivs

    integer :: ncols, nrows, i
    double precision :: xllcorner, yllcorner, cellsize, nodata, catch_area
    character(len=1024) :: gauge, temp_filename, folder, arg
    character(len=1024) :: gaugelist, outputfolder, hru_file
    double precision, allocatable, dimension(:,:) :: gauges_all
    integer, allocatable, dimension(:,:) :: hru_class_meta, rivs_meta
    logical :: input_is_valid
    integer :: fn_log, ioerr, k

    input_is_valid = .true.

    ! Set these variables to blank otherwise the error checking won't pick them up!
    gauge = ''
    folder = ''
    hru_file = ''

    ! This gets the command argument from the command line - specifies the gauge ID you want to build the HRUs for
    temp_filename = 'DTA_calculate_hrus.log'
    fn_log = 999
    open(fn_log, file = temp_filename, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(temp_filename)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(fn_log,*) '--- calculate_hrus.f90 ---'
    write(fn_log,*) ''
    print *, '--- Starting calculate_hrus ---'

    i=0
    outputfolder=''

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-gaugelist')) then
            CALL get_command_argument(i+1, gaugelist)
        else if (are_equal(arg, '-hru_class_file')) then
            CALL get_command_argument(i+1, hru_file)
        else if (are_equal(arg, '-output_folder')) then
            CALL get_command_argument(i+1, outputfolder)
        endif
        i = i + 1
    enddo


    if (len_trim(gaugelist) == 0) then
        print *, '-gaugelist not specified'
        input_is_valid = .false.
    else if (len_trim(hru_file) == 0) then
        print *, '-hru_class_file not specified'
        input_is_valid = .false.
    else if (len_trim(outputfolder) == 0) then
        print *, '-output_folder not specified'
        input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'calculate_hrus command options '
        print *, '-gaugelist <gauge list>   input text file list of gauges'
        print *, '-hru_class_file <hru_class_file>   HRU Classifiers File'
        print *, '-outputfolder <output folder>   select output folder pathway'
        stop
    endif

    ! Step 1.  Read the gauge list and the HRU classifier file
    call read_hruclass_file(hru_class, hru_file, fn_log)
    gaugelist = trim(gaugelist)
    call read_numeric_list(gaugelist, 1, 1, gauges_all)
    write(fn_log,*) 'Caclulating HRUs for ', size(gauges_all, 1), ' gauges'

    ! Step 2. Loop through each catchment
    do i = 1, size(gauges_all, 1)

        write(gauge, *) int(gauges_all(i,1))
        gauge = adjustl(gauge)
        PRINT *, trim(gauge)

        write(fn_log,*) ''
        write(fn_log,*) gauge

        ! Step 4 - Determine HRUs - loop through all the different classifiers and each time update HRUs

        ! Always split up the gauge by sub-catchments and always read rivers layer
        ! If no other hru classifiers are specified then this is the model running in a 'lumped' mode for
        ! a single HRU per catchment

        write(fn_log,*) 'Classifying HRUs'

        temp_filename = trim(outputfolder)//trim(gauge)//'_mask.asc'
        call read_ascii_grid(temp_filename, catch_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
        call reclass_grid(int(catch_data), int(nodata), hru_class(8)%class_array,  hru_class(8)%meta)

        allocate(hru_class_meta(size(hru_class(8)%meta, 1), 2))
        allocate(hru_class_array(size(hru_class(8)%class_array, 1), size(hru_class(8)%class_array, 2)))
        hru_class_meta(:, 1) = hru_class(8)%meta(:, 1)
        hru_class_meta(:, 2) = hru_class(8)%meta(:, 1)

        hru_class_array = hru_class(8)%class_array
        deallocate(catch_data)

        temp_filename = trim(outputfolder)//'/'//trim(gauge)//'_riv.asc'
        call read_ascii_grid(temp_filename, catch_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        ! Reclass rivers the same as the gauges - can't just reclass grid as can mis-classify the river
        allocate(rivs(size(catch_data, 1), size(catch_data, 2)))
        allocate(rivs_meta(size(hru_class(8)%meta, 1), size(hru_class(8)%meta, 2)))

        rivs = 0
        do k = 1, size(hru_class(8)%meta, 1)
            where (int(catch_data).eq.hru_class(8)%meta(k, 2))
                rivs = hru_class(8)%meta(k,1)
            end where
        end do
        rivs_meta = hru_class(8)%meta
        deallocate(catch_data)

        catch_area = (count(hru_class_array.gt.0) * cellsize**2)/1000000

        do k = 1, 7

            if (len_trim(hru_class(k)%fn_end).gt.0) then

                ! This hru class has data and needs to be added as a classifier
                ! Different functions for topographic data which are split up by fractional units

                if (hru_class(k)%class_type(1:4).eq.'TOPO') then
                    ! Read in the data and make class array by reclassifying data by fractional values
                    temp_filename = trim(outputfolder)//trim(gauge)//trim(hru_class(k)%fn_end)//'.asc'
                    call read_ascii_grid(temp_filename, hru_class(k)%data, ncols, nrows, xllcorner, &
                        yllcorner, cellsize, nodata)
                    ! Reclassify data by fractional values and then add to HRU class array
                    call calc_topo_class_frac(hru_class, k, nodata, fn_log)
                    call add_class_hru_array(hru_class_array, hru_class_meta, hru_class(k)%class_array, rivs)
                    
                else
                    ! Read in the data, reclass and add it to hru layer
                    temp_filename = trim(outputfolder)//trim(gauge)//trim(hru_class(k)%fn_end)//'.asc'
                    call read_ascii_grid(temp_filename, hru_class(k)%data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)
                    call reclass_grid(int(hru_class(k)%data), int(nodata), hru_class(k)%class_array, hru_class(k)%meta)
                    call add_class_hru_array(hru_class_array, hru_class_meta, hru_class(k)%class_array, rivs)

                end if

                deallocate(hru_class(k)%data, hru_class(k)%class_array)

            end if

        end do

    ! STEP 4.  Calculate sharing between HRUs

    temp_filename = trim(outputfolder)//'/'//trim(gauge)//'_atb.asc'
    call read_ascii_grid(temp_filename, atb, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

    temp_filename = trim(outputfolder)//'/'//trim(gauge)//'_area.asc'
    call read_ascii_grid(temp_filename, area, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

    temp_filename = trim(outputfolder)//'/'//trim(gauge)//'_dem.asc'
    call read_ascii_grid(temp_filename, dem, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

    write(fn_log,*) 'Calculating fluxes between HRUs'

    call sharing_bw_hrus(hru_class_array, hru_share_percent, riv_share_percent, &
        hru_properties, atb, rivs, dem, area, nodata)

    temp_filename = trim(outputfolder)//'/'//trim(gauge)//'_dyna_hru.dat'
    !  Write out each HRU and write out the class array

    write(fn_log,*) 'Writing dyna_hru file'
    call write_dta_file_hru_share(temp_filename, hru_properties, hru_share_percent, riv_share_percent, &
        hru_class_array, hru_class_array_atb, hru_class_meta, hru_class, gauge, catch_area)
    ! Write out HRU class array and the meta table

    write(fn_log,*) 'Writing HRU ascii file'
    temp_filename = trim(outputfolder)//trim(gauge)//'_hru_array.asc'
    call write_ascii_grid_int(temp_filename, hru_class_array_atb, ncols, nrows, xllcorner, yllcorner, cellsize, 0)

    write(fn_log,*) 'Writing HRU meta file'
    temp_filename = trim(outputfolder)//trim(gauge)//'_hru_meta.dat'
    call write_hru_meta(hru_class_meta, hru_properties, hru_class, temp_filename, gauge)

    deallocate(hru_class_meta, hru_class_array, rivs, rivs_meta, area)
    deallocate(dem, atb, hru_properties, hru_class_array_atb, hru_share_percent, riv_share_percent)

    end do

    print *, '--- Finished calculate_hrus ---'
    write(fn_log, *) ''
    write(fn_log,*) 'Successfully finished calculate_hrus.f90'
    close(fn_log)


end program calculate_hrus
