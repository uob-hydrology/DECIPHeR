program preprocess

    ! GC 12th May 2016

    ! This code preprocesses all the necessary data for a model run for a particular catchment

    use dta_preprocess
    use dta_utility

    implicit none

    ! Declare variables

    character(len=1024) :: gauge, catch_fn, arg, temp_fn, gaugelist, dta_data_fn
    character(len=1024) :: root_fn, hru_file
    character(len=1024) :: flowpoint_fn, rivdata_fn, flowconn_fn
    character(len=1024) :: catchmaskfolder, outputfolder, catchmaskfile
    double precision, allocatable, dimension(:,:) :: catch_mask
    double precision, allocatable, dimension(:,:) :: flow_conn, riv_point_data, flow_conn_catch, riv_point_data_catch
    double precision, allocatable, dimension(:,:) :: flow_point, flow_point_catch, catch_data
    double precision, allocatable, dimension(:,:) :: dta_data
    double precision, allocatable, dimension(:,:) :: gauges_all
    double precision :: xllcorner, yllcorner, cellsize, nodata, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata
    integer :: ncols, nrows, i, catch_ncols, catch_nrows, fn_log, ioerr
    logical :: input_is_valid, writecatchfiles

    type(hru_class_struct) :: hru_class(8)


    ! Set these variables to blank otherwise the error checking won't pick them up!
    gaugelist = ''
    root_fn = ''
    outputfolder = ''
    catchmaskfolder = ''

    !Initialise variables
    writecatchfiles = .true.
    input_is_valid = .true.

    temp_fn = 'DTA_preprocess.log'
    fn_log = 999
    open(fn_log, file = temp_fn, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(temp_fn)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(fn_log,*) '--- preprocess.f90 ---'
    write(fn_log,*) ''
    print *, '--- Starting preprocess ---'

    i=0

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-gaugelist')) then
            CALL get_command_argument(i+1, gaugelist)
        else if (are_equal(arg, '-hru_class_file')) then
            CALL get_command_argument(i+1, hru_file)
        else if (are_equal(arg, '-root_fn')) then
            CALL get_command_argument(i+1, root_fn)
        else if (are_equal(arg, '-catchmask_folder')) then
            CALL get_command_argument(i+1, catchmaskfolder)
        else if (are_equal(arg, '-output_folder')) then
            CALL get_command_argument(i+1, outputfolder)
        else if (are_equal(arg, '-disable_writecatchfiles')) then
            writecatchfiles = .false.
        endif
        i = i + 1
    enddo

    if (len_trim(gaugelist) == 0) then
        print *, '-gaugelist not specified'
        input_is_valid = .false.
    else if (len_trim(hru_file) == 0) then
        print *, '-hru_class_file not specified'
        input_is_valid = .false.
    else if (len_trim(root_fn) == 0) then
        print *, '-root_fn not specified'
        input_is_valid = .false.
    else if (len_trim(outputfolder) == 0) then
        print *, '-output_folder not specified'
        input_is_valid = .false.
    else if (len_trim(catchmaskfolder) == 0) then
        print *, '-catchmask_folder not specified'
        input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'dta_preprocess command options '
        print *, '-gaugelist <gauge list>   input text file list of gauges'
        print *, '-hru_class_file <hru_class_file>   HRU Classifiers File'
        print *, '-catchmask_folder <catchment mask folder>  select catchment mask folder'
        print *, '-root_fn <root_fn file>   select dem'
        print *, '-output_folder <output folder>   select output folder pathway'
        print *, '-disable_writecatchfiles only include flag if you do not want catchment dem,atb,area,slope written out'
        stop
    endif

    ! STEP 1 - read in the HRU classifier file and the gauges you need to process
    call read_hruclass_file(hru_class, hru_file, fn_log)

    gaugelist = trim(gaugelist)
    call read_numeric_list(gaugelist, 1, 1, gauges_all)

    write(fn_log,*) 'Pre-processing ', size(gauges_all, 1), ' gauges'

    ! STEP 2 - subset all the necessary data - can be switched off if you have already done this

    if (writecatchfiles.eqv..true.) then

        ! Subset the DEM, slope, rivers, accumulated area and topographic index and catch mask
        ! Subset the flow routing files

        ! DEM
        PRINT *, 'Subsetting DEMs'
        write(fn_log,*) 'Subsetting DEMS for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = trim(root_fn)//'.asc'
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_dem.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 3)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

        ! Rivers
        PRINT *, 'Subsetting Rivers'

        write(fn_log,*) 'Subsetting Rivers for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = trim(root_fn)//'_riv_id_check.asc'
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_riv.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

                deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

        ! Gauge mask and routing files
        PRINT *, 'Subsetting Gauge Mask and Routing Files'
        write(fn_log,*) 'Subsetting Gauge Mask and Routing Files for ', size(gauges_all, 1), ' gauges'

        dta_data_fn = trim(root_fn)//'_mask.asc'
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)


        temp_fn = trim(root_fn)//'_flow_conn.txt'
        call read_numeric_list(temp_fn, 10, 1, flow_conn)

        ! Read in the river data file (need distance)
        temp_fn = trim(root_fn)//'_river_data.txt'
        call read_numeric_list(temp_fn, 8, 1, riv_point_data)

        ! Read in the flow point data file
        temp_fn = trim(root_fn)//'_flow_point.txt'
        call read_numeric_list(temp_fn, 6, 1, flow_point)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_mask.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

            call subset_route_files(catch_data, flow_conn, flow_point, riv_point_data, &
                flow_conn_catch, flow_point_catch, riv_point_data_catch)

            ! Write the flow connectivity file
            flowconn_fn = trim(outputfolder)//'/'//trim(gauge)//'_flow_conn.dat'
            rivdata_fn = trim(outputfolder)//'/'//trim(gauge)//'_riv_data.dat'
            flowpoint_fn = trim(outputfolder)//'/'//trim(gauge)//'_flow_point.dat'

            call write_routingfiles(flowconn_fn, flow_conn_catch, rivdata_fn, riv_point_data_catch, &
                flowpoint_fn, flow_point_catch)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

        ! Topographic Index - maske sure to use the one that DOES NOT accumulate downstream
        PRINT *, 'Subsetting Topographic Index'
        write(fn_log,*) 'Subsetting Topographic Index for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = trim(root_fn)//'_riv_mask_atb.asc'
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'/gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_atb.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

        ! Slope
        PRINT *, 'Subsetting Slope'
        write(fn_log,*) 'Subsetting Slope for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = trim(root_fn)//'_riv_mask_mfd_slope.asc'
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_slope.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 6)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

        ! Accumulated area
        PRINT *, 'Subsetting Accumulated Area'
        write(fn_log,*) 'Subsetting Accumulated Area for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = trim(root_fn)//'_riv_mask_area.asc'
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_area.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 3)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

    end if

    ! Next look at HRU Classifiers file what else needs to be subsetted to classify your HRUs

    ! Rainfall
    if (len_trim(hru_class(4)%filename).gt.0) then

        PRINT *, 'Subsetting Rainfall Grid'
        write(fn_log,*) 'Subsetting Rainfall Grid for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = hru_class(4)%filename
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//trim(hru_class(4)%fn_end)//'.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

    end if

    ! PET
    if (len_trim(hru_class(5)%filename).gt.0) then

        PRINT *, 'Subsetting PET Grid'
        write(fn_log,*) 'Subsetting PET Grid for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = hru_class(5)%filename
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//trim(hru_class(5)%fn_end)//'.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

    end if

    ! Model Structure
    if (len_trim(hru_class(6)%filename).gt.0) then

        PRINT *, 'Subsetting Model Structure'
        write(fn_log,*) 'Subsetting Model Structure Grid for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = hru_class(6)%filename
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//trim(hru_class(6)%fn_end)//'.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

    end if

    ! Parameters
    if (len_trim(hru_class(7)%filename).gt.0) then

        PRINT *, 'Subsetting Parameter Grid'
        write(fn_log,*) 'Subsetting Parameter Grid for ', size(gauges_all, 1), ' gauges'
        dta_data_fn = hru_class(7)%filename
        call read_ascii_grid(dta_data_fn, dta_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        do i = 1, size(gauges_all, 1)

            write(gauge, *) int(gauges_all(i,1))
            gauge = adjustl(gauge)

            catchmaskfile = trim(catchmaskfolder)//'gauge_'//trim(gauge)//'.asc'

            ! Read in the catchment mask
            call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
                catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            call subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
                catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

            catch_fn = trim(outputfolder)//'/'//trim(gauge)//trim(hru_class(7)%fn_end)//'.asc'
            call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
                catch_cellsize, real(-9999, 8), 4)

            deallocate(catch_mask, catch_data)

        end do

        deallocate(dta_data)

    end if

    print *, '--- Finished preprocess ---'
    write(fn_log,*) ''
    write(fn_log,*) 'Successfully finished preprocess.f90'
    close(fn_log)

end program preprocess
