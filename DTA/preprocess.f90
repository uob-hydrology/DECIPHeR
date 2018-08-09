program preprocess

    ! GC 12th May 2016

    ! This code preprocesses all the necessary data for a model run for a particular catchment

    ! NB : Need to change all these to the projects directory - do this once finished
    ! use subroutines in module
    use dta_preprocess
    use dta_utility

    implicit none

    ! Declare variables

    character(len=1024) :: gauge, uk_fn, catch_fn, arg, temp_fn
    character(len=1024) :: catchmaskfile, demfile, outputfolder, raingridfile
    double precision, allocatable, dimension(:,:) :: catch_mask, gauges, catch_area
    double precision, allocatable, dimension(:,:) :: flow_conn, riv_point_data, flow_conn_catch, riv_point_data_catch
    double precision, allocatable, dimension(:,:) :: flow_point, flow_point_catch
    double precision, allocatable, dimension(:,:) :: catch_data, rain_grid, river_data
    integer, allocatable, dimension(:,:) :: rain_grid_dta, catch_data_dta, river_data_dta, fu_grid, fu_table
    double precision :: catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata, nodata_w
    integer :: catch_ncols, catch_nrows, i
    logical :: input_is_valid, writecatchfiles
    character(64) :: col_headers(6)

    ! Set these variables to blank otherwise the error checking won't pick them up!
    gauge = ''
    demfile = ''
    outputfolder = ''
    catchmaskfile = ''
    demfile = ''
    raingridfile = ''

    !Initialise variables
    nodata_w = -9999
    writecatchfiles = .true.
    input_is_valid = .true.

    ! This gets the command argument from the command line - specifies the gauge ID you want to build the HRUs for
    i = 0
    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-gauge')) then
            CALL get_command_argument(i+1, gauge)
        else if (are_equal(arg, '-catchmaskfile')) then
            CALL get_command_argument(i+1, catchmaskfile)
        else if (are_equal(arg, '-disable_writecatchfiles')) then
            writecatchfiles = .false.
            !CALL get_command_argument(i+1, writecatchfiles)
        else if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, demfile)
        else if (are_equal(arg, '-raingridfile')) then
            CALL get_command_argument(i+1, raingridfile)
        else if (are_equal(arg, '-outputfolder')) then
            CALL get_command_argument(i+1, outputfolder)
        endif
        i = i + 1
    enddo

    if (len_trim(gauge) == 0) then
        print *, '-gauge not specified'
        input_is_valid = .false.
    else if (len_trim(catchmaskfile) == 0) then
        print *, '-catchmaskfile not specified'
        input_is_valid = .false.
    else if (len_trim(demfile) == 0) then
        print *, '-dem not specified'
        input_is_valid = .false.
        else if (len_trim(raingridfile) == 0) then
        print *, '-raingridfile not specified'
        input_is_valid = .false.
    else if (len_trim(outputfolder) == 0) then
        print *, '-outputfolder not specified'
        input_is_valid = .false.
    endif

    if(input_is_valid .eqv. .false.) then
        print *, 'dta_preprocess command options '
        print *, '-gauge <gauge id>   select gauge ID'
        print *, '-catchmaskfile <catchment mask file>  select catchment mask file'
        print *, '-demfile <dem file>   select dem'
        print *, '-outputfolder <output folder>   select output folder pathway'
        print *, '-raingridfile <rain grid file>   select the rain grid you want to use'
        print *, '-disable_writecatchfiles only include flag if you do not want catchment dem,atb,area,slope written out'
        stop
    endif

    !  Step 1. Read in the catchment mask

    catchmaskfile = trim(catchmaskfile)//'.asc'

    call read_ascii_grid(catchmaskfile, catch_mask, catch_ncols, catch_nrows, &
        catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

    !  Step 2.  Cut out DEM, slope, topographic index and accumulated area for the specified catchment and write out the masked variables to file

    !  This takes around 8mins to process so included an option not to wite this out if you just want to re-run the functional units
    if (writecatchfiles.eqv..true.) then

        PRINT *, 'Subsetting DEM, slope, accumulated area and topographic index files'

        ! DEM
        uk_fn = trim(demfile)//'.asc'
        catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_dem.asc'
        call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, &
            catch_xllcorner, catch_yllcorner, catch_nodata, catch_data)
        call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
            catch_cellsize, nodata_w, 3)

        ! Slope
        uk_fn = trim(demfile)//'_mfd_slope.asc'
        catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_slope.asc'
        call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, &
            catch_xllcorner, catch_yllcorner, catch_nodata, catch_data)
        call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
            catch_cellsize, nodata_w, 5)

        ! Topographic index
        uk_fn = trim(demfile)//'_atb.asc'
        catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_atb.asc'
        call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, &
            catch_xllcorner, catch_yllcorner, catch_nodata, catch_data)
        call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
            catch_cellsize, nodata_w, 4)

        ! Accumulated area
        uk_fn = trim(demfile)//'_area.asc'
        catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_area.asc'
        call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_nodata, catch_data)
        call write_ascii_grid(catch_fn, catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
            catch_cellsize, nodata_w, 3)

            PRINT *, 'Created catchment dem, accumulated area, slope, topographic index ascii grids'

    end if

    ! Step 3.  Create functional units file
    !  Options currently to read in rainfall data and catchment data (need to add geology, soils, etc.)
    !  Writes out table that lists each functional unit, UK_Rainfall_ID, Catchment Rainfall_ID, Catchment ID

    PRINT *, 'Getting rain grid, catchment IDs and river data'
    PRINT *, 'Creating functional units'

    ! Get rain data
    uk_fn = trim(raingridfile)//'.asc'
    call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_nodata, rain_grid)

    ! Get Catchment mask data
    uk_fn = trim(demfile)//'_mask.asc'
    call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_nodata, catch_data)

    ! Get rivers data - Add the rivers to table_fu, reclassify the rivers so they go from 1 -> n
    uk_fn = trim(demfile)//'_riv_500m_100m_riv_id_check.asc'
    call subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_nodata, river_data)

    ! Create functional units table and grid and write out to file
    call create_funcunits(rain_grid, catch_data, river_data, fu_table, fu_grid, rain_grid_dta, catch_data_dta, &
        river_data_dta, gauges)

    ! Write out all the results to file

    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_raingrid.asc'
    call write_ascii_grid_int(catch_fn, rain_grid_dta, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
        catch_cellsize, 0)
    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_catch.asc'
    call write_ascii_grid_int(catch_fn, catch_data_dta, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
        catch_cellsize, 0)
    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_funcunits.asc'
    call write_ascii_grid_int(catch_fn, fu_grid, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
        catch_cellsize, 0)
    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_riv.asc'
    call write_ascii_grid_int(catch_fn, river_data_dta, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, &
        catch_cellsize, 0)

    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_funcunit_table.txt'
    open(99, file = catch_fn, status = 'unknown')

    write (99, '(19A)') 'Funcunit_ID', tab, &
        'Catch_ID', tab, &
        'Catch_Num', tab, &
        'Rain_ID', tab, &
        'Rain_Num', tab, &
        'River_ID', tab, &
        'River_Num'

    do i = 1, size(fu_table, 1)
        write (99, 71) fu_table(i, :)
    end do

71  FORMAT(i0, 1x, i0, 1x, i0, 1x, i0, 1x, i0, 1x, i0, 1x, i0)

    PRINT *, 'Created functional units file'

    !  Step 4.  Create the initialisation file
    ! Specifies all the gauges, their downstream gauges and the area of each gauge

    ! Read in the flow connectivity file and catchment areas - this needs to be changed but at the moment we are
    ! not using the initialisation file...
    temp_fn = '/projects/The_Env_Virtual_observatory/DynaTOP_data/MERGED_rivers_uk_50m/catch3/station_river_gauge.txt'
    call read_numeric_list(temp_fn, 5, 1, catch_area)

    temp_fn = trim(demfile)//'_flow_conn.txt'
    call read_numeric_list(temp_fn, 10, 1, flow_conn)

    call create_init_file(gauges, gauge, flow_conn, catch_area, outputfolder)

    ! Read in the river data file (need distance and slope)
    temp_fn = trim(demfile)//'_river_data.txt'
    call read_numeric_list(temp_fn, 6, 1, riv_point_data)

    ! Read in the flow point data file
    temp_fn = trim(demfile)//'_flow_point.txt'
    call read_numeric_list(temp_fn, 6, 1, flow_point)

    call create_flowroute_files(gauges, river_data, flow_conn, flow_point, riv_point_data, &
        flow_conn_catch, flow_point_catch, riv_point_data_catch)

    ! Write the flow connectivity file

    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_flow_conn.dat'

    open(100, file = catch_fn, status = 'unknown')

    write (100, '(19A)') 'Node_ID', tab, &
        'x_easting', tab, &
        'y_northing', tab, &
        'Type', tab, &
        'ID_down', tab, &
        'x_easting2', tab, &
        'y_northing2', tab, &
        'delta_dist', tab, &
        'delta_h', tab, &
        'slope'

    do i = 1 , size(flow_conn_catch, 1)
        write(100, '(I0,A,F0.1,A,F0.1,A,I0,A,I0,A,F0.1,A,F0.1,A,F0.3,A,F0.3,A,F0.3)') &
            int(flow_conn_catch(i, 1)), tab, &
            flow_conn_catch(i, 2), tab, &
            flow_conn_catch(i, 3), tab, &
            int(flow_conn_catch(i, 4)), tab, &
            int(flow_conn_catch(i, 5)), tab, &
            flow_conn_catch(i, 6), tab, &
            flow_conn_catch(i, 7), tab, &
            flow_conn_catch(i, 8), tab, &
            flow_conn_catch(i, 9), tab, &
            flow_conn_catch(i, 10)
    end do

    close (100)

    ! Write the river data file
    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_riv_data.dat'

    col_headers(1) = 'riv_id'
    col_headers(2) = 'area'
    col_headers(3) = 'dist'
    col_headers(4) = 'secion_dist'
    col_headers(5) = 'slope'
    col_headers(6) = 'elevation'

    call write_numeric_list(catch_fn, col_headers, riv_point_data_catch, 3)

    catch_fn = trim(outputfolder)//'/'//trim(gauge)//'_flow_point.dat'

    open(100, file = catch_fn, status = 'unknown')

    write (100, '(19A)') 'Node_ID', tab, &
        'x_easting', tab, &
        'y_northing', tab, &
        'Type', tab, &
        'dist', tab, &
        'h_elevation'

    do i = 1 , size(flow_point_catch, 1)
        write(100, '(I0,A,F0.1,A,F0.1,A,I0,A,F0.3,A,F0.3)') &
            int(flow_point_catch(i, 1)), tab, &
            flow_point_catch(i, 2), tab, &
            flow_point_catch(i, 3), tab, &
            int(flow_point_catch(i, 4)), tab, &
            flow_point_catch(i, 5), tab, &
            flow_point_catch(i, 6)
    end do

    close (100)

end program preprocess
