module dta_preprocess

    implicit none

    type hru_class_struct
        character(len=1024) :: class_type
        character(len=1024) :: filename
        character(len=1024) :: fn_end
        double precision, allocatable, dimension (:,:) :: data
        integer, allocatable, dimension(:,:) :: class_array
        integer :: num_frac
        double precision, allocatable, dimension(:) :: frac
        double precision, allocatable, dimension(:) :: frac_val
        integer, allocatable, dimension(:, :) :: meta
    end type hru_class_struct

contains

    subroutine read_hruclass_file(hru_class, hru_file, fn_log)

        implicit none

        ! Input related variables
        character(len=1024) :: buffer, label, hru_file, comment_check, name
        integer :: pos, num_frac
        integer, parameter :: fh = 15, fh_topo = 16
        integer :: fn_log
        integer :: ios = 0, ios_topo = 0
        integer :: line = 0, line_topo = 0
        type(hru_class_struct), intent(in out):: hru_class(8)

        open(fh, file=hru_file)

        hru_class(1)%class_type = 'TOPO_AREA_CLASS'
        hru_class(2)%class_type = 'TOPO_SLOPE_CLASS'
        hru_class(3)%class_type = 'TOPO_ELEV_CLASS'
        hru_class(4)%class_type = 'INPUT_RAIN_CLASS'
        hru_class(5)%class_type = 'INPUT_PET_CLASS'
        hru_class(6)%class_type = 'MODELSTRUCT_CLASS'
        hru_class(7)%class_type = 'PARAM_CLASS'

        hru_class%filename = ''
        hru_class%fn_end = ''

        ios = 0
        ios_topo = 0
        ! ios is negative if an end of record condition is encountered or if
        ! an endfile condition was detected.  It is positive if an error was
        ! detected.  ios is zero otherwise.

        do while (ios == 0)
            read(fh, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1
                comment_check = buffer(1:1)
                select case (comment_check)
                    case('!')
                    ! Comment line so can skip
                    case default
                        !Split label and data.
                        pos = scan(buffer, "=")
                        label = trim(buffer(1:pos-1))
                        buffer = buffer(pos+1:)
                        select case (label)
                            case ('TOPO_CLASS')

                                ! Split the catchnment up by topographic classes
                                ! Read in different classes from classinc.dat file
                                open(fh_topo, file = trim(buffer))

                                do while (ios_topo == 0)
                                    read(fh_topo, '(A)', iostat=ios_topo) buffer
                                    if (ios_topo == 0) then
                                        line_topo = line_topo + 1
                                        comment_check = buffer(1:1)
                                        select case (comment_check)
                                            case('!')
                                            ! Comment line so can skip
                                            case default
                                                !Get label and classes
                                                read(buffer, *) name, num_frac
                                                select case (name)
                                                    case ('AREA')
                                                        hru_class(1)%fn_end = '_area'
                                                        hru_class(1)%num_frac = num_frac
                                                        allocate(hru_class(1)%frac(num_frac))
                                                        allocate(hru_class(1)%frac_val(num_frac))
                                                        read(fh_topo, *) hru_class(1)%frac(:)
                                                        write(fn_log, *) 'Read TOPO_AREA_CLASS: ', num_frac, ' classes'
                                                    case ('SLOPE')
                                                        hru_class(2)%fn_end = '_slope'
                                                        hru_class(2)%num_frac = num_frac
                                                        allocate(hru_class(2)%frac(num_frac))
                                                        allocate(hru_class(2)%frac_val(num_frac))
                                                        read(fh_topo, *) hru_class(2)%frac(:)
                                                        write(fn_log, *) 'Read TOPO_SLOPE_CLASS: ', num_frac, ' classes'
                                                    case ('ELEV')
                                                        hru_class(3)%fn_end = '_dem'
                                                        hru_class(3)%num_frac = num_frac
                                                        allocate(hru_class(3)%frac(num_frac))
                                                        allocate(hru_class(3)%frac_val(num_frac))
                                                        read(fh_topo, *) hru_class(3)%frac(:)
                                                        write(fn_log, *) 'Read TOPO_ELEV_CLASS: ', num_frac, ' classes'
                                                    case default
                                                        print *, 'Skipping invalid TOPO_CLASS ', trim(name), ' on line ', line
                                                        write(fn_log, *) 'Skipping invalid TOPO_CLASS ', trim(name), &
                                                            ' on line ', line_topo
                                                end select
                                        end select
                                    end if
                                end do
                                close(fh_topo)

                            case ('INPUT_RAIN_CLASS')
                                hru_class(4)%filename = trim(buffer)
                                hru_class(4)%fn_end = '_rain'
                                write(fn_log, *) 'Read INPUT_RAIN_CLASS: ', trim(hru_class(4)%filename)
                            case ('INPUT_PET_CLASS')
                                hru_class(5)%filename = trim(buffer)
                                hru_class(5)%fn_end = '_pet'
                                write(fn_log, *) 'Read INPUT_PET_CLASS: ', trim(hru_class(5)%filename)
                            case ('MODELSTRUCT_CLASS')
                                hru_class(6)%filename = trim(buffer)
                                hru_class(6)%fn_end = '_modelstruct'
                                write(fn_log, *) 'Read : MODELSTRUCT_CLASS', trim(hru_class(6)%filename)
                            case ('PARAM_CLASS')
                                hru_class(7)%filename = trim(buffer)
                                hru_class(7)%fn_end = '_param'
                                write(fn_log, *) 'Read PARAM_CLASS: ', trim(hru_class(7)%filename)
                            case default
                                print *, 'Skipping invalid label ', trim(label), ' on line ', line
                        end select
                end select
            end if
        end do

        write(fn_log, *) 'Successfully read in HRU Classifier File'
        write(fn_log, *) ''
        close(fh)
        close(fh_topo)

    end subroutine read_hruclass_file

    !  This subroutine cuts out data for a catchment from a uk dataset.

    subroutine subset_catch(dta_data, catch_mask, ncols, nrows, xllcorner, yllcorner, cellsize, nodata, &
        catch_data, catch_ncols, catch_nrows, catch_xllcorner, catch_yllcorner, catch_cellsize, catch_nodata)

        use dta_utility

        implicit none

        double precision, intent(in) :: dta_data(:,:), catch_mask(:,:)
        double precision, intent(in) :: xllcorner, yllcorner, nodata, cellsize
        integer, intent(in) :: ncols, nrows

        ! Locals
        double precision :: catch_xllcorner, catch_yllcorner, catch_nodata
        double precision :: catch_cellsize, ix(2), iy(2)
        integer :: catch_ncols, catch_nrows
        double precision, allocatable :: catch_data(:,:)

        allocate (catch_data(size(catch_mask, 1), size(catch_mask, 2)))

        iy(1) = nrows - ((catch_yllcorner - yllcorner)/cellsize)
        iy(2) = iy(1) - catch_nrows+1
        ix(1) = ((catch_xllcorner - xllcorner)/cellsize) +1
        ix(2) = ix(1) + catch_ncols -1

        catch_data = dta_data(int(iy(2)):int(iy(1)), int(ix(1)):int(ix(2)))

        where (catch_mask.eq.catch_nodata)
            catch_data = -9999
        end where

    end subroutine subset_catch

    ! This subroutine subsets the flow routing files for the principal basin

    subroutine subset_route_files(mask_data, flow_conn, flow_point, riv_point_data, &
        flow_conn_catch, flow_point_catch, riv_point_data_catch)

        double precision, intent(in) :: flow_conn(:,:), riv_point_data(:,:), mask_data(:,:), flow_point(:,:)
        double precision, allocatable, intent(out) :: flow_conn_catch(:,:), riv_point_data_catch(:,:), flow_point_catch(:,:)

        !Locals
        integer :: i, ix1, j, riv_point_count, ncatch
        integer, allocatable :: gauges(:)
        integer, dimension(size(mask_data, 1), size(mask_data, 2)) :: mask

        ncatch = 0
        mask = 0

        ! Get all the gauge numbers you are interested in from the mask

        do i = 1, size(mask_data, 1)
            do j = 1, size(mask_data, 2)
                if ((mask_data(i,j).gt.0).and.(mask(i,j).lt.1)) then
                    ncatch = ncatch + 1
                    where(mask_data.eq.mask_data(i,j))
                        mask=1
                    end where
                end if
            end do
        end do

        allocate(gauges(ncatch))

        mask = 0
        ncatch = 0

        do i = 1, size(mask_data, 1)
            do j = 1, size(mask_data, 2)
                if ((mask_data(i,j).gt.0).and.(mask(i,j).lt.1)) then
                    ncatch = ncatch + 1
                    gauges(ncatch) = int(mask_data(i,j))
                    where(mask_data.eq.mask_data(i,j))
                        mask=1
                    end where
                end if
            end do
        end do

        ! Subset the river point data file to just the gauges of interest
        riv_point_count = 0

        do j=1,size(gauges,1)
            do i = 1, size(riv_point_data, 1)
                if (riv_point_data(i, 1).eq.gauges(j)) then
                    riv_point_count = riv_point_count+1
                end if
            end do
        end do
        allocate(riv_point_data_catch(riv_point_count, 6))

        ix1 = 0

        do j=1,size(gauges,1)
            do i = 1, size(riv_point_data, 1)
                if (riv_point_data(i, 1).eq.gauges(j)) then
                    ix1 = ix1+1
                    riv_point_data_catch(ix1, 1:6) =riv_point_data(i,1:6)
                end if
            end do
        end do

        ! Subset the flow connectivity file for just the gauges of interest
        allocate(flow_conn_catch(size(gauges, 1), 10))

        do i = 1, size(gauges, 1)
            ix1 = minloc(abs(flow_conn(:, 1) - gauges(i)), 1)
            flow_conn_catch(i, 1:10) = flow_conn(ix1, 1:10)

        end do

        allocate(flow_point_catch(size(gauges, 1), 6))

        do i = 1, size(gauges, 1)

            ix1 = minloc(abs(flow_point(:, 1) - gauges(i)), 1)
            flow_point_catch(i, 1:6) = flow_point(ix1, 1:6)

        end do

    end subroutine subset_route_files

    subroutine write_routingfiles(flowconn_fn, flow_conn_catch, rivdata_fn, riv_point_data_catch, &
        flowpoint_fn, flow_point_catch)

        use dta_utility

        double precision, intent(in), dimension(:,:) :: flow_conn_catch, riv_point_data_catch, flow_point_catch
        character(len=1024), intent(in) :: flowconn_fn, rivdata_fn, flowpoint_fn
        integer :: i

        ! Locals
        character(64) :: col_headers(6)

        open(100, file = flowconn_fn, status = 'unknown')

        write(100, '(19A)') &
            'Node_ID', tab, &
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

        col_headers(1) = 'riv_id'
        col_headers(2) = 'area'
        col_headers(3) = 'dist'
        col_headers(4) = 'section_dist'
        col_headers(5) = 'elevation'
        col_headers(6) = 'slope'

        call write_numeric_list(rivdata_fn, col_headers, riv_point_data_catch, 3)

        open(100, file = flowpoint_fn, status = 'unknown')

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

    end subroutine write_routingfiles

end module dta_preprocess
