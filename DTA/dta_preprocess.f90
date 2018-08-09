!% Toby Dunne
!% Apr 2016
module dta_preprocess

    implicit none

contains

    !  This subroutine cuts out data for a catchment from a uk dataset.  Catchments were created by Toby!

    subroutine subset_catch(uk_fn, catch_mask, catch_ncols, catch_nrows, &
        catch_xllcorner, catch_yllcorner, catch_nodata, catch_data)

        use dta_utility

        implicit none

        double precision, intent(in) :: catch_mask(:,:)
        double precision, intent(in) :: catch_xllcorner, catch_yllcorner, catch_nodata
        integer, intent(in) :: catch_ncols, catch_nrows
        character(len=1024), intent(in) :: uk_fn
        double precision, allocatable, intent(out) :: catch_data(:,:)

        ! Locals
        double precision, allocatable, dimension(:,:) :: uk_data
        double precision :: xllcorner, yllcorner, nodata, cellsize
        double precision :: ix(2), iy(2)
        integer :: ncols, nrows

        allocate (catch_data(size(catch_mask, 1), size(catch_mask, 2)))

        !  Read in the UK dataset

        call read_ascii_grid(uk_fn, uk_data, ncols, nrows, xllcorner, yllcorner, cellsize, nodata)

        iy(1) = nrows - ((catch_yllcorner - yllcorner)/cellsize)
        iy(2) = iy(1) - catch_nrows+1
        ix(1) = ((catch_xllcorner - xllcorner)/cellsize) +1
        ix(2) = ix(1) + catch_ncols -1

        catch_data = uk_data(int(iy(2)):int(iy(1)), int(ix(1)):int(ix(2)))

        where (catch_mask.eq.catch_nodata)
            catch_data = -9999
        end where

    end subroutine subset_catch

    subroutine create_funcunits(rain_grid, catch_data, river_data, fu_table, fu_grid, rain_grid_dta, catch_data_dta, &
        river_data_dta, gauges)

        use dta_utility

        implicit none

        double precision, intent(in) :: catch_data(:,:), rain_grid(:,:), river_data(:,:)
        integer, allocatable, intent(in out) :: fu_grid(:,:), fu_table(:,:)
        double precision, allocatable, intent(out) :: gauges(:,:)

        ! Locals
        integer :: x, y, ix(2), krain, kcatch, num_fu, fu_count
        integer, allocatable, dimension(:, :) :: class_all
        integer, allocatable, dimension(:, :), intent(out) :: rain_grid_dta, catch_data_dta, river_data_dta

        ! Find all the unique elements in each array

        allocate (rain_grid_dta(size(rain_grid, 1), size(rain_grid, 2)))
        rain_grid_dta = 0
        krain=1

        do y = 1, size(rain_grid, 1)
            do x = 1, size(rain_grid, 2)
                if (rain_grid(y,x).eq.-9999) then
                    cycle
                else if (rain_grid_dta(y,x).ne.0) then
                    cycle
                else
                    where (rain_grid.eq.rain_grid(y,x))
                        rain_grid_dta = krain
                    end where
                    krain=krain+1
                end if
            end do
        end do

        allocate (catch_data_dta(size(catch_data, 1), size(catch_data, 2)))
        catch_data_dta = 0

        ! Check whether kcatch should start off as 1 rather than 0
        kcatch = 1

        allocate (river_data_dta(size(catch_data, 1), size(river_data, 2)))
        river_data_dta = 0

        do y = 1, size(catch_data, 1)
            do x = 1, size(catch_data, 2)
                if (catch_data(y,x).eq.-9999) then
                    cycle
                else if (catch_data_dta(y,x).ne.0) then
                    cycle
                else
                    where (catch_data.eq.catch_data(y,x))
                        catch_data_dta = kcatch
                    end where
                    where ((catch_data.eq.catch_data(y,x)).and.(river_data.eq.catch_data(y,x)))
                        river_data_dta = kcatch
                    end where
                    kcatch=kcatch+1
                end if
            end do
        end do

        ! Put in a warning here
        if (count(river_data > 0).ne.count(river_data_dta>0)) then

            print *, 'ERROR : River IDs need to match catchment IDs'
            stop

        end if

        allocate(gauges(kcatch-1, 7))
        gauges = 0
        kcatch =1

        ! Populate a list of gauges
        do y = 1, size(catch_data, 1)
            do x = 1, size(catch_data, 2)
                if (catch_data(y,x).eq.-9999) then
                    cycle
                else if (ANY(gauges(:, 2).eq.catch_data(y,x))) then
                    cycle
                else
                    gauges(kcatch, 1) = kcatch
                    gauges(kcatch, 2) = catch_data(y,x)
                    kcatch = kcatch+1
                end if
            end do
        end do

        gauges(:, 2) = gauges(:, 2) / 1000

        allocate (class_all(kcatch-1, krain-1))
        class_all=0

        do y = 1, size(catch_data, 1)
            do x = 1, size(catch_data, 2)

                ix(1) = int(catch_data_dta(y,x))
                ix(2) = int(rain_grid_dta(y,x))

                ! Cells outside the catchment
                if ((ix(1).eq.0).and.(ix(2).eq.0)) then
                    cycle
                else
                    class_all(ix(1), ix(2)) = class_all(ix(1), ix(2)) + 1
                end if

            end do
        end do

        num_fu = count(class_all.gt.0)

        allocate (fu_table(num_fu, 7))
        allocate (fu_grid(size(catch_data, 1), size(catch_data, 2)))

        fu_count = 0
        fu_grid = 0
        fu_table = 0

        do y = 1, size(class_all, 1)
            do x = 1, size(class_all, 2)
                if (class_all(y, x).gt.0) then
                    fu_count = fu_count + 1;
                    fu_table(fu_count, 1) = fu_count
                    where ((catch_data_dta.eq.y).and.(rain_grid_dta.eq.x))
                        fu_grid = fu_count
                    end where
                end if
            end do
        end do

        do y = 1, size(fu_grid, 1)
            do x = 1, size(fu_grid, 2)
                if (fu_grid(y,x).ne.0) then
                    fu_table(int(fu_grid(y,x)), 2) = int(catch_data_dta(y,x))
                    fu_table(int(fu_grid(y,x)), 3) = int(catch_data(y,x))
                    fu_table(int(fu_grid(y,x)), 4) = int(rain_grid_dta(y,x))
                    fu_table(int(fu_grid(y,x)), 5) = int(rain_grid(y,x))
                    if (river_data(y, x).gt.0) then
                        fu_table(int(fu_grid(y,x)), 6) = int(river_data_dta(y,x))
                        fu_table(int(fu_grid(y,x)), 7) = int(river_data(y,x))
                    end if
                end if
            end do
        end do

    end subroutine create_funcunits

    ! Subroutine to create the initialisation file

    subroutine create_init_file(gauges, gauge, flow_conn, catch_area, outputfolder)

        double precision, intent(in out) :: gauges(:,:)
        character(len=1024), intent(in) :: gauge, outputfolder
        double precision, intent(in) :: catch_area(:,:), flow_conn(:,:)

        ! Locals
        character(len=1024) :: temp_fn
        integer :: i, ds_count, k, ix, j
        double precision :: downstream, upstream, outlet

        ! For each gauge get its downstream gauge and catchment area

        do i = 1, size(gauges, 1)

            ix = minloc(abs(flow_conn(:,1) - (gauges(i,2)*1000)), 1)
            gauges(i, 3) = flow_conn(ix, 5) / 1000
            ix = minloc(abs(catch_area(:,1) - gauges(i,2)), 1)
            ! Divide catchment area by 100000 to get to km2
            gauges(i, 4) = catch_area(ix, 5)*1e-6

        end do

        read(gauge, *) outlet
        ix = minloc(abs(gauges(:, 2) - outlet), 1)
        outlet = gauges(ix,3)

        temp_fn = trim(outputfolder)//'/'//trim(gauge)//'_init.dat'
        open(99, file = temp_fn, status = 'unknown')

        ! Write out the total number of gauges

        write(99, '(i0)') size(gauges, 1)

        !  Keep on looping through all the gauges until they are all processed

        do while (sum(gauges(:, 6)) < size(gauges, 1))

            ! Flag all the headwater cells
            ! If it is a headwater then it won't have any downstream matches
            ! Only check for cells that haven't been processed (iterate around all the gauges)

            outer : do i = 1, size(gauges, 1)

                ! Skip the gauge if it has already been processed,
                if (gauges(i, 6).eq.0) then

                    do j = 1, size(gauges, 1)

                        ! Found a downstream gauge
                        if (gauges(i,2).eq.gauges(j,3)) then
                            ! Check if it has already been processed - if it hasn't skip
                            if (gauges(j,6).eq.0) then
                                cycle outer
                            end if
                        end if

                    end do

                    ! It is a headwater so mark as 1
                    gauges(i,5) = 1

                end if

            end do outer

            ! Loop through all the flagged headwater streams

            do i = 1, size(gauges, 1)

                if (gauges(i, 5).eq.1) then

                    ! Fifth column tells you whether this gauge is a 'current' headwater
                    ! Sixth column tells you if a gauge is processed (0 - not processed, 1- processed)
                    ! Seventh column tells you the downstreams in order

                    gauges(i, 6) = 1
                    gauges(i, 7) = 1
                    ds_count = 1

                    ! It is a headwater so find all the downstreams

                    do while (downstream.ne.outlet)

                        ix = minloc(abs(gauges(:, 7) - ds_count), 1)
                        upstream = gauges(ix,2)
                        downstream = gauges(ix,3)
                        ! Label the downstreams in the 6th column in order
                        ix = minloc(abs(gauges(:, 2) - downstream), 1)
                        gauges(ix, 7) = ds_count+1
                        ds_count = ds_count+1

                    end do

                    ! Write out to file
                    ! First line

                    ix = minloc(abs(gauges(:,2) - gauges(i,2)), 1)
                    write(99, '(i0, 1x, i0, 1x, f0.2)') int(gauges(ix, 1)), ds_count-2, gauges(ix, 4)

                    do k = 2, ds_count-1

                        ix = minloc(abs(gauges(:, 7) - k), 1)

                        if (k.eq.ds_count-1) then
                            write(99, '(i0, 1x)') int(gauges(ix, 1))
                        else
                            write(99, '(i0, 1x)', advance = 'no') int(gauges(ix, 1))
                        end if

                    end do

                    do k = 2, ds_count-1

                        ix = minloc(abs(gauges(:, 7) - k), 1)
                        if (k.eq.ds_count-1) then
                            write(99, '(f0.1, 1x)') gauges(ix, 4)
                        else
                            write(99, '(f0.1, 1x)', advance = 'no') gauges(ix, 4)
                        end if

                    end do

                    ! Reset the downstream ordering to zero
                    gauges(:, 7) = 0
                    downstream = 0
                    upstream = 0

                end if

            end do

            ! Reset the current headwater cells
            gauges(:, 5) = 0

        end do

    end subroutine create_init_file

    ! This subroutine creates the files needed for the flow routing
    ! Flow connectivity file and river_data file

    subroutine create_flowroute_files(gauges, river_data, flow_conn, flow_point, riv_point_data, &
    flow_conn_catch, flow_point_catch, riv_point_data_catch)

        double precision, intent(in out) :: gauges(:,:)
        double precision, intent(in) :: flow_conn(:,:), riv_point_data(:,:), river_data(:,:), flow_point(:,:)
        double precision, allocatable, intent(out) :: flow_conn_catch(:,:), riv_point_data_catch(:,:), flow_point_catch(:,:)

        !Locals
        integer :: i, ix1, j, riv_point_count

        ! Subset the flow connectivity file for just the gauges of interest
        allocate(flow_conn_catch(size(gauges, 1), 10))

        do i = 1, size(gauges, 1)

            ix1 = minloc(abs(flow_conn(:, 1) - (gauges(i, 2)*1000)), 1)
            flow_conn_catch(i, 1:10) = flow_conn(ix1, 1:10)

        end do

        ! Subset the river point data file to just the gauges of interest
        riv_point_count = count(river_data > 0)
        allocate(riv_point_data_catch(riv_point_count, 6))

        ix1 = 1

        do i = 1, size(gauges, 1)

            do j = 1, size(riv_point_data, 1)

                if (riv_point_data(j, 1).eq.(gauges(i, 2)*1000)) then

                    riv_point_data_catch(ix1, 1:6) = riv_point_data(j, 1:6)
                    ix1 = ix1+1

                end if

            end do

        end do

        allocate(flow_point_catch(size(gauges, 1), 6))

        do i = 1, size(gauges, 1)

            ix1 = minloc(abs(flow_point(:, 1) - (gauges(i, 2)*1000)), 1)
            flow_point_catch(i, 1:6) = flow_point(ix1, 1:6)

        end do


    end subroutine create_flowroute_files

end module dta_preprocess
