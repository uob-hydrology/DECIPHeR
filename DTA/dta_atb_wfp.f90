!% Toby Dunne
!% Mar 2016
!% Modified by Gemma Coxon to calculate slope/topographic index properly when no downslope (i.e. coastal cells or edges of DEMs)
!% Oct 2018
module dta_atb_wfp
    implicit none
    double precision :: routefrac_cardinal
    double precision :: routefrac_ordinal
    double precision :: routedist_cardinal
    double precision :: routedist_ordinal
contains

    subroutine find_peaks(nrows, ncols, dem, peak_list)
        use dta_stack
        use dta_utility
        implicit none
        integer :: ncols
        integer :: nrows
        double precision :: dem(nrows,ncols)
        type(point_type), allocatable, dimension(:) :: peak_list

        type(stack_type) :: peak_q

        integer :: x
        integer :: y
        double precision :: dem_centre

        type(point_type) :: newpeak
        !print *, 'find peaks'

        peak_q = new_stack()

        do y=2,nrows-1
            do x=2,ncols-1
                dem_centre = dem(y, x)
                if ((dem_centre > -9000) .and. &
                    (dem(y-1, x-1) <= dem_centre) .and. &
                    (dem(y-1, x  ) <= dem_centre) .and. &
                    (dem(y-1, x+1) <= dem_centre) .and. &
                    (dem(y,   x-1) <= dem_centre) .and. &
                    (dem(y,   x+1) <= dem_centre) .and. &
                    (dem(y+1, x-1) <= dem_centre) .and. &
                    (dem(y+1, x  ) <= dem_centre) .and. &
                    (dem(y+1, x+1) <= dem_centre)) then

                    newpeak%x = x
                    newpeak%y = y
                    call peak_q%push(newpeak)
                endif
            enddo
        enddo


        !        do y=1,nrows
        !            newpeak%x = 1
        !            newpeak%y = y
        !            if (dem(newpeak%y, newpeak%x) > -9000) then
        !                call peak_q%enqueue(newpeak)
        !            endif
        !            newpeak%x = ncols
        !            newpeak%y = y
        !            if (dem(newpeak%y, newpeak%x) > -9000) then
        !                call peak_q%enqueue(newpeak)
        !            endif
        !        enddo
        !        do x=1,ncols
        !            newpeak%x = x
        !            newpeak%y = 1
        !            if (dem(newpeak%y, newpeak%x) > -9000) then
        !                call peak_q%enqueue(newpeak)
        !            endif
        !            newpeak%x = x
        !            newpeak%y = nrows
        !            if (dem(newpeak%y, newpeak%x) > -9000) then
        !                call peak_q%enqueue(newpeak)
        !            endif
        !        enddo

        allocate(peak_list(peak_q%item_count))
        call peak_q%to_array(peak_list)
        call peak_q%free()

    end subroutine


    subroutine find_mask_bounds(nrows, ncols, catch_mask, cell_list, peak_list)
        ! if a mask is specified, ensure the cells on the mask boundary are added as start points
        ! since flow stops at the boundary - a peak may not actually exist inside a masked section
        ! this will ensure that the search will re-start at boundaries
        ! start points are sorted by height to ensure they are processed in correct order
        use dta_utility
        use dta_stack
        implicit none
        integer :: ncols
        integer :: nrows
        integer :: catch_mask(nrows,ncols)
        type(point_type), allocatable, dimension(:) :: cell_list
        ! existing peaks
        type(point_type), allocatable, dimension(:) :: peak_list

        type(stack_type) :: cell_q

        integer :: x
        integer :: y
        integer :: catch_id_centre

        type(point_type) :: newpeak
        print *, 'find peaks'

        cell_q = new_stack()

        do y=1,size(peak_list)
            call cell_q%push(peak_list(y))
        end do

        do y=2,nrows-1
            do x=2,ncols-1
                catch_id_centre = catch_mask(y, x)
                if (catch_id_centre < -9000) then
                    cycle
                endif
                if ((catch_mask(y-1, x-1) /= catch_id_centre) .or. &
                    (catch_mask(y-1, x  ) /= catch_id_centre) .or. &
                    (catch_mask(y-1, x+1) /= catch_id_centre) .or. &
                    (catch_mask(y,   x-1) /= catch_id_centre) .or. &
                    (catch_mask(y,   x+1) /= catch_id_centre) .or. &
                    (catch_mask(y+1, x-1) /= catch_id_centre) .or. &
                    (catch_mask(y+1, x  ) /= catch_id_centre) .or. &
                    (catch_mask(y+1, x+1) /= catch_id_centre)) then

                    newpeak%x = x
                    newpeak%y = y
                    call cell_q%push(newpeak)
                endif
            enddo
        enddo

        allocate(cell_list(cell_q%item_count))
        call cell_q%to_array(cell_list)
        call cell_q%free()

    end subroutine

    subroutine calc_atb(nrows, ncols, &
        riv_nrows, riv_ncols, &
        catch_mask_nrows, catch_mask_ncols, &
        dem, riv, catch_mask, &
        peak_list, a, c, slope, cellsize)
        use dta_stack
        use dta_utility
        implicit none
        integer :: ncols
        integer :: nrows
        integer :: riv_ncols ! 0 if disabled
        integer :: riv_nrows ! 0 if disabled
        integer :: catch_mask_ncols ! 0 if disabled
        integer :: catch_mask_nrows ! 0 if disabled
        double precision :: dem(nrows,ncols)
        double precision :: riv(riv_nrows,riv_ncols) !optional (size 0 if disabled)
        integer :: catch_mask(catch_mask_nrows,catch_mask_ncols) !optional (size 0 if disabled)
        double precision :: c(nrows,ncols)
        double precision :: a(nrows,ncols)
        double precision :: slope(nrows,ncols)
        double precision :: cellsize
        logical :: mask(nrows,ncols)

        integer :: cellState(nrows, ncols)
        type(point_type), dimension(:) :: peak_list
        type(stack_type) :: openCells
        integer :: ipeak
        type(point_type) :: peak
        type(point_type) :: node
        type(time_type) :: start_time, now_time, last_print_time

        integer :: processed_count
        integer :: total

        routefrac_cardinal = 0.5
        routefrac_ordinal = 0.25

        routedist_cardinal = 1 * cellsize
        routedist_ordinal = sqrt(2.0) * cellsize

        c(:,:) = -9999.0d0
        a(:,:) = (cellsize * cellsize)
        slope(:,:) = -9999.0d0
        cellState = 0
        processed_count = 0

        openCells = new_stack()

        mask = dem > -9000

        total = sum(count(mask, 2))
        !print *, 'atb start', total

        CALL timer_get(start_time)
        last_print_time = start_time

        do ipeak=1,size(peak_list)
            ! read from peak_list in reverse order
            peak = peak_list(size(peak_list) - (ipeak-1))
            if(cellState(peak%Y, peak%X) == 0) then
                cellState(peak%Y, peak%X) = 1

                call openCells%push(peak)

                do
                    if(openCells%item_count == 0) then
                        exit
                    endif
                    node = openCells%pop()

                    call processNode(nrows, ncols, &
                        riv_nrows, riv_ncols, &
                        catch_mask_nrows, catch_mask_ncols, &
                        dem, riv, catch_mask, cellState, c, a, slope, node)

                    call addCandidates(nrows, ncols, catch_mask_nrows, catch_mask_ncols, &
                        dem, catch_mask, cellState, openCells, node)
                    processed_count = processed_count + 1

                    !if(mod(processed_count,10000) == 0) then
                    !    CALL timer_get(now_time)
                    !    call timer_estimate(processed_count, total, start_time, now_time, last_print_time)
                    !endif
                enddo
            endif
        enddo
        call openCells%free()
    end subroutine calc_atb

    function isCandidate(nrows, ncols, &
        catch_mask_nrows, catch_mask_ncols, &
        dem, catch_mask, cellState, nodeY, nodeX) result(is_candidate)
        implicit none
        integer :: ncols
        integer :: nrows
        integer :: catch_mask_ncols ! 0 if disabled
        integer :: catch_mask_nrows ! 0 if disabled
        double precision :: dem(nrows,ncols)
        integer :: catch_mask(catch_mask_nrows,catch_mask_ncols) !optional (size 0 if disabled)
        integer :: cellState(nrows,ncols)
        integer :: nodeY
        integer :: nodeX
        logical :: is_candidate

        integer :: j
        integer :: k
        integer :: x
        integer :: y

        integer :: processed_count

        processed_count = 0

        do j=-1,1
            do k=-1,1
                if ((j==0 .and. k==0)) then
                    cycle
                endif
                !neighbourIndex = neighbourIndex + 1
                y = nodeY + j
                x = nodeX + k
                if (x < 1 .or. y < 1 .or. x > ncols .or. y > nrows ) then
                    processed_count = processed_count + 1
                    cycle
                endif

                ! candidate when all upstream cells are processed

                ! neighbour is processed (cellState == 2)
                ! neighbour dem cell is not upstream
                ! single check takes into account: 3 cases
                ! - neighbour less than
                ! - neighbour equal to
                ! - neighbour NaN (-9999) so less than match
                if (cellState(y,x) == 2 &
                    .or. dem(y, x) <= dem(nodeY, nodeX) ) then
                    processed_count = processed_count + 1
                else
                    ! if neighbour is in a different catchment mask mark neighbour as processed
                    if(catch_mask_nrows > 0)then
                        if(catch_mask(y, x) /= catch_mask(nodeY, nodeX))then
                            processed_count = processed_count + 1
                        endif
                    else
                        exit
                    endif
                endif
            enddo
        enddo

        !print *, 'isCandidate', processed_count

        is_candidate = (processed_count == 8)
    end function

    subroutine addCandidates(nrows, ncols, &
        catch_mask_nrows, catch_mask_ncols, &
        dem, catch_mask, cellState, openCells, node)
        use dta_stack
        use dta_utility
        implicit none
        !argument declares
        integer :: ncols
        integer :: nrows
        integer :: catch_mask_ncols ! 0 if disabled
        integer :: catch_mask_nrows ! 0 if disabled
        double precision :: dem(nrows,ncols)
        integer :: catch_mask(catch_mask_nrows,catch_mask_ncols) !optional (size 0 if disabled)
        integer :: cellState(nrows,ncols)
        type(stack_type) :: openCells
        type(point_type) :: node

        ! local declares
        type(point_type) :: newnode
        integer :: j
        integer :: k
        integer :: x
        integer :: y

        !print *, 'addCandidates'

        do j=-1,1
            do k=-1,1
                if (j==0.and.k==0) then
                    !%                     fprintf(' X ')
                    cycle
                endif
                y = node%Y + j
                x = node%X + k
                if (x < 1 .or. y < 1 .or. x > ncols .or. y > nrows ) then
                    !                     fprintf(' X ')
                    cycle
                endif

                ! don't check already queued or NaN cells
                if (cellState(y,x) /= 0 .or. dem(y,x) < -9000 ) then
                    !                     fprintf(' - ')
                    cycle
                endif
                if(isCandidate(nrows, ncols, catch_mask_nrows, catch_mask_ncols,&
                    dem, catch_mask, cellState, y, x)) then

                    !                     fprintf(' 1 ')
                    cellState(y,x) = 1

                    newnode%y = y
                    newnode%x = x

                    call openCells%push(newnode)
                else
                    !%                     fprintf(' 0 ')
                endif
            enddo
            !%             fprintf(' \n')
        enddo
    end subroutine

    subroutine processNode(nrows, ncols, &
        riv_nrows, riv_ncols, &
        catch_mask_nrows, catch_mask_ncols, &
        dem, riv, catch_mask, cellState, c, a, slope, node)
        use dta_stack
        use dta_utility
        implicit none
        integer :: ncols
        integer :: nrows
        integer :: riv_ncols ! 0 if disabled
        integer :: riv_nrows ! 0 if disabled
        integer :: catch_mask_ncols ! 0 if disabled
        integer :: catch_mask_nrows ! 0 if disabled
        double precision :: dem(nrows,ncols)
        double precision :: riv(riv_nrows,riv_ncols) !optional (size 0 if disabled)
        integer :: catch_mask(catch_mask_nrows,catch_mask_ncols) !optional (size 0 if disabled)
        integer :: cellState(nrows,ncols)
        double precision :: c(nrows,ncols)
        double precision :: a(nrows,ncols)
        double precision :: slope(nrows,ncols)
        type(point_type) :: node

        ! local declares
        double precision :: wfp(3,3)
        double precision :: slopes(3,3)
        integer :: j
        integer :: k
        integer :: x
        integer :: y
        ! sum of the weighted flow proportions
        double precision :: tanb
        double precision :: routefrac
        double precision :: routedist
        ! number of downstream cells
        integer :: nroute
        logical :: stop_flag

        ! mark as processed
        cellState(node%Y, node%X) = 2

        !% tanB is zero on cells that flow out of the dem (no downstream)
        ! tanB is the sum of the wfp for downstream cells

        wfp(:,:) = 0
        slopes(:,:) = 0

        nroute = 0
        ! calculate flow proportions for downstream cells
        do j=-1,1
            do k=-1,1
                if((j==0.and.k==0)) then
                    cycle
                endif
                y = node%Y + j
                x = node%X + k
                if (x < 1 .or. y < 1 .or. x > ncols .or. y > nrows ) then
                    cycle
                endif

                stop_flag = .false.

                ! stop flow across catch_mask boundaries (if catch_mask enabled)
                if (catch_mask_nrows > 0) then
                    if (catch_mask(y,x) /= catch_mask(node%Y, node%X)) then
                        stop_flag = .true. !
                    endif
                endif

                if (stop_flag .eqv. .false.) then
                    ! calculate wfp for all downstream
                    if(dem(y,x) > -9000 .and. dem(y, x) < dem(node%Y, node%X)) then
                        if(j == 0 .or. k == 0) then
                            routefrac = routefrac_cardinal
                            routedist = routedist_cardinal
                        else
                            routefrac = routefrac_ordinal
                            routedist = routedist_ordinal
                        endif

                        wfp(j+2, k+2) = (dem(node%Y, node%X) - dem(y, x)) * routefrac

                        ! actual slopes in each direction
                        slopes(j+2, k+2) = (dem(node%Y, node%X) - dem(y, x)) / routedist
                        nroute = nroute + 1
                    endif
                endif
            enddo
        enddo

        if(nroute>0) then
            slope(node%Y, node%X) = sum(slopes) / nroute

            ! tanB is the sum of outflow weights as per Jim's original code
            tanB = sum(wfp)

            ! Toby - should tanb be slope since that is what ??
            !tanB = slope(node%Y, node%X)

            if(tanB > 0) then

                c(node%Y, node%X) = a(node%Y, node%X) / tanB

                stop_flag = .false.
                ! stop accumulation if cell is a river (if rivers enabled)
                if(riv_nrows > 0) then
                    if(riv(node%Y, node%X) > 0) then
                        stop_flag = .true.
                    endif
                endif
                if(stop_flag.eqv..false.)then
                    !% all downstream cells have their area updated in proportion to wfp
                    !% upstream cells will have wpf=0, so no change
                    do j=-1,1
                        do k=-1,1
                            if((j==0.and.k==0)) then
                                cycle
                            endif
                            !neighbourIndex = neighbourIndex + 1
                            y = node%Y + j
                            x = node%X + k
                            if (x < 1 .or. y < 1 .or. x > ncols .or. y > nrows ) then
                                cycle
                            endif

                            a(y, x) = a(y, x) + c(node%Y, node%X) * wfp(j+2, k+2)
                        enddo
                    enddo
                endif
            endif

        else if ((nroute.eq.0) .and. (dem(node%Y, node%X) > - 9000)) then

        ! Fix for coastal cells with no downslope - use average of inflow slopes

            ! Check to see if has a no data cell around it (i.e. the sea!)
            stop_flag = .true.
            do j=-1,1
            do k=-1,1
                if((j==0.and.k==0)) then
                    cycle
                endif
                y = node%Y + j
                x = node%X + k
                if (x < 1 .or. y < 1 .or. x > ncols .or. y > nrows ) then
                    cycle
                endif

            if(dem(y,x) < -9000) then
                stop_flag = .false.
            end if

            end do
            end do


            if (stop_flag .eqv. .false.) then
                 ! slope and wfp is calculated from upslope cells

               do j=-1,1
               do k=-1,1

                if((j==0.and.k==0)) then
                    cycle
                endif
                y = node%Y + j
                x = node%X + k
                if (x < 1 .or. y < 1 .or. x > ncols .or. y > nrows ) then
                    cycle
                endif

                if(dem(y,x) > -9000) then
                     if(j == 0 .or. k == 0) then
                            routefrac = routefrac_cardinal
                            routedist = routedist_cardinal
                        else
                            routefrac = routefrac_ordinal
                            routedist = routedist_ordinal
                        endif

                      wfp(j+2, k+2) = (dem(y, x) - dem(node%Y, node%X)) * routefrac
                      slopes(j+2, k+2) = (dem(y, x) - dem(node%Y, node%X)) / routedist
                      nroute = nroute + 1
                 end if
                end do
                end do

            slope(node%Y, node%X) = sum(slopes) / nroute

            ! tanB is the sum of outflow weights as per Jim's original code
            tanB = sum(wfp)

            if(tanB > 0) then

                c(node%Y, node%X) = a(node%Y, node%X) / tanB

            endif
            end if
        endif

    end subroutine

end module dta_atb_wfp
