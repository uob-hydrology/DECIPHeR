!% Toby Dunne
!% Mar 2016
module dta_catch_cut
    implicit none

contains

    subroutine find_candidate_river_cells(nrows, ncols, riv_mask_grid, candidates, &
        candidate_count, row, col, &
        search_radius, &
        search_radius_done)
        implicit none
        integer :: nrows
        integer :: ncols
        logical :: riv_mask_grid(nrows, ncols)
        ! rows:max_candidates, cols:3 (1=1, 2=x, 3=y, 4=processed)
        integer :: candidates(:,:)
        ! end index of candidates - find_candidate_river_cells can be called multiple times
        ! each time adding more
        integer :: candidate_count
        ! search start point
        integer :: row
        integer :: col
        integer, intent(in) :: search_radius
        integer, intent(in) :: search_radius_done

        integer :: bound_top
        integer :: bound_bot
        integer :: bound_lef
        integer :: bound_rig
        integer :: x_s
        integer :: y_s
        integer :: i_c
        integer :: index_c
        double precision :: dist
        logical is_new
        ! rows:max_candidates
        double precision :: candidate_dist(size(candidates,1))

        integer max_candidate
        ! how many candidates found in this search
        integer found_candidates

        max_candidate = size(candidates,1)

        found_candidates = 0

        bound_top = max(1, row - search_radius)
        bound_bot = min(nrows, row + search_radius)
        bound_lef = max(1, col - search_radius)
        bound_rig = min(ncols, col + search_radius)

        !%search the local region for the closest river cell
        do y_s=bound_top,bound_bot
            do x_s=bound_lef,bound_rig
                if(riv_mask_grid(y_s,x_s)) then
                    dist = sqrt(real((row - y_s)*(row - y_s) + (col-x_s)*(col-x_s)))

                    ! search in circle < search_radius
                    if (ceiling(dist) >= search_radius_done &
                        .and. ceiling(dist) <= search_radius) then
                        index_c = -1
                        do i_c=1,candidate_count
                            ! match on cell x,y - every cell on river checked
                            if(candidates(i_c,2) == x_s .and. &
                                candidates(i_c,3) == y_s ) then
                                index_c = i_c
                                exit
                            endif
                        enddo
                        is_new = .false.

                        if(index_c == -1 .and. candidate_count < max_candidate) then
                            index_c = candidate_count + 1
                            candidate_count = candidate_count + 1
                            found_candidates = found_candidates + 1
                            is_new = .true.
                            candidates(index_c,1) = 1
                        endif
                        if(index_c == -1) then
                            !print *, 'max candidates excceded', search_radius,'cells'
                            exit
                        else
                            if(is_new .or. dist < candidates(index_c,2)) then
                                candidate_dist(index_c) = dist
                                candidates(index_c,2) = x_s
                                candidates(index_c,3) = y_s
                            endif
                        endif
                    endif
                endif
            enddo
            if(candidate_count >= max_candidate) then
                exit
            endif
        enddo



    end subroutine


    subroutine calc_catch_cut_mfd( nrows, ncols, dem, starty, startx, mask, &
        miny, maxy, minx, maxx, area_count)
        !%CATCH_CUT Find all upstream cells from the given point, labells with flow
        !%length
        !%   mask - 0 for non catchment cell
        !%        > 0 = distance to outlet
        !%   area_count = number of cells in the catchment (multiply by cellsize^2
        !%   to get area)
        use dta_stack
        use dta_utility
        implicit none
        integer, intent(in) :: nrows
        integer, intent(in) :: ncols
        double precision, intent(in) :: dem (nrows, ncols)
        integer, intent(in) :: starty
        integer, intent(in) :: startx

        double precision, intent(inout) :: mask (nrows, ncols)
        ! bounds
        integer, intent(out) :: miny, maxy, minx, maxx
        integer, intent(out) :: area_count

        type(stack_type) openCells
        type(point_type) node

        double precision :: dist_cardinal
        double precision :: dist_ordinal

        integer :: xl, yl, x, y, x_n, y_n
        double precision :: elevation
        double precision :: counter
        double precision :: upCounter

        area_count = 0
        ! clear the entire mask
        mask(:,:) = 0.0d0

        openCells = new_stack()

        node%x = startx
        node%y = starty

        call openCells%push(node)
        mask(starty, startx) = 1

        dist_cardinal = 1
        dist_ordinal = sqrt(2.)

        miny = starty
        maxy = starty
        minx = startx
        maxx = startx

        do while(openCells%item_count > 0)
            node = openCells%pop()
            x = node%x
            y = node%y

            ! expand bounds
            if(y < miny) then
                miny = y
            endif
            if(y > maxy) then
                maxy = y
            endif

            if(x < minx) then
                minx = x
            endif
            if(x > maxx) then
                maxx = x
            endif

            elevation = dem(y, x)
            counter = mask(y, x)

            !print *, 'dequeue', node%x, node%y, counter

            do yl = -1,1
                y_n = y + yl
                if(y_n < 1 .or. y_n > nrows) then
                    cycle
                endif
                do xl = -1,1
                    x_n = x + xl
                    if (xl == 0 .and. yl == 0) then
                        cycle
                    endif
                    if(x_n<1 .or. x_n > ncols) then
                        cycle
                    endif

                    if (xl == 0 .or. yl == 0) then
                        upCounter = counter + dist_cardinal
                    else
                        upCounter = counter + dist_ordinal
                    endif

                    !% label and push to stack any cell that is higher than this cell
                    !% label only unlebelled cells (mask..=0)
                    !% re-label if the distance is less
                    if(dem(y_n,x_n) >= elevation .and. &
                        (mask(y_n,x_n) < 0.5 .or. &
                        mask(y_n,x_n) > upCounter)) then
                        node%x = x_n
                        node%y = y_n

                        !print *, 'enqueue', x, y, '=>', upCounter, '--', xl, yl

                        call openCells%push(node)
                        if(mask(y_n,x_n) < 0.5) then
                            area_count = area_count + 1
                        endif
                        mask(y_n,x_n) = upCounter
                    endif
                enddo
            enddo
        enddo

        call openCells%free()

    end subroutine



    subroutine calc_catch_cut_flow_dir( nrows, ncols, flow_dir_grid, starty, startx, &
        mask, miny, maxy, minx, maxx, area_count)
        !%CATCH_CUT Find all upstream cells from the given point, labells with flow
        !%length
        !%   mask - 0 for non catchment cell
        !%        > 0 = distance to outlet
        !%   area_count = number of cells in the catchment (multiply by cellsize^2
        !%   to get area)
        use dta_queue
        use dta_utility
        implicit none
        integer, intent(in) :: nrows
        integer, intent(in) :: ncols
        integer(1), intent(in) :: flow_dir_grid (nrows, ncols)
        integer, intent(in) :: starty
        integer, intent(in) :: startx

        double precision, intent(inout) :: mask (nrows, ncols)
        ! bounds
        integer, intent(out) :: miny, maxy, minx, maxx
        integer, intent(out) :: area_count

        type(queue_type) openCells
        type(point_type) node

        double precision :: dist_cardinal
        double precision :: dist_ordinal

        integer :: xl, yl, x, y, x_n, y_n
        double precision :: counter
        double precision :: upCounter
        integer neighbour_index
        integer :: loop_count_check
        integer :: max_loop_count_check

        max_loop_count_check = nrows * ncols

        area_count = 0
        loop_count_check = 0
        ! clear the entire mask
        mask(:,:) = 0.0d0

        openCells = new_queue()

        node%x = startx
        node%y = starty

        !call openCells%push(node)
        call openCells%enqueue(node)
        mask(starty, startx) = 1

        dist_cardinal = 1
        dist_ordinal = sqrt(2.)

        miny = starty
        maxy = starty
        minx = startx
        maxx = startx

        do while(openCells%item_count > 0)
            loop_count_check = loop_count_check + 1
            if(loop_count_check > max_loop_count_check) then
                print *, 'error catch cut stuck'
                exit
            endif

            !node = openCells%pop()
            node = openCells%dequeue()
            x = node%x
            y = node%y

            ! expand bounds
            if(y < miny) then
                miny = y
            endif
            if(y > maxy) then
                maxy = y
            endif

            if(x < minx) then
                minx = x
            endif
            if(x > maxx) then
                maxx = x
            endif

            counter = mask(y, x)

            !print *, 'dequeue', node%x, node%y, counter

            neighbour_index = 0
            do yl = -1,1
                y_n = y + yl
                do xl = -1,1
                    x_n = x + xl
                    if ( xl == 0 .and. yl == 0) then
                        cycle
                    endif
                    neighbour_index = neighbour_index + 1
                    if(y_n < 1 .or. y_n > nrows &
                        .or. x_n<1 .or. x_n > ncols) then
                        cycle
                    endif

                    if (xl == 0 .or. yl == 0) then
                        upCounter = counter + dist_cardinal
                    else
                        upCounter = counter + dist_ordinal
                    endif

                    if (flow_dir_grid(y_n,x_n) == dir_in(neighbour_index) &
                        .and. ((mask(y_n,x_n) < 0.5) .or. (mask(y_n,x_n) > upCounter))) then
                        node%x = x_n
                        node%y = y_n

                        !print *, 'enqueue', x, y, '=>', upCounter, '--', xl, yl

                        !call openCells%push(node)
                        call openCells%enqueue(node)
                        if(mask(y_n,x_n) < 0.5) then
                            area_count = area_count + 1
                        endif
                        mask(y_n,x_n) = upCounter
                    endif
                enddo
            enddo
        enddo

        call openCells%free()

    end subroutine

    subroutine calc_catch_cut_flow_dir_mask( nrows, ncols, flow_dir_grid, starty, startx, &
        mask, miny, maxy, minx, maxx, area_count)
        !%CATCH_CUT Find all upstream cells from the given point
        !% logical mask returned
        !%   area_count = number of cells in the catchment (multiply by cellsize^2
        !%   to get area)
        use dta_queue
        use dta_utility
        implicit none
        integer, intent(in) :: nrows
        integer, intent(in) :: ncols
        integer(1), intent(in) :: flow_dir_grid (nrows, ncols)
        integer, intent(in) :: starty
        integer, intent(in) :: startx

        logical, intent(inout) :: mask (nrows, ncols)
        ! bounds
        integer, intent(out) :: miny, maxy, minx, maxx
        integer, intent(out) :: area_count

        type(queue_type) openCells
        type(point_type) node

        integer :: xl, yl, x, y, x_n, y_n
        integer neighbour_index
        integer :: loop_count_check
        integer :: max_loop_count_check

        max_loop_count_check = nrows * ncols

        area_count = 0
        loop_count_check = 0

        openCells = new_queue()

        node%x = startx
        node%y = starty

        !call openCells%push(node)
        call openCells%enqueue(node)
        mask(starty, startx) = .true.

        miny = starty
        maxy = starty
        minx = startx
        maxx = startx

        do while(openCells%item_count > 0)
            loop_count_check = loop_count_check + 1
            if(loop_count_check > max_loop_count_check) then
                print *, 'error catch cut stuck', loop_count_check , max_loop_count_check
                exit
            endif

            !node = openCells%pop()
            node = openCells%dequeue()
            x = node%x
            y = node%y

            ! expand bounds
            if(y < miny) then
                miny = y
            endif
            if(y > maxy) then
                maxy = y
            endif

            if(x < minx) then
                minx = x
            endif
            if(x > maxx) then
                maxx = x
            endif

            neighbour_index = 0
            do yl = -1,1
                y_n = y + yl
                do xl = -1,1
                    x_n = x + xl
                    if ( xl == 0 .and. yl == 0) then
                        cycle
                    endif
                    neighbour_index = neighbour_index + 1
                    if(y_n < 1 .or. y_n > nrows &
                        .or. x_n<1 .or. x_n > ncols) then
                        cycle
                    endif

                    if (flow_dir_grid(y_n,x_n) == dir_in(neighbour_index) &
                        .and. (mask(y_n,x_n).eqv..false.) ) then
                        node%x = x_n
                        node%y = y_n

                        call openCells%enqueue(node)
                        area_count = area_count + 1

                        mask(y_n,x_n) = .true.
                    endif
                enddo
            enddo
        enddo

        call openCells%free()

    end subroutine

end module dta_catch_cut
