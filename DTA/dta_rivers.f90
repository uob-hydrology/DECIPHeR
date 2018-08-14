module dta_rivers
    implicit none
contains

    ! distance calculated in d8 single flow direction
    ! measuring moves upstream in single flow direction only
    ! ensuring that neighbouring but non linked river cells are treated correctly
    subroutine river_find_dist(ncols, nrows, &
        riv_mask_grid, &
        flow_dir_grid, &
        riv_dist_grid, &
        is_from_sea, outlet_points)
        use dta_utility
        use dta_stack
        use dta_catch_cut
        implicit none
        integer :: nrows
        integer :: ncols
        logical :: riv_mask_grid(nrows,ncols)
        integer(1) :: flow_dir_grid(:,:)
        logical :: is_from_sea
        type(point_type), allocatable :: outlet_points(:)
        type(point_type) :: point

        double precision :: riv_dist_grid(nrows,ncols)

        ! locals
        integer :: i
        double precision :: dist_cardinal, dist_ordinal
        double precision :: dist
        type(stack_type) :: stack
        integer :: x,y,x_n,y_n,xl,yl
        integer :: neighbour_index
        integer :: loop_count_check, max_loop_count_check

        max_loop_count_check = nrows * ncols

        stack = new_stack()

        dist_cardinal = 1
        dist_ordinal = sqrt(2.0d0)

        do i=1,size(outlet_points)

            !print*,'outlet', i,'----------------------', outlet_points(i)%y, outlet_points(i)%x
            loop_count_check = 0

            ! only start if river is not labelled
            ! could already be labelled when:
            ! * sea cell flows into other sea cell
            ! * when trying to label outlets not linked to the sea
            if(riv_dist_grid(outlet_points(i)%y, outlet_points(i)%x) < 0.001) then
                riv_dist_grid(outlet_points(i)%y, outlet_points(i)%x) = 1.0d0
            else
                cycle
            endif

            call stack%push(outlet_points(i))

            do while(stack%item_count > 0)
                loop_count_check = loop_count_check + 1
                if(loop_count_check > max_loop_count_check) then
                    print *, 'error river stuck'
                    exit
                endif

                point = stack%pop()
                x = point%x
                y = point%y

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

                        if(riv_mask_grid(y_n,x_n)) then
                            ! if the neighbour river cell flows into current
                            if (flow_dir_grid(y_n,x_n) == dir_in(neighbour_index)) then

                                if(is_from_sea .and. riv_dist_grid(y_n,x_n) > 0.001) then
                                    ! sometimes the flow direction of adjacent sea cells causes a loop
                                    ! need to check if already processed
                                    ! in the case of checking outlets,
                                    !   already processed cells can be overwritten when further from the sea
                                    !
                                    !print*,'already processed',riv_dist_grid(y_n,x_n), loop_count_check,y_n,x_n
                                    cycle
                                endif

                                if (xl == 0 .or. yl == 0) then
                                    dist = dist_cardinal
                                else
                                    dist = dist_ordinal
                                endif
                                riv_dist_grid(y_n,x_n) = riv_dist_grid(y,x) + dist

                                point%x = x_n
                                point%y = y_n
                                !print*,'push',y_n,x_n
                                call stack%push(point)

                            endif

                        endif
                    enddo
                enddo
            enddo
        end do

        call stack%free()
    end subroutine


    subroutine river_find_outlets(nrows, ncols, &
        riv_mask_grid, sea_mask_grid, sea_outlets)
        !%RIVER_FIND_OUTLETS Find all river cells with a nodata neighbour (sea)
        !%
        use dta_stack
        use dta_utility
        implicit none
        integer :: nrows
        integer :: ncols
        ! .true. for river cell
        ! .false. for not river cell
        logical :: riv_mask_grid(nrows, ncols)
        ! .true.for the sea
        ! .false. for not the sea
        logical :: sea_mask_grid(nrows, ncols)
        ! output, list of cell indexes
        type(point_type), allocatable, dimension(:) :: sea_outlets

        type(stack_type) :: queue
        type(point_type) :: item
        integer :: row, col
        integer :: yl, xl, x, y
        integer :: neighbour_count

        queue = new_stack()

        !% starting cells are all river cells that are next to nodata values
        !print *, 'find starting cells (any river cell next to a sea=nodata)'

        do row=1,nrows
            do col=1,ncols
                !% ignore river cells on sea cells
                if (riv_mask_grid(row, col) &
                    .and.(sea_mask_grid(row, col).eqv..false.)) then
                    !% count nan neighbours
                    neighbour_count = 0

                    do yl = -1,1
                        y = row + yl
                        if(y < 1 .or. y > nrows) then
                            cycle
                        endif
                        do xl = -1,1
                            x = col + xl
                            if (xl == 0 .and. yl == 0) then
                                cycle
                            endif
                            if(x<1 .or. x > ncols) then
                                cycle
                            endif
                            if (sea_mask_grid(y,x).eqv..true. ) then
                                neighbour_count = neighbour_count + 1
                                exit
                            endif
                        end do
                    end do
                    !% river cell that is next to a nodata
                    !% add to queue
                    if(neighbour_count > 0) then
                        item%x = col
                        item%y = row
                        call queue%push(item)
                    endif
                endif
            end do
        end do

        !print *, 'found outlets', queue%item_count

        allocate(sea_outlets(queue%item_count))
        call queue%to_array(sea_outlets)
        call queue%free()
    end subroutine

    !
    subroutine river_find_headwater_pass1(nrows, ncols, riv_mask_grid, sea_outlets,&
        riv_dist_grid, headwater_pass1)
        !function [ riv_dist_grid, headwater_pass1_yx ] = headwater_pass1( riv, sea_outlets )
        !%HEADWATER_PASS1 measure river distances from outlets
        !%   riv: 0 not river cell 1: river cell
        !%   Use this function with the EA river network to find the candidate
        !%   headwaters
        !%
        !% grid to hold the shortest distance to the sea for each river cell
        !% distance is calculated without regard to flow direction, and can work on d8 or d4 river file
        !% distance suitable for locating headwaters
        !% initially -99 for non river cells and -1 for river cells
        !% First the algorithm finds all rivercells with a nodata neighbour and
        !% labels them with '1' and adds them to FIFO queue
        !% Second, each item in the queue is processed adding a distance
        !%
        !% headwater_pass1_yx potential headwaters - this pass1 version will greatly
        !% estimate, but makes a good starting point for headwaters_pass2
        use dta_queue
        use dta_utility
        implicit none
        integer :: nrows
        integer :: ncols
        logical :: riv_mask_grid(nrows,ncols)
        type(point_type) ::sea_outlets(:)
        double precision :: riv_dist_grid(nrows,ncols)
        type(point_type), allocatable, dimension(:) :: headwater_pass1

        integer :: total_river_cells
        type(queue_type) :: openCells
        type(queue_type) :: hw_pass1_queue
        integer :: ii
        integer :: yl, xl, y, x
        type(point_type) :: node
        type(point_type) :: node_tmp
        double precision :: counter, upCounter
        double precision :: dist_cardinal, dist_ordinal
        integer :: neighbour_count
        !integer :: progress_count

        total_river_cells = count(riv_mask_grid)

        ! set riv_dist_grid -1 on all river cells
        where (riv_mask_grid) riv_dist_grid = -1

        hw_pass1_queue = new_queue()
        openCells = new_queue()

        do ii = 1,size(sea_outlets)
            call openCells%enqueue(sea_outlets(ii))
            !% mark starting points
            riv_dist_grid(sea_outlets(ii)%y, sea_outlets(ii)%x) = 1
        end do

        !progress_count = 0

        !% estimated total add 10% some items are requeued
        !progress = progress_display(total_river_cells * 1.1);

        !print *, 'calculate pass 1 headwaters and river cell distance to sea'

        dist_cardinal = 1
        dist_ordinal = sqrt(2.0d0)

        do while(openCells%item_count > 0)
            node = openCells%dequeue()

            !progress_count = progress_count+1;
            !if(mod(progress_count,100) == 0)
            !    progress.show_progress(progress_count);
            !end

            counter = riv_dist_grid(node%y, node%x)
            neighbour_count = 0

            do yl = -1,1
                y = node%y + yl
                if(y<1 .or. y > nrows) then
                    cycle
                endif
                do xl = -1,1
                    x = node%x + xl
                    if (xl == 0 .and. yl == 0) then
                        cycle
                    end if
                    if(x<1 .or. x > ncols) then
                        cycle
                    end if

                    if (xl == 0 .or. yl == 0) then
                        upCounter = counter + dist_cardinal
                    else
                        upCounter = counter + dist_ordinal
                    end if

                    !% if the neighbour has not been marked
                    !% or it has already been marked, but this flow is shorter
                    !% (this may result in duplicate queueing, so need to check if already in queue)
                    if ((nint(riv_dist_grid(y, x)) == -1) &
                        .or. (riv_dist_grid(y, x) > upCounter)) then
                        neighbour_count = neighbour_count + 1
                        riv_dist_grid(y, x) = upCounter

                        node_tmp%y = y
                        node_tmp%x = x

                        call openCells%enqueue(node_tmp)
                    elseif (riv_dist_grid(y, x) > counter) then
                        neighbour_count = neighbour_count + 1
                    end if
                end do
            end do

            !% while counting the cells, if there are no directly connected upstream
            !% neighbours and more than 2 cells from the sea, add to headwaters list
            !% there will be mislabelled headwaters in this list, headwater_pass2
            !% should be used to cut down this list further
            if(neighbour_count == 0 .and. counter > 2) then
                call hw_pass1_queue%enqueue(node)
            endif

        end do

        call openCells%free()

        allocate(headwater_pass1(hw_pass1_queue%item_count))
        call hw_pass1_queue%to_array(headwater_pass1)
        call hw_pass1_queue%free()

        !print *, 'finished pass 1 headwaters and river cell distance to sea'

    end subroutine river_find_headwater_pass1


    subroutine river_find_headwater_pass2(nrows, ncols, riv_dist_grid, headwater_pass1, cellsize, &
        search_dist, move_downstream_dist, &
        headwater_pass2 )
        use dta_queue
        use dta_utility
        implicit none
        integer :: nrows
        integer :: ncols
        double precision :: riv_dist_grid(nrows,ncols)
        type(point_type) :: headwater_pass1(:)
        double precision :: cellsize
        double precision :: search_dist
        double precision :: move_downstream_dist
        type(point_type), allocatable, dimension(:) :: headwater_pass2

        !locals
        double precision :: searchLength
        double precision :: minCellLength
        double precision :: moveDownstreamLength
        logical :: headwater_grid(nrows,ncols)
        type(queue_type) :: headwater_queue
        type(queue_type) :: openCells
        integer :: state_mask(nrows,ncols)
        type(point_type) :: node
        type(point_type) :: node_check
        type(point_type) :: node_tmp
        integer :: i
        double precision :: headValue
        logical :: isHeadwater
        double precision :: downstreamTargetValue
        double precision :: downstreamFoundValue
        double precision :: neighbourValue
        integer :: downstream_y, downstream_x
        integer :: yl, xl, y, x


        !%HEADWATER_PASS2 Check each hw from pass1 by checking for any nearby upstream
        !%   the smallest headwater stream is approximatly half of search distance

        if(move_downstream_dist > search_dist) then
            print *,'moveDownstreamDistance must be >= searchDistance'
            stop
        endif

        !% how far to search for another upstream
        !% another upstream is determined by another cell that is higher than the current
        searchLength = search_dist / cellsize ! dist stored as cell unit length
        !% minimum distance from outlet to start checking headwaters
        minCellLength = searchLength / 2

        !% how far to move down from headwater
        !% avoid having headwaters on ridges
        !% This will help to match the river networks more closely as the
        !% very high headwater cells are more likely to flow the wrong
        !% direction than cells further down into the valley
        moveDownstreamLength = move_downstream_dist / cellsize


        !% flag 0 not headwater, 1 headwater
        headwater_grid(:,:) = .false.

        headwater_queue = new_queue()

        !% pre-allocate river_find_headwater_pass2
        !headwater_pass2_yx = int32(zeros(headwater_i1,2));
        !headwater_i2 = 0;

        !progress = progress_display(headwater_i1);

        openCells = new_queue()
        state_mask(:,:) = 0

        !% downstreamMask is the ends of the headwaters that have been cut off
        !% used for checking
        !if(nargout > 2)
        !    downstreamMask = int32(zeros(size(riv_dist_grid)));
        !else
        !    downstreamMask = [];
        !end

        !progress_count=0;

        !print *, 'headwater from pass 1', size(headwater_pass1)

        do i=1,size(headwater_pass1)
            !progress_count = progress_count+1;

            !if(mod(i,100) == 0)
            !    progress.show_progress(progress_count);
            !end

            node = headwater_pass1(i)
            headValue = riv_dist_grid(node%y,node%x)

            !% don't process matches so close to the outlet nodata
            if(headValue > minCellLength) then
                !% keep track of which cells have been queues from this cell
                !% prevent re-queueing
                state_mask(node%y,node%x) = i

                call openCells%enqueue(node)

                !% check all the cells surrounding the pass1 headwater
                !% if there are any cells with a higher count, then it is not a
                !% headwater
                isHeadwater = .true.

                !% for moving the headwater downstream, keep track of the location
                !% of the target downstream cell

                !% when searching to move downstream, keep above this value
                downstreamTargetValue = riv_dist_grid(node%y,node%x) - moveDownstreamLength
                !% the actual value found
                downstreamFoundValue = riv_dist_grid(node%y,node%x)
                downstream_y = -1
                downstream_x = -1

                do while (openCells%item_count > 0 .and. isHeadwater)
                    node_check = openCells%dequeue()
                    do yl = -1, 1
                        y = node_check%y + yl
                        if(y<1 .or. y > nrows) then
                            cycle
                        endif
                        do xl = -1,1
                            x = node_check%x + xl
                            if (xl == 0 .and. yl == 0) then
                                cycle
                            endif
                            if(x<1 .or. x > ncols) then
                                cycle
                            endif

                            neighbourValue = riv_dist_grid(y, x)
                            !%find lowest downstream which is still above the downstreamTargetValue
                            if(neighbourValue < downstreamFoundValue .and. &
                                neighbourValue >= downstreamTargetValue) then
                                downstreamFoundValue = neighbourValue
                                downstream_y = y
                                downstream_x = x
                            endif

                            if(neighbourValue > headValue .or. &
                                headwater_grid(y,x)) then
                                !% NOT A HEADWATER
                                !% - linked to another headwater
                                !% - linked to cell with higher value
                                isHeadwater = .false.
                            else if (neighbourValue > minCellLength .and. &
                                neighbourValue > headValue - searchLength .and. &
                                state_mask(y,x) /= i) then
                                state_mask(y,x) = i

                                node_tmp%y = y
                                node_tmp%x = x

                                call openCells%enqueue(node_tmp)

                            !                        if(display >= 2) then
                            !                            subplot(2,2,1)
                            !                            show_debug_grid(riv_dist_grid, node(1),node(2), 10, 'plain');
                            !                            subplot(2,2,4)
                            !                            show_debug_grid(stateMask, node(1),node(2), 10, 'filter', i);
                            !                            drawnow();
                            !                        endif
                            endif
                        end do
                    end do
                end do
                call openCells%clear()

                !% in order to remove the very tiny headwaters.
                if(moveDownstreamLength > 0 .and. &
                    isHeadwater .and. downstream_y /= -1) then

                    node%y = downstream_y
                    node%x = downstream_x
                    !if(~isempty(downstreamMask))
                    !    downstreamMask(downstream_y, downstream_x) = 3;
                    !end
                endif

                !        if(display >= 1)
                !            subplot(2,2,1)
                !            show_debug_grid(riv_dist_grid, node(1),node(2), 10, 'plain');
                !            subplot(2,2,4)
                !            show_debug_grid(stateMask, node(1),node(2), 10, 'filter', i);
                !            drawnow();
                !        end

                if(isHeadwater) then
                    !headwater_i2 = headwater_i2 + 1
                    !headwater_pass2_yx(headwater_i2,:) = node
                    call headwater_queue%enqueue(node)
                    headwater_grid(node%y,node%x) = .true.

                endif
            end if
        end do

        call openCells%free()

        !print *, 'pass 2 complete'

        !% trim river_find_headwater_pass2 to the actual number found
        !headwater_pass2_yx = headwater_pass2_yx(1:headwater_i2,:);

        allocate(headwater_pass2(headwater_queue%item_count))
        call headwater_queue%to_array(headwater_pass2)
        call headwater_queue%free()


    !%a = double(reshape(countMask, size(countMask,1) * size(countMask,2), 1));
    !%a(a<2) = NaN;
    !%hist(a );

    !% number of river cells not counted (not connected to sea/nodata)
    !non_counted = sum(sum(riv_dist_grid<-0.5));
    !fprintf('river cells not processed (no link to sea) %d\n', non_counted);
    !% number of headwater cells
    !fprintf('headwater cells %d\n', headwater_i2);


    end subroutine river_find_headwater_pass2


    subroutine river_single_flow(nrows, ncols, dem_grid, headwater, riv_label_grid, cellsize)
        use dta_utility
        implicit none
        integer :: nrows
        integer :: ncols
        double precision :: dem_grid(nrows, ncols)
        type(point_type) :: headwater(:)
        integer :: riv_label_grid(nrows, ncols)
        double precision :: cellsize
        integer :: nhw
        integer :: x, y, label

        !% Function to calculate UK rivers
        !% GC 23.02
        !% TD 2016/03/04 - error checking, display
        !% input: sink filled dem_grid with minimum slope
        !% input: headwaters stored in y, x

        riv_label_grid(:,:) = 0

        !print *, 'river_single_flow from ',size(headwater),'headwaters'

        !% Set all NaNs to -90
        where(dem_grid < -90) riv_label_grid = -99

        !progress = progress_display(size(hw_YX, 1));

        !skip_list = zeros(size(hw_YX, 1), 1);

        !if(display == 1)
        !    scratch = dem_grid;
        !    scratch(scratch<-90) = 500;
        !end

        do nhw = 1, size(headwater)
            y = headwater(nhw)%y
            x = headwater(nhw)%x

            if(riv_label_grid(y, x) == 0) then
                label = nhw * 1000
                call river_single_flow_point(nrows, ncols, dem_grid, cellsize, riv_label_grid, &
                    x, y, label, &
                    .false.)
            else
                ! label headwater as skipped
                headwater(nhw)%y = -1 * y
                headwater(nhw)%x = -1 * x
            endif
        end do
    end subroutine river_single_flow

    subroutine river_single_flow_point(nrows, ncols, dem_grid, cellsize, riv_label_grid, &
        x_start, y_start, label, &
        enable_fast_join_neighbour)
        implicit none
        integer, intent(in) :: ncols
        integer, intent(in) :: nrows
        double precision :: dem_grid(nrows, ncols)
        double precision, intent(in) :: cellsize
        integer :: riv_label_grid(nrows, ncols)
        integer, intent(in) :: x_start
        integer, intent(in) :: y_start
        integer, intent(in) :: label
        logical, intent(in) :: enable_fast_join_neighbour

        double precision :: dist
        double precision :: dist_cardinal
        double precision :: dist_ordinal
        logical :: isdone
        double precision :: max_slope
        integer :: riv_value
        double precision :: slope

        integer :: x
        integer :: y
        integer :: x_max
        integer :: x_n
        integer :: xl

        integer :: y_max
        integer :: y_n
        integer :: yl

        x = x_start
        y = y_start

        dist_cardinal = cellsize
        dist_ordinal = sqrt(2.0) * cellsize

        isdone = .false.
        !        if(riv_label_grid(y, x) /= 1) then
        !            isdone = .true.
        !            print*,'skip', riv_label_grid(y, x)
        !        endif

        do while (isdone.eqv..false.)
            riv_value = riv_label_grid(y, x)

            riv_label_grid(y, x) = label

            !% find max slope
            max_slope = -1
            x_max = 0
            y_max = 0

            !% Do local search and find maximum gradient downstream
            do yl = -1,1
                y_n = y + yl
                if(y_n < 1 .or. y_n > nrows) then
                    !% hit boundary
                    isdone = .true.
                    exit
                endif
                do xl = -1,1
                    if (xl == 0 .and. yl == 0) then
                        cycle
                    endif
                    x_n = x + xl
                    if(x_n<1 .or. x_n > ncols) then
                        isdone = .true.
                        exit
                    endif

                    riv_value = riv_label_grid(y_n, x_n)

                    !if(enable_find_link.eqv..false.) then
                    !% already labelled as this river segment
                    !% handle flat areas (should not be in sink fill noflat dem_grid)
                    if(riv_value == label) then
                        cycle
                    endif
                    !endif

                    if(enable_fast_join_neighbour) then
                        !% if any neighbour is nodata or joins to an exising river
                        !% cell, then join to this regardless of slope
                        !% this ensures no clumping and clean junctions in D8
                        ! this check is no longer required
                        if(riv_value < -90 .or. &
                            ((riv_value > 0) .and. (riv_value /= label))) then
                            ! x_end, y_end set to the joined river cell
                            isdone = .true.
                            exit
                        endif
                    endif

                    if (xl == 0 .or. yl == 0) then
                        dist = dist_cardinal
                    else
                        dist = dist_ordinal
                    endif
                    slope = (dem_grid(y, x) - dem_grid(y_n, x_n)) / dist
                    if (slope > max_slope) then
                        max_slope = slope
                        x_max = x_n
                        y_max = y_n
                    endif
                enddo
            enddo
            if ((isdone.eqv..false.).and. max_slope<0) then
                ! when at edge
                !print *, 'max gradient ', max_slope
                isdone = .true.
            endif

            x = x_max
            y = y_max
        end do
    end subroutine river_single_flow_point


    subroutine river_find_lowest(nrows, ncols, dem_grid, riv_dist_grid, x_start, y_start, x_out, y_out)
        use dta_stack
        use dta_utility
        implicit none
        integer, intent(in) :: ncols
        integer, intent(in) :: nrows
        double precision :: dem_grid(nrows, ncols)
        double precision :: riv_dist_grid(nrows, ncols)
        integer, intent(in) :: x_start
        integer, intent(in) :: y_start
        integer, intent(out) :: x_out
        integer, intent(out) :: y_out

        double precision :: lowest

        integer :: x
        integer :: y
        integer :: x_n
        integer :: xl

        integer :: y_n
        integer :: yl
        type(stack_type) :: stack
        type(point_type) :: item


        stack = new_stack()

        lowest = riv_dist_grid(y_start, x_start)
        x_out = x_start
        y_out = y_start

        item%x = x_start
        item%y = y_start
        call stack%push(item)

        do while (stack%item_count > 0)

            item = stack%pop()
            x = item%x
            y = item%y

            riv_dist_grid(y, x) = 1

            if(dem_grid(y,x) < lowest) then
                lowest = dem_grid(y,x)
                x_out = x
                y_out = y
            endif


            !% Do local search and find maximum gradient downstream
            do yl = -1,1
                y_n = y + yl
                if(y_n < 1 .or. y_n > nrows) then
                    !% hit boundary
                    cycle
                endif
                do xl = -1,1
                    if (xl == 0 .and. yl == 0) then
                        cycle
                    endif
                    x_n = x + xl
                    if(x_n<1 .or. x_n > ncols) then
                        cycle
                    endif

                    ! scan all connected river cells where dist == -1
                    if(abs(riv_dist_grid(y_n, x_n) + 1) < 0.0001) then
                        item%x = x_n
                        item%y = y_n
                        call stack%push(item)
                    endif
                enddo
            enddo
        end do

        call stack%free()

    end subroutine river_find_lowest

end module dta_rivers
