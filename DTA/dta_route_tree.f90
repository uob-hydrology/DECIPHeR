module dta_route_tree
    implicit none



contains


    subroutine routing_dta(nrows, ncols, point_row_list, point_list, sea_outlets, dem_grid, &
        flow_dir_grid, &
        riv_mask, &
        riv_dist_grid, &
        river_reach_increment_m, river_reach_increment_m_min, &
        point_col_id, &
        node_list, riv_label_grid, fp)
        !%ROUTING_DTA Create routing tree structure
        !%   by linking linkages based on where they sit on the river
        !%   dem_grid - should be masked to catchment boundary - to ensure river cells are
        !%   next to a nodata cell (unless force_sea_outlets is set)
        !%   force_sea_outlets - set the start points for the downstream endpoint of
        !%   the river
        use dta_queue
        use dta_rivers
        use dta_utility
        use dta_riv_tree_node
        implicit none
        integer, intent(in) :: nrows
        integer, intent(in) :: ncols
        double precision, intent(in) :: point_row_list(:,:)
        type(point_type), intent(in) :: point_list(:)
        type(point_type), intent(in) :: sea_outlets(:)
        double precision, intent(in) :: dem_grid(nrows, ncols)
        integer(1), intent(in) :: flow_dir_grid(:,:)
        logical, intent(in) :: riv_mask(nrows, ncols)
        double precision, intent(in) :: riv_dist_grid(nrows, ncols)
        double precision, intent(in) :: river_reach_increment_m, river_reach_increment_m_min
        integer, intent(in) :: point_col_id
        ! size is size(point_list) + size(sea_outlets)
        type(riv_tree_node), allocatable :: node_list(:)
        integer :: riv_label_grid(nrows, ncols)

        integer :: tmp

        !locals
        integer :: outlet_count
        integer :: i
        integer :: point_row, point_col
        integer :: existing_gauge_id
        type(riv_tree_node), allocatable :: node_list_copy(:)
        type(point_type), allocatable :: reach_intervals(:)
        integer :: down
        integer :: remove_count
        integer :: node_index
        integer :: gauge_end_index
        integer :: new_length
        integer :: x, y, fp

        double precision :: river_reach_increment_cells

        outlet_count = size(sea_outlets)

        !print *, 'point_list size', size(point_row_list,1)
        !print *, 'sea_outlets size', outlet_count
        !% flat list of all gauge nodes and outlets
        new_length = outlet_count + size(point_row_list,1)
        allocate(node_list(new_length))


        riv_label_grid(:,:) = 0
        !% update all rivers to be -1
        !% this is so that gauge cells can be labelled with gauge index, which could
        !% be 1
        where (riv_mask) riv_label_grid = -1
        where (dem_grid < -90) riv_label_grid = -99

        !% label river cells with gauge locations
        !% gauge will be on top of the river cell
        do i=1,size(point_row_list,1)

            point_row = point_list(i)%y
            point_col = point_list(i)%x

            !call NorthingEastingToRowCol(point_row_list(i,point_col_north), point_row_list(i,point_col_east),&
            !    nrows, xllcorner, yllcorner, cellsize, point_row, point_col)

            if(riv_label_grid(point_row, point_col) /= -1) then
                if(riv_mask(point_row, point_col) .eqv..false.) then
                    write(fp,*) 'gauge not on river cell:', point_row_list(i,point_col_id)

                    !old_value = riv(point_row(i), point_col(i));
                    !riv(point_row(i), point_col(i)) = 1;
                    !show_debug_grid(riv, point_row(i),point_col(i), 10, 'none', 0);
                    !riv(point_row(i), point_col(i)) = old_value;
                else
                    existing_gauge_id = nint(point_row_list(riv_label_grid(point_row, point_col),point_col_id))
                    !% note this shouldn't occurr as pre-processor filters these out
                    write(fp,*) 'gauge located on existing gauge: this gauge: ', &
                        nint(point_row_list(i,point_col_id)), existing_gauge_id

                    !old_value = riv(point_row(i), point_col(i));
                    !riv(point_row(i), point_col(i)) = 1;
                    !show_debug_grid(riv, point_row(i),point_col(i), 10, 'none', 0);
                    !riv(point_row(i), point_col(i)) = old_value;
                endif
            endif

            !% first items in node_list reserved for the sea outlet nodes
            node_index = outlet_count + i

            node_list(node_index) = new_riv_tree_node()
            node_list(node_index)%node_type = NODE_TYPE_GAUGE
            node_list(node_index)%gauge_id = nint(point_row_list(i,point_col_id))
            node_list(node_index)%row = point_row
            node_list(node_index)%col = point_col

            !% label river cells with node_indexes
            riv_label_grid(point_row, point_col) = node_index
        enddo

        !% process reach lengths
        if(river_reach_increment_m > 0.001) then
            river_reach_increment_cells = river_reach_increment_m;

            call routing_reach_intervals(nrows, ncols, sea_outlets, &
                riv_label_grid, riv_dist_grid, &
                river_reach_increment_cells, &
                reach_intervals)

            gauge_end_index = size(node_list)
            !% allocate new space for all the reach increments

            new_length = gauge_end_index + size(reach_intervals)

            ! moves the node_list to node_list_copy
            ! node_list becomes deallocated
            call move_alloc(node_list, node_list_copy)
            allocate(node_list(new_length))
            node_list(1:gauge_end_index) = node_list_copy

            do i=1,size(reach_intervals)

                if (riv_label_grid(reach_intervals(i)%y, reach_intervals(i)%x) /= -1) then
                    write(fp,*) 'calculated reach increment not on river cell'
                    !%riv(point_row(i), point_col(i)) = 1;
                    !show_debug_grid(riv, point_row(i), point_col(i), 10, 'none', 0);
                endif

                node_index = gauge_end_index + i
                node_list(node_index) = new_riv_tree_node()
                node_list(node_index)%node_type = NODE_TYPE_RIVER
                node_list(node_index)%gauge_id = -1*i
                node_list(node_index)%row = reach_intervals(i)%y
                node_list(node_index)%col = reach_intervals(i)%x

                !% label river cells with node_indexes
                riv_label_grid(reach_intervals(i)%y, reach_intervals(i)%x) = node_index
            end do


        endif

        write(fp,*) 'routing_link_rivers node_list:', size(node_list)
                !call routing_link_rivers(nrows, ncols, sea_outlets, node_list, riv_label_grid )

        do i=1,size(node_list)
            node_list(i)%downstream_index = 0
        end do

        call routing_link_rivers_sfd(nrows, ncols, flow_dir_grid, sea_outlets, node_list, riv_label_grid, fp)

                !% note riv grid is now labelled with the stations
                !% could output here if it is useful

        call build_upstream_from_downstream(node_list)

        !        print *, 'create upstream links:', size(node_list)
        !        !% create the downstream neighbour links
        !        do i=1,size(node_list)
        !            print*,i,'upstream:',node_list(i)%upstream_count
        !            do jj=1,node_list(i)%upstream_count
        !                upstream_index = node_list(i)%upstream_indexes(jj)
        !                node_list(upstream_index)%downstream_index = i
        !            end do
        !        end do

        do i=1,size(node_list)
            if(node_list(i)%downstream_index == 0) then
                call riv_node_recurse_label(i, node_list, node_list(i)%gauge_id)
            endif
        end do

        write(fp,*) 'routing_set_node_properties:', size(node_list)
        call routing_set_node_properties(nrows, ncols, node_list, dem_grid, riv_dist_grid)

        remove_count = 0

        !% remove all river reach from the tree where the distance to the next gauge
        !% is less than threshold
        ! there is some problem with removing - disabled for now
        ! result of disabling is there could be some reach close to next upstream
        do i=1,size(node_list)
            if(node_list(i)%downstream_index > 0) then
                down = node_list(i)%downstream_index
                !% if current is gauge and downstream is a reach interval
                if(node_list(i)%node_type == NODE_TYPE_GAUGE .and. node_list(down)%node_type == NODE_TYPE_RIVER) then
                    if(node_list(i)%downstream_dist < river_reach_increment_m_min) then
                        !% river reach too short, remove
                        print *, 'remove short downstream', &
                            node_list(i)%gauge_id, node_list(i)%downstream_dist
                        call riv_node_remove(node_list, down)
                        remove_count = remove_count + 1
                    end if
                end if
            end if
        end do

        !% re-calculate the distances if any nodes removed
        if(remove_count > 0) then
            call routing_set_node_properties(nrows, ncols, node_list, dem_grid, riv_dist_grid)
        endif

        !% in order to give the reach lengths
        !% update all gaugeid's to leave 3 digits for reach length id's
        !% reach length will be gauge_id001,gauge_id002...gauge_id00n
        !% sea outlets have negative id ensuring they don't clash with gauges if
        !% over 1000 outlets (lowest gauge id is 1001)
        !%
        do i=1,size(node_list)
            if(node_list(i)%node_type == NODE_TYPE_SEA ) then
                node_list(i)%gauge_id = node_list(i)%gauge_id * 1000
            endif
            if(node_list(i)%node_type == NODE_TYPE_GAUGE ) then
                node_list(i)%gauge_id = node_list(i)%gauge_id * 1000
            endif
        end do

        do i=1,size(node_list)
            if(node_list(i)%node_type == NODE_TYPE_SEA) then
                !% sea outlets have negative values
                tmp = node_list(i)%gauge_id - 1
                call riv_node_recurse_label(i, node_list, tmp)
            endif
        end do


        do x=1,ncols
            do y=1,nrows
                if(riv_label_grid(y,x) > 0) then
                    riv_label_grid(y,x) = node_list(riv_label_grid(y,x))%gauge_id
                endif
            end do
        end do


        ! slow way!!!
        !do i=1,size(node_list)
         !   where (riv_label_grid==i) riv_label_grid = node_list(i)%gauge_id
        !end do


    end subroutine routing_dta

    subroutine routing_link_rivers_sfd(nrows, ncols, flow_dir_grid, sea_outlets, node_list, riv_label_grid, fp)
        !%routing_link_rivers Link river points based on where they sit on rivers
        !%   Input: sea_outlets - list of downstream points [y1,x1; y2, x2; ...]
        !%   Input: riv - river mask 0 no data, -1: normal river cell, >0 river cell
        !%   labelled with with indexes into node_list
        !%   node_list
        !%   Input: node_list - list of nodes (with unknown linkages)
        !%
        !%   Output: node_list - list of nodes with downstream


        use dta_stack
        use dta_utility
        use dta_rivers
        use dta_riv_tree_node
        use dta_catch_cut
        implicit none
        integer :: nrows
        integer :: ncols
        integer(1) :: flow_dir_grid(:,:)
        type(point_type) :: sea_outlets(:)
        type(riv_tree_node) :: node_list(:)
        integer :: riv_label_grid(nrows, ncols)

        !locals
        type(stack_type) :: stack
        integer :: outlet_index
        type(point_type) :: point
        integer :: node_index
        integer :: label
        integer :: x, y, fp
        integer :: xl, yl, x_n, y_n
        integer :: neighbour_index
        integer :: neighbour_value

        stack = new_stack()

        !% label all outlets
        do outlet_index = 1,size(sea_outlets)
            point = sea_outlets(outlet_index)

            node_list(outlet_index) = new_riv_tree_node()
            node_list(outlet_index)%node_type = NODE_TYPE_SEA
            !% no gauge id as this is the outlet to sea
            node_list(outlet_index)%gauge_id = -1*outlet_index
            node_list(outlet_index)%row = point%y;
            node_list(outlet_index)%col = point%x

            riv_label_grid(point%y,point%x) = outlet_index
        end do

        do node_index=1,size(node_list)
            if(node_list(node_index)%enabled .eqv..false.) then
                cycle
            endif

            label = node_index

            point%x = node_list(node_index)%col
            point%y = node_list(node_index)%row

            call stack%push(point)

            do while(stack%item_count > 0)
                !loop_count_check = loop_count_check + 1
                !if(loop_count_check > max_loop_count_check) then
                !    print *, 'error catch cut stuck'
                !    exit
                !endif

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


                        ! if the neighbour river cell flows into current
                        if (flow_dir_grid(y_n,x_n) == dir_in(neighbour_index)) then
                            ! this will be an index of a node that flows down into the current
                            neighbour_value = riv_label_grid(y_n,x_n)

                            ! not a river, so ignore it
                            if(neighbour_value == 0 .or. neighbour_value < -90) then
                                cycle
                            endif

                            ! this is a labelled river - a node which flows into here
                            if(riv_label_grid(y_n,x_n) > 0) then
                                if (riv_label_grid(y_n,x_n) /= label) then
                                    node_list(riv_label_grid(y_n,x_n))%downstream_index = label
                                else
                                    write(fp,*) 'node flows to self', node_list(label)%gauge_id
                                endif
                            else if(riv_label_grid(y_n,x_n) > -90) then
                                riv_label_grid(y_n,x_n) = label

                                point%x = x_n
                                point%y = y_n

                                call stack%push(point)
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

        ! sea cells cannot have downstream cells
        ! the above flow direction calculation may produce these where a sea cell flows from another sea cell
        ! reset them all to 0 here
        do node_index=1,size(node_list)
            if(node_list(node_index)%node_type == NODE_TYPE_SEA) then
                node_list(node_index)%downstream_index = 0
            endif
        end do

        call stack%free()
    end subroutine routing_link_rivers_sfd


    subroutine routing_set_node_properties(nrows, ncols, node_list, dem_grid, riv_dist_grid )
        !%routing_set_node_properties calculate the node physical properties
        !%   height, river distance
        !%   delta_height, delta distance and slope
        use dta_riv_tree_node
        implicit none
        integer :: nrows
        integer :: ncols
        type(riv_tree_node) :: node_list(:)
        double precision :: dem_grid(nrows,ncols)
        double precision :: riv_dist_grid(nrows,ncols)

        !locals
        integer :: i
        type(riv_tree_node) :: downstream

        !% calculate the neighbour physical properties
        do i=1,size(node_list)
            if(node_list(i)%enabled) then
                node_list(i)%point_h = dem_grid(node_list(i)%row, node_list(i)%col)
                node_list(i)%point_dist = riv_dist_grid(node_list(i)%row, node_list(i)%col)
            endif
        end do
        do i=1,size(node_list)
            if((node_list(i)%enabled) .and. (node_list(i)%downstream_index > 0)) then
                downstream = node_list(node_list(i)%downstream_index)
                node_list(i)%downstream_dist = (node_list(i)%point_dist - downstream%point_dist)
                node_list(i)%downstream_delta_h = node_list(i)%point_h - downstream%point_h
                node_list(i)%downstream_slope = node_list(i)%downstream_delta_h / node_list(i)%downstream_dist
            endif
        end do

    end subroutine

    subroutine routing_reach_intervals( nrows, ncols, sea_outlets, &
        riv_label_grid, riv_dist_grid, increment, reach_points )
        !%routing_reach_intervals Move up from sea creating points atincrements up
        !%the river network
        !%   Input: sea_outlets - list of downstream points [y1,x1; y2,x2; ...]
        !%   Input: riv - river mask 0 no data, -1: normal river cell, >0 river cell
        !%          labelled with with indexes into node_list (used to restart
        !%          counting at known points)
        !%   Input: riv_dist_grid - flow length distance for each river cell
        !%   Increment: distance increment in units is to match riv_dist_grid
        !%
        !%   Output: row, col - list of points on river network at regular intervals
        use dta_utility
        use dta_queue
        implicit none
        integer :: nrows
        integer :: ncols
        type(point_type), intent(in) :: sea_outlets(:)
        integer, intent(in) :: riv_label_grid(nrows,ncols)
        double precision, intent(in) :: riv_dist_grid(nrows,ncols)
        double precision, intent(in) :: increment
        type(point_type), allocatable :: reach_points(:)

        ! locals
        type(queue_type) :: result_queue
        type(queue_type) :: search_queue
        integer :: riv_scratch_grid(nrows,ncols)
        logical :: state_grid(nrows,ncols)
        double precision :: target_grid(nrows,ncols)
        integer :: outlet_index
        type(point_type) :: cell_yx
        type(point_type) :: node_tmp
        integer :: riv_value
        integer :: next_riv_value
        double precision :: target_value
        logical :: increment_dist
        logical :: add_marker
        integer :: yl, xl, y_n, x_n

        !% add points to list
        result_queue = new_queue()
        search_queue = new_queue()

        ! make a copy of the river grid
        riv_scratch_grid(:,:) = riv_label_grid(:,:)

        state_grid(:,:) = .false.
        target_grid(:,:) = 0

        !% start at each outlet and move upstream
        do outlet_index = 1,size(sea_outlets)
            cell_yx = sea_outlets(outlet_index)

            riv_scratch_grid(cell_yx%y,cell_yx%x) = outlet_index
            call search_queue%enqueue(cell_yx)
            state_grid(cell_yx%y,cell_yx%x) = .true.

            !% state grid labelled with target
            target_grid(cell_yx%y,cell_yx%x) = increment

            do while(search_queue%item_count > 0)
                cell_yx = search_queue%Dequeue()

                riv_value = riv_scratch_grid(cell_yx%y,cell_yx%x)
                next_riv_value = riv_value

                target_value = target_grid(cell_yx%y,cell_yx%x)
                increment_dist = .false.
                add_marker = .false.
                !% search for any gauges, or any increment in target
                do yl = -1,1
                    y_n = cell_yx%y + yl
                    if (y_n < 1 .or. y_n > nrows) then
                        cycle
                    endif
                    do xl = -1,1
                        if (xl == 0 .and. yl == 0) then
                            cycle
                        endif
                        x_n = cell_yx%x + xl
                        if(x_n<1 .or. x_n > ncols) then
                            cycle
                        endif
                        if(state_grid(y_n, x_n) .eqv..false.) then
                            if( riv_scratch_grid(y_n, x_n) > 0 ) then
                                !% this an already marked station
                                !% reset the increment counter
                                !% don't add a new mark
                                next_riv_value = riv_scratch_grid(y_n,x_n)
                                increment_dist = .true.
                                add_marker = .false.
                            elseif(riv_dist_grid(y_n, x_n) > target_value) then
                                !% this cell triggers a new interval
                                !% reset the increment counter
                                !% add a new mark
                                increment_dist = .true.
                                add_marker = .true.
                            endif
                        endif
                    end do
                end do

                if (increment_dist) then
                    target_value = riv_dist_grid(cell_yx%y,cell_yx%x) + increment
                    if (add_marker) then
                        call result_queue%enqueue(cell_yx)
                    end if
                end if

                !%search and queue neighbours
                do yl = -1,1
                    y_n = cell_yx%y + yl
                    if (y_n < 1 .or. y_n > nrows) then
                        cycle
                    endif
                    do xl = -1,1
                        if (xl == 0 .and. yl == 0) then
                            cycle
                        endif
                        x_n = cell_yx%x + xl
                        if (x_n<1 .or. x_n > ncols) then
                            cycle
                        endif

                        if(riv_scratch_grid(y_n, x_n) /= 0 &
                            .and. (state_grid(y_n, x_n) .eqv..false.)) then
                            !% move the label upstream
                            riv_scratch_grid(y_n,x_n) = next_riv_value
                            state_grid(y_n, x_n) = .true.
                            target_grid(y_n, x_n) = target_value;
                            node_tmp%y = y_n
                            node_tmp%x = x_n
                            call search_queue%Enqueue(node_tmp)
                        endif
                    end do
                end do
            end do
        end do


        call search_queue%free()

        !% copy points back to a n by 2 matrix
        allocate(reach_points(result_queue%item_count))
        call result_queue%to_array(reach_points)
        call result_queue%free()
    end subroutine routing_reach_intervals








end module dta_route_tree
