module dta_routing_file
contains

    subroutine routing_file_find_rivers(nrows, ncols, node_list, riv_label_grid, &
        river_point_lists, total_river_cells)
        use dta_utility
        use dta_riv_tree_node
        use dta_stack
        ! builds lists of points for each river id
        !
        ! used for reading river information from grids
        ! finds all river points, and puts them in lists matching their respective riv_id
        ! each index matching node_list contains a list of points where all the river cells are
        implicit none
        integer :: nrows
        integer :: ncols
        type(riv_tree_node) :: node_list(:)
        integer :: riv_label_grid(nrows,ncols)

        type(point_list_type), dimension(:) :: river_point_lists(:)
        integer :: total_river_cells


        ! locals
        integer :: x, y
        integer :: riv_label
        ! for each node/riv_id create a stack to store the found points (don't know in advance how many
        ! these stacks get copied to arrays of the correct size and returned
        type(stack_type), allocatable, dimension(:) :: river_point_stacks
        integer :: i
        type(point_type) :: point

        allocate(river_point_stacks(size(node_list)))

        do i = 1,size(node_list)
            river_point_stacks(i) = new_stack()
        end do

        do x=1,ncols
            do y=1,nrows
                ! look through all cells and make list of the point locations
                ! equivilent to matlab:
                ! for ii=1:length(node_list)
                !    area_values = area(riv_labelled==node_list(ii).gauge_id);
                !    dist_values = riv_dist(riv_labelled==node_list(ii).gauge_id);

                if(riv_label_grid(y,x) /= 0) then
                    riv_label = riv_label_grid(y,x)
                    ! lookup the index that matches the node value
                    do i = 1,size(node_list)
                        if(node_list(i)%gauge_id == riv_label) then
                            point%y = y
                            point%x = x
                            call river_point_stacks(i)%push(point)
                            exit
                        endif
                    end do
                endif
            end do
        end do


        total_river_cells = 0
        !convert all the stacks to arrays
        do i = 1,size(node_list)
            total_river_cells = total_river_cells + river_point_stacks(i)%item_count
            allocate(river_point_lists(i)%list(river_point_stacks(i)%item_count))
            call river_point_stacks(i)%to_array(river_point_lists(i)%list)
            call river_point_stacks(i)%free()
        enddo
        !print *, 'total_river_cells', total_river_cells




    end subroutine routing_file_find_rivers


    subroutine routing_file_grids_to_list(nrows, ncols, cellsize, &
        node_list, &
        river_point_lists, total_river_cells, &
        riv_dist_grid, area_grid, dem_grid, &
        river_data)
        use dta_utility
        use dta_riv_tree_node
        implicit none
        ! used by route_river_file
        ! converts the grid data to a list for each river cell
        ! needs to be run after atb - as it requires the river area contributions
        integer :: nrows
        integer :: ncols
        double precision :: cellsize
        type(riv_tree_node) :: node_list(:)
        type(point_list_type), allocatable, dimension(:) :: river_point_lists
        double precision :: riv_dist_grid(nrows,ncols)
        double precision :: area_grid(nrows,ncols), dem_grid(nrows, ncols)
        !% riv_id, area, dist, section_dist, slope, height
        double precision, allocatable, dimension(:,:) :: river_data

        !locals
        integer :: riv_cell_count
        integer :: total_river_cells
        integer :: i, j
        type(point_type) :: point, min_point
        integer :: max_slope_index
        double precision :: max_slope, dist_cardinal, dist_ordinal
        integer :: river_data_start
        integer :: river_data_end
        double precision, allocatable, dimension(:) :: river_col_tmp
        double precision, allocatable, dimension(:) :: river_col_tmp2

        double precision :: min_dist

        ! total number of river cells
        !total_river_cells = count(riv_dist_grid > 0.001)

        allocate(river_data(total_river_cells, 8))
        river_data(:,:) = 0
        river_data_start = 1
        !%figure

        dist_cardinal = 1 * cellsize
        dist_ordinal = sqrt(2.0) * cellsize

        do i=1,size(node_list)
            !area_values = area(riv_labelled==node_list(ii).gauge_id);
            !dist_values = riv_dist(riv_labelled==node_list(ii).gauge_id);

            !% each cell can have:
            !%  1    id
            !%  2    area
            !%  3    dist (total to final outlet)
            !%  4    section_dist (to bottom of this river_id)
            !%
            !%  5    elevation
            !%  6    slope (single flow direction slope to next river cell)
            !%  7    ds_dist
            !%  8    ds_elevation

            riv_cell_count = size(river_point_lists(i)%list)
            river_data_end = river_data_start + riv_cell_count-1

            river_data(river_data_start:river_data_end,1) = node_list(i)%gauge_id

            allocate(river_col_tmp(riv_cell_count))
            allocate(river_col_tmp2(riv_cell_count))

            do j= 1,riv_cell_count
                point = river_point_lists(i)%list(j)
                river_col_tmp(j) = riv_dist_grid(point%y, point%x)
            end do

            ! write the distance to outlet values
            river_data(river_data_start:river_data_end,3) = river_col_tmp

            ! update distance values to be distance to reach outlet
            min_dist = minval(river_col_tmp)
            river_col_tmp(:) = river_col_tmp(:) - min_dist
            river_data(river_data_start:river_data_end,4) = river_col_tmp

            ! read the area accumulations from the grid
            do j= 1,riv_cell_count
                point = river_point_lists(i)%list(j)
                river_col_tmp(j) = area_grid(point%y, point%x)
            end do

            river_data(river_data_start:river_data_end,2) = river_col_tmp

            ! read the elevations of each river cell from grid
            do j = 1, riv_cell_count
                point = river_point_lists(i)%list(j)
                river_col_tmp(j) = dem_grid(point%y, point%x)
            end do
            river_data(river_data_start:river_data_end,5) = river_col_tmp

            ! find the downstram values of elevation and distance
            do j = 1, riv_cell_count

                point = river_point_lists(i)%list(j)

                ! find downstream point (greatest slope)
                call cell_flow_dir(ncols, nrows, dem_grid, &
                    dist_cardinal, dist_ordinal, point%x, point%y, &
                    max_slope_index, min_point, max_slope)
                if(point%x == min_point%x .and. point%y == min_point%y) then
                    print*, 'No slope'
                endif
                if(dem_grid(min_point%y, min_point%x) < -90) then
                    print*, 'flow out of dem'
                endif

                if(riv_dist_grid(min_point%y, min_point%x) >= riv_dist_grid(point%y, point%x)) then

                    ! sea outlet flow direction may not take no-data into account
                    ! if centre cell distance is next to the boundary, use 0 as downstream
                    if(riv_dist_grid(point%y, point%x) < dist_ordinal+0.001) then
                        river_col_tmp(j) = 0
                        ! approximate slope (neighbour is going to be higher)
                        ! to get a reasonable slope, reduce height by the rise to the neighbour

                        !print*, dem_grid(point%y, point%x)
                        !print*, dem_grid(min_point%y, min_point%x)
                        river_col_tmp2(j) = dem_grid(point%y, point%x) - &
                            (dem_grid(min_point%y, min_point%x) - dem_grid(point%y, point%x))
                    else
                        print*, 'invalid downstream distance (not river cell)', riv_dist_grid(min_point%y, min_point%x)
                        stop
                    endif
                else
                    river_col_tmp(j) = riv_dist_grid(min_point%y, min_point%x)
                    river_col_tmp2(j) = dem_grid(min_point%y, min_point%x)
                endif



            end do

            ! downstream dist
            river_data(river_data_start:river_data_end,7) = river_col_tmp
            !downstream elevation
            river_data(river_data_start:river_data_end,8) = river_col_tmp2

            deallocate(river_col_tmp)
            deallocate(river_col_tmp2)

            river_data_start = river_data_end + 1

        end do

        ! slope         =      (elevation  - ds_elevation) /    (dist            - ds_dist)
        river_data(:,6) = (river_data(:,5) - river_data(:,8)) / (river_data(:,3) - river_data(:,7))

    end subroutine routing_file_grids_to_list

end module dta_routing_file
