module dta_riv_tree_node
    implicit none
    integer, parameter :: V_MODE_CONSTANT = 1
    integer, parameter :: V_MODE_SLOPE = 2

    integer, parameter :: NODE_TYPE_SEA = 1
    integer, parameter :: NODE_TYPE_GAUGE = 2
    integer, parameter :: NODE_TYPE_RIVER = 3 ! classed as 'r' river or 'e' (eturary=river before any gauge)



    type riv_tree_node
        !% 1=sea outlet
        !% 2=gauge
        !% 3=river reach increment
        integer :: node_type
        integer :: gauge_id
        !%list of cells directly connected upstream
        integer, allocatable, dimension(:) :: upstream_indexes
        integer upstream_count

        !%list of all upstream indexes from the entire tree
        integer, allocatable, dimension(:) :: upstream_tree_indexes
        integer :: upstream_tree_count

        integer :: row
        integer :: col

        logical :: enabled

        !% single index of the directly connected downstream
        integer :: downstream_index

        !% extra physical properties
        double precision :: point_h !% height of point on dem
        double precision :: point_dist !% distance of point on river to outlet

        double precision :: downstream_dist
        double precision :: downstream_delta_h
        double precision :: downstream_slope

        double precision :: reach_delay
        double precision :: total_downstream_delay
        double precision :: frac_sub_area
        double precision :: catch_area

    end type riv_tree_node

contains

    function new_riv_tree_node()
        implicit none
        type(riv_tree_node) :: new_riv_tree_node

        new_riv_tree_node%node_type = 0
        new_riv_tree_node%gauge_id = 0
        new_riv_tree_node%enabled = .true.

        ! when counts 0 index arrays are not allocated
        new_riv_tree_node%upstream_count = 0
        new_riv_tree_node%upstream_tree_count = 0

        new_riv_tree_node%downstream_index = 0

        new_riv_tree_node%row = 0
        new_riv_tree_node%col = 0

        new_riv_tree_node%point_h = 0
        new_riv_tree_node%point_dist = 0
        new_riv_tree_node%downstream_dist = 0
        new_riv_tree_node%downstream_delta_h = 0
        new_riv_tree_node%downstream_slope = 0
        new_riv_tree_node%total_downstream_delay = 0

    end function new_riv_tree_node

    !% build the full upstream linkages for each node
    !% uses the upstream_indexes to calculate the full upstream tree for each node
    ! every node with have a list of every upstream node
    recursive subroutine riv_node_full_upstream(target_index, node_index, node_list)
        use dta_utility
        implicit none
         ! index of the node who's tree is being computed
        integer, intent(in) :: target_index
        ! recurse node index
        integer, intent(in) :: node_index
        type(riv_tree_node) :: node_list(:)
        ! locals
        integer :: i
        integer :: new_count
        integer :: up_node_index


        new_count = node_list(node_index)%upstream_tree_count + node_list(node_index)%upstream_count

        !% copy the direct links
        do i=1,node_list(node_index)%upstream_count
            up_node_index = node_list(node_index)%upstream_indexes(i)
            call list_add_item_reserve(node_list(target_index)%upstream_tree_indexes, &
                up_node_index, &
                node_list(target_index)%upstream_tree_count, &
                new_count)
        enddo

        do i=1,node_list(node_index)%upstream_count
            up_node_index = node_list(node_index)%upstream_indexes(i)
            call riv_node_full_upstream(target_index, up_node_index, node_list)
        end do

    end subroutine riv_node_full_upstream


    !% uses the upstream_indexes to calculate the full upstream tree
    !% each node will have the full list up nodes that are upstream
    recursive subroutine riv_node_river_delay(node_list, node_index, downstream_delay)
        use dta_utility
        implicit none
        ! recurse node index
        integer, intent(in) :: node_index
        type(riv_tree_node) :: node_list(:)
        ! summed up delay of downstream nodes
        double precision , intent(in) :: downstream_delay
        ! locals
        integer :: i
        integer :: up_node_index

        node_list(node_index)%total_downstream_delay = downstream_delay

        ! sum this delay upstream
        do i=1,node_list(node_index)%upstream_count
            up_node_index = node_list(node_index)%upstream_indexes(i)
            call riv_node_river_delay(node_list, up_node_index, &
                downstream_delay + node_list(node_index)%reach_delay)
        end do

    end subroutine riv_node_river_delay





    recursive subroutine riv_node_recurse_label(node_index, node_list, counter)
        implicit none
        integer, intent(in) :: node_index
        type(riv_tree_node), intent(inout) :: node_list(:)
        integer, intent(inout) :: counter

        !locals
        integer :: i
        integer :: up
        integer :: tmp

        do i=1,node_list(node_index)%upstream_count
            up = node_list(node_index)%upstream_indexes(i)

            if(node_list(up)%node_type == NODE_TYPE_RIVER) then !% label river reaches
                node_list(up)%gauge_id = counter
                if(counter > 0) then
                    counter = counter + 1;
                else ! sea outlet indexes are negative
                    counter = counter - 1;
                endif
            endif
        enddo
        do i=1,node_list(node_index)%upstream_count
            up = node_list(node_index)%upstream_indexes(i)
            if (node_list(up)%node_type == NODE_TYPE_GAUGE) then !% keep gauges, restart numbering with gauge_id
                tmp = node_list(up)%gauge_id + 1
                call riv_node_recurse_label(up, node_list, tmp)
            else
                ! note counter modified by call
                call riv_node_recurse_label(up, node_list, counter)
            endif
        enddo

    end subroutine riv_node_recurse_label

    subroutine riv_node_remove(node_list, remove_index)
        use dta_utility
        implicit none
        type(riv_tree_node), intent(inout) :: node_list(:)
        integer :: remove_index
        !locals
        integer :: new_down_index
        integer :: i
        integer :: check_index

        node_list(remove_index)%enabled = .false.

        !% the removed node might be downstream of many other nodes
        !% rebuild the upstream indexes of the new downstream node
        new_down_index = node_list(remove_index)%downstream_index
        if(node_list(new_down_index)%enabled.eqv..false.) then
            print*,'WARNING update link to removed node'
        endif

        !% the updated upstreams will be the new_downstream nodes
        !% upstreams (with the removed node taken out) plus all the
        !% upstreams of the removed node

        ! remove the index thant points up to the removed node from the downstream node
        !print*,'remove', remove_index
        !print*,'before',node_list(new_down_index)%upstream_count
        !print*, node_list(new_down_index)%upstream_indexes(1:node_list(new_down_index)%upstream_count)
        do i=1,node_list(new_down_index)%upstream_count
            check_index = node_list(new_down_index)%upstream_indexes(i)
            if (check_index == remove_index) then
                !print*, 'update upstream', i, remove_index
                !print*, node_list(new_down_index)%upstream_indexes(1:node_list(new_down_index)%upstream_count)
                call list_remove_at(node_list(new_down_index)%upstream_indexes, &
                    i, &
                    node_list(new_down_index)%upstream_count)
                !print*, node_list(new_down_index)%upstream_indexes(1:node_list(new_down_index)%upstream_count)
            endif
        end do
        !print*,'after',node_list(new_down_index)%upstream_count
        !print*, node_list(new_down_index)%upstream_indexes(1:node_list(new_down_index)%upstream_count)

        ! add all the upstreams from the removed node
        do i=1,node_list(remove_index)%upstream_count
            check_index = node_list(remove_index)%upstream_indexes(i)
            node_list(check_index)%downstream_index = new_down_index

            call list_add_item(node_list(new_down_index)%upstream_indexes, &
                check_index, &
                node_list(new_down_index)%upstream_count)
        enddo
        node_list(remove_index)%upstream_count = 0



    end subroutine riv_node_remove

    function riv_node_get_label(node)
        implicit none
        type(riv_tree_node) :: node
        character(32) :: riv_node_get_label
        integer gauge_id

        !% use absolute value of gauge_id as graphvis doesn't like the
        gauge_id = abs(node%gauge_id)

        if(node%node_type == NODE_TYPE_SEA) then


            !prefix = 's'; % sea
            write (riv_node_get_label,'(A,I0)') 's', gauge_id
        else if (node%node_type == NODE_TYPE_GAUGE) then
            !prefix = 'g'; % gauge
            write (riv_node_get_label,'(A,I0)') 'g', gauge_id
        else if (node%node_type == NODE_TYPE_RIVER) then !% river reach
            if(node%gauge_id > 0) then
                !prefix = 'r'; % river after catchment
                write (riv_node_get_label,'(A,I0)') 'r', gauge_id

            else
                    !prefix = 'e'; % esturary (between sea and first gauge)
                write (riv_node_get_label,'(A,I0)') 'e', gauge_id

            endif
        else

            !prefix = 'u'; % unknown
            write (riv_node_get_label,'(A,I0)') 'u', gauge_id

        endif


            !% minus sign, type is given with the prefix
            !label = sprintf('%s%d', prefix, abs(int32(obj.gauge_id)));


    end function

    subroutine riv_tree_write_graph(filename, node_list)
        implicit none
        character(1024) :: filename
        type(riv_tree_node), intent(inout) :: node_list(:)
        ! locals
        character(1024) :: graph_image_file
        integer :: ii, jj
        character(32) :: label
        character(32) :: label_next
        integer :: upstream
        character(1024) :: command

        graph_image_file = trim(filename) // '.png'
        !% create a tree network file that can be graphed with graphvis 'dot'

        open (10, file = filename, status = 'unknown')

        write (10,'(A)') 'digraph G {'
        write (10,'(A)') 'rankdir="LR";'

        do ii=1,size(node_list)

            if(node_list(ii)%enabled.eqv..false.) then
                cycle
            endif

            if(node_list(ii)%node_type == NODE_TYPE_GAUGE) then
                label = riv_node_get_label(node_list(ii))
                write (10,'(A,A)') trim(label),'[color="red"];'

            elseif(node_list(ii)%node_type == NODE_TYPE_SEA) then
                if(node_list(ii)%upstream_count > 0) then
                    label = riv_node_get_label(node_list(ii))
                    write (10,'(A,A)') trim(label),'[color="blue"];'
                endif
            else
                ! default node color
            endif
        end do

        do ii=1,size(node_list)

            if(node_list(ii)%enabled.eqv..false.) then
                cycle
            endif

            do jj=1,node_list(ii)%upstream_count
                upstream = node_list(ii)%upstream_indexes(jj);


                label = riv_node_get_label(node_list(ii))
                label_next = riv_node_get_label(node_list(upstream))

                write (10,'(A,'' -> '',A,'' ;'')') trim(label), &
                    trim(label_next)

                !fprintf(fid, '%s -> %s ;\n', ...
                !    current.GetLabel(), upstream.GetLabel());
            end do
        end do
        write (10,'(A)') '}'

        close(10)

        !% graphviz-2.38\release\bin\dot -Tpng -Gsize=9,15\! -Gdpi=1200 -oriv.png riv_graph.dot

        write (command, '(A,A,A,A,A,A)') &
            'graphviz-2.38\release\bin\dot.exe', &
            ' -Tpng ', &
            ' -o', trim(graph_image_file), &
            ' ', trim(filename)

        print *, trim(command)
            !system(command, '-echo');


    end subroutine riv_tree_write_graph


    subroutine riv_tree_write(node_list, nrows, xllcorner, yllcorner, cellsize, file_prefix)
        use dta_utility
        implicit none
        type(riv_tree_node) :: node_list(:)
        integer :: nrows
        double precision :: xllcorner, yllcorner, cellsize
        character(1024) :: file_prefix

        ! locals
        character(1024) :: flow_conn_file
        character(1024) :: flow_point_file
        integer :: i
        integer :: down
        double precision :: curr_northing, curr_easting
        double precision :: down_northing, down_easting

        integer :: flow_fid, point_fid
        flow_fid = 100
        point_fid = 101

        flow_conn_file = trim(file_prefix)//'_flow_conn.txt'
        flow_point_file = trim(file_prefix)//'_flow_point.txt'

        print *, 'write ', trim(flow_conn_file)
        print *, 'write ', trim(flow_point_file)

        open (flow_fid, file = flow_conn_file, status = 'unknown')
        open (point_fid, file = flow_point_file, status = 'unknown')

        write (flow_fid, '(19A)') 'Node_ID', tab, &
            'x_easting', tab, &
            'y_northing', tab, &
            'Type', tab, &
            'ID_down', tab, &
            'x_easting2', tab, &
            'y_northing2', tab, &
            'delta_dist', tab, &
            'delta_h', tab, &
            'slope'

        write (point_fid, '(11A)') 'Node_ID', tab, &
            'x_easting', tab, &
            'y_northing', tab, &
            'Type', tab, &
            'dist', tab, &
            'h_elevation'


        do i=1,size(node_list)
            call RowColToNorthingEasting( node_list(i)%row, node_list(i)%col, &
                nrows, xllcorner, yllcorner, cellsize, &
                curr_northing, curr_easting, .true.)

            if (node_list(i)%enabled) then
                write(point_fid, '(I0,A,F0.1,A,F0.1,A,I0,A,F0.3,A,F0.3)') &
                    node_list(i)%gauge_id, tab, &
                    curr_easting,  tab, &
                    curr_northing, tab, &
                    node_list(i)%node_type, tab, &
                    node_list(i)%point_dist, tab, &
                    node_list(i)%point_h

                if (node_list(i)%downstream_index > 0) then

                    down = node_list(i)%downstream_index
                    if(node_list(down)%enabled.eqv..false.) then
                        print *, 'link to removed node still exists'
                    endif

                    call RowColToNorthingEasting( node_list(down)%row, node_list(down)%col, &
                        nrows, xllcorner, yllcorner, cellsize, &
                        down_northing, down_easting, .true. )

                    write(flow_fid, '(I0,A,F0.1,A,F0.1,A,I0,A,I0,A,F0.1,A,F0.1,A,F0.3,A,F0.3,A,F0.3)') &
                        node_list(i)%gauge_id, tab, &
                        curr_easting,  tab, &
                        curr_northing, tab, &
                        node_list(i)%node_type, tab, &
                        node_list(down)%gauge_id, tab, &
                        down_easting, tab, &
                        down_northing, tab, &
                        node_list(i)%downstream_dist, tab, &
                        node_list(i)%downstream_delta_h, tab, &
                        node_list(i)%downstream_slope
                endif
            endif
        end do


    end subroutine riv_tree_write

    subroutine riv_tree_read(node_list, &
        nrows, xllcorner, yllcorner, cellsize, &
        flow_conn_file, flow_point_file)
        use dta_utility
        implicit none
        type(riv_tree_node), allocatable :: node_list(:)
        integer :: nrows
        double precision :: xllcorner, yllcorner, cellsize
        character(1024) :: flow_conn_file
        character(1024) :: flow_point_file

        ! locals
        double precision, allocatable, dimension(:,:) :: flow_conn_data
        double precision, allocatable, dimension(:,:) :: flow_point_data

        call read_numeric_list(flow_point_file, 6, 1, flow_point_data)
        call read_numeric_list(flow_conn_file, 10, 1, flow_conn_data)

        call riv_tree_read_impl(node_list, &
            nrows, xllcorner, yllcorner, cellsize, &
            flow_conn_data, flow_point_data)

    end subroutine

    subroutine riv_tree_read_impl(node_list, &
        nrows, xllcorner, yllcorner, cellsize, &
        flow_conn_data, flow_point_data)
        use dta_utility
        implicit none
        type(riv_tree_node), allocatable :: node_list(:)
        integer :: nrows
        double precision :: xllcorner, yllcorner, cellsize
        double precision, allocatable, dimension(:,:) :: flow_conn_data
        double precision, allocatable, dimension(:,:) :: flow_point_data

        ! locals
        integer :: i, j
        integer :: row, col
        integer :: node_id, down_id
        integer :: node_index

        allocate(node_list(size(flow_point_data,1)))

        do i=1,size(node_list)

            !flow_point_data(2)!; % easting
            !flow_point_data(3)!; % northing
            if(nrows>0) then
                call NorthingEastingToRowCol(flow_point_data(i,3), flow_point_data(i,2), &
                    nrows, xllcorner, yllcorner, cellsize, row, col)
            else
                !if nrows is not set - just set the node row/col to the easting/northing
                row = nint(flow_point_data(i,3))
                col = nint(flow_point_data(i,2))
            endif

            node_list(i) = new_riv_tree_node()
            node_list(i)%gauge_id = nint(flow_point_data(i,1))
            node_list(i)%col = col
            node_list(i)%row = row
            node_list(i)%node_type = nint(flow_point_data(i,4))
            node_list(i)%point_dist = flow_point_data(i,5)
            node_list(i)%point_h = flow_point_data(i,6)
        end do

        do i=1,size(flow_conn_data,1)

            node_id = nint(flow_conn_data(i,1))
            down_id = nint(flow_conn_data(i,5))

            !% find matching node from the node_list
            node_index = 0
            do j=1,size(node_list)
                if(node_list(j)%gauge_id == node_id) then
                    node_index = j
                    exit
                endif
            enddo
            if(node_index == 0) then
                print *, 'Point file and flow file do not match', node_id
                stop
            endif

            !% find the matching downstream node from node_list
            do j=1,size(node_list)
                if(node_list(j)%gauge_id == down_id) then

                    !print*,'down',j, node_list(node_index)%gauge_id,'->',node_list(j)%gauge_id
                    node_list(node_index)%downstream_index = j
                    exit
                end if
            end do


            node_list(node_index)%downstream_dist = flow_conn_data(i,8)
            node_list(node_index)%downstream_delta_h = flow_conn_data(i,9)
            node_list(node_index)%downstream_slope = flow_conn_data(i,10)

        end do

        call build_upstream_from_downstream(node_list)

        !% build the full upstream linkages for each node
        !% uses the upstream_indexes to calculate the full upstream tree for each node

        ! every node with have a list of every upstream node
        !% these lists will determine which rows in the routing table to read from to sum up the flow at each point
        do i = 1,size(node_list)
            if (allocated(node_list(i)%upstream_tree_indexes)) then
                deallocate(node_list(i)%upstream_tree_indexes)
            endif
            node_list(i)%upstream_tree_count = 0

            call riv_node_full_upstream(i, i, node_list)
        enddo

    end subroutine riv_tree_read_impl

    subroutine build_upstream_from_downstream(node_list)
        use dta_utility
        implicit none
        type(riv_tree_node) :: node_list(:)

        integer :: i
        integer :: down_index

        do i=1,size(node_list)
            if (allocated(node_list(i)%upstream_indexes)) then
                deallocate(node_list(i)%upstream_indexes)
            endif
            node_list(i)%upstream_count = 0
        end do

        !% re-build the upstream links from the downstream info
        do i=1,size(node_list)
            down_index = node_list(i)%downstream_index
            if (down_index > 0) then
                ! add the upstream link to this node
                ! on the downstream node
                call list_add_item(node_list(down_index)%upstream_indexes, i, &
                    node_list(down_index)%upstream_count)

                !print*, node_list(i)%gauge_id,'<-',node_list(down_index)%gauge_id
            endif
        enddo
    end subroutine build_upstream_from_downstream

    subroutine riv_tree_read_dyna(node_list, &
        nrows, xllcorner, yllcorner, cellsize)

        use dta_utility
        implicit none
        type(riv_tree_node), allocatable :: node_list(:)
        integer :: nrows
        double precision :: xllcorner, yllcorner, cellsize

        ! locals
        double precision, allocatable, dimension(:,:) :: flow_conn_data
        double precision, allocatable, dimension(:,:) :: flow_point_data

        call read_numeric_list_fid(503, 6, 1, flow_point_data)
        call read_numeric_list_fid(501, 10, 1, flow_conn_data)

        call riv_tree_read_impl(node_list, &
            nrows, xllcorner, yllcorner, cellsize, &
            flow_conn_data, flow_point_data)

    end subroutine riv_tree_read_dyna


end module dta_riv_tree_node
