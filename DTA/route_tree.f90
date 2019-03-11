!% Toby Dunne
!% Mar 2016

program route_tree
    use dta_utility
    use dta_route_tree
    use dta_riv_tree_node
    use dta_rivers
    use dta_catch_cut
    implicit none

    character(1024) in_dem_file
    character(1024) in_point_file
    character(1024) in_riv_file

    character(1024) tmp_char

    character(1024) dem_output_prefix
    character(1024) riv_output_prefix
    character(1024) riv_output_prefix2

    type(time_type) :: start_time, end_time, run_start_time

    CHARACTER(len=1024) :: arg

    double precision, allocatable, dimension(:,:) :: riv_grid
    double precision, allocatable, dimension(:,:) :: dem_grid
    ! point list from file
    double precision, allocatable, dimension(:,:) :: read_point_row_list
    ! point enabled/disabled flag for filtering
    logical, allocatable, dimension(:) :: read_point_row_enabled
    type(point_type), allocatable, dimension(:) :: read_point_list

    ! point list after filtering
    double precision, allocatable, dimension(:,:) :: point_row_list
    type(point_type), allocatable, dimension(:) :: point_list


    logical, allocatable, dimension(:,:) :: riv_mask_grid
    logical, allocatable, dimension(:,:)  :: sea_mask_grid
    double precision, allocatable, dimension(:,:) :: riv_dist_grid

    integer, allocatable, dimension(:,:) :: riv_label_grid

    integer :: point_col_id
    integer :: point_col_east
    integer :: point_col_north
    integer :: point_col_node_type
    integer :: point_col_ref_area
    integer :: point_col_area
    integer :: point_col_score
    integer :: point_col_dist
    integer :: point_col_out_bound

    integer i, j, ioerr
    integer point_count
    integer :: fp

    integer :: ncols, nrows, riv_ncols, riv_nrows
    double precision :: xllcorner, yllcorner, cellsize, riv_xllcorner, riv_yllcorner, riv_cellsize
    double precision :: double_nodata
    double precision :: filter_area
    double precision :: filter_score
    double precision :: riv_reach_increment_m
    double precision :: riv_reach_increment_m_min

    type(point_type), allocatable :: sea_outlets(:)
    type(riv_tree_node), allocatable :: node_list(:)

    logical :: input_is_valid

    integer :: filter_stat_na
    integer :: filter_valid ! valid but might be filtered (not na)
    integer :: filter_stat_area
    integer :: filter_stat_score
    integer :: filter_stat_area_and_score

    integer :: filter_stat_duplicate

    integer(1), allocatable :: flow_dir_grid (:,:)
    double precision, allocatable :: slope_grid (:,:)


    CALL timer_get(run_start_time)

    !call utility_test()
    !stop

    filter_score = 20
    filter_area = 1
    riv_reach_increment_m = 0
    ! default to half of riv_reach_increment_m
    ! when searching for reach lengths distance is
    ! measured from the downstream node and works up
    ! if a reach ends up being too close to the upstream gauge it will be removed
    riv_reach_increment_m_min = 0

    in_dem_file = ''
    in_riv_file = ''
    in_point_file = ''

    tmp_char = ''

    input_is_valid = .true.

    tmp_char = 'DTA_route_tree.log'
    fp=fp
    open(fp, file = tmp_char, status="unknown", action="write", iostat=ioerr)

    if(ioerr/=0) then
        print*,'error opening output file: ', trim(tmp_char)
        print*,'ensure the directory exists and correct write permissions are set'
        stop
    endif

    write(fp,*) '--- route_tree.f90 ---'
    write(fp,*) ''
    print *, '--- Starting route_tree ---'

    i=0

    do
        CALL get_command_argument(i, arg)
        if(len_trim(arg) == 0) exit
        if (are_equal(arg, '-dem')) then
            CALL get_command_argument(i+1, in_dem_file)
        elseif (are_equal(arg, '-river')) then
            CALL get_command_argument(i+1, in_riv_file)
        elseif (are_equal(arg, '-points')) then
            CALL get_command_argument(i+1, in_point_file)
        elseif (are_equal(arg, '-filter_area')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) filter_area
        elseif (are_equal(arg, '-filter_score')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) filter_score
        elseif (are_equal(arg, '-reach_length')) then
            CALL get_command_argument(i+1, arg)
            read (arg,*) riv_reach_increment_m
            riv_reach_increment_m_min = riv_reach_increment_m / 2
        endif
        i = i + 1
    enddo

    if(check_file_arg(in_dem_file,'-dem').eqv..false.) then
        input_is_valid = .false.
    endif
    if(check_file_arg(in_riv_file,'-river').eqv..false.) then
        input_is_valid = .false.
    endif
    !if(check_file_arg(in_point_file,'-points').eqv..false.) then
    !    input_is_valid = .false.
    !endif

    if(input_is_valid .eqv. .false.) then
        print *, 'command options'
        print *, '-dem <file.asc>   select dem ascii grid file'
        print *, '-river <file.asc>   select river dist ascii grid file'
        print *, '-points <file.txt>  points on the river mask ascii grid file'
        print *, '      (e.g. station_river_gauge.txt from catch_cut.e)'
        print *, '      if points not specified, only sea outlets will be used'
        print *, ''
        print *, 'optional filter points (can be combined)'
        print *, '-filter_area  <min>  filter points by area (default 1, 0 disable)'
        print *, '-filter_score <min>  filter points by score (default 50, 0 disable)'
        print *, ''
        print *, 'optional automatic reach length between gauges'
        print *, '-reach_length <dist> distance in metres (or units of cellsize)'
        stop
    endif

    write(fp,*) 'Command Options Read In'
    write(fp,*) 'dem: ', trim(in_dem_file)
    write(fp,*) 'river file: ', trim(in_riv_file)
    write(fp,*) 'points on river:   ', trim(in_point_file)
    write(fp,*) 'filter area:   ', filter_area
    write(fp,*) 'filter score:  ', filter_score
    write(fp,*) 'reach length:   ', riv_reach_increment_m
    write(fp,*) ''

    print *, 'dem ', trim(in_dem_file)
    print *, 'river file ', trim(in_riv_file)
    print *, 'points on river', trim(in_point_file)
    print *, 'filter score', filter_score
    print *, 'filter area', filter_area
    print *, 'reach length', riv_reach_increment_m

    if(riv_reach_increment_m > 0.01d0)then
        write (dem_output_prefix, '(A,A,I0)') in_dem_file(1:len_trim(in_dem_file)-4), '_rl',int(riv_reach_increment_m)
        write (riv_output_prefix, '(A,A,I0)') in_riv_file(1:len_trim(in_riv_file)-4), '_rl',int(riv_reach_increment_m)
    else
        ! output based on dem filename
        dem_output_prefix = in_dem_file(1:len_trim(in_dem_file)-4)
        ! output based on riv filename
        riv_output_prefix = in_riv_file(1:len_trim(in_riv_file)-4)
    endif
    riv_output_prefix2 = in_riv_file(1:len_trim(in_riv_file)-4)

    print *, 'read river grid:', trim(in_riv_file)
    write(fp,*) 'read river grid:', trim(in_riv_file)
    CALL timer_get(start_time)
    call read_ascii_grid(in_riv_file, riv_grid, &
        riv_ncols, riv_nrows, riv_xllcorner, riv_yllcorner, riv_cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read river grid', start_time, end_time)

    allocate(riv_mask_grid(riv_nrows, riv_ncols))
    riv_mask_grid(:,:) = .false.
    where (riv_grid>0.001) riv_mask_grid = .true.
    deallocate(riv_grid)



    point_col_id = 1
    point_col_east = 2
    point_col_north = 3
    point_col_node_type = 4
    point_col_ref_area = 5
    point_col_area = 6
    point_col_score = 7
    point_col_dist = 8
    point_col_out_bound = 9
    if(len_trim(in_point_file) > 0) then
        print *, 'read point list: ', trim(in_point_file)
        write(fp,*) 'read point list: ', trim(in_point_file)
        call read_numeric_list(in_point_file, 8, 1, read_point_row_list)
    else
        ! no gauge points, generate sea outlet only
        allocate(read_point_row_list(0,8))
    endif

    filter_stat_area = 0
    filter_stat_score = 0
    filter_stat_area_and_score = 0
    filter_stat_duplicate = 0
    filter_stat_na = 0

    ! Filter Points
    !% * remove all with low match score
    !% * remove tiny catchments

    point_count = size(read_point_row_list,1)
    allocate(read_point_row_enabled(point_count))
    allocate(read_point_list(point_count))
    read_point_row_enabled(:) = .true.

    do i = 1,size(read_point_row_list,1)

        call NorthingEastingToRowCol(&
            read_point_row_list(i,point_col_north), read_point_row_list(i,point_col_east), &
            riv_nrows, riv_xllcorner, riv_yllcorner, riv_cellsize,&
            read_point_list(i)%y, read_point_list(i)%x)

        if(nint(read_point_row_list(i, point_col_node_type)) /= 2) then
            write(fp,*) 'custom point: ', (read_point_row_list(i, point_col_id)), (read_point_row_list(i, point_col_node_type))
        else
            if(read_point_row_list(i, point_col_area) < 0.1) then
                read_point_row_enabled(i) = .false.
                filter_stat_na = filter_stat_na + 1
            elseif((filter_score > 0.001 .and. read_point_row_list(i, point_col_score) > filter_score) &
                .or. (read_point_row_list(i, point_col_area) < filter_area)) then
                read_point_row_enabled(i) = .false.

                ! filter on both conditions
                if((filter_score > 0.001 .and. read_point_row_list(i, point_col_score) > filter_score) &
                    .and. (read_point_row_list(i, point_col_area) < filter_area)) then
                    filter_stat_area_and_score = filter_stat_area_and_score + 1

                    write(fp,*) 'point removed area/score', &
                        nint(read_point_row_list(i,point_col_id)), &
                        read_point_row_list(i,point_col_ref_area), &
                        read_point_row_list(i,point_col_score)

                elseif  (filter_score > 0.001 .and. read_point_row_list(i, point_col_score) > filter_score) then ! filter on just score
                    filter_stat_score = filter_stat_score + 1

                    write(fp,*) 'point removed score', &
                        nint(read_point_row_list(i,point_col_id)), &
                        read_point_row_list(i,point_col_ref_area), &
                        read_point_row_list(i,point_col_score)

                else ! filter on just area
                    filter_stat_area = filter_stat_area + 1

                    write(fp,*) 'point removed area', &
                        nint(read_point_row_list(i,point_col_id)), &
                        read_point_row_list(i,point_col_ref_area), &
                        read_point_row_list(i,point_col_score)

                endif
            endif
        endif
    end do

    !% * remove any catchments that are on duplicate xy locations (both removed)
    do i = 1,size(read_point_row_list,1)
        if(read_point_row_enabled(i)) then
            ! point area and score is ok - check if duplicated
            do j=1,size(read_point_row_list,1)
                if ((i==j) .or. (read_point_row_enabled(j).eqv..false.)) then
                    cycle
                endif
                ! make sure there isn't another catchment on the same point
                ! remove all duplicates
                ! could do something more clever by keeping the one with the best score
                if( (read_point_list(i)%x == read_point_list(j)%x) &
                    .and. (read_point_list(i)%y == read_point_list(j)%y)) then

                    ! i may be already disabled in the case that there are more than 2 duplicates
                    if(read_point_row_enabled(i)) then
                        read_point_row_enabled(i) = .false.
                        filter_stat_duplicate = filter_stat_duplicate + 1
                        write(fp,*) 'point removed duplicate location', nint(read_point_row_list(i,point_col_id))
                    endif

                    read_point_row_enabled(j) = .false.
                    filter_stat_duplicate = filter_stat_duplicate + 1
                    write(fp,*) '      removed duplicate location', nint(read_point_row_list(j,point_col_id))

                endif
            end do
        endif
    end do

    point_count = count(read_point_row_enabled)


    print *, 'Filters catchments'
    print *, 'Initial Points        : ', size(read_point_row_enabled)
    filter_valid = size(read_point_row_enabled) - filter_stat_na
    print *, 'Valid Points         : ', filter_valid
    print *, 'Remaining Points     : ', point_count

    write(fp,*) ''
    write(fp,*) 'Filters catchments'
    write(fp,*) 'Removed Area          : ', filter_stat_area
    write(fp,*) 'Removed Score         : ', filter_stat_score
    write(fp,*) 'Removed Area and Score: ', filter_stat_area_and_score
    write(fp,*) 'Removed Duplicate     : ', filter_stat_duplicate
    write(fp,*) 'Removed NA            : ', filter_stat_na
    write(fp,*) 'Initial Points        : ', size(read_point_row_enabled)
    filter_valid = size(read_point_row_enabled) - filter_stat_na
    write(fp,*) 'Valid Points         : ', filter_valid
    write(fp,*) 'Remaining Points     : ', point_count
    write(fp,*) '% of valid removed   : ',100- real(point_count)/real(filter_valid)*100

    if(point_count == 0) then
        print *, 'No remaining points, try removing filters'
    !    stop
    endif

    allocate(point_row_list(point_count,size(read_point_row_list,2)))
    allocate(point_list(point_count))

    ! copy all the enabled points
    j = 0
    do i = 1,size(read_point_row_list,1)
        if (read_point_row_enabled(i)) then
            j = j+1
            point_row_list(j,:) = read_point_row_list(i,:)
            point_list(j) = read_point_list(i)
        endif
    end do
    deallocate(read_point_row_list)
    deallocate(read_point_row_enabled)
    deallocate(read_point_list)


    print *, 'read dem grid:', trim(in_dem_file)
    write(fp,*)''
    write(fp,*) 'read dem grid', trim(in_dem_file)
    CALL timer_get(start_time)
    call read_ascii_grid(in_dem_file, dem_grid, &
        ncols, nrows, xllcorner, yllcorner, cellsize, double_nodata)
    CALL timer_get(end_time)
    !call timer_print('read dem grid', start_time, end_time)

    if (nrows /= riv_nrows .or. ncols /= riv_ncols) then
        print *, 'river file does not match dem'
        stop
    endif

!    if (len_trim(sea_outlet_input_file) > 0) then
!        call read_numeric_list(sea_outlet_input_file, 2, 1, read_manual_sea_outlets)
!
!        print *, 'outlets from file', size(read_manual_sea_outlets,1)
!        allocate(sea_outlets(size(read_manual_sea_outlets,1)))
!
!        do i=1,size(read_manual_sea_outlets,1)
!            !                            northing                      easting
!            call NorthingEastingToRowCol(read_manual_sea_outlets(i,2), read_manual_sea_outlets(i,1), &
!                nrows, xllcorner, yllcorner, cellsize,&
!                sea_outlets(i)%y, sea_outlets(i)%x)
!        end do
!
!        deallocate(read_manual_sea_outlets)
!    else
        write(fp,*) 'find outlets joining sea'
        allocate(sea_mask_grid(nrows,ncols))
        sea_mask_grid = .false.
        where(dem_grid<-90) sea_mask_grid = .true.

        call river_find_outlets(nrows, ncols, &
            riv_mask_grid, sea_mask_grid, sea_outlets)

        deallocate(sea_mask_grid)
        !tmp_char = in_dem_file(1:len_trim(in_dem_file)-4)//'_outlet.txt'
        !call write_point_list(tmp_char, sea_outlets, nrows, xllcorner, yllcorner, cellsize)
!    endif

    allocate (flow_dir_grid(nrows,ncols))
    allocate (slope_grid(nrows,ncols))

    call calc_flow_direction( nrows, ncols, dem_grid, cellsize, flow_dir_grid, slope_grid)
    deallocate(slope_grid)

    allocate(riv_dist_grid(nrows,ncols))

    allocate(riv_label_grid(nrows,ncols))

    riv_dist_grid(:,:) = 0

    call river_find_dist(ncols, nrows, &
        riv_mask_grid, &
        flow_dir_grid, &
        riv_dist_grid, &
        .true., sea_outlets)

    ! calculate distance for all cells that are not linked to sea_outlets
    ! updated riv_dist_grid where distance hasn't been done yet
    call river_find_dist(ncols, nrows, &
        riv_mask_grid, &
        flow_dir_grid, &
        riv_dist_grid, &
        .false., point_list)

    riv_dist_grid = riv_dist_grid * cellsize

    call routing_dta(nrows, ncols, point_row_list, point_list, sea_outlets, dem_grid, &
        flow_dir_grid, &
        riv_mask_grid, &
        riv_dist_grid, &
        riv_reach_increment_m, riv_reach_increment_m_min, &
        point_col_id, &
        node_list, riv_label_grid, fp)

    deallocate(riv_mask_grid)

    tmp_char = trim(riv_output_prefix)//'_dist.asc'

    print*,'Write: ', trim(tmp_char)
    write(fp,*) ''
    write(fp,*) 'Write: ', trim(tmp_char)
    call write_ascii_grid(tmp_char, riv_dist_grid, &
        ncols, nrows, &
        xllcorner, yllcorner, cellsize, 0.0d0, 3)


    !    tmp_char = trim(riv_output_prefix2)//'_dist.asc'
    !
    !    print *, 'Write river dist grid: ', trim(tmp_char)
    !    call write_ascii_grid(tmp_char, riv_dist_grid, ncols, nrows, xllcorner, yllcorner, cellsize, 0.0d0, 1)

    tmp_char = trim(riv_output_prefix)//'_id.asc'
    print *, 'Write river label grid: ', trim(tmp_char)
    write(fp,*) 'Write river label grid: ', trim(tmp_char)
    call write_ascii_grid_int(tmp_char, riv_label_grid, ncols, nrows, xllcorner, yllcorner, cellsize, 0)

    !tmp_char = trim(dem_output_prefix)//'_tree.gdot'
    !print *, 'Write graph: ', trim(tmp_char)
    !call riv_tree_write_graph(tmp_char, node_list)

    call riv_tree_write(node_list, nrows, xllcorner, yllcorner, cellsize, dem_output_prefix)

    CALL timer_get(end_time)
    !call timer_print('route_tree', run_start_time, end_time)
        write(fp,*) ''
    write(fp,*) 'Successfully finished route_tree.f90'
    print *, '--- Finished route_tree ---'
    close(fp)
    stop

end program route_tree





