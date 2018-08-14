!% Toby Dunne
!% Mar 2016
module dta_catch_mask
    implicit none

contains

    recursive subroutine add_catch_mask(nrows, ncols, mask_dirs, &
        node_list, this_node, &
        mask_grid, xllcorner, yllcorner, cellsize, fp)
        use dta_riv_tree_node
        use dta_utility
        implicit none
        integer :: nrows
        integer :: ncols
        character(1024), intent(in) :: mask_dirs(:)
        type(riv_tree_node), intent(inout) :: node_list(:)
        integer, intent(in)::this_node, fp
        integer, intent(inout) :: mask_grid(nrows,ncols)
        double precision, intent(in) :: xllcorner, yllcorner, cellsize

        ! locals
        integer::i,j
        integer :: gid
        character(1024)::mask_filename
        logical::found
        double precision, allocatable:: this_mask(:,:)
        integer, allocatable:: section(:,:)
        integer :: this_ncols, this_nrows
        double precision :: this_xllcorner, this_yllcorner, this_cellsize, this_nodata
        integer :: up, up_index
        integer :: b_start_y, b_start_x, b_end_y, b_end_x
        integer :: a_start_y, a_start_x, a_end_y, a_end_x
        character(16) :: prefixes(4)

        prefixes(1) = 'gauge_'
        prefixes(2) = 'reach_'
        prefixes(3) = 'sea_'
        prefixes(4) = 'other_'


        gid = node_list(this_node)%gauge_id
        ! catchment files of node type 2 (gauges) don't have the 000 suffix
        ! - new mode masks all have 000 suffix
!        if(node_list(this_node)%node_type == 2) then
!            gid = gid / 1000
!        endif

        found = .false.

        do i=1,size(mask_dirs)
            do j=1,size(prefixes)
                write(mask_filename,'(A,A,A,I0,A)') trim(mask_dirs(i)), '/', trim(prefixes(j)), gid, '.asc'
                !print*,trim(mask_filename)
                if(file_exists(mask_filename)) then
                    found = .true.
                    exit
                endif
            end do
            if(found) then
                exit
            endif
        end do

        if(found.eqv..false.) then
            write(fp,*) 'WARNING missing catchment mask file: ', &
                node_list(this_node)%gauge_id, 'type: ',node_list(this_node)%node_type
        else
            write(fp,*) 'USE:', trim(mask_filename)

            call read_ascii_grid(mask_filename, this_mask, &
                this_ncols, this_nrows, this_xllcorner, this_yllcorner, this_cellsize, this_nodata)

            call match_coordinates(cellsize, &
                this_ncols, this_nrows, this_xllcorner, this_yllcorner, &
                ncols, nrows, xllcorner, yllcorner, &
                b_start_y, b_start_x, b_end_y, b_end_x, &
                a_start_y, a_start_x, a_end_y, a_end_x  )

            allocate(section(this_nrows,this_ncols))

            ! take section from the full mask
            section = mask_grid(b_start_y:b_end_y, b_start_x:b_end_x)
            ! update the section where this_mask has a value
            where (this_mask > 0) section = int(node_list(this_node)%gauge_id)

            ! write back to full mask
            mask_grid(b_start_y:b_end_y, b_start_x:b_end_x) = section

            deallocate(this_mask)
            deallocate(section)

            !imagesc(mask);
            !drawnow();
            !%pause(1);



        endif
        do up=1,node_list(this_node)%upstream_count

            up_index = node_list(this_node)%upstream_indexes(up)

            !print*,node_list(this_node)%gauge_id,'->',node_list(up_index)%gauge_id

            call add_catch_mask(nrows, ncols, mask_dirs, &
                node_list, up_index, &
                mask_grid, xllcorner, yllcorner, cellsize, fp)
        end do

    end subroutine

end module dta_catch_mask
