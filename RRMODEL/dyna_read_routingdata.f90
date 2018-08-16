module dyna_read_routingdata
    implicit none
contains

    !
    !==========================================================
    !  SUBROUTINE FOR READING IN THE FLOW ROUTING DATA
    !==========================================================
    !
    subroutine read_routingdata(cat_route_vmode, &
        dyna_hru, &
        nac, &
        node_to_flow_mapping, &
        num_rivers, &
        route_riv, &
        route_tdh, &
        routing_mode)

        use dta_route_processing
        use dta_riv_tree_node
        use dyna_common_types

        implicit none

        ! Argument Declares

        integer :: cat_route_vmode, num_rivers
        integer :: routing_mode, nac
        integer, dimension(:), allocatable :: node_to_flow_mapping
        type(route_river_info_type) :: route_riv
        type(route_time_delay_hist_type) :: route_tdh
        type(dyna_hru_type) :: dyna_hru(nac)

        ! Local Declares
        integer :: i, k, z

        ! End Declares

        ! Toby Dunne and Gemma Coxon

        ! default to the original routing mode
        routing_mode = ROUTE_MODE_NON_ROUTED
        route_riv%is_enabled = .false.

        ! Read in the routing files

        call route_processing_read_info_dyna(route_riv)

        route_riv%node_list%catch_area = 0
        route_riv%node_list%frac_sub_area = 0

        do i = 1, size(route_riv%node_list)

            ! Calculate the catchment areas
            route_riv%node_list(i)%catch_area = route_riv%node_list(i)%catch_area + &
                sum(route_riv%river_data(:, 2), int(route_riv%river_data(:, 1)) &
                == route_riv%node_list(i)%gauge_id)

            do z = 1, size(dyna_hru%ipriv, 1)

                if (dyna_hru(z)%ipriv == i) then

                    route_riv%node_list(i)%frac_sub_area = dyna_hru(z)%ac &
                        + route_riv%node_list(i)%frac_sub_area

                end if

            end do

            do k = 1, route_riv%node_list(i)%upstream_tree_count
                ! Calculate the fractional areas of each sub-catchment

                do z = 1, size(dyna_hru%ipriv, 1)

                    if (route_riv%node_list(i)%upstream_tree_indexes(k) &
                        == dyna_hru(z)%ipriv) then

                        route_riv%node_list(i)%frac_sub_area = dyna_hru(z)%ac + &
                            route_riv%node_list(i)%frac_sub_area

                    end if

                end do

                route_riv%node_list(i)%catch_area = route_riv%node_list(i)%catch_area + &
                    sum(route_riv%river_data(:, 2), int(route_riv%river_data(:, 1)) &
                    == route_riv%node_list(route_riv%node_list(i)%upstream_tree_indexes(k))%gauge_id)
            end do
        end do

        ! Set the velocity mode for the time delay histogram - only one option at the moment so set to 1
        cat_route_vmode = 1
        route_tdh%v_mode = cat_route_vmode

        ! Allocate the velocity parameter
        allocate (route_tdh%v_param(1))
        allocate(node_to_flow_mapping(num_rivers))

        ! mapping required as there is no input q for the sea river reach
        do i = 1, num_rivers
            node_to_flow_mapping(i) = i
        end do

    end subroutine

end module dyna_read_routingdata
