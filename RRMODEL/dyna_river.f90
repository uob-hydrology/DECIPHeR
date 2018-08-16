module dyna_river
contains
    !
    !===============================================================
    !  ROUTINE TO CALCULATE THE QB CALCULATIONS AND ROUTE THE FLOW
    !  TO THE RIVERS
    !===============================================================
    !
    subroutine river (dyna_riv, &
        nstep, &
        num_rivers, &
        it, &
        node_to_flow_mapping, &
        q, &
        route_riv, &
        route_tdh)

        use dyna_common_types
        use dta_route_processing
        implicit none
        ! Argument Declares

        integer :: nstep
        integer :: num_rivers
        integer :: it
        double precision :: q(num_rivers, nstep)

        ! GC added for routing 14th June
        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        integer :: node_to_flow_mapping(:)

        ! Local Declares
        integer :: iriv

        ! End declares

        !  SUM UP THE TOTAL DISCHARGES FOR THE OUTFALL
        do iriv = 1, num_rivers
            dyna_riv(iriv)%qout = dyna_riv(iriv)%qb + dyna_riv(iriv)%qof
            dyna_riv(iriv)%sumqout = dyna_riv(iriv)%sumqout + dyna_riv(iriv)%qout
            dyna_riv(iriv)%sumqof = dyna_riv(iriv)%sumqof + dyna_riv(iriv)%qof
        end do

        ! Use Toby's new code to route the flow
        call route_process_run_step(route_riv, route_tdh, &
            node_to_flow_mapping, &
            dyna_riv%qout, &
            q(:, it))

        ! Convert the discharge to m/ts
        ! Divide by the fractional area of the catchment

        do iriv = 1, num_rivers
            q(iriv, it) = q(iriv, it) / route_riv%node_list(node_to_flow_mapping(iriv))%frac_sub_area
        end do

    end subroutine

end module dyna_river
