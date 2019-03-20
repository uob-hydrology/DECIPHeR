module dyna_topmod
contains
    !
    !=================================================================
    !  DynaTOP 1-2 (more modular structure for development)
    !=================================================================
    !  Keith Beven, Lausanne, 1997: Jim Freer Lancaster 12/1/00
    !  Updated by Gemma Coxon and Toby Dunne - update to Fortran 2003, made more modular
    !  Works with new routing routines and routes flows through a regional river netowrk
    !

    subroutine Topmod (nac, &
        nstep, &
        num_rivers, &
        acc, &
        dt, &
        dtt, &
        dyna_hru, &
        dyna_riv, &
        node_to_flow_mapping, &
        ntt, &
        pe_step, &
        q, &
        r_gau_step, &
        rivers, &
        route_riv, &
        route_tdh, &
        smax, &
        srmax, &
        szm, &
        t0dt, &
        td, &
        wt)

        use dyna_common_types
        use dyna_rain
        use dyna_dynamic_dist
        use dyna_root_zone1
        use dyna_unsat_zone1
        use dyna_results_ac
        use dyna_river
        use dyna_results_ts
        use dta_route_processing
        use dyna_satzone_analyticalsolver

        implicit none

        ! Argument Declares
        integer :: nac
        integer :: nstep
        integer :: num_rivers
        doubleprecision :: acc
        doubleprecision :: dt
        doubleprecision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        integer :: ntt
        integer :: node_to_flow_mapping(:)

        double precision, allocatable, dimension(:,:) :: pe_step
        double precision, allocatable, dimension(:,:) :: r_gau_step

        double precision :: q(num_rivers, nstep)
        doubleprecision :: rivers(nac, num_rivers)

        ! Parameters
        doubleprecision, dimension(:) :: smax
        doubleprecision, dimension(:) :: srmax
        doubleprecision, dimension(:) :: szm
        doubleprecision, dimension(:) :: t0dt
        doubleprecision, dimension(:) :: td
        doubleprecision :: wt

        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh

        ! Local Declares
        integer :: ia
        integer :: it
        doubleprecision :: sat_ac

        ! End declares

        dtt = dt / dble (ntt)

        !==========================================================
        !  START LOOP ON TIME STEPS (MAIN MODEL LOOP)
        !=========================================================

        do it = 1, nstep

            !  INITIALISE QB FOR EVERY RIVER REACH FOR THIS TIMESTEP
            !
            dyna_riv%qof = 0
            dyna_riv%qb = 0
            dyna_hru%qin = 0

            !
            !  CALL THE RAIN ROUTINE TO SPLIT UP THE RAIN INPUTS
            !  AND DISTRIBUTE THE EVAPORATION DATA
            !
            call rain (nac, &
                dyna_hru, &
                it, &
                pe_step, &
                r_gau_step)
            !
            !================================================================
            !
            !  FIRST PART OF THE DYNAMIC TOPMODEL FORMULATION
            !
            !  Distribute current QBF values to other elements
            !  Old time step QBF used to form new time step QIN values
            !
            !  CALL DYNAMIC_DIST() TO DO THESE CALCULATIONS
            !===========================================================
            !
            call dynamic_dist (nac, &
                num_rivers, &
                dtt, &
                dyna_hru, &
                dyna_riv, &
                rivers)

            !=================================================================
            !  START LOOP ON Fractional Area Cells
            !=================================================================

            sat_ac = 0.D0

            do ia = 1, nac

                dyna_hru(ia)%quz = 0.D0
                dyna_hru(ia)%uz = 0.D0
                dyna_hru(ia)%ex = 0.D0
                dyna_hru(ia)%exs = 0.D0
                dyna_hru(ia)%exus = 0.D0

                !
                !=============================================================
                !  CALL THE APPROPRIATE ROOT ZONE FORMULATION (1 is standard)
                !  THIS CALCULATES THE MAIN P OR EVAP INPUTS/OUTPUTS TO SRZ
                !=============================================================
                if (dyna_hru(ia)%ims_rz.eq.1) then

                    call root_zone1 (nac, &
                        dyna_hru, &
                        ia, &
                        it, &
                        srmax)

                end if
                !
                !=============================================================
                !  CALL THE APPROPRIATE SUZ ZONE FOMULATION (1 is standard)
                !  THIS CALCULATES QUZ AND CHANGES IN THE SUZ STORE
                !=============================================================

                if (dyna_hru(ia)%ims_uz.eq.1) then

                    call unsat_zone1 (nac, &
                        dt, &
                        dtt, &
                        dyna_hru, &
                        ia, &
                        it, &
                        td)

                end if

                !=============================================================
                !  CALL THE KINEMATIC SOLUTION
                !  Implicit time stepping. Iterative solution of kinematic equation
                !  for downslope flows at this time step
                !=============================================================
                !

                    call satzone_analyticalsolver (nac, &
                        dtt, &
                        dyna_hru, &
                        ia, &
                        it, &
                        ntt, &
                        smax, &
                        szm)

                !
                !=============================================================
                !  CALL THE FRACTIONAL AREA RESULTS SUBROUTINE
                !=============================================================
                !
                call results_ac (nac, &
                    num_rivers, &
                    dyna_hru, &
                    dyna_riv, &
                    ia)

                !==========================================================
                !  END OF Fractional area LOOP
                !==========================================================

                end do

                !=============================================================
                !  CALL THE RIVER ROUTING SUBROUTINE
                !=============================================================

                call river (dyna_riv, &
                    nstep, &
                    num_rivers, &
                    it, &
                    node_to_flow_mapping, &
                    q, &
                    route_riv, &
                    route_tdh)

                !=============================================================
                !  CALL THE TIMESTEP RESULTS SUBROUTINE - CHECKS WATER BALANCE
                !=============================================================

                call results_ts (dtt, &
                    nac, &
                    num_rivers, &
                    dyna_hru, &
                    dyna_riv)

                !==========================================================
                !    END OF TIME STEP LOOP
                !==========================================================

                end do

                return

            end subroutine

        end module dyna_topmod
