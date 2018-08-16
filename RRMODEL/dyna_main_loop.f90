module dyna_main_loop
    implicit none
contains
    !
    !===========================================================
    !  MainLoop Subroutine to set parameters, initialise model and run Model
    !===========================================================
    !
    subroutine mainloop (all_pm_names, &
        nac, &
        nstep, &
        num_rivers, &
        acc, &
        dt, &
        dtt, &
        dyna_hru, &
        mcpar, &
        ntt, &
        node_to_flow_mapping, &
        num_par_types, &
        pe_step, &
        qobs_riv_step_start, &
        r_gau_step, &
        rivers, &
        route_riv, &
        route_tdh, &
        sum_ac_riv, &
        wt)

        use dyna_common_types
        use dyna_param_setup
        use dyna_initialise_run
        use dyna_init_new_route
        use dyna_topmod
        use dyna_write_output
        use dta_route_processing

        implicit none

        ! Argument Declares

        integer :: nac
        integer :: nstep
        integer :: num_rivers
        doubleprecision :: dt
        doubleprecision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)

        ! Parameters
        character(len=20), dimension(7) :: all_pm_names !names of all possible parameters
        doubleprecision, allocatable, dimension(:,:) :: mcpar
        integer :: num_par_types
        integer :: ntt
        doubleprecision :: wt
        doubleprecision :: acc
        doubleprecision, allocatable, dimension(:) :: srmax
        doubleprecision, allocatable, dimension(:) :: chv
        doubleprecision, allocatable, dimension(:) :: chvdt
        doubleprecision, allocatable, dimension(:) :: lnto
        doubleprecision, allocatable, dimension(:) :: smax
        doubleprecision, allocatable, dimension(:) :: srinit
        doubleprecision, allocatable, dimension(:) :: szm
        doubleprecision, allocatable, dimension(:) :: t0dt
        doubleprecision, allocatable, dimension(:) :: td

        ! Input and flow data
        double precision, allocatable, dimension(:,:) :: pe_step
        double precision, allocatable, dimension(:,:) :: r_gau_step
        double precision :: qobs_riv_step_start(num_rivers)

        ! River Data
        doubleprecision :: rivers(nac, num_rivers)
        double precision, dimension(:,:), allocatable :: q
        ! river tree and river cell information for entire river network
        type(route_river_info_type) :: route_riv
        ! histogram for entire river network
        type(route_time_delay_hist_type) :: route_tdh
        integer :: node_to_flow_mapping(:)
        doubleprecision :: sum_ac_riv (num_rivers)

        ! End declares

        !
        !=========================================================
        !  PUT THE PARAM VALUES INTO THE VARIABLE NAMES FOR THE MODEL TYPE
        !=========================================================
        !
        call param_setup (chv, &
            chvdt, &
            dt, &
            lnto, &
            num_par_types, &
            mcpar, &
            smax, &
            srinit, &
            srmax, &
            szm, &
            t0dt, &
            td)
        !
        !=========================================================
        !  INITIALISE VARIABLES, ROUTING AND STORES
        !=========================================================
        !
        call initialise_run (dyna_hru, &
            dyna_riv, &
            nac, &
            num_rivers)

        call init_new_route (nac, &
            nstep, &
            num_rivers, &
            chvdt, &
            dt, &
            dyna_hru, &
            mcpar, &
            q, &
            qobs_riv_step_start, &
            route_riv, &
            route_tdh, &
            smax, &
            srinit, &
            srmax, &
            sum_ac_riv, &
            szm, &
            t0dt)

        !
        !=========================================================
        !  CALL THE MAIN TOPMODEL STRUCTURE
        !=========================================================
        !
        call Topmod (nac, &
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

        !
        !=========================================================
        !  WRITE OUTPUTS
        !=========================================================

        call write_output(dyna_hru, &
            mcpar, &
            nac, &
            all_pm_names, &
            num_par_types, &
            num_rivers, &
            q)

        deallocate(mcpar)

        return
                                                                        
    end subroutine mainloop

end module dyna_main_loop
