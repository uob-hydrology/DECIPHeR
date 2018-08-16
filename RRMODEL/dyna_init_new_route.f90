module dyna_init_new_route

contains
    !
    !==========================================================
    !
    !  JEF 03/07/16
    !  SUBROUTINE FOR ONLY THE NEW ROUTING SCHEME FOR THE NEW
    !  MODEL, INCLUDES:
    !  1) Calculate all the starting flows (if some missing)
    !     a) Do this by the mean of the flows
    !  2) TYPE 1 Do a super simple homogeneous initialisation by:
    !     a) All subcatchments are set at the main outfall flow
    !     b) So ALL SD's calculated to that one value
    !     c) No need to consider downslope contributions for each
    !     subcatchment
    !
    !==========================================================
    !
    subroutine init_new_route (nac, &
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

        use dyna_common_types
        use dta_route_processing

        implicit none

        ! Argument Declares
        integer :: nac
        integer :: nstep
        integer :: num_rivers

        double precision :: dt
        type(dyna_hru_type) :: dyna_hru(nac)

        double precision, dimension(:,:), allocatable :: q
        double precision :: qobs_riv_step_start(num_rivers) !(n_riv)
        double precision :: mean_q_ac_weighted
        double precision :: sum_ac_riv(num_rivers)

        ! Parameters
        double precision, dimension(:) :: smax
        double precision, dimension(:) :: srinit
        double precision, dimension(:) :: srmax
        double precision, dimension(:) :: szm
        double precision, dimension(:) :: t0dt
        double precision, dimension(:,:) :: mcpar
        double precision, dimension(:) :: chvdt

        ! MULTI POINT RIVER ROUTING
        ! Toby Dunne April 2016
        ! these type are defined in modules as part of the dta files
        type(route_river_info_type) :: route_riv
        type(route_time_delay_hist_type) :: route_tdh

        ! Local Declares
        integer :: i
        integer :: ia
        integer :: iriv

        ! End declares
        !
        !========================
        !  DYNAMIC INITIALISATION

        ! At the moment the velocity parameter is set to the first value of channel velocity
        ! Routing codes can't support multiple values of the v_param
        route_tdh%v_param(1) = chvdt(1)
        route_tdh%timestep = dt

        ! Initialise the time delay histogram with velocity and timestep
        call route_processing_build_hist(route_riv, route_tdh)

        ! JEF 03/07/16 - First determine the mean area weighted starting flow
        !
        mean_q_ac_weighted = 0.D0
        do iriv = 1, num_rivers
            mean_q_ac_weighted = mean_q_ac_weighted + qobs_riv_step_start(iriv) * sum_ac_riv(iriv)
        enddo

        do ia = 1, nac
                !
                ! JEF 03/07/16 Only initialise to mean_q_ac_weighted for the whole catchment
                ! and set all dyna_hru(ia)%qbf to this
                !
                dyna_hru(ia)%sd = - SZM (dyna_hru(ia)%ipar ) * (Log (mean_q_ac_weighted / dyna_hru(ia)%st)  &
                    - T0DT (dyna_hru(ia)%ipar))
                dyna_hru(ia)%qbf = mean_q_ac_weighted

                !
                ! JEF 03/07/16 Now check and make sure nothing is going beyond certain limits
                !
                if (dyna_hru(ia)%sd .lt.0.D0) then
                    dyna_hru(ia)%sd = 0.D0
                elseif (dyna_hru(ia)%sd .ge.Smax (dyna_hru(ia)%ipar ) ) then
                    dyna_hru(ia)%sd = Smax (dyna_hru(ia)%ipar )
                    dyna_hru(ia)%qbf = 0.D0
                elseif (dyna_hru(ia)%sd .lt.Smax (dyna_hru(ia)%ipar ) .and.dyna_hru(ia)%qbf .gt.  &
                    ( Smax (dyna_hru(ia)%ipar ) - dyna_hru(ia)%sd ) ) then
                    dyna_hru(ia)%qbf = (Smax (dyna_hru(ia)%ipar ) - dyna_hru(ia)%sd )
                endif

            dyna_hru(ia)%qbft1 = dyna_hru(ia)%qbf
            dyna_hru(ia)%qint1 = dyna_hru(ia)%qbft1 * dyna_hru(ia)%ac
            dyna_hru(ia)%qbfini = dyna_hru(ia)%qbf

        enddo
        !
        !  INITIALISE SRZ VALUES HERE
        !  SRinit IS INITIAL ROOT ZONE STORAGE DEFICIT BELOW FIELD CAPACITY
        do ia = 1, nac

            dyna_hru(ia)%suz = 0.D0

            if (SRinit (dyna_hru(ia)%ipar ) .gt.SRmax (dyna_hru(ia)%ipar ) ) then
                srinit (dyna_hru(ia)%ipar ) = SRmax (dyna_hru(ia)%ipar )
                mcpar (dyna_hru(ia)%ipar, 4) = SRinit (dyna_hru(ia)%ipar )
            endif
            dyna_hru(ia)%srz = (SRmax (dyna_hru(ia)%ipar ) - SRinit (dyna_hru(ia)%ipar ) )
            dyna_hru(ia)%sdini = dyna_hru(ia)%sd
            dyna_hru(ia)%srzini = dyna_hru(ia)%srz
            dyna_hru(ia)%suzini = dyna_hru(ia)%suz
        end do
        !
        !  INITIALISE DISCHARGE ARRAY
        !
        call checked_allocate(q, num_rivers, nstep)
        do i = 1, nstep
            do iriv = 1, num_rivers
                q (iriv, i) = 0.D0
            end do
        end do

    end subroutine
end module dyna_init_new_route


