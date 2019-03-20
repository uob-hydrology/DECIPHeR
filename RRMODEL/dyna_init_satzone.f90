module dyna_init_satzone

contains
    !
    !==========================================================
    !
    !  GC 17/02/19
    !  SUBROUTINE FOR INITIALISING THE SAT STORE OF THE MODEL
    !
    !==========================================================
    !
    subroutine init_satzone(nac, &
        nstep, &
        num_rivers, &
        chvdt, &
        dt, &
        dyna_hru, &
        mcpar, &
        q, &
        qobs_riv_step_start, &
        rivers, &
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
        double precision :: rivers(nac, num_rivers)

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
        integer :: sol_flag

        double precision :: D
        double precision :: q0
        double precision :: meanatb, meanT0DT, meanSZM
        double precision :: qb_tmp
        double precision :: mean_q_ac_weighted_guess
        double precision :: ndiff, qb_diff, qb_diff_prev

        ! End declares
        !
        !========================
        !  DYNAMIC INITIALISATION

        ! At the moment the velocity parameter is set to the first value of channel velocity
        ! Routing codes can't support multiple values of the v_param - this needs to be improved!
        route_tdh%v_param(1) = chvdt(1)
        route_tdh%timestep = dt

        ! Initialise the time delay histogram with velocity and timestep
        call route_processing_build_hist(route_riv, route_tdh)

        ! First determine the mean area weighted starting flow
        mean_q_ac_weighted = 0
        do iriv = 1, num_rivers
            mean_q_ac_weighted = mean_q_ac_weighted + qobs_riv_step_start(iriv) * sum_ac_riv(iriv)
        enddo

        ! Then determine the mean topographic index and T0DT value so you can calculate the mean deficit
        meanatb = 0
        meanT0DT = 0
        meanSZM = 0

        do ia = 1, nac
            meanatb = meanatb + (dyna_hru(ia)%st * dyna_hru(ia)%ac)
            meanT0DT = meanT0DT + (T0DT(dyna_hru(ia)%ipar) * dyna_hru(ia)%ac)
            meanSZM = meanSZM + (szm(dyna_hru(ia)%ipar) * dyna_hru(ia)%ac)

            ! Initialise q0 for each HRU (used throughout the code)
            dyna_hru(ia)%q0 = exp(T0DT(dyna_hru(ia)%ipar) - dyna_hru(ia)%st)
        enddo

        ! Take an initial guess at what your sat store should be based on the mean deficit
        q0 = exp(meanT0DT-meanatb)
        D = -meanSZM* log(mean_q_ac_weighted/q0)

        mean_q_ac_weighted_guess = mean_q_ac_weighted
        ndiff = 0.0001;
        qb_diff = 1;

        ! Set sol_flag = 1
        sol_flag = 1

        ! Check that this flow is physically possible (i.e. can you gain this flow from qbf contributions only)

        call calculate_qbmax(dyna_hru, qb_tmp, nac, num_rivers, dt, rivers)
        
        if (qb_tmp.le.mean_q_ac_weighted) then
            
           ! This will be the best solution - all deficit stores are set to zero and qbf to q0
           qb_diff = 0
           sol_flag = 2
        end if
        

        do while (abs(qb_diff).gt.0.01)

            call calculate_qb(D, dyna_hru, qb_tmp, szm, T0DT, smax, meanatb, &
                            meanT0DT, nac, num_rivers, dt, rivers)

            qb_diff_prev = qb_diff;
            qb_diff = (qb_tmp-mean_q_ac_weighted)*1000

            if (qb_diff.gt.0.01) then

                mean_q_ac_weighted_guess = mean_q_ac_weighted_guess - ndiff

            elseif (qb_diff.lt.-0.01) then

                mean_q_ac_weighted_guess = mean_q_ac_weighted_guess + ndiff

            end if

            ! Calculate new starting deficit
            D = -meanSZM* log(mean_q_ac_weighted_guess/q0)

            ! Checks to make sure things aren't oscillating
            if ((qb_diff_prev.gt.0).and.(qb_diff.lt.0)) then
                ! Code is oscillating between over and underestimating - try a
                ! smaller value of ndiff
                if (ndiff.gt.1e-9) then
                    ndiff = ndiff/10
                else
                    sol_flag = 0
                    if (abs(qb_diff).lt.abs(qb_diff_prev)) then
                        exit
                    else
                        ! Re-calculate for previous solution
                            call calculate_qb(D, dyna_hru, qb_tmp, szm, T0DT, smax, meanatb, &
                            meanT0DT, nac, num_rivers, dt, rivers)
                            exit
                        endif
                    endif
                endif

            enddo

            ! Write out results to .res file
            write(40, *) ''
            write(40, *) '! INITIALISATION'

            if (sol_flag.eq.1) then
                write(40, *) 'Solution found!'
            else if (sol_flag.eq.0) then
                write(40, *) 'No solution found - going for the minimum solution'
            else if (sol_flag.eq.2) then
                write(40, *) 'Maximum QBF is less than mean_q_ac_weighted - deficits set to 0 and max QBFs'
            end if

            write(40, *) 'Starting qobs ', mean_q_ac_weighted * 1000, 'mm/day'
            write(40, *) 'Best estimate of QB riv', qb_tmp * 1000, 'mm/day'
            write(40, *) 'Finishing qt0 ', mean_q_ac_weighted_guess * 1000, 'mm/day'
            write(40, *) ''

            ! Initialise a lot of soil root zone and sat zone variables
            do ia = 1, nac

                ! Sat Zone Variables
                dyna_hru(ia)%qbft1 = dyna_hru(ia)%qbf
                dyna_hru(ia)%qint1 = dyna_hru(ia)%qbft1 * dyna_hru(ia)%ac
                dyna_hru(ia)%qbfini = dyna_hru(ia)%qbf
                dyna_hru(ia)%suz = 0

                !Root Zone Variables
                if (SRinit (dyna_hru(ia)%ipar) .gt.SRmax (dyna_hru(ia)%ipar ) ) then
                    srinit (dyna_hru(ia)%ipar) = SRmax (dyna_hru(ia)%ipar )
                    mcpar (dyna_hru(ia)%ipar, 4) = SRinit (dyna_hru(ia)%ipar )
                endif

                dyna_hru(ia)%srz = (SRmax (dyna_hru(ia)%ipar ) - SRinit (dyna_hru(ia)%ipar ) )
                dyna_hru(ia)%sdini = dyna_hru(ia)%sd
                dyna_hru(ia)%srzini = dyna_hru(ia)%srz
                dyna_hru(ia)%suzini = dyna_hru(ia)%suz

            enddo
            !
            !  INITIALISE DISCHARGE ARRAY
            !
            call checked_allocate(q, num_rivers, nstep)
            do i = 1, nstep
                do iriv = 1, num_rivers
                    q (iriv, i) = 0
                end do
            end do

        end subroutine

        subroutine calculate_qb(D, &
            dyna_hru, &
            qb_tmp, &
            szm, &
            T0DT, &
            smax, &
            meanatb, &
            meanT0DT, &
            nac, &
            num_rivers, &
            dt, &
            rivers)

            use dyna_common_types

            implicit none

            ! Argument Declares
            integer :: nac
            integer :: num_rivers

            double precision :: dt
            type(dyna_hru_type) :: dyna_hru(nac)

            ! Parameters
            double precision, dimension(:) :: smax
            double precision, dimension(:) :: szm
            double precision, dimension(:) :: t0dt

            double precision :: D
            double precision :: meanatb, meanT0DT
            double precision :: qb_tmp
            double precision :: rivers(nac, num_rivers)

            ! Local Declares
            integer :: ia
            integer :: iriv

            ! End declares

            do ia = 1 , nac

                dyna_hru(ia)%sd = D + szm(dyna_hru(ia)%ipar) * ((meanatb - dyna_hru(ia)%st) + (T0DT(dyna_hru(ia)%ipar) - meanT0DT))
                dyna_hru(ia)%qbf = dyna_hru(ia)%q0 * exp(-dyna_hru(ia)%sd/szm(dyna_hru(ia)%ipar))

                if (dyna_hru(ia)%sd.lt.0) then
                    dyna_hru(ia)%sd = 0
                    dyna_hru(ia)%qbf = exp(T0DT(dyna_hru(ia)%ipar)-dyna_hru(ia)%st)
                elseif (dyna_hru(ia)%sd.ge.smax(dyna_hru(ia)%ipar)) then
                    dyna_hru(ia)%sd = smax(dyna_hru(ia)%ipar)
                    dyna_hru(ia)%qbf = 0
                endif

            enddo

            ! Calculate contibutions to river - we do not try and optimise for each river,
            ! we treat it as a single river stretch and optimise to that
            qb_tmp = 0;

            do ia = 1 , nac

                do iriv = 1 , num_rivers

                    qb_tmp = qb_tmp + (dt * dyna_hru(ia)%qbf * dyna_hru(ia)%wtrans(dyna_hru(ia)%ntrans+1) &
                    * rivers(ia, iriv) * dyna_hru(ia)%ac);

                enddo

            enddo

        end subroutine calculate_qb

        subroutine calculate_qbmax(dyna_hru, &
            qb_tmp, &
            nac, &
            num_rivers, &
            dt, &
            rivers)

            use dyna_common_types

            implicit none

            ! Argument Declares
            integer :: nac
            integer :: num_rivers

            double precision :: dt
            type(dyna_hru_type) :: dyna_hru(nac)
            double precision :: qb_tmp
            double precision :: rivers(nac, num_rivers)

            ! Local Declares
            integer :: ia
            integer :: iriv

            ! End declares

            do ia = 1 , nac

                ! Assume that all HRUs are at max saturation
                dyna_hru(ia)%sd = 0
                dyna_hru(ia)%qbf = dyna_hru(ia)%q0 

            enddo

            ! Calculate contibutions to river - we do not try and optimise for each river,
            ! we treat it as a single river stretch and optimise to that
            qb_tmp = 0;

            do ia = 1 , nac

                do iriv = 1 , num_rivers

                    qb_tmp = qb_tmp + (dt * dyna_hru(ia)%qbf * dyna_hru(ia)%wtrans(dyna_hru(ia)%ntrans+1) &
                    * rivers(ia, iriv) * dyna_hru(ia)%ac);

                enddo

            enddo

        end subroutine calculate_qbmax

    end module dyna_init_satzone


