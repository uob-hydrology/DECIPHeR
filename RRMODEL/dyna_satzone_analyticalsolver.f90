module dyna_satzone_analyticalsolver

contains
    !
    !===============================================================
    !  ROUTINE FOR CALCULATING THE KINEMATIC SOLUTION
    !  EXPONENTIAL FORMULATION BY ROSS WOODS, WOUTER KNOBEN AND GEMMA
    !  COXON
    !===============================================================
    !
    subroutine satzone_analyticalsolver (nac, &
        dtt, &
        dyna_hru, &
        ia, &
        it, &
        ntt, &
        smax, &
        szm)

        use dyna_common_types

        implicit none
        ! Argument Declares
        integer :: nac
        double precision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        integer :: ia
        integer :: it
        integer :: ntt
        double precision, dimension(:) :: smax
        double precision, dimension(:) :: szm

        ! Local Declares
        integer :: itt

        double precision :: q1, q2, qin2
        double precision :: u, u2
        double precision :: m2, tc
        double precision :: S_temp, t_temp

        double precision :: storebal
        double precision :: storeend
        double precision :: storein
        double precision :: storeini
        double precision :: storeout
        ! End declares

        ! Define constants - Q0 is defined in the initalisation, mean slope is in the hru_meta file
        q1 = dyna_hru(ia)%q0 * cos(dyna_hru(ia)%mslp)
        q2 = dyna_hru(ia)%q0 *cos(dyna_hru(ia)%mslp) * exp(-1*cos(dyna_hru(ia)%mslp) &
            * smax(dyna_hru(ia)%ipar)/szm(dyna_hru(ia)%ipar))
        m2 = szm(dyna_hru(ia)%ipar)/cos(dyna_hru(ia)%mslp)

        ! Water balance calcs
        storeini = dyna_hru(ia)%sd
        storein = 0
        storeout = 0

        ! Solve downslope kinematic wave equation - main loop
        do itt = 1, ntt

            ! Linear interpolation for inflows for inner time loop
            qin2 = (dyna_hru(ia)%qint1 + itt * (dyna_hru(ia)%qin - dyna_hru(ia)%qint1) / ntt) / dyna_hru(ia)%ac;

            ! Inputs
            u  = qin2 + dyna_hru(ia)%uz;
            u2 = q2 + u;

            ! Select the proper case
            if (dyna_hru(ia)%sd.le.Smax(dyna_hru(ia)%ipar)) then
                if (u2.gt.0) then
                    if (u2*dtt/m2.lt.1e-10) then
                        ! Very small u2 values can lead to undesirable sd values
                        ! when the terms in the equation have very different orders
                        ! of magnitude. This bit detects and prevents that
                        dyna_hru(ia)%sd = m2*log( (1/exp(-dyna_hru(ia)%sd/m2) + q1*dtt/m2) - &
                            u2*dtt/m2*(1/exp(-dyna_hru(ia)%sd/m2)-0.5*q1*dtt/m2))
                    else
                        ! Regular case for u2 values that are not very small
                        dyna_hru(ia)%sd = m2*log( q1/u2 + (1/exp(-dyna_hru(ia)%sd/m2) - q1/u2)*exp(-u2*dtt/m2))
                    end if
                else
                    ! Case where u2 = 0
                    dyna_hru(ia)%sd = m2*log( exp(dyna_hru(ia)%sd/m2) + q1*dtt/m2)
                end if
            else  ! S > Smax
                tc = (dyna_hru(ia)%sd-smax(dyna_hru(ia)%ipar))/u
                if (tc.gt.dtt) then
                    dyna_hru(ia)%sd = dyna_hru(ia)%sd - u*dtt
                else
                    t_temp = dtt-tc
                    S_temp = smax(dyna_hru(ia)%ipar)
                    dyna_hru(ia)%sd = m2*log( q1/u2 + (1/exp(-S_temp/m2) - q1/u2)*exp(-u2*t_temp/m2) )
                end if
            end if

            ! Find Qbf
            dyna_hru(ia)%qbf = (dyna_hru(ia)%sd - storeini)/dtt + &
                dyna_hru(ia)%uz + dyna_hru(ia)%qin / dyna_hru(ia)%ac;

            ! WB check
            storein = storein + dtt * ((-1 * dyna_hru(ia)%uz) - dyna_hru(ia)%qin / dyna_hru(ia)%ac)

        end do

        ! If SD is less than zero calculate the saturation excess EXS
        if (dyna_hru(ia)%sd.lt.0) then
            dyna_hru(ia)%exs = (dyna_hru(ia)%qbf - dyna_hru(ia)%q0)*dtt;
            dyna_hru(ia)%exs = 0 - dyna_hru(ia)%sd;
            dyna_hru(ia)%sumexs = dyna_hru(ia)%sumexs + (dyna_hru(ia)%exs * dyna_hru(ia)%ac);
            dyna_hru(ia)%sd = 0;
        end if

        ! Ensure QBF is not greater than Q0
        if (dyna_hru(ia)%qbf.gt.dyna_hru(ia)%q0) then
            dyna_hru(ia)%exs = dyna_hru(ia)%exs + (dyna_hru(ia)%qbf - dyna_hru(ia)%q0)*dtt;
            dyna_hru(ia)%sumexs = dyna_hru(ia)%sumexs + (((dyna_hru(ia)%qbf - dyna_hru(ia)%q0)*dtt) * dyna_hru(ia)%ac);
            dyna_hru(ia)%qbf = dyna_hru(ia)%q0;
        end if

        ! WB check
        storeout = storeout + dyna_hru(ia)%exs
        storeout = storeout + dyna_hru(ia)%qbf * dtt;
        storeend = dyna_hru(ia)%sd
        storebal = storeini + storein + storeout - storeend

        if ((storebal.gt.0.0000001).or.(storebal.lt.- 0.0000001)) then
            print * , ' Storebal SD wrong KINEMATIC ', it, ia, storebal,   &
                storeini, storein, storeout, storeend, dyna_hru(ia)%uz , dyna_hru(ia)%qin , &
                dyna_hru(ia)%ac , Smax (dyna_hru(ia)%ipar ) , dtt, ntt
        end if

    end subroutine

end module dyna_satzone_analyticalsolver



