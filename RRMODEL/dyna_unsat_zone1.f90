module dyna_unsat_zone1
contains
    !
    !===============================================================
    !  ROUTINES FOR CALCULATING THE UNSAT ZONE FUNCTIONS
    !  UNSAT_ZONE1 = BASIC UNSAT ZONE FORMULATION
    !===============================================================

    subroutine unsat_zone1 (nac, &
    dt, &
    dtt, &
    dyna_hru, &
    ia, &
    it, &
    td)

        use dyna_common_types
        implicit none

! Argument Declares
        integer :: nac
        double precision :: dt
        double precision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        integer :: ia
        integer :: it
        double precision, dimension(:) :: td

! Local Declares
        double precision :: storebal
        double precision :: storeend
        double precision :: storein
        double precision :: storeini
        double precision :: storeout
        double precision :: uz2
! End declares

        !  UZ CALCULATIONS (EX IS USED FOR SATURATION EXCESS)
        !  ALL SUZ BECOMES EX IF dyna_hru(IA)%sd IS ZERO
        uz2 = 0.D0
        dyna_hru(ia)%uz = 0.D0
        storeini = dyna_hru(ia)%suz
        storein = 0
        storeout = 0

        !
        !  ADD PEX TO SUZ FROM THE ROOT ZONE CALCULATIONS
        !
        dyna_hru(ia)%suz = dyna_hru(ia)%suz + dyna_hru(ia)%pex
        storein = dyna_hru(ia)%pex

        if (dyna_hru(ia)%sd .le.0.D0.and.dyna_hru(ia)%suz .gt.0.D0) then
            dyna_hru(ia)%exus = dyna_hru(ia)%suz
            dyna_hru(ia)%sumexus = dyna_hru(ia)%sumexus + (dyna_hru(ia)%exus * dyna_hru(ia)%ac)
            storeout = dyna_hru(ia)%exus
            dyna_hru(ia)%suz = 0.D0
        elseif (dyna_hru(ia)%sd .gt.0.D0.and.dyna_hru(ia)%suz .gt.dyna_hru(ia)%sd ) then
            dyna_hru(ia)%exus = (dyna_hru(ia)%suz - dyna_hru(ia)%sd )
            dyna_hru(ia)%sumexus = dyna_hru(ia)%sumexus + (dyna_hru(ia)%exus * dyna_hru(ia)%ac)
            storeout = storeout + dyna_hru(ia)%exus
            dyna_hru(ia)%suz = dyna_hru(ia)%sd
        endif
        !
        !  CALCULATE DRAINAGE FROM SUZ USING dyna_hru(IA)%sd AT START OF TIME STEP
        !
        if (dyna_hru(ia)%sd .gt.0.D0.and.dyna_hru(ia)%suz .gt.0.D0) then
            uz2 = dt * (dyna_hru(ia)%suz / (dyna_hru(ia)%sd * TD (dyna_hru(ia)%ipar ) ) )
            if (uz2 .gt.dyna_hru(ia)%suz ) uz2 = dyna_hru(ia)%suz
            dyna_hru(ia)%suz = dyna_hru(ia)%suz - uz2
            if (dyna_hru(ia)%suz .lt.0.0000001D0) then
                uz2 = uz2 + dyna_hru(ia)%suz
                dyna_hru(ia)%suz = 0.0D0
            endif
            dyna_hru(ia)%quz = uz2                                               
            storeout = storeout + dyna_hru(ia)%quz                                   
        endif
                                                                        
        if(uz2.gt.0.D0) dyna_hru(ia)%uz = uz2 / dtt
                                                                        
        if(dyna_hru(ia)%suz.lt.0.D0) print*,' SUZ less than zero help'
        !
        storeend = dyna_hru(ia)%suz
        storebal = storeini + storein - storeout - storeend
        if(storebal.gt.0.0000000001.or.storebal.lt.-0000000001) then
            print*,' Storebal SUZ wrong ', it,ia,storebal                   
        endif
    end subroutine

end module dyna_unsat_zone1
