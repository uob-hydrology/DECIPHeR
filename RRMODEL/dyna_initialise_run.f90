module dyna_initialise_run
contains
    !
    !===============================================================
    !  ROUTINES FOR INITIALISING THE MODEL RUN VARIABLES
    !===============================================================
    !

    subroutine initialise_run (dyna_hru, &
            dyna_riv, &
            nac, &
            num_rivers)

        use dyna_common_types

        implicit none
        ! Argument Declares

        integer :: nac, num_rivers
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)

        ! End declares

        !  INITIALISE MODEL VARIABLES

        dyna_hru%sum_aeout = 0
        dyna_hru%sum_pptin = 0
        dyna_hru%sum_pein = 0

        dyna_hru%sumexs = 0
        dyna_hru%sumexus = 0

        dyna_riv%sumqb = 0
        dyna_riv%sumqout = 0
        dyna_riv%sumqof = 0

    end subroutine

end module dyna_initialise_run
