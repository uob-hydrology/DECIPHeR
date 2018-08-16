module dyna_results_ac
contains
    !
    !===============================================================
    !  ROUTINE TO REDISTRIBUTE QIN AND QBF, ADD OVERLAND FLOW TO QOF
    !===============================================================
    !
    subroutine results_ac (nac, &
        num_rivers, &
        dyna_hru, &
        dyna_riv, &
        ia)

        use dyna_common_types

        implicit none
        ! Argument Declares
        integer :: nac
        integer :: num_rivers
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        integer :: ia

        ! End declares

        !  PUT CURRENT QIN AND QBF VALUES INTO QINt1 and QBFt1
        dyna_hru(ia)%qint1 = dyna_hru(ia)%qin
        dyna_hru(ia)%qbft1 = dyna_hru(ia)%qbf

        if (dyna_hru(ia)%exs .lt.0.D0.or.dyna_hru(ia)%exus .lt.0.D0) then
            print  * ,' Problem with EX ', dyna_hru(ia)%exs , dyna_hru(ia)%exus
        endif

        dyna_hru(ia)%ex = dyna_hru(ia)%exs + dyna_hru(ia)%exus
        if (dyna_hru(ia)%ex .gt.0.D0) then
            dyna_riv(dyna_hru(ia)%ipriv)%qof = dyna_riv(dyna_hru(ia)%ipriv)%qof + &
               dyna_hru(ia)%ex * dyna_hru(ia)%ac
        endif

    end subroutine

end module dyna_results_ac
