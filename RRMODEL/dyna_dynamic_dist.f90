module dyna_dynamic_dist
contains
    !
    !===============================================================
    !  ROUTINE FOR CALCULATING DYNAMIC REDISTRIBUTION OF QBF(ts-1)
    !===============================================================
    !
    !  FIRST PART OF THE DYNAMIC TOPMODEL FORMULATION
    !
    !  Distribute current QBF values to other elements
    !  Old time step QBF used to form new time step QIN values
    !  NB doing this before the new time step calculations allows for
    !  calculation of redistribution to the same element, which
    !  otherwise reduces inputs into river.
    !  Otherwise fully iterative solution required over all increments
    !
    !  CALCULATE QIN(ts) VALUES FROM PREVIOUS QBF(ts-1) VALUES
    !  ROUTE QBF(ts-1) VALUES TO THE RIVER INTO QB
    !
    subroutine dynamic_dist (nac, &
        num_rivers, &
        dtt, &
        dyna_hru, &
        dyna_riv, &
        rivers)

        use dyna_common_types

        implicit none
        ! Argument Declares
        integer :: nac
        integer :: num_rivers
        double precision :: dtt
        type(dyna_hru_type) :: dyna_hru(nac)
        type(dyna_riv_type) :: dyna_riv(num_rivers)
        double precision :: rivers(nac, num_rivers)

        ! Local Declares
        integer :: ia
        integer :: iriv
        integer :: iw

        ! End declares

        !  REDISTRIBUTE QBF(ts-1) VALUES TO QIN(ts) VALUES
                                                                        
        do ia = 1, nac

            do iw = 1, dyna_hru(ia)%ntrans
                dyna_hru(dyna_hru(ia)%itrans(iw))%qin = dyna_hru(dyna_hru(ia)%itrans(iw))%qin + dyna_hru(ia)%qbf  &
                    * dyna_hru(ia)%wtrans(iw) * dyna_hru(ia)%ac
            end do

            !  ADD QBF INPUTS DIRECTLY INTO THE RIVER NETWORK
            !  ADD QBF INPUTS TO THE RIVERS FOR EACH RIVER REACH
            !
            do iriv = 1, num_rivers

                dyna_riv(iriv)%qb = dyna_riv(iriv)%qb + dtt * dyna_hru(ia)%qbf &
                    * dyna_hru(ia)%wtrans(dyna_hru(ia)%ntrans + 1) &
                    * rivers (ia, iriv) * dyna_hru(ia)%ac
                                                                        
                dyna_riv(iriv)%sumqb = dyna_riv(iriv)%sumqb + dtt * dyna_hru(ia)%qbf &
                    * dyna_hru(ia)%wtrans(dyna_hru(ia)%ntrans + 1) &
                    * rivers (ia, iriv) * dyna_hru(ia)%ac
            end do
        end do

    end subroutine

end module dyna_dynamic_dist
