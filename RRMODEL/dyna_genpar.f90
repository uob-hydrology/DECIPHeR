module dyna_genpar
contains
    !
    !=======================================================================
    !     GETS THE RANDOM NUMBER VALUES FOR MC RUNS
    !=======================================================================
    !
    subroutine genpar (mcpar, &
        mcpar_ll, &
        mcpar_ul, &
        num_mcpar, &
        num_par_types)

        use dyna_common_types
        use dyna_random

        implicit none

        ! Argument Declares
        double precision, allocatable, dimension(:,:) :: mcpar
        double precision, dimension(:,:) :: mcpar_ll
        double precision, dimension(:,:) :: mcpar_ul
        integer, dimension(:) :: num_mcpar
        integer :: num_par_types

        ! Local Declares
        integer :: ig
        integer :: igen
        double precision :: rand(1)

        ! End declares

        allocate(mcpar(num_par_types, num_mcpar(1)))

        do ig = 1, num_par_types
            do igen = 1, num_mcpar (ig)
                call ranaru (rand, 1)
                mcpar (ig, igen) = mcpar_ll (ig, igen) + rand(1) *       &
                    (mcpar_ul (ig, igen) - mcpar_ll (ig, igen) )
            end do
        end do

    end subroutine

end module dyna_genpar
