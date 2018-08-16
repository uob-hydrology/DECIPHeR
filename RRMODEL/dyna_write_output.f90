module dyna_write_output
contains
!
!===============================================================
!  ROUTINE TO WRITE DYMOND OUTPUT
!===============================================================
!
    subroutine write_output(dyna_hru, &
        mcpar, &
        nac, &
        all_pm_names, &
        num_par_types, &
        num_rivers, &
        q)

        use dyna_common_types

        ! Function to write output for dynamic topmodel

        implicit none


        double precision :: q(:,:)
        double precision :: mcpar(:,:)
        integer :: i, ia
        integer :: num_par_types
        integer :: num_rivers
        integer :: nac
        type(dyna_hru_type) :: dyna_hru(nac)
        character(len=1024) :: line_format
        character(len=20), dimension(:) :: all_pm_names !names of all possible parameters
        character(len=20) :: pm_class = 'Param_class'

        double precision, dimension(num_rivers) :: tot_reach_ppt, tot_reach_pet
        double precision, dimension(num_rivers) :: tot_reach_aet, tot_reach_frac

        !format
        !CHARACTER(LEN=200), PARAMETER :: FMT0 = "(T5,T5,T5,T5,T5,T5,T5)"
        CHARACTER(LEN=200), PARAMETER :: FMT0 = "(A16,A10,  A12,  A10 , A11,  A12,  A11,  A10)"
        CHARACTER(LEN=200), PARAMETER :: FMT1 = "(I12,F10.4,F12.4,F10.4,F10.4,F13.4,F11.4,F10.4)"

        !writing to the new .res file - Parameters first

        write(40, '(A)') '! PARAMETERS'
        write(40,FMT0) [pm_class,all_pm_names(1:7)]

        do i = 1, num_par_types
            write(40,FMT1) i, mcpar(i,1:7)
        end do

        write(40, '(A)') ''

        ! writing to the .res file - summary stats

        tot_reach_pet = 0
        tot_reach_aet = 0
        tot_reach_frac = 0

        ! Calculate summary stats

        do ia = 1, nac

            tot_reach_ppt(dyna_hru(ia)%ipriv) = tot_reach_ppt(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sum_pptin*dyna_hru(ia)%ac)
            tot_reach_pet(dyna_hru(ia)%ipriv) = tot_reach_pet(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sum_pein*dyna_hru(ia)%ac)
            tot_reach_aet(dyna_hru(ia)%ipriv) = tot_reach_aet(dyna_hru(ia)%ipriv) + (dyna_hru(ia)%sum_aeout*dyna_hru(ia)%ac)

            tot_reach_frac(dyna_hru(ia)%ipriv) = tot_reach_frac(dyna_hru(ia)%ipriv) + dyna_hru(ia)%ac

        end do

        do ia = 1, num_rivers

            tot_reach_ppt(ia) = (tot_reach_ppt(ia) / tot_reach_frac(ia))*1000
            tot_reach_pet(ia) = (tot_reach_pet(ia) / tot_reach_frac(ia))*1000
            tot_reach_aet(ia) = (tot_reach_aet(ia) / tot_reach_frac(ia))*1000

        end do

        write (line_format, '(A,I0, A)') '(',num_rivers,'(1X, f0.2), A)'

        write(40, '(A)') '! SUMMARY STATS PER RIVER REACH'

        write(40, '(A)') ''
        write(40, '(A)') '! PRECIP TOTALS (mm)'
        write(40, line_format) tot_reach_ppt

        write(40, '(A)') ''
        write(40, '(A)') '! PET TOTALS (mm)'
        write(40, line_format) tot_reach_pet

        write(40, '(A)') ''
        write(40, '(A)') '! ET TOTALS (mm)'
        write(40, line_format) tot_reach_aet


        write (line_format, '(A,I0, A)') '(',size(q, 1),'(1X, f0.8), A)'

        do i = 1, size(q, 2)
            write(41, line_format) q(:, i)
        end do

    end subroutine

end module dyna_write_output
