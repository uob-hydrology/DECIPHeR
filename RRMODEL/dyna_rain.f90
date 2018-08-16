module dyna_rain
contains
!
!===============================================================
!  ROUTINE TO DISTRBITUE RAIN AND PET TO EACH HRU
!===============================================================
!
      subroutine rain (nac, &
      dyna_hru, &
      it, &
      pe_step, &
      r_gau_step)

          use dyna_common_types
          implicit none
! Argument Declares
          integer :: nac
          type(dyna_hru_type) :: dyna_hru(nac)
          integer :: it
          double precision, dimension(:,:) :: pe_step
          double precision, dimension(:,:) :: r_gau_step

! Local Declares
          integer :: ia
! End declares

          do ia = 1, nac
             dyna_hru(ia)%p = r_gau_step (dyna_hru(ia)%ippt, it)
             dyna_hru(ia)%ep = pe_step (dyna_hru(ia)%ipet, it)
          end do

      end subroutine

end module dyna_rain
