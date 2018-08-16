module dyna_root_zone1
contains
!                                                                       
!===============================================================        
!  ROUTINES FOR CALCULATING THE ROOT ZONE FUNCTIONS                     
!  ROOT_ZONE1 = BASIC ROOT ZONE FORMULATION                             
!===============================================================        
!                                                                       
      subroutine root_zone1 (nac, &
      dyna_hru, &
      ia, &
      it, &
      srmax)
          use dyna_common_types
          implicit none

! Argument Declares
          integer :: nac
          type(dyna_hru_type) :: dyna_hru(nac)
          integer :: ia
          integer :: it
          double precision, dimension(:) :: srmax

! Local Declares
          double precision :: storebal
          double precision :: storeend
          double precision :: storein
          double precision :: storeini
          double precision :: storeout
! End declares

          dyna_hru(ia)%pex = 0.D0 
          storein = dyna_hru(ia)%p
          storeini = dyna_hru(ia)%srz 
          storeout = 0.D0 

          !  Explicit ROOT ZONE CALCULATIONS given P in this time step

          dyna_hru(ia)%srz = dyna_hru(ia)%srz + dyna_hru(ia)%p
          dyna_hru(ia)%sum_pptin = dyna_hru(ia)%sum_pptin + dyna_hru(ia)%p

          if (dyna_hru(ia)%srz .gt.SRmax (dyna_hru(ia)%ipar ) ) then

             dyna_hru(ia)%pex = (dyna_hru(ia)%srz - SRmax (dyna_hru(ia)%ipar ) )
             dyna_hru(ia)%srz = SRmax (dyna_hru(ia)%ipar)

          endif 

          storeout = dyna_hru(ia)%pex 

          !  CALCULATE EVAPOTRANSPIRATION FROM ROOT ZONE DEFICIT                  

          dyna_hru(ia)%ea = 0

          if (dyna_hru(ia)%ep .gt.0) then 

              dyna_hru(ia)%ea = dyna_hru(ia)%ep * (dyna_hru(ia)%srz / SRMAX (dyna_hru(ia)%ipar ) )

              if (dyna_hru(ia)%ea.gt.dyna_hru(ia)%ep ) dyna_hru(ia)%ea = dyna_hru(ia)%ep
              if (dyna_hru(ia)%ea.ge.dyna_hru(ia)%srz ) dyna_hru(ia)%ea = dyna_hru(ia)%srz

              dyna_hru(ia)%srz = dyna_hru(ia)%srz - dyna_hru(ia)%ea
              dyna_hru(ia)%sum_pein = dyna_hru(ia)%sum_pein + dyna_hru(ia)%ep
              dyna_hru(ia)%sum_aeout = dyna_hru(ia)%sum_aeout + dyna_hru(ia)%ea

              storeout = storeout + dyna_hru(ia)%ea

          endif 

          ! Water balance calculations

          storeend = dyna_hru(ia)%srz 
          storebal = storeini + storein - storeout - storeend 
          if (storebal.gt.0.0000000001.or.storebal.lt. - 0000000001) then 
             print * , ' Storebal SRZ wrong ', it, ia, storebal 
             STOP
          endif 

      end subroutine

end module
