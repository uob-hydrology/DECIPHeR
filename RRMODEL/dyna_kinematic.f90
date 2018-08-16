module dyna_kinematic
contains
!                                                                       
!===============================================================        
!  ROUTINE FOR CALCULATING THE KINEMATIC SOLUTION                       
!  EXPONENTIAL FORMULATION                                              
!===============================================================        
!                                                                       
      subroutine kinematic (nac, &
      acc, &
      dtt, &
      dyna_hru, &
      ia, &
      it, &
      ntt, &
      smax, &
      szm, &
      t0dt, &
      wt)
          use dyna_common_types
          implicit none
! Argument Declares
integer :: nac
          double precision :: acc
          double precision :: dtt
          type(dyna_hru_type) :: dyna_hru(nac)
          integer :: ia
          integer :: it
          integer :: ntt
          double precision, dimension(:) :: smax
          double precision, dimension(:) :: szm
          double precision, dimension(:) :: t0dt
          double precision :: wt

! Local Declares
          double precision :: c
          double precision :: const1
          double precision :: ct1
          integer :: iter
          integer :: itt
          double precision :: qbf1
          double precision :: qin1
          double precision :: qin2
          double precision :: qit
          double precision :: rhs
          double precision :: storebal
          double precision :: storeend
          double precision :: storein
          double precision :: storeini
          double precision :: storeout
! End declares

          !  QINt1(IA) = PREVIOUS TIME STEP INPUTS QIN(IA) - or QIN1              
          !  QIN(IA) = ALL DOWNSLOPE QBF(IA) FLUXES * dyna_hru(ia)%ac summed for all increm
          !  to all downslope TotanB increments from previous time step outputs   
          !  QBFt1(IA) = PREVIOUS TIMESTEPS OUTPUT FROM EACH ELEMENT QBF(IA) - or 
          !  QBF(IA) = INITIALLY PREVIOUS TIMESTEPS OUTFLOW BUT IS ITERATIVELY SOL
          !  IN THE INNER LOOP                                                    
          qin1 = dyna_hru(ia)%qint1 / dyna_hru(ia)%ac 
          qbf1 = dyna_hru(ia)%qbft1 
                                                                        
          storeini = dyna_hru(ia)%sd 

          if (dyna_hru(ia)%sd .lt.0.D0) print * , ' SD lt zero ', it, ia, dyna_hru(ia)%sd 
          storein = 0.D0 
          storeout = 0.D0 
          !                                                                       
          !  SOLUTION OF THE DOWNSLOPE KINEMATIC WAVE EQUATIONS - MAIN LOOP       
          !

          do 30 itt = 1, ntt 
          !                                                                       
          !  Initialise implicit solution for QBF. This form of wave speed
          !  EQUATION IS FOR AN EXPONENTIAL TRANSMISSIVITY PROFILE                
          !                                                                       
          !  NB DX IS TREATED IMPLICITLY BY USING Q AS FLUX PER UNIT AREA         
          !                                                                       
             ct1 = - 0.5D0 * (qbf1 + qin1) / SZM (dyna_hru(ia)%ipar )
             const1 = (1.D0 - wt) * ct1 * dtt * (qbf1 - qin1 - dyna_hru(ia)%uz ) 

          !                                                                       
          !  LINEAR INTERPOLATION FOR INFLOWS FOR INNER TIME LOOP                 
          !                                                                       
             qin2 = (dyna_hru(ia)%qint1 + dble (itt) * (dyna_hru(ia)%qin - dyna_hru(ia)%qint1 )     &
             / dble (ntt) ) / dyna_hru(ia)%ac                                       
          !                                                                       
          !  START ITERATION LOOP FOR SOLUTION FOR QBF                            
          !  NB THIS IS EQUIVALENT TO THE NEWTON-RAPHSON (LINEAR IN QBF)          
          !                                                                       
          !  ITERATE UNTIL SOLUTION CONVERGES ACCORDING TO ACC CRITERIA           
          !                                                                       
             do 10 iter = 1, 30 
                                                                        
                c = - 0.5D0 * (dyna_hru(ia)%qbf + qin2) / SZM (dyna_hru(ia)%ipar )
                rhs = qbf1 + wt * c * dtt * ( ( - 1.D0 * qin2) - dyna_hru(ia)%uz )  &
                + const1                                                    
                qit = dyna_hru(ia)%qbf 
                dyna_hru(ia)%qbf = rhs / (1.D0 - wt * c * dtt) 
                if (dabs (dyna_hru(ia)%qbf - qit) .lt.acc) goto 20
          !                                                                       
          !  CHECK -1 FOR NEGATIVE - JEF                                          
          !                                                                       
                if (dyna_hru(ia)%qbf .lt.0.D0) dyna_hru(ia)%qbf = 0.D0 

          10    end do 
          !                                                                       
          !  THE SOLUTION IS FOUND: CONTINUE                                      
          !                                                                       
20        continue 
          !                                                                       
          !  CHECK THAT THE MAX DEFICIT HAS NOT BEEN REACHED, IF SO NO DOWNSLOPE  
          !  FLOW (DEPENDENT ON SMAX)                                             
          !                                                                       
             if (dyna_hru(ia)%sd .gt.Smax (dyna_hru(ia)%ipar ) ) then
                dyna_hru(ia)%qbf = 0.D0 
             endif 
          !                                                                       
          !  CHECK FOR QBF GT THAN POSSIBLE MAX SATURATED FLOW (LN(A/TANB)*EXP(T0)
          !  IN THIS CASE THE EXTRA FLOW IS ADDED INTO EX AND QBF IS SET TO THE   
          !  MAXIMUM SATURATED FLOW. THEREFORE SD(IA) IS SET TO ZERO (SATURATED)  
          !                                                                       
             dyna_hru(ia)%sd_change = 0.D0 
                                                                        
             if (dyna_hru(ia)%qbf .gt. (dyna_hru(ia)%st * dexp (T0DT (dyna_hru(ia)%ipar ) ) ) ) then
          !                                                                       
          !  NOTE IF THIS WAS TO RUN LIKE THE OLDE TOPMOD, EX WOULD NOT RESULT IN 
          !  EXCESS FLOW AND A NEGATIVE SD WOULD BE KEPT TO MAINTAIN A SATURATED  
          !  AREA (comment out ex line and SD(IA) = 0.0 line to change this)      
          !                                                                       
             dyna_hru(ia)%qbf = dyna_hru(ia)%st * dexp (T0DT (dyna_hru(ia)%ipar ) )
             endif 
          !                                                                       
          !  CALCULATE THE CHANGES TO SD WITH THE INFLOWS AND OUTFLOWS            
          !                                                                       
             dyna_hru(ia)%sd = dyna_hru(ia)%sd + dtt * ( ( - 1.D0 * dyna_hru(ia)%uz ) + dyna_hru(ia)%qbf    &
             - dyna_hru(ia)%qin / dyna_hru(ia)%ac )                                         
                                                                        
             storein = storein + dtt * ( ( - 1.D0 * dyna_hru(ia)%uz ) - dyna_hru(ia)%qin    &
             / dyna_hru(ia)%ac )                                                    
                                                                        
             dyna_hru(ia)%sd_change = dyna_hru(ia)%sd_change + ( ( - 1.D0 * dyna_hru(ia)%uz )       &
             + (dtt * (dyna_hru(ia)%qbf - dyna_hru(ia)%qin ) ) / dyna_hru(ia)%ac )                  
          !                                                                       
          !  END OF THE MAIN DOWNSLOPE KINEMATIC WAVE SOLUTION                    
          !                                                                       
          30 end do 
          !                                                                       
          !  DO NOT ALLOW QBF TO REDUCE SD BELOW SMAX                             
          !                                                                       
          if (dyna_hru(ia)%sd .gt.Smax (dyna_hru(ia)%ipar ) ) then
             dyna_hru(ia)%qbf = dyna_hru(ia)%qbf - (dyna_hru(ia)%sd - Smax (dyna_hru(ia)%ipar ) ) / dtt
             !if (dyna_hru(ia)%qbf .lt.0.D0) print * , ' WOW dyna_hru(IA)%qbf lt zero '
             dyna_hru(ia)%sd = Smax (dyna_hru(ia)%ipar )
          endif 
          storeout = storeout + dyna_hru(ia)%qbf * dtt 
          !                                                                       
          !  IF SD IS LESS THAN ZERO CALCULATE THE SATURATION EXCESS EXS          
          !                                                                       
          if (dyna_hru(ia)%sd .lt.0.D0) then 
             dyna_hru(ia)%exs = 0.D0 - dyna_hru(ia)%sd 
             dyna_hru(ia)%sumexs = dyna_hru(ia)%sumexs + (dyna_hru(ia)%exs * dyna_hru(ia)%ac)
             dyna_hru(ia)%sd = 0.D0 
             storeout = storeout + dyna_hru(ia)%exs 
          endif 
                                                                        
          qin1 = qin2 
          qbf1 = dyna_hru(ia)%qbf 
          !                                                                       
          storeend = dyna_hru(ia)%sd 

          !                                                                       
          !  NOTE STOREBAL IS A DEFICIT CHANGE, STOREIN IS SUMMED AS A NEGATIVE   
          !                                                                       
          storebal = storeini + storein + storeout - storeend 
          if (storebal.gt.0.0000000001D0.or.storebal.lt. - 0.000000001D0)   &
          then                                                              

             print * , ' Storebal SD wrong KINEMATIC ', it, ia, storebal,   &
             storeini, storein, storeout, storeend, dyna_hru(ia)%uz , dyna_hru(ia)%qin , &
             dyna_hru(ia)%ac , Smax (dyna_hru(ia)%ipar ) , dtt, ntt
          endif 
      end subroutine

end module dyna_kinematic
