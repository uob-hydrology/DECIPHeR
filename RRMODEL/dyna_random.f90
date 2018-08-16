module dyna_random
! Toby Dunne 28/01/2016 Random number generator module merged
! from ranaru and ranin
! Module created removes need for common block when using these functions
      implicit none
!     ..
!     .. SCALARS PRIVATE TO MODULE ..
!
      doubleprecision, save :: cc, cd, cm
      integer, save :: i97, j97
!     ..
!     .. ARRAYS PRIVATE TO MODULE ..
!
      double precision, save :: cu (97)

      private :: cc, cd, cm, i97, j97, cu

contains
!                                                                       
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     *** RANARU ***                                                   c
!                                                                      c
!     CREATES A REAL ARRAY ARU OF LENGTH LEN.                          c
!     THE VALUES IN ARU ARE UNIFORMLY DISTRIBUTED IN THE RANGE         c
!     0.0 <= ARU(I) < 1.0.                                             c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                       
      subroutine ranaru (aru, length)
      implicit none
!                                                                       
!     .. SCALAR ARGUMENTS ..                                            
!                                                                       
      integer length
!     ..                                                                
!     .. ARRAY ARGUMENTS ..                                             
!                                                                       
      doubleprecision aru ( * ) 

!     ..                                                                
!     .. LOCAL SCALARS ..                                               
!                                                                       
      doubleprecision uni 
      integer ivec 
!     ..                                                                
!     .. COMMON BLOCKS ..                                               
!                                                                       
      !common / ranset / cu, cc, cd, cm, i97, j97
!     ..                                                                
      do 100 ivec = 1, length
         uni = cu (i97) - cu (j97) 
         if (uni.lt.0.0d0) uni = uni + 1.0d0 
         cu (i97) = uni 
         i97 = i97 - 1 
         if (i97.eq.0) i97 = 97 
         j97 = j97 - 1 
         if (j97.eq.0) j97 = 97 
         cc = cc - cd 
         if (cc.lt.0.0d0) cc = cc + cm 
         uni = uni - cc 
         if (uni.lt.0.0d0) uni = uni + 1.0d0 
         aru (ivec) = uni 
  100 end do 
!                                                                       
      return 
      end subroutine ranaru                         

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!     *** RANIN ***                                                    c
!                                                                      c
!     INITIALISES THE RANDOM NUMBER GENERATOR USING TWO INTEGER SEEDS  c
!     IJ AND KL,  WHERE 0 <= IJ <= 31328  AND  0 <= KL <= 30081.       c
!     IT IS ONLY CALLED ONCE IN A NORMAL RUN, BUT IS NOT CALLED IF     c
!     A RUN IS RESTARTED USING RANRST.                                 c
!                                                                      c
!     PERMUTATIONS OF THESE SEEDS GIVE 900 MILLION SEQUENCES OF        c
!     PSEUDO-RANDOM VALUES, EACH OF LENGTH 10**30, AND EACH OF WHICH   c
!     IS SAID TO BE INDEPENDENT OF THE OTHERS.                         c
!                                                                      c
!     THIS RANDOM NUMBER GENERATOR IS DUE TO MARSAGLIA AND ZAMAN(1987) c
!     WITH SLIGHT MODIFICATIONS BY JAMES(1988).                        c
!     IT IS PORTABLE ACROSS MOST 32-BIT MACHINES.                      c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     IT IS RECOMMENDED THAT THE COMMON /RANSET/ .... STATEMENT SHOULD
!     BE INCLUDED IN THE USER'S MAIN PROGRAM UNIT.
!
      subroutine ranin (ij, kl)
!
!     .. SCALAR ARGUMENTS ..
!
      integer ij, kl
!     ..
!     .. LOCAL SCALARS ..
!
      doubleprecision s, t
      integer i, ii, j, jj, k, l, m
!     ..
!     .. INTRINSIC FUNCTIONS ..
!
      intrinsic mod
!     ..
!     .. COMMON BLOCKS ..
!
!      common / ranset / cu, cc, cd, cm, i97, j97
!     ..
      i = mod (ij / 177, 177) + 2
      j = mod (ij, 177) + 2
      k = mod (kl / 169, 178) + 1
      l = mod (kl, 169)
      do 2 ii = 1, 97
         s = 0.0d0
         t = 0.5d0
         do 3 jj = 1, 24
            m = mod (mod (i * j, 179) * k, 179)
            i = j
            j = k
            k = m
            l = mod (53 * l + 1, 169)
            if (mod (l * m, 64) .ge.32) s = s + t
            t = 0.5d0 * t
    3    continue
         cu (ii) = s
    2 continue
      cc = 362436.d0 / 16777216.d0
      cd = 7654321.d0 / 16777216.d0
      cm = 16777213.d0 / 16777216.d0
      i97 = 97
      j97 = 33
!
      return
      end subroutine ranin

end module dyna_random
