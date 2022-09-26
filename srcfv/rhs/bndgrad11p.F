       !  coef are obtained with sympy for FD scheme
       j = 5
       !{a8: -3/56, a1: -1/63, a5: 6/5, a7: 4/21, a2: 3/28, a9: 1/105, a10: -1/1260, a6: -1/2, a3: -4/7, a0: 1/840, a4: -11/30}
       a0  =   1.d0/840.d0
       a1  = - 1.d0/63.d0
       a2  =   3.d0/28.d0
       a3  = - 4.d0/7.d0
       a4  = -11.d0/30.d0
       a5  =   6.d0/5.d0
       a6  = - 0.5d0
       a7  =   4.d0/21.d0
       a8  = - 3.d0/56.d0
       a9  =   1.d0/105.d0
       a10 = - 1.d0/1260.d0

!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_11pi.F"
#include "rhs/gradop_prim_11pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"
       
       j = 4
       ! {a8: 3/40, a1: 1/24, a5: -21/20, a7: -1/4, a2: -3/8, a9: -1/72, a10: 1/840, a6: 7/12, a3: -319/420, a0: -1/360, a4: 7/4}
       ! FD37
       a0  = -  1.d0/360.d0
       a1  =    1.d0/24.d0
       a2  = -  3.d0/8.d0
       a3  = -319.d0/420.d0
       a4  =    1.75d0
       a5  = - 21.d0/20.d0
       a6  =    7.d0/12.d0
       a7  = -  0.25d0
       a8  =    0.075d0
       a9  = -  1.d0/72.d0
       a10 =    1.d0/840.d0
       
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_11pi.F"
#include "rhs/gradop_prim_11pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"


       j = 3
       ! {a8: -1/6, a1: -2/9, a5: 28/15, a7: 8/15, a2: -341/280, a9: 2/63, a10: -1/360, a6: -7/6, a3: 8/3, a0: 1/90, a4: -7/3}
       ! FD28
       a0=    1.d0/90.d0
       a1= -  2.d0/9.d0
       a2= -341.d0/280.d0
       a3=    8.d0/3.d0
       a4= -  7.d0/3.d0
       a5=   28.d0/15.d0
       a6= -  7.d0/6.d0
       a7=    8.d0/15.d0
       a8= -  1.d0/6.d0
       a9=    2.d0/63.d0
       a10= - 1.d0/360.d0
       
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_11pi.F"
#include "rhs/gradop_prim_11pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"
       
       j = 2
       ! {a8: 9/14, a1: -4609/2520, a5: -63/10, a7: -2, a2: 9/2, a9: -1/8, a10: 1/90, a6: 21/5, a3: -6, a0: -1/10, a4: 7}
       ! FD19
       a0= -   0.1d0
       a1= -4609.d0/2520.d0
       a2=     4.5d0
       a3= -   6.d0
       a4=     7.d0
       a5= -  6.3d0
       a6=    21.d0*0.2d0
       a7= -   2.d0
       a8=     9.d0/14.d0
       a9= -   0.125d0
       a10=    1.d0/90.d0
       
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_11pi.F"
#include "rhs/gradop_prim_11pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"
       
       j = 1 
       
       !  !off-centered grad at first cell center FD010
       ! using off centered fluxes expressions of DNC9 for velo (because vel_wall=0) 
       !{a0: -7381/2520, a1: 10, a2: -45/2, a3: 40, a4: -105/2, a5: 252/5, a6: -35, a7: 120/7, a8: -45/8, a9: 10/9, a10: -1/10}
     
       !a0  =- 7381.d0/2520.d0
       !a1  =  10.d0
       !a2  =- 22.5d0
       !a3  =  40.d0
       !a4  =- 52.5d0
       !a5  =  50.4d0
       !a6  =- 35.d0
       !a7  = 120.d0/7.d0
       !a8  =- 11.25d0
       !a9  =  10.d0/9.d0
       !a10 =-  0.1d0
         ! from Berland JCP 2009  (drp scheme o4 optimon 11p)
       a0  = - 2.391602219538d0
       a1  =   5.832490322294d0
       a2  = - 7.650218001182d0
       a3  =   7.907810563576d0
       a4  = - 5.922599052629d0
       a5  =   3.071037015445d0
       a6  = - 1.014956769726d0
       a7  =   0.170022256519d0
       a8  =   0.002819958377d0
       a9  = - 0.004791009708d0
       a10 = - 0.000013063429d0
       
       denom = 1.d0/2520.d0
       c3d0  =    230.d0 * denom
       c3d1  =   5081.d0 * denom
       c3d2  = - 7469.d0 * denom
       c3d3  =  11281.d0 * denom
       c3d4  = -13619.d0 * denom
       c3d5  =  12421.d0 * denom
       c3d6  = - 8327.d0 * denom
       c3d7  =   3973.d0 * denom
       c3d8  = - 1277.d0 * denom
       c3d9  =    248.d0 * denom
       c3d10 = -   26.d0 * denom
           
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1,im
#include "rhs/gradop_prim_11pi.F"
#include "rhs/gradop_prim_11pj_offcentered1stcell.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"
       
      
  
