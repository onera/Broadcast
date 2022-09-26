       !  coef are obtained with sympy for FD scheme
       j = 4
       ! {a8: 1/280, a1: 1/14, a5: -1/2, a7: -1/28, a2: -1/2, a6: 1/6, a3: -9/20, a0: -1/168, a4: 5/4}
       ! FD35
       a0  = - 1.d0/168.d0
       a1  =   1.d0/14.d0
       a2  = - 0.5d0
       a3  = - 9.d0*0.05d0
       a4  =   1.25d0
       a5  = - 0.5d0
       a6  =   1.d0/6.d0
       a7  = - 1.d0/28.d0
       a8  =   1.d0/280.d0
       
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_9pi.F"
#include "rhs/gradop_prim_9pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"
       
       j = 3
       ! {a8: -1/168, a1: -2/7, a5: 2/3, a7: 2/35, a2: -19/20, a6: -1/4, a3: 2, a0: 1/56, a4: -5/4}
       ! FD26
       a0=    1.d0/56.d0
       a1= -  2.d0/7.d0
       a2= - 19.d0*0.05d0
       a3=    2.d0
       a4= -  1.25d0
       a5=    2.d0/3.d0
       a6= -  0.25d0
       a7=    2.d0/35.d0
       a8= -  1.d0/168.d0

!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_9pi.F"
#include "rhs/gradop_prim_9pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"
       
       j = 2
       ! {a8: 1/56, a1: -223/140, a5: -7/4, a7: -1/6, a2: 7/2, a6: 7/10, a3: -7/2, a0: -1/8, a4: 35/12}
       ! FD17
       a0= -   0.125d0
       a1= - 223.d0/140.d0
       a2=     3.5d0
       a3= -   3.5d0
       a4=    35.d0/12.d0
       a5= -   1.75d0
       a6=     0.7d0
       a7= -   1.d0/6.d0
       a8=     1.d0/56.d0
       
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1, im
#include "rhs/gradop_prim_9pi.F"
#include "rhs/gradop_prim_9pj_offcentered.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh_i.F"

       j = 1 ! FD08 for grad thermo
      ! using off centered fluxes expressions of DNC7 for velo (because vel_wall=0) 
      !{a0: -761/280, a1: 8, a2: -14, a3: 56/3, a4: -35/2, a5: 56/5, a6: -14/3, a7: 8/7, a8: -1/8}
    
    
      ! a0= -761.d0/280.d0
      ! a1= 8.d0
      ! a2= -14.d0
      ! a3= 56.d0/3.d0
      ! a4= -17.5d0
      ! a5=  11.2d0
      ! a6= -14.d0/3.d0
      ! a7= 8.d0/7.d0
      ! a8= -0.125d0
      ! from Berland JCP 2009  (drp scheme o4 optimon 11p)
     
       a0= - 2.225833963270d0
       a1=   4.827779580575d0
       a2= - 5.001388453836d0
       a3=   3.911103941646d0
       a4= - 2.115267458633d0
       a5=   0.718882784412d0
       a6= - 0.115276430895d0
       a7= ZERO
       a8= ZERO
       
       
       c3d0  =    92.d0 * denom
       c3d1  =  1547.d0 * denom
       c3d2  = -1861.d0 * denom
       c3d3  =  2171.d0 * denom
       c3d4  = -1917.d0 * denom
       c3d5  =  1191.d0 * denom
       c3d6  = - 489.d0 * denom
       c3d7  =   119.d0 * denom
       c3d8  = -  13.d0 * denom
           
!$AD II-LOOP
!DIR$ IVDEP      
       do i = 1,im
#include "rhs/gradop_prim_9pi.F"
#include "rhs/gradop_prim_9pj_offcentered1stcell.F"
#include "rhs/gradient_prim.F"
       enddo
#include "rhs/gradingh.F"  
