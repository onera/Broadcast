       !  coef are obtained with sympy for FD scheme
        j = 3
        ! {a1: -2/5, a5: 2/15, a2: -7/12, a6: -1/60, a3: 4/3, a0: 1/30, a4: -1/2}
        ! FD24
        a0=   1.d0/30.d0
        a1= - 2.d0/5.d0
        a2= - 7.d0/12.d0
        a3=   4.d0/3.d0
        a4= - 0.5d0
        a5=   2.d0/15.d0
        a6= - 1.d0/60.d0
        
!$AD II-LOOP
!DIR$ IVDEP      
        do i = 1, im
#include "rhs/gradop_prim_7pi.F"
#include "rhs/gradop_prim_7pj_offcentered.F"
#include "rhs/gradient_prim.F"
        enddo
#include "rhs/gradingh_i.F"

        j = 2       
        ! {a1: -77/60, a5: -1/4, a2: 5/2, a6: 1/30, a3: -5/3, a0: -1/6, a4: 5/6}
        ! FD15
        a0= - 1.d0/6.d0
        a1= -77.d0/60.d0
        a2=   2.5d0
        a3= - 5.d0/3.d0
        a4=   5.d0/6.d0
        a5= - 0.25d0
        a6=   1.d0/30.d0
        
!$AD II-LOOP
!DIR$ IVDEP      
        do i = 1, im
#include "rhs/gradop_prim_7pi.F"
#include "rhs/gradop_prim_7pj_offcentered.F"
#include "rhs/gradient_prim.F"
        enddo
#include "rhs/gradingh_i.F"

        j = 1
        ! FD06 for grad on thermo
        ! Grad on Velocities is deducted from vel = 0 at wall 
        !{a0: -49/20, a1: 6, a2: -15/2, a3: 20/3, a4: -15/4, a5: 6/5, a6: -1/6}
        !a0 = -  49.d0*0.05d0
        !a1 =     6.d0
        !a2 = -  7.5d0
        !a3 =    20.d0/3.d0
        !a4 = - 3.75d0
        !a5 =    1.2d0
        !a6 = -  1.d0/6.d0
        ! from Berland JCP 2009  (drp scheme o4 optimon 7p)
        a0= - 2.225833963270d0
        a1=   4.827779580575d0
        a2= - 5.001388453836d0
        a3=   3.911103941646d0
        a4= - 2.115267458633d0
        a5=   0.718882784412d0
        a6= - 0.115276430895d0
        
        c3d0  =  9.d0 * denom
        c3d1  = 93.d0 * denom
        c3d2  =-78.d0 * denom
        c3d3  = 57.d0 * denom
        c3d4  =-28.d0 * denom
        c3d5  =  8.d0 * denom
        c3d6  =- 1.d0 * denom
      
!$AD II-LOOP
!DIR$ IVDEP      
        do i = 1,im
#include "rhs/gradop_prim_7pi.F"
#include "rhs/gradop_prim_7pj_offcentered1stcell.F"
#include "rhs/gradient_prim.F"
#include "rhs/gradingh_j.F"
        enddo
#include "rhs/gradingh_i.F"
!
