    
! 2nd direction

    k_sensor1 = ABS(p(i,j-1) - TWO*p(i,j) + p(i,j+1)) / &
                ABS(p(i,j-1) + TWO*p(i,j) + p(i,j+1))
                
    k_sensor2 = ABS(p(i,j-2) - TWO*p(i,j-1) + p(i,j)) / &
                ABS(p(i,j-2) + TWO*p(i,j-1) + p(i,j))
                
!     ducros1 is done in dissipation_ducros_x 

    divu      = gradu(i,j-1,1)+gradv(i,j-1,2)                           
    divu2     = divu * divu
    vort2     = (gradv(i,j-1,1)-gradu(i,j-1,2)) * (gradv(i,j-1,1)-gradu(i,j-1,2))
    ducros2   = divu2/(divu2+vort2+1d-15)
    
      !sq = sound/dx
    dxm2       = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sq+1d-15)*divu) )              
    
    coef      = MAX( k_sensor1*ducros1*dxm1, &
                     k_sensor2*ducros2*dxm2)
    
     ! extension of shock capturing
    
!    k_sensor1 = ABS(p(i,j) - TWO*p(i,j+1) + p(i,j+2)) / &
!                ABS(p(i,j) + TWO*p(i,j+1) + p(i,j+2))
!                
!    k_sensor2 = ABS(p(i,j-3) - TWO*p(i,j-2) + p(i,j-1)) / &
!                ABS(p(i,j-3) + TWO*p(i,j-2) + p(i,j-1))
!                
!    divu2     = (gradu(i,j+1,1)+gradv(i,j+1,2)) * (gradu(i,j+1,1)+gradv(i,j+1,2))
!    vort2     = (gradv(i,j+1,1)-gradu(i,j+1,2)) * (gradv(i,j+1,1)-gradu(i,j+1,2))
!    ducros1   = divu2/(divu2+vort2+1d-15)
!    
!    divu2     = (gradu(i,j-2,1)+gradv(i,j-2,2)) * (gradu(i,j-2,1)+gradv(i,j-2,2))
!    vort2     = (gradv(i,j-2,1)-gradu(i,j-2,2)) * (gradv(i,j-2,1)-gradu(i,j-2,2))
!    ducros2   = divu2/(divu2+vort2+1d-15)
!      
!    coef      = MAX( k_sensor1*ducros1, &
!                     k_sensor2*ducros2, coef)
    
    eps2      = k2*coef
    eps4      = MAX(ZERO,k4 - 12.d0 * eps2)   
    

    diffro   = HALF * (w(i,j,1) - w(i,j-1,1))
    diffrou  = HALF * (w(i,j,2) - w(i,j-1,2))
    diffrov  = HALF * (w(i,j,3) - w(i,j-1,3))
    diffrow  = HALF * (w(i,j,4) - w(i,j-1,4))
    diffroe  = HALF * (w(i,j,5) - w(i,j-1,5))
    !
