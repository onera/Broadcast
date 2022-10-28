    !1st direction
    
    k_sensor1 = ABS(p(i-1,j) - TWO*p(i,j) + p(i+1,j)) / &
                ABS(p(i-1,j) + TWO*p(i,j) + p(i+1,j))
                
    k_sensor2 = ABS(p(i-2,j) - TWO*p(i-1,j) + p(i,j)) / &
                ABS(p(i-2,j) + TWO*p(i-1,j) + p(i,j))
    
                
    divu      = (gradu(i,j,1)+gradv(i,j,2))            
    divu2     = divu * divu
    vort2     = (gradv(i,j,1)-gradu(i,j,2)) * (gradv(i,j,1)-gradu(i,j,2))
    ducros1   = divu2/(divu2+vort2+1d-15)
    ! dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sq+1.d-15)*divu) )
    dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sqrt(c2r*nx2)+1.d-15)*divu) )
 
     
    divu      = (gradu(i-1,j,1)+gradv(i-1,j,2))
    divu2     = divu * divu
    vort2     = (gradv(i-1,j,1)-gradu(i-1,j,2)) * (gradv(i-1,j,1)-gradu(i-1,j,2))
    ducros2   = divu2/(divu2+vort2+1d-15)
    ! dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i-1,j)/(sq+1.d-15)*divu) )                            
    dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i-1,j)/(sqrt(c2l*nx2)+1.d-15)*divu) )                            
    
    !coef      = MAX( k_sensor1*ducros1*dxm1, &
    !                 k_sensor2*ducros2*dxm2)
    
    coef      = max(k_sensor1, k_sensor2) * max(ducros1, ducros2) * max(dxm1, dxm2)
   
!     ! extension of shock capturing

    
!   
    !k_sensor1 = ABS(p(i,j) - TWO*p(i+1,j) + p(i+2,j)) / &
    !            ABS(p(i,j) + TWO*p(i+1,j) + p(i+2,j))
    !          
    !k_sensor2 = ABS(p(i-3,j) - TWO*p(i-2,j) + p(i-1,j)) / &
    !            ABS(p(i-3,j) + TWO*p(i-2,j) + p(i-1,j))
    !          
    !divu      = gradu(i+1,j,1)+gradv(i+1,j,2)
    !divu2     = divu*divu
    !vort2     = (gradv(i+1,j,1)-gradu(i+1,j,2)) * (gradv(i+1,j,1)-gradu(i+1,j,2))
    !ducros1   = divu2/(divu2+vort2+1d-15)
    !dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i+1,j)/(sq+1.d-15)*divu) )
!
  !
    !divu      = gradu(i-2,j,1)+gradv(i-2,j,2)
    !divu2     = divu*divu
    !vort2     = (gradv(i-2,j,1)-gradu(i-2,j,2)) * (gradv(i-2,j,1)-gradu(i-2,j,2))
    !ducros2   = divu2/(divu2+vort2+1d-15)
    !dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i-2,j)/(sq+1.d-15)*divu) )
    
    !coef      = MAX( k_sensor1*ducros1*dxm1, &
    !                 k_sensor2*ducros2*dxm2, coef)
    
    !coef      = MAX(sensor(i-1,j  ,1), sensor(i,j  ,1), &
    !                sensor(i-1,j-1,1), sensor(i,j-1,1), &
    !                sensor(i-1,j+1,1), sensor(i,j+1,1)  )
    !
    !coef      = MAX(sensor(i-3,j,1), sensor(i-2,j,1), &
    !                sensor(i-1,j,1), sensor(i  ,j,1), &
    !                sensor(i+1,j,1), sensor(i+2,j,1)  )
    !
    !coef      = MAX(sensor(i-2,j,1), &
    !                sensor(i-1,j,1), sensor(i  ,j,1), &
    !                sensor(i+1,j,1))
    !
    !rspec      = rconv(i,j,1)
    
!    
    eps2      = k2*coef
    eps4      = MAX(ZERO,k4-eps2 * 12.d0) ! to follow Sciacovelli CF 2021

    diffro   = HALF * (w(i,j,1) - w(i-1,j,1))
    diffrou  = HALF * (w(i,j,2) - w(i-1,j,2))
    diffrov  = HALF * (w(i,j,3) - w(i-1,j,3))
    diffrow  = HALF * (w(i,j,4) - w(i-1,j,4))
    diffroe  = HALF * (w(i,j,5) - w(i-1,j,5))
    !
