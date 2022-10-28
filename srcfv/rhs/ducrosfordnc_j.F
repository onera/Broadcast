    
! 2nd direction

    k_sensor1 = ABS(p(i,j-1) - TWO*p(i,j) + p(i,j+1)) / &
                ABS(p(i,j-1) + TWO*p(i,j) + p(i,j+1))
                
    k_sensor2 = ABS(p(i,j-2) - TWO*p(i,j-1) + p(i,j)) / &
                ABS(p(i,j-2) + TWO*p(i,j-1) + p(i,j))
                
!     ducros1 is done in dissipation_ducros_x 

    divu      = (gradu(i,j,1)+gradv(i,j,2))            
    divu2     = divu * divu
    vort2     = (gradv(i,j,1)-gradu(i,j,2)) * (gradv(i,j,1)-gradu(i,j,2))
    ducros1   = divu2/(divu2+vort2+1d-15)
    ! dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sq+1.d-15)*divu) )
    dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sqrt(c2r*nx2)+1.d-15)*divu) )

    divu      = gradu(i,j-1,1)+gradv(i,j-1,2)                           
    divu2     = divu * divu
    vort2     = (gradv(i,j-1,1)-gradu(i,j-1,2)) * (gradv(i,j-1,1)-gradu(i,j-1,2))
    ducros2   = divu2/(divu2+vort2+1.d-15)
    
      !sq = sound/dx
    ! dxm2       = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j-1)/(sq+1.d-15)*divu) )              
    dxm2       = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j-1)/(sqrt(c2l*nx2)+1.d-15)*divu) )              
    
    !coef      = MAX( k_sensor1*ducros1*dxm1, &
    !                 k_sensor2*ducros2*dxm2)
    
    coef      = max(k_sensor1, k_sensor2) * max(ducros1, ducros2) * max(dxm1, dxm2)

    
    ! extension of shock capturing
    
    !k_sensor1 = ABS(p(i,j) - TWO*p(i,j+1) + p(i,j+2)) / &
    !            ABS(p(i,j) + TWO*p(i,j+1) + p(i,j+2))
    !          
    !k_sensor2 = ABS(p(i,j-3) - TWO*p(i,j-2) + p(i,j-1)) / &
    !            ABS(p(i,j-3) + TWO*p(i,j-2) + p(i,j-1))
    !          
    !divu      = gradu(i,j+1,1)+gradv(i,j+1,2)
    !divu2     = divu*divu
    !vort2     = (gradv(i,j+1,1)-gradu(i,j+1,2)) * (gradv(i,j+1,1)-gradu(i,j+1,2))
    !ducros1   = divu2/(divu2+vort2+1d-15)
    !dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j+1)/(sq+1.d-15)*divu) )
!
  !
    !divu      = gradu(i,j-2,1)+gradv(i,j-2,2)
    !divu2     = divu*divu
    !vort2     = (gradv(i,j-2,1)-gradu(i,j-2,2)) * (gradv(i,j-2,1)-gradu(i,j-2,2))
    !ducros2   = divu2/(divu2+vort2+1d-15)
    !dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j-2)/(sq+1.d-15)*divu) )
    !
    !coef      = MAX( k_sensor1*ducros1*dxm1, &
    !                 k_sensor2*ducros2*dxm2, coef)

    !coef      = MAX(sensor(i  ,j-1,2), sensor(i  ,j,2), &
    !                sensor(i-1,j-1,2), sensor(i-1,j,2), &
    !                sensor(i+1,j-1,2), sensor(i+1,j,2)  )
    !
    !coef      = MAX(sensor(i,j-3,2), sensor(i,j-2,2), &
    !                sensor(i,j-1,2), sensor(i,j  ,2), &
    !                sensor(i,j+1,2), sensor(i,j+2,2)  )
    !
    !coef      = MAX(sensor(i,j-2,2), &
    !                sensor(i,j-1,2), sensor(i,j,2), &
    !                sensor(i,j+1,2))
    !
    !rspec      = rconv(i,j,2)

    
    eps2      = k2*coef
    eps4      = MAX(ZERO, k4 - 12.d0 * eps2)  ! 12.d0 to follow Sciacovelli CF 2021 
    

    diffro   = HALF * (w(i,j,1) - w(i,j-1,1))
    diffrou  = HALF * (w(i,j,2) - w(i,j-1,2))
    diffrov  = HALF * (w(i,j,3) - w(i,j-1,3))
    diffrow  = HALF * (w(i,j,4) - w(i,j-1,4))
    diffroe  = HALF * (w(i,j,5) - w(i,j-1,5))
    !
