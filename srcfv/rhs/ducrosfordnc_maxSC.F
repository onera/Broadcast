    !1st direction

    !$AD II-LOOP
    do l = j-gh, j+gh-1
    !$AD II-LOOP
    !DIR$ IVDEP            
    do k = i-gh, i+gh-1
    
    k_sensor1 = ABS(p(k-1,l) - TWO*p(k,l) + p(k+1,l)) / &
                ABS(p(k-1,l) + TWO*p(k,l) + p(k+1,l))
                
    ! k_sensor2 = ABS(p(k-2,l) - TWO*p(k-1,l) + p(k,l)) / &
    !             ABS(p(k-2,l) + TWO*p(k-1,l) + p(k,l))
                
    divu2     = (gradu(k,l,1)+gradv(k,l,2)) * (gradu(k,l,1)+gradv(k,l,2))
    vort2     = (gradv(k,l,1)-gradu(k,l,2)) * (gradv(k,l,1)-gradu(k,l,2))
    ducros1   = divu2/(divu2+vort2+1d-15)
    
    dxm1       = HALF*(ONE-tanh(2.5d0+10.d0*vol(k,l)/(sq+1d-15)*(gradu(k,l,1)+gradv(k,l,2))) )

    ! divu2     = (gradu(k-1,l,1)+gradv(k-1,l,2)) * (gradu(k-1,l,1)+gradv(k-1,l,2))
    ! vort2     = (gradv(k-1,l,1)-gradu(k-1,l,2)) * (gradv(k-1,l,1)-gradu(k-1,l,2))
    ! ducros2   = divu2/(divu2+vort2+1d-15)
    
    ! dxm2       = HALF*(ONE-tanh(2.5d0+10.d0*vol(k,l)/(sq+1d-15)*(gradu(k-1,l,1)+gradv(k-1,l,2))) ) 

    
    ! coef      = MAX( k_sensor1*ducros1, &
    !                  k_sensor2*ducros2)

    ! coef      = MAX( k_sensor1*ducros1*dxm1, &
    !                  k_sensor2*ducros2*dxm2)

    coef      = k_sensor1*ducros1*dxm1
        
    eps2arra(k,l)      = k2*coef
    ! eps2      = ONE

    enddo
    enddo

    eps2 = MAX(eps2arra(i,j), MAX(eps2arra(i-1,j), MAX(eps2arra(i,j-1), MAX(eps2arra(i+1,j), MAX(eps2arra(i,j+1), &
        MAX(eps2arra(i+1,j+1), MAX(eps2arra(i+1,j-1), MAX(eps2arra(i-1,j+1), MAX(eps2arra(i-1,j-1), & 
        MAX(eps2arra(i+2,j), MAX(eps2arra(i-2,j), MAX(eps2arra(i-3,j), &
        MAX(eps2arra(i,j+2), MAX(eps2arra(i,j-2), eps2arra(i,j-3)))))))))))))) )

    ! eps4      = MAX(ZERO,k4-eps2) 
    eps4      = MAX(ZERO,k4 - 12.d0 * eps2) 
    ! eps4      = MAX(ZERO,k4 - 3.d0/70.d0 * eps2)

    diffro   = ( w(i,j,1) - w(i-1,j,1) )
    diffrou  = ( w(i,j,2) - w(i-1,j,2) )
    diffrov  = ( w(i,j,3) - w(i-1,j,3) )
    diffrow  = ( w(i,j,4) - w(i-1,j,4) )
    diffroe  = ( w(i,j,5) - w(i-1,j,5) )
    !
