      !shock capturing via ducros
      k_sensor1 = ABS(p(i,j-1) - TWO*p(i,j) + p(i,j+1)) / &
                  ABS(p(i,j-1) + TWO*p(i,j) + p(i,j+1))
                
      k_sensor2 = ABS(p(i,j-2) - TWO*p(i,j-1) + p(i,j)) / &
                  ABS(p(i,j-2) + TWO*p(i,j-1) + p(i,j))
                  
      !k_sensor1 = min(10.d0*k_sensor1, ONE)            
      !k_sensor2 = min(10.d0*k_sensor2, ONE) 
      
      soundm = SQRT(gam*rgaz*tloc(i  ,j-1))
      soundp = SQRT(gam*rgaz*tloc(i  ,j  ))
                
      divu2     = (gradu(i,j,1)+gradv(i,j,2)) * (gradu(i,j,1)+gradv(i,j,2))
      vort2     = (gradv(i,j,1)-gradu(i,j,2)) * (gradv(i,j,1)-gradu(i,j,2))
      dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*invsn/soundp*sqrt(divu2)) )
      ducrosr   = ONE - divu2/(divu2+vort2+1d-15)*k_sensor1                
      
      divu2     = (gradu(i,j-1,1)+gradv(i,j-1,2)) * (gradu(i,j-1,1)+gradv(i,j-1,2))
      vort2     = (gradv(i,j-1,1)-gradu(i,j-1,2)) * (gradv(i,j-1,1)-gradu(i,j-1,2))
      dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*invsn/soundp*sqrt(divu2))) 
      ducrosl   = ONE - divu2/(divu2+vort2+1d-15)*k_sensor2
