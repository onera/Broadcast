      rol   = w(i-i0,j-j0,1)
      ror   = w(i   ,j   ,1)
      rom1l = ONE/rol
      rom1r = ONE/ror
      iro   = 2.d0/(rol+ror)
      
      uu    = HALF*(w(i-i0,j-j0,2)*rom1l + w(i,j,2)*rom1r)
      vv    = HALF*(w(i-i0,j-j0,3)*rom1l + w(i,j,3)*rom1r)
        
      cc    = gam*rgaz*HALF*(tloc(i-i0,j-j0)+tloc(i,j)) 
      norm  = nx(i,j,kdir)*nx(i,j,kdir)+ny(i,j,kdir)*ny(i,j,kdir)
      
      mutot =     HALF * gampr  * ( mu  (i-i0,j-j0) + mu  (i,j) )      &
               +  HALF * gamprt * ( mut (i-i0,j-j0) + mut (i,j) ) 
      
      

      specdiff = mutot * norm * volf(i,j,kdir) * iro
      
      coef(i,j,1,kdir) = 0.5d0 * ( abs(  uu*nx(i,j,kdir)  &
                                       + vv*ny(i,j,kdir)) &
                                       + sqrt(cc*norm)    &
                                       + TWO*specdiff     )


      
      !mutot = HALF * ( mu(i-i0,j-j0) + mu(i,j) + w(i-i0,j-j0,6) + w(i,j,6) )
      mutot = HALF * ( mu(i-i0,j-j0) + mu(i,j) + mut(i-i0,j-j0) + mut(i,j) )

      difftur = mutot * sigmainv * norm * volf(i,j,kdir) * iro
      
      coef(i,j,2,kdir) = 0.5d0 * ( abs(  uu*nx(i,j,kdir)  &
                                       + vv*ny(i,j,kdir)) &
                                       + TWO*difftur      )
