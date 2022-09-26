      xa  = x0(i  ,j  ) ; ya  = y0(i  ,j  )
      xb  = x0(i+1,j  ) ; yb  = y0(i+1,j  )
      xc1 = x0(i  ,j+1) ; yc1 = y0(i  ,j+1)
      xd  = x0(i+1,j+1) ; yd  = y0(i+1,j+1)
      !
      abx = xb -xa ; aby = yb -ya
      acx = xc1-xa ; acy = yc1-ya
      q1 = HALF*abs( abx*acy  - acx*aby )
      !
      dcx = xc1-xd ; dcy = yc1-yd
      dbx = xb -xd ; dby = yb -yd
      q2 = HALF*abs( dcx*dby  - dbx*dcy )
      !
      vol(i,j) = q1 + q2
