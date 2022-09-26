      sc1 = nx(i,j,1)
      sc2 = ny(i,j,1)
      sn  = sqrt(sc1*sc1 + sc2*sc2)
      invsn = ONE/sn
      nxloc = sc1*invsn
      nyloc = sc2*invsn
      
      
      dxl = HALF*(x0(i,j) + x0(i,j+1)) - xc(i-1,j)
      dxr = HALF*(x0(i,j) + x0(i,j+1)) - xc(i  ,j)
      
      
      dyl = HALF*(y0(i,j) + y0(i,j+1)) - yc(i-1,j)
      dyr = HALF*(y0(i,j) + y0(i,j+1)) - yc(i  ,j)
