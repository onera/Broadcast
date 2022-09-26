      sc1 = nx(i,j,2)
      sc2 = ny(i,j,2)
      sn  = sqrt(sc1*sc1 + sc2*sc2)
      invsn = ONE/sn
      nxloc = sc1*invsn
      nyloc = sc2*invsn
      !face1
      dxl = HALF*(x0(i,j) + x0(i+1,j)) - xc(i,j-1)
      dxr = HALF*(x0(i,j) + x0(i+1,j)) - xc(i,j  )
      
      !face2
      dyl = HALF*(y0(i,j) + y0(i+1,j)) - yc(i,j-1)
      dyr = HALF*(y0(i,j) + y0(i+1,j)) - yc(i,j  )
