      volm1 = ONE/vol(i,j)
      dxm1 = HALF*(nx(i,j,1)+nx(i+1,j  ,1))*volm1
      dxm2 = HALF*(nx(i,j,2)+nx(i  ,j+1,2))*volm1
      
      dym1 = HALF*(ny(i,j,1)+ny(i+1,j  ,1))*volm1
      dym2 = HALF*(ny(i,j,2)+ny(i  ,j+1,2))*volm1
