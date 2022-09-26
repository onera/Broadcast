       sc1 = nx(i,j,1)
       sc2 = ny(i,j,1)
       sn  = sqrt(sc1*sc1 + sc2*sc2)
       invsn = ONE/sn
       nxloc = sc1*invsn
       nyloc = sc2*invsn
       
       hn(i,j,1,1) = fxro1  -  dissro1  
       !  
       hn(i,j,2,1) = fxrou1 - dissrou1 - (fvrou1 * nxloc + gvrou1 * nyloc)*sn 
       !  
       hn(i,j,3,1) = fxrov1 - dissrov1 - (fvrov1 * nxloc + gvrov1 * nyloc)*sn
       !  
       hn(i,j,4,1) = fxrow1 - dissrow1 - (fvrow1 * nxloc + gvrow1 * nyloc)*sn
       ! hn(i,j,4,1) = ZERO
       !  
       hn(i,j,5,1) = fxroe1 - dissroe1 - (fvroe1 * nxloc + gvroe1 * nyloc)*sn

