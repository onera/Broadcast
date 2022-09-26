       sc1 = nx(i,j,2)
       sc2 = ny(i,j,2)
       sn  = sqrt(sc1*sc1 + sc2*sc2)
       invsn = ONE/sn
       nxloc = sc1*invsn
       nyloc = sc2*invsn
       
       hn(i,j,1,2) = fxro2  -  dissro2 - source(1)
       !  
       hn(i,j,2,2) = fxrou2 - dissrou2 - (fvrou2 * nxloc + gvrou2 * nyloc)*sn - source(2)
       !  
       hn(i,j,3,2) = fxrov2 - dissrov2 - (fvrov2 * nxloc + gvrov2 * nyloc)*sn - source(3)
       !  
       hn(i,j,4,2) = fxrow2 - dissrow2 - (fvrow2 * nxloc + gvrow2 * nyloc)*sn - source(4)
       !  
       hn(i,j,5,2) = fxroe2 - dissroe2 - (fvroe2 * nxloc + gvroe2 * nyloc)*sn - source(5)

