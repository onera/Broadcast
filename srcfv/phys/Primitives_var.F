
    ro = w(i,j,1)
    rom1 = ONE/ro
    velx(i,j) = w(i,j,2) * rom1
    vely(i,j) = w(i,j,3) * rom1
    velz(i,j) = w(i,j,4) * rom1
    !
    ec  = HALF*( velx(i,j)*velx(i,j) &
               + vely(i,j)*vely(i,j) &
               + velz(i,j)*velz(i,j))
    !
    eloc = (w(i,j,5) - ec*ro)*rom1
    !
    tloc(i,j) = eloc*cvm1
    !
    p(i,j)  = ro*rgaz*tloc(i,j)
    
    
  
