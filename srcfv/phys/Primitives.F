
    ro = w(i,j,1)
    rom1 = ONE/ro
    velx(i,j) = w(i,j,2) * rom1
    vely(i,j) = w(i,j,3) * rom1
    velz(i,j) = w(i,j,4) * rom1
    !
    ec  = HALF*( velx(i,j)*velx(i,j) &
               + vely(i,j)*vely(i,j) &
               + velz(i,j)*velz(i,j))

    ! ec  = HALF*( velx(i,j)*velx(i,j) &
    !            + vely(i,j)*vely(i,j))
    !
    eloc = (w(i,j,5) - ec*ro)*rom1
    !
    tloc(i,j) = eloc*cvm1
    !
    p(i,j)  = (gam-ONE)*ro*eloc
    !  p(i,j)  = ro*rgaz*tloc(i,j)
    
    htot= (w(i,j,5) + p(i,j))*rom1
    
    f(i,j,1)  = w(i,j,2) 
    f(i,j,2)  = w(i,j,2) * velx(i,j) + p(i,j)
    f(i,j,3)  = w(i,j,2) * vely(i,j)
    f(i,j,4)  = w(i,j,2) * velz(i,j)
    f(i,j,5)  = w(i,j,2) * htot
    
    g(i,j,1)  = w(i,j,3) 
    g(i,j,2)  = w(i,j,3) * velx(i,j) 
    g(i,j,3)  = w(i,j,3) * vely(i,j) + p(i,j)
    g(i,j,4)  = w(i,j,3) * velz(i,j)
    g(i,j,5)  = w(i,j,3) * htot
    
    !h(i,j,1)  = w(i,j,4) 
    !h(i,j,2)  = w(i,j,4) * velx(i,j) 
    !h(i,j,3)  = w(i,j,4) * vely(i,j) 
    !h(i,j,4)  = w(i,j,4) * velz(i,j) + p(i,j)
    !h(i,j,5)  = w(i,j,4) * htot
    
