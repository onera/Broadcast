
    ror = w(i,j,1)
    rom1 = ONE/ror
    velxi = w(i,j,2) * rom1
    velyi = w(i,j,3) * rom1
    velzi = w(i,j,4) * rom1
    !
    ec  = HALF*( velxi*velxi &
               + velyi*velyi &
               + velzi*velzi)
    !
    eloc = (w(i,j,5) - ec*ror)*rom1
    !
    tloc(i,j) = eloc*cvm1
    !
    pint  = ror*rgaz*tloc(i,j)
    
    htot= (w(i,j,5) + pint)*rom1
    
    f(i,j,1)  = w(i,j,2) 
    f(i,j,2)  = w(i,j,2) * velxi + pint
    f(i,j,3)  = w(i,j,2) * velyi
    f(i,j,4)  = w(i,j,2) * velzi
    f(i,j,5)  = w(i,j,2) * htot
    
    g(i,j,1)  = w(i,j,3) 
    g(i,j,2)  = w(i,j,3) * velxi 
    g(i,j,3)  = w(i,j,3) * velyi + pint
    g(i,j,4)  = w(i,j,3) * velzi
    g(i,j,5)  = w(i,j,3) * htot
    
    !h(i,j,1)  = w(i,j,4) 
    !h(i,j,2)  = w(i,j,4) * velx(i,j) 
    !h(i,j,3)  = w(i,j,4) * vely(i,j) 
    !h(i,j,4)  = w(i,j,4) * velz(i,j) + pint
    !h(i,j,5)  = w(i,j,4) * htot
    
