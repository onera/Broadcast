      !idir
       grhoi = ( b1*(   w(i+1,j  , 1) -    w(i-1,j  , 1)) &
               + b2*(   w(i+2,j  , 1) -    w(i-2,j  , 1)) &
               + b3*(   w(i+3,j  , 1) -    w(i-3,j  , 1)) )
                            
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) &
               + b2*(velx(i+2,j     ) - velx(i-2,j     )) &
               + b3*(velx(i+3,j     ) - velx(i-3,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) &
               + b2*(vely(i+2,j     ) - vely(i-2,j     )) &
               + b3*(vely(i+3,j     ) - vely(i-3,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) &
               + b2*(velz(i+2,j     ) - velz(i-2,j     )) &
               + b3*(velz(i+3,j     ) - velz(i-3,j     )) )
                               
       gti   = ( b1*(tloc(i+1,j     ) - tloc(i-1,j     )) &
               + b2*(tloc(i+2,j     ) - tloc(i-2,j     )) &
               + b3*(tloc(i+3,j     ) - tloc(i-3,j     )) )
                                
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) &
               + b2*(  mu(i+2,j     ) -   mu(i-2,j     )) &
               + b3*(  mu(i+3,j     ) -   mu(i-3,j     )) )
       
       gf1i  = ( b1*(   f(i+1,j  , 1) -    f(i-1,j  , 1)) &
               + b2*(   f(i+2,j  , 1) -    f(i-2,j  , 1)) &
               + b3*(   f(i+3,j  , 1) -    f(i-3,j  , 1)) )
              
       gf2i  = ( b1*(   f(i+1,j  , 2) -    f(i-1,j  , 2)) &
               + b2*(   f(i+2,j  , 2) -    f(i-2,j  , 2)) &
               + b3*(   f(i+3,j  , 2) -    f(i-3,j  , 2)) )
               
       gf3i  = ( b1*(   f(i+1,j  , 3) -    f(i-1,j  , 3)) &
               + b2*(   f(i+2,j  , 3) -    f(i-2,j  , 3)) &
               + b3*(   f(i+3,j  , 3) -    f(i-3,j  , 3)) )
               
       gf4i  = ( b1*(   f(i+1,j  , 4) -    f(i-1,j  , 4)) &
               + b2*(   f(i+2,j  , 4) -    f(i-2,j  , 4)) &
               + b3*(   f(i+3,j  , 4) -    f(i-3,j  , 4)) )
               
       gf5i  = ( b1*(   f(i+1,j  , 5) -    f(i-1,j  , 5)) &
               + b2*(   f(i+2,j  , 5) -    f(i-2,j  , 5)) &
               + b3*(   f(i+3,j  , 5) -    f(i-3,j  , 5)) )   
               
       gg1i  = ( b1*(   g(i+1,j  , 1) -    g(i-1,j  , 1)) &
               + b2*(   g(i+2,j  , 1) -    g(i-2,j  , 1)) &
               + b3*(   g(i+3,j  , 1) -    g(i-3,j  , 1)) )
              
       gg2i  = ( b1*(   g(i+1,j  , 2) -    g(i-1,j  , 2)) &
               + b2*(   g(i+2,j  , 2) -    g(i-2,j  , 2)) &
               + b3*(   g(i+3,j  , 2) -    g(i-3,j  , 2)) )
               
       gg3i  = ( b1*(   g(i+1,j  , 3) -    g(i-1,j  , 3)) &
               + b2*(   g(i+2,j  , 3) -    g(i-2,j  , 3)) &
               + b3*(   g(i+3,j  , 3) -    g(i-3,j  , 3)) )
               
       gg4i  = ( b1*(   g(i+1,j  , 4) -    g(i-1,j  , 4)) &
               + b2*(   g(i+2,j  , 4) -    g(i-2,j  , 4)) &
               + b3*(   g(i+3,j  , 4) -    g(i-3,j  , 4)) )
               
       gg5i  = ( b1*(   g(i+1,j  , 5) -    g(i-1,j  , 5)) &
               + b2*(   g(i+2,j  , 5) -    g(i-2,j  , 5)) &
               + b3*(   g(i+3,j  , 5) -    g(i-3,j  , 5)) )                
       
