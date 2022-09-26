      !idir
       grhoi = ( b1*(   w(i+1,j  , 1) -    w(i-1,j  , 1)) &
               + b2*(   w(i+2,j  , 1) -    w(i-2,j  , 1)) &
               + b3*(   w(i+3,j  , 1) -    w(i-3,j  , 1)) &
               + b4*(   w(i+4,j  , 1) -    w(i-4,j  , 1)) &
               + b5*(   w(i+5,j  , 1) -    w(i-5,j  , 1)) )
                            
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) &
               + b2*(velx(i+2,j     ) - velx(i-2,j     )) &
               + b3*(velx(i+3,j     ) - velx(i-3,j     )) &
               + b4*(velx(i+4,j     ) - velx(i-4,j     )) &
               + b5*(velx(i+5,j     ) - velx(i-5,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) &
               + b2*(vely(i+2,j     ) - vely(i-2,j     )) &
               + b3*(vely(i+3,j     ) - vely(i-3,j     )) &
               + b4*(vely(i+4,j     ) - vely(i-4,j     )) &
               + b5*(vely(i+5,j     ) - vely(i-5,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) &
               + b2*(velz(i+2,j     ) - velz(i-2,j     )) &
               + b3*(velz(i+3,j     ) - velz(i-3,j     )) &
               + b4*(velz(i+4,j     ) - velz(i-4,j     )) &
               + b5*(velz(i+5,j     ) - velz(i-5,j     )) )
                               
       gti   = ( b1*(tloc(i+1,j     ) - tloc(i-1,j     )) &
               + b2*(tloc(i+2,j     ) - tloc(i-2,j     )) &
               + b3*(tloc(i+3,j     ) - tloc(i-3,j     )) &
               + b4*(tloc(i+4,j     ) - tloc(i-4,j     )) &
               + b5*(tloc(i+5,j     ) - tloc(i-5,j     )) )
                                
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) &
               + b2*(  mu(i+2,j     ) -   mu(i-2,j     )) &
               + b3*(  mu(i+3,j     ) -   mu(i-3,j     )) &
               + b4*(  mu(i+4,j     ) -   mu(i-4,j     )) &
               + b5*(  mu(i+5,j     ) -   mu(i-5,j     )) )
       
                           
       
