      !idir
       grhoi = ( b1*(   w(i+1,j  , 1) -    w(i-1,j  , 1)) &
               + b2*(   w(i+2,j  , 1) -    w(i-2,j  , 1)) )
                            
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) &
               + b2*(velx(i+2,j     ) - velx(i-2,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) &
               + b2*(vely(i+2,j     ) - vely(i-2,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) &
               + b2*(velz(i+2,j     ) - velz(i-2,j     )) )
                               
       gti   = ( b1*(tloc(i+1,j     ) - tloc(i-1,j     )) &
               + b2*(tloc(i+2,j     ) - tloc(i-2,j     )) )
                                
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) &
               + b2*(  mu(i+2,j     ) -   mu(i-2,j     )) )
       
                           
       
