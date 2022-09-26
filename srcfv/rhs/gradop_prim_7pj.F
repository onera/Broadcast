       
      !jdir
       grhoj = ( b1*(   w(i  ,j+1, 1) -    w(i  ,j-1, 1)) &
               + b2*(   w(i  ,j+2, 1) -    w(i  ,j-2, 1)) &
               + b3*(   w(i  ,j+3, 1) -    w(i  ,j-3, 1)) )
                                                         
       guj   = ( b1*(velx(i  ,j+1   ) - velx(i  ,j-1   )) &
               + b2*(velx(i  ,j+2   ) - velx(i  ,j-2   )) &
               + b3*(velx(i  ,j+3   ) - velx(i  ,j-3   )) )
                                                            
       gvj   = ( b1*(vely(i  ,j+1   ) - vely(i  ,j-1   )) &
               + b2*(vely(i  ,j+2   ) - vely(i  ,j-2   )) &
               + b3*(vely(i  ,j+3   ) - vely(i  ,j-3   )) )
                                                           
       gwj   = ( b1*(velz(i  ,j+1   ) - velz(i  ,j-1   )) &
               + b2*(velz(i  ,j+2   ) - velz(i  ,j-2   )) &
               + b3*(velz(i  ,j+3   ) - velz(i  ,j-3   )) )
                                                                               
       gtj   = ( b1*(tloc(i  ,j+1   ) - tloc(i  ,j-1   )) &
               + b2*(tloc(i  ,j+2   ) - tloc(i  ,j-2   )) &
               + b3*(tloc(i  ,j+3   ) - tloc(i  ,j-3   )) )
                                                                                  
       gmuj  = ( b1*(  mu(i  ,j+1   ) -   mu(i  ,j-1   )) &
               + b2*(  mu(i  ,j+2   ) -   mu(i  ,j-2   )) &
               + b3*(  mu(i  ,j+3   ) -   mu(i  ,j-3   )) )
                           
       
