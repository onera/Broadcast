      !jdir
       grhoj = ( b1*(   w(i,j+1, 1) -    w(i,j-1, 1)) &
               + b2*(   w(i,j+2, 1) -    w(i,j-2, 1)) &
               + b3*(   w(i,j+3, 1) -    w(i,j-3, 1)) &
               + b4*(   w(i,j+4, 1) -    w(i,j-4, 1)) &
               + b5*(   w(i,j+5, 1) -    w(i,j-5, 1)) )
                 
       guj   = ( b1*(velx(i,j+1   ) - velx(i,j-1   )) &
               + b2*(velx(i,j+2   ) - velx(i,j-2   )) &
               + b3*(velx(i,j+3   ) - velx(i,j-3   )) &
               + b4*(velx(i,j+4   ) - velx(i,j-4   )) &
               + b5*(velx(i,j+5   ) - velx(i,j-5   )) )
       
       gvj   = ( b1*(vely(i,j+1   ) - vely(i,j-1   )) &
               + b2*(vely(i,j+2   ) - vely(i,j-2   )) &
               + b3*(vely(i,j+3   ) - vely(i,j-3   )) &
               + b4*(vely(i,j+4   ) - vely(i,j-4   )) &
               + b5*(vely(i,j+5   ) - vely(i,j-5   )) )
       
       gwj   = ( b1*(velz(i,j+1   ) - velz(i,j-1   )) &
               + b2*(velz(i,j+2   ) - velz(i,j-2   )) &
               + b3*(velz(i,j+3   ) - velz(i,j-3   )) &
               + b4*(velz(i,j+4   ) - velz(i,j-4   )) &
               + b5*(velz(i,j+5   ) - velz(i,j-5   )) )

       gtj   = ( b1*(tloc(i,j+1   ) - tloc(i,j-1   )) &
               + b2*(tloc(i,j+2   ) - tloc(i,j-2   )) &
               + b3*(tloc(i,j+3   ) - tloc(i,j-3   )) &
               + b4*(tloc(i,j+4   ) - tloc(i,j-4   )) &
               + b5*(tloc(i,j+5   ) - tloc(i,j-5   )) )

       gmuj  = ( b1*(  mu(i,j+1   ) -   mu(i,j-1   )) &
               + b2*(  mu(i,j+2   ) -   mu(i,j-2   )) &
               + b3*(  mu(i,j+3   ) -   mu(i,j-3   )) &
               + b4*(  mu(i,j+4   ) -   mu(i,j-4   )) &
               + b5*(  mu(i,j+5   ) -   mu(i,j-5   )) )
       
                           
       
