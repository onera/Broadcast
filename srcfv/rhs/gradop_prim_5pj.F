      !jdir
       grhoj = ( b1*(   w(i,j+1, 1) -    w(i,j-1, 1)) &
               + b2*(   w(i,j+2, 1) -    w(i,j-2, 1)) )
                 
       guj   = ( b1*(velx(i,j+1   ) - velx(i,j-1   )) &
               + b2*(velx(i,j+2   ) - velx(i,j-2   )) )
       
       gvj   = ( b1*(vely(i,j+1   ) - vely(i,j-1   )) &
               + b2*(vely(i,j+2   ) - vely(i,j-2   )) )
       
       gwj   = ( b1*(velz(i,j+1   ) - velz(i,j-1   )) &
               + b2*(velz(i,j+2   ) - velz(i,j-2   )) )

       gtj   = ( b1*(tloc(i,j+1   ) - tloc(i,j-1   )) &
               + b2*(tloc(i,j+2   ) - tloc(i,j-2   )) )

       gmuj  = ( b1*(  mu(i,j+1   ) -   mu(i,j-1   )) &
               + b2*(  mu(i,j+2   ) -   mu(i,j-2   )) )
       
                           
       
