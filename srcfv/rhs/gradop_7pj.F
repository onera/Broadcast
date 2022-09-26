      !jdir
       guj   = ( b1*(velx(i,j+1   ) - velx(i,j-1   )) &
               + b2*(velx(i,j+2   ) - velx(i,j-2   )) &
               + b3*(velx(i,j+3   ) - velx(i,j-3   )) )
       
       gvj   = ( b1*(vely(i,j+1   ) - vely(i,j-1   )) &
               + b2*(vely(i,j+2   ) - vely(i,j-2   )) &
               + b3*(vely(i,j+3   ) - vely(i,j-3   )) )
       
       gwj   = ( b1*(velz(i,j+1   ) - velz(i,j-1   )) &
               + b2*(velz(i,j+2   ) - velz(i,j-2   )) &
               + b3*(velz(i,j+3   ) - velz(i,j-3   )) )

       gmuj  = ( b1*(  mu(i,j+1   ) -   mu(i,j-1   )) &
               + b2*(  mu(i,j+2   ) -   mu(i,j-2   )) &
               + b3*(  mu(i,j+3   ) -   mu(i,j-3   )) )
       
                           
       
