      !idir
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) &
               + b2*(velx(i+2,j     ) - velx(i-2,j     )) &
               + b3*(velx(i+3,j     ) - velx(i-3,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) &
               + b2*(vely(i+2,j     ) - vely(i-2,j     )) &
               + b3*(vely(i+3,j     ) - vely(i-3,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) &
               + b2*(velz(i+2,j     ) - velz(i-2,j     )) &
               + b3*(velz(i+3,j     ) - velz(i-3,j     )) )
                               
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) &
               + b2*(  mu(i+2,j     ) -   mu(i-2,j     )) &
               + b3*(  mu(i+3,j     ) -   mu(i-3,j     )) )
       
                           
       
