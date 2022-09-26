      !idir
       grhoj = ( b1*(   w(i  ,j+1, 1) -    w(i  ,j-1, 1)) &
               + b2*(   w(i  ,j+2, 1) -    w(i  ,j-2, 1)) &
               + b3*(   w(i  ,j+3, 1) -    w(i  ,j-3, 1)) &
               + b4*(   w(i  ,j+4, 1) -    w(i  ,j-4, 1)) )
              
       guj   = ( b1*(velx(i  ,j+1   ) - velx(i  ,j-1   )) &
               + b2*(velx(i  ,j+2   ) - velx(i  ,j-2   )) &
               + b3*(velx(i  ,j+3   ) - velx(i  ,j-3   )) &
               + b4*(velx(i  ,j+4   ) - velx(i  ,j-4   )) )
       
       gvj   = ( b1*(vely(i  ,j+1   ) - vely(i  ,j-1   )) &
               + b2*(vely(i  ,j+2   ) - vely(i  ,j-2   )) &
               + b3*(vely(i  ,j+3   ) - vely(i  ,j-3   )) &
               + b4*(vely(i  ,j+4   ) - vely(i  ,j-4   )) )
       
       gwj   = ( b1*(velz(i  ,j+1   ) - velz(i  ,j-1   )) &
               + b2*(velz(i  ,j+2   ) - velz(i  ,j-2   )) &
               + b3*(velz(i  ,j+3   ) - velz(i  ,j-3   )) &
               + b4*(velz(i  ,j+4   ) - velz(i  ,j-4   )) )
          
       gtj   = ( b1*(tloc(i  ,j+1   ) - tloc(i  ,j-1   )) &
               + b2*(tloc(i  ,j+2   ) - tloc(i  ,j-2   )) &
               + b3*(tloc(i  ,j+3   ) - tloc(i  ,j-3   )) &
               + b4*(tloc(i  ,j+4   ) - tloc(i  ,j-4   )) )
           
       gmuj  = ( b1*(  mu(i  ,j+1   ) -   mu(i  ,j-1   )) &
               + b2*(  mu(i  ,j+2   ) -   mu(i  ,j-2   )) &
               + b3*(  mu(i  ,j+3   ) -   mu(i  ,j-3   )) &
               + b4*(  mu(i  ,j+4   ) -   mu(i  ,j-4   )) )
       
       gf1j  = ( b1*(   f(i  ,j+1, 1) -    f(i  ,j-1, 1)) &
               + b2*(   f(i  ,j+2, 1) -    f(i  ,j-2, 1)) &
               + b3*(   f(i  ,j+3, 1) -    f(i  ,j-3, 1)) &
               + b4*(   f(i  ,j+4, 1) -    f(i  ,j-4, 1)) )
              
       gf2j  = ( b1*(   f(i  ,j+1, 2) -    f(i  ,j-1, 2)) &
               + b2*(   f(i  ,j+2, 2) -    f(i  ,j-2, 2)) &
               + b3*(   f(i  ,j+3, 2) -    f(i  ,j-3, 2)) &
               + b4*(   f(i  ,j+4, 2) -    f(i  ,j-4, 2)) )
               
       gf3j  = ( b1*(   f(i  ,j+1, 3) -    f(i  ,j-1, 3)) &
               + b2*(   f(i  ,j+2, 3) -    f(i  ,j-2, 3)) &
               + b3*(   f(i  ,j+3, 3) -    f(i  ,j-3, 3)) &
               + b4*(   f(i  ,j+4, 3) -    f(i  ,j-4, 3)) )
               
       gf4j  = ( b1*(   f(i  ,j+1, 4) -    f(i  ,j-1, 4)) &
               + b2*(   f(i  ,j+2, 4) -    f(i  ,j-2, 4)) &
               + b3*(   f(i  ,j+3, 4) -    f(i  ,j-3, 4)) &
               + b4*(   f(i  ,j+4, 4) -    f(i  ,j-4, 4)) )
               
       gf5j  = ( b1*(   f(i  ,j+1, 5) -    f(i  ,j-1, 5)) &
               + b2*(   f(i  ,j+2, 5) -    f(i  ,j-2, 5)) &
               + b3*(   f(i  ,j+3, 5) -    f(i  ,j-3, 5)) &
               + b4*(   f(i  ,j+4, 5) -    f(i  ,j-4, 5)) )   
               
       gg1j  = ( b1*(   g(i  ,j+1, 1) -    g(i  ,j-1, 1)) &
               + b2*(   g(i  ,j+2, 1) -    g(i  ,j-2, 1)) &
               + b3*(   g(i  ,j+3, 1) -    g(i  ,j-3, 1)) &
               + b4*(   g(i  ,j+4, 1) -    g(i  ,j-4, 1)) )
              
       gg2j  = ( b1*(   g(i  ,j+1, 2) -    g(i  ,j-1, 2)) &
               + b2*(   g(i  ,j+2, 2) -    g(i  ,j-2, 2)) &
               + b3*(   g(i  ,j+3, 2) -    g(i  ,j-3, 2)) &
               + b4*(   g(i  ,j+4, 4) -    g(i  ,j-4, 4)) )
               
       gg3j  = ( b1*(   g(i  ,j+1, 3) -    g(i  ,j-1, 3)) &
               + b2*(   g(i  ,j+2, 3) -    g(i  ,j-2, 3)) &
               + b3*(   g(i  ,j+3, 3) -    g(i  ,j-3, 3)) &
               + b4*(   g(i  ,j+4, 4) -    g(i  ,j-4, 4)) )
               
       gg4j  = ( b1*(   g(i  ,j+1, 4) -    g(i  ,j-1, 4)) &
               + b2*(   g(i  ,j+2, 4) -    g(i  ,j-2, 4)) &
               + b3*(   g(i  ,j+3, 4) -    g(i  ,j-3, 4)) &
               + b4*(   g(i  ,j+4, 4) -    g(i  ,j-4, 4)) )
               
       gg5j  = ( b1*(   g(i  ,j+1, 5) -    g(i  ,j-1, 5)) &
               + b2*(   g(i  ,j+2, 5) -    g(i  ,j-2, 5)) &
               + b3*(   g(i  ,j+3, 5) -    g(i  ,j-3, 5)) &
               + b4*(   g(i  ,j+4, 4) -    g(i  ,j-4, 4)) )                
       
