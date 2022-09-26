      
       grhoj   = ( a0   *    w(i  ,1 , 1) &
                 + a1   *    w(i  ,2 , 1) &
                 + a2   *    w(i  ,3 , 1) &
                 + a3   *    w(i  ,4 , 1) &
                 + a4   *    w(i  ,5 , 1) )
       
       guj     = ( c3d0 * velx(i  ,1    ) &
                 + c3d1 * velx(i  ,2    ) &
                 + c3d2 * velx(i  ,3    ) &
                 + c3d3 * velx(i  ,4    ) &
                 + c3d4 * velx(i  ,5    ) )
       
       gvj     = ( c3d0 * vely(i  ,1    ) &
                 + c3d1 * vely(i  ,2    ) &
                 + c3d2 * vely(i  ,3    ) &
                 + c3d3 * vely(i  ,4    ) &
                 + c3d4 * vely(i  ,5    ) )
       
       gwj     = ( c3d0 * velz(i  ,1    ) &
                 + c3d1 * velz(i  ,2    ) &
                 + c3d2 * velz(i  ,3    ) &
                 + c3d3 * velz(i  ,4    ) &
                 + c3d4 * velz(i  ,5    ) )
       
                    
       gtj  = ( a0      * tloc(i  ,1    ) &
              + a1      * tloc(i  ,2    ) &
              + a2      * tloc(i  ,3    ) &
              + a3      * tloc(i  ,4    ) &
              + a4      * tloc(i  ,5    ) )
                                          
                                    
       gmuj = ( a0      *   mu(i  ,1    ) &
              + a1      *   mu(i  ,2    ) &
              + a2      *   mu(i  ,3    ) &
              + a3      *   mu(i  ,4    ) &
              + a4      *   mu(i  ,5    ) )
                                                     
       
