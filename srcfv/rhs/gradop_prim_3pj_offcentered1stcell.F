      
       grhoj   = ( a0   *    w(i  ,1 , 1) &
                 + a1   *    w(i  ,2 , 1) &
                 + a2   *    w(i  ,3 , 1) )
       
       guj     = ( c3d0 * velx(i  ,1    ) &
                 + c3d1 * velx(i  ,2    ) &
                 + c3d2 * velx(i  ,3    ) )
       
       gvj     = ( c3d0 * vely(i  ,1    ) &
                 + c3d1 * vely(i  ,2    ) &
                 + c3d2 * vely(i  ,3    ) )
       
       gwj     = ( c3d0 * velz(i  ,1    ) &
                 + c3d1 * velz(i  ,2    ) &
                 + c3d2 * velz(i  ,3    ) )
       
                    
       gtj  = ( a0      * tloc(i  ,1    ) &
              + a1      * tloc(i  ,2    ) &
              + a2      * tloc(i  ,3    ) )
                                          
                                    
       gmuj = ( a0      *   mu(i  ,1    ) &
              + a1      *   mu(i  ,2    ) &
              + a2      *   mu(i  ,3    ) )
                                                     
       
