      
       grhoj   = ( a0   *    w(i  ,1 , 1) &
                 + a1   *    w(i  ,2 , 1) &
                 + a2   *    w(i  ,3 , 1) &
                 + a3   *    w(i  ,4 , 1) &
                 + a4   *    w(i  ,5 , 1) &
                 + a5   *    w(i  ,6 , 1) &
                 + a6   *    w(i  ,7 , 1) &
                 + a7   *    w(i  ,8 , 1) &
                 + a8   *    w(i  ,9 , 1) )
       
       guj     = ( c3d0 * velx(i  ,1    ) &
                 + c3d1 * velx(i  ,2    ) &
                 + c3d2 * velx(i  ,3    ) &
                 + c3d3 * velx(i  ,4    ) &
                 + c3d4 * velx(i  ,5    ) &
                 + c3d5 * velx(i  ,6    ) &
                 + c3d6 * velx(i  ,7    ) &
                 + c3d7 * velx(i  ,8    ) &
                 + c3d8 * velx(i  ,9    ) )
       
       gvj     = ( c3d0 * vely(i  ,1    ) &
                 + c3d1 * vely(i  ,2    ) &
                 + c3d2 * vely(i  ,3    ) &
                 + c3d3 * vely(i  ,4    ) &
                 + c3d4 * vely(i  ,5    ) &
                 + c3d5 * vely(i  ,6    ) &
                 + c3d6 * vely(i  ,7    ) &
                 + c3d7 * vely(i  ,8    ) &
                 + c3d8 * vely(i  ,9    ) )
       
       gwj     = ( c3d0 * velz(i  ,1    ) &
                 + c3d1 * velz(i  ,2    ) &
                 + c3d2 * velz(i  ,3    ) &
                 + c3d3 * velz(i  ,4    ) &
                 + c3d4 * velz(i  ,5    ) &
                 + c3d5 * velz(i  ,6    ) &
                 + c3d6 * velz(i  ,7    ) &
                 + c3d7 * velz(i  ,8    ) &
                 + c3d8 * velz(i  ,9    ) )
       
                    
       gtj  = ( a0      * tloc(i  ,1    ) &
              + a1      * tloc(i  ,2    ) &
              + a2      * tloc(i  ,3    ) &
              + a3      * tloc(i  ,4    ) &
              + a4      * tloc(i  ,5    ) &
              + a5      * tloc(i  ,6    ) &
              + a6      * tloc(i  ,7    ) &
              + a7      * tloc(i  ,8    ) &
              + a8      * tloc(i  ,9    ) )
                                          
                                    
       gmuj = ( a0      *   mu(i  ,1    ) &
              + a1      *   mu(i  ,2    ) &
              + a2      *   mu(i  ,3    ) &
              + a3      *   mu(i  ,4    ) &
              + a4      *   mu(i  ,5    ) &
              + a5      *   mu(i  ,6    ) &
              + a6      *   mu(i  ,7    ) &
              + a7      *   mu(i  ,8    ) &
              + a8      *   mu(i  ,9    ) )
                                                     
       
