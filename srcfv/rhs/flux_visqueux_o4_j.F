                 
    ! grad for viscous fluxes o4 - 5p
    !
    ! for coef ref Zhing et al, JCP2000 ou Shen et al AIAAP 2008
    ! 1/16 = 0.0625 ccross = 1/12*1/16
    ! TWENTYFOURTH = ONE/24.d0
    ! ccross = TWELFTH*0.0625d0

    volm1 = volf(i,j,2)

    nx_N =   HALF*( nx(i+1,j-1,1) + nx(i+1,j,1) )
    nx_S = - HALF*( nx(i  ,j-1,1) + nx(i  ,j,1) )
    nx_O =   HALF*( nx(i  ,j+1,2) + nx(i  ,j,2) )
    nx_E = - HALF*( nx(i  ,j-1,2) + nx(i  ,j,2) )
    !
    ny_N =   HALF*( ny(i+1,j-1,1) + ny(i+1,j,1) )
    ny_S = - HALF*( ny(i  ,j-1,1) + ny(i  ,j,1) )
    ny_O =   HALF*( ny(i  ,j+1,2) + ny(i  ,j,2) )
    ny_E = - HALF*( ny(i  ,j-1,2) + ny(i  ,j,2) )            
    
    val_N = ccross  * ( -        (-velx(i-1,j-2)+ 9.d0 * velx(i-1,j-1) + 9.d0 * velx(i-1,j  ) - velx(i-1,j+1)) &
                        + 7.d0 * (-velx(i  ,j-2)+ 9.d0 * velx(i  ,j-1) + 9.d0 * velx(i  ,j  ) - velx(i  ,j+1)) &
                        + 7.d0 * (-velx(i+1,j-2)+ 9.d0 * velx(i+1,j-1) + 9.d0 * velx(i+1,j  ) - velx(i+1,j+1)) &
                        -        (-velx(i+2,j-2)+ 9.d0 * velx(i+2,j-1) + 9.d0 * velx(i+2,j  ) - velx(i+2,j+1)) )
                     
    val_S = ccross  * ( -        (-velx(i-2,j-2)+ 9.d0 * velx(i-2,j-1) + 9.d0 * velx(i-2,j  ) - velx(i-2,j+1)) &
                        + 7.d0 * (-velx(i-1,j-2)+ 9.d0 * velx(i-1,j-1) + 9.d0 * velx(i-1,j  ) - velx(i-1,j+1)) &
                        + 7.d0 * (-velx(i  ,j-2)+ 9.d0 * velx(i  ,j-1) + 9.d0 * velx(i  ,j  ) - velx(i  ,j+1)) &
                        -        (-velx(i+1,j-2)+ 9.d0 * velx(i+1,j-1) + 9.d0 * velx(i+1,j  ) - velx(i+1,j+1)) )
                     
    val_E = TWENTYFOURTH * (- velx(i  ,j  ) + 26.d0 * velx(i  ,j-1) - velx(i  ,j-2))
    
    val_O = TWENTYFOURTH * (- velx(i  ,j+1) + 26.d0 * velx(i  ,j  ) - velx(i  ,j-1))
    
    
    ux = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    uy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = ccross  * ( -        (-vely(i-1,j-2)+ 9.d0 * vely(i-1,j-1) + 9.d0 * vely(i-1,j  ) - vely(i-1,j+1)) &
                        + 7.d0 * (-vely(i  ,j-2)+ 9.d0 * vely(i  ,j-1) + 9.d0 * vely(i  ,j  ) - vely(i  ,j+1)) &
                        + 7.d0 * (-vely(i+1,j-2)+ 9.d0 * vely(i+1,j-1) + 9.d0 * vely(i+1,j  ) - vely(i+1,j+1)) &
                        -        (-vely(i+2,j-2)+ 9.d0 * vely(i+2,j-1) + 9.d0 * vely(i+2,j  ) - vely(i+2,j+1)) )
                     
    val_S = ccross  * ( -        (-vely(i-2,j-2)+ 9.d0 * vely(i-2,j-1) + 9.d0 * vely(i-2,j  ) - vely(i-2,j+1)) &
                        + 7.d0 * (-vely(i-1,j-2)+ 9.d0 * vely(i-1,j-1) + 9.d0 * vely(i-1,j  ) - vely(i-1,j+1)) &
                        + 7.d0 * (-vely(i  ,j-2)+ 9.d0 * vely(i  ,j-1) + 9.d0 * vely(i  ,j  ) - vely(i  ,j+1)) &
                        -        (-vely(i+1,j-2)+ 9.d0 * vely(i+1,j-1) + 9.d0 * vely(i+1,j  ) - vely(i+1,j+1)) )
                     
    val_E = TWENTYFOURTH * (- vely(i  ,j  ) + 26.d0 * vely(i  ,j-1) - vely(i  ,j-2))
    
    val_O = TWENTYFOURTH * (- vely(i  ,j+1) + 26.d0 * vely(i  ,j  ) - vely(i  ,j-1))
    
    vx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    vy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = ccross  * ( -        (-velz(i-1,j-2)+ 9.d0 * velz(i-1,j-1) + 9.d0 * velz(i-1,j  ) - velz(i-1,j+1)) &
                        + 7.d0 * (-velz(i  ,j-2)+ 9.d0 * velz(i  ,j-1) + 9.d0 * velz(i  ,j  ) - velz(i  ,j+1)) &
                        + 7.d0 * (-velz(i+1,j-2)+ 9.d0 * velz(i+1,j-1) + 9.d0 * velz(i+1,j  ) - velz(i+1,j+1)) &
                        -        (-velz(i+2,j-2)+ 9.d0 * velz(i+2,j-1) + 9.d0 * velz(i+2,j  ) - velz(i+2,j+1)) )
                     
    val_S = ccross  * ( -        (-velz(i-2,j-2)+ 9.d0 * velz(i-2,j-1) + 9.d0 * velz(i-2,j  ) - velz(i-2,j+1)) &
                        + 7.d0 * (-velz(i-1,j-2)+ 9.d0 * velz(i-1,j-1) + 9.d0 * velz(i-1,j  ) - velz(i-1,j+1)) &
                        + 7.d0 * (-velz(i  ,j-2)+ 9.d0 * velz(i  ,j-1) + 9.d0 * velz(i  ,j  ) - velz(i  ,j+1)) &
                        -        (-velz(i+1,j-2)+ 9.d0 * velz(i+1,j-1) + 9.d0 * velz(i+1,j  ) - velz(i+1,j+1)) )
                     
    val_E = TWENTYFOURTH * (- velz(i  ,j  ) + 26.d0 * velz(i  ,j-1) - velz(i  ,j-2))
    
    val_O = TWENTYFOURTH * (- velz(i  ,j+1) + 26.d0 * velz(i  ,j  ) - velz(i  ,j-1))
    
    wx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    wy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = ccross  * ( -        (-tloc(i-1,j-2)+ 9.d0 * tloc(i-1,j-1) + 9.d0 * tloc(i-1,j  ) - tloc(i-1,j+1)) &
                        + 7.d0 * (-tloc(i  ,j-2)+ 9.d0 * tloc(i  ,j-1) + 9.d0 * tloc(i  ,j  ) - tloc(i  ,j+1)) &
                        + 7.d0 * (-tloc(i+1,j-2)+ 9.d0 * tloc(i+1,j-1) + 9.d0 * tloc(i+1,j  ) - tloc(i+1,j+1)) &
                        -        (-tloc(i+2,j-2)+ 9.d0 * tloc(i+2,j-1) + 9.d0 * tloc(i+2,j  ) - tloc(i+2,j+1)) )
                     
    val_S = ccross  * ( -        (-tloc(i-2,j-2)+ 9.d0 * tloc(i-2,j-1) + 9.d0 * tloc(i-2,j  ) - tloc(i-2,j+1)) &
                        + 7.d0 * (-tloc(i-1,j-2)+ 9.d0 * tloc(i-1,j-1) + 9.d0 * tloc(i-1,j  ) - tloc(i-1,j+1)) &
                        + 7.d0 * (-tloc(i  ,j-2)+ 9.d0 * tloc(i  ,j-1) + 9.d0 * tloc(i  ,j  ) - tloc(i  ,j+1)) &
                        -        (-tloc(i+1,j-2)+ 9.d0 * tloc(i+1,j-1) + 9.d0 * tloc(i+1,j  ) - tloc(i+1,j+1)) )
                     
    val_E = TWENTYFOURTH * (- tloc(i  ,j  ) + 26.d0 * tloc(i  ,j-1) - tloc(i  ,j-2))
    
    val_O = TWENTYFOURTH * (- tloc(i  ,j+1) + 26.d0 * tloc(i  ,j  ) - tloc(i  ,j-1))
    
    tx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    ty = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    
    uu     = 0.0625d0 * (-velx(i  ,j-2) + 9.d0 * velx(i  ,j-1) + 9.d0 * velx(i  ,j  ) - velx(i  ,j+1))
    vv     = 0.0625d0 * (-vely(i  ,j-2) + 9.d0 * vely(i  ,j-1) + 9.d0 * vely(i  ,j  ) - vely(i  ,j+1))
    ww     = 0.0625d0 * (-velz(i  ,j-2) + 9.d0 * velz(i  ,j-1) + 9.d0 * velz(i  ,j  ) - velz(i  ,j+1))
    !
    mmu    = 0.0625d0 * (  -mu(i  ,j-2) + 9.d0 *   mu(i  ,j-1) + 9.d0 *   mu(i  ,j  ) -   mu(i  ,j+1))
    lambda = mmu *cpprandtl
    !
    fvrou2 = TWOTHIRD*mmu*( TWO*ux -     vy        )    
    fvrov2 =          mmu*(     uy +     vx        )    
    fvrow2 =          mmu*(                    wx  )    
    fvroe2 = (lambda*tx + uu*fvrou2 + vv*fvrov2 + ww * fvrow2)
    ! fvroe2 = (lambda*tx + uu*fvrou2 + vv*fvrov2)
    
    gvrou2 =          mmu*(     uy +     vx      )  
    gvrov2 = TWOTHIRD*mmu*(    -ux + TWO*vy      )  
    gvrow2 =          mmu*(                   wy )  
    gvroe2 = (lambda*ty + uu*gvrou2 + vv*gvrov2 + ww * gvrow2)
    ! gvroe2 = (lambda*ty + uu*gvrou2 + vv*gvrov2)
    
    !!  fvrou2 = TWOTHIRD*mmu*( TWO*ux -     vy  -  wz  )    
    !!  fvrov2 =          mmu*(     uy +     vx         )    
    !!  fvrow2 =          mmu*(     uz           +  wx  )    
    !!  fvroe2 = (lambda*tx + uu*fvrou2 + vv*fvrov2 + ww * fvrow2)
    !!  
    !!  gvrou2 =          mmu*(     uy +     vx       )  
    !!  gvrov2 = TWOTHIRD*mmu*(    -ux + TWO*vy -  wz )  
    !!  gvrow2 =          mmu*(              vz +  wy )  
    !!  gvroe2 = (lambda*ty + uu*gvrou2 + vv*gvrov2 + ww*gvrow2)
    
  
