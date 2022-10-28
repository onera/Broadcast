                 
    ! grad for viscous fluxes o4 - 5p
    ! for coef ref Zhing et al, JCP2000 ou Shen et al AIAAP 2008
    ! 1/16 = 0.0625 ccross = 1/12*1/16
    ! TWENTYFOURTH = ONE/24.d0
    ! ccros = TWELFTH*0.0625d0

    volm1 = volf(i,j,1)
    
    nx_N =   HALF*( nx(i+1,j  ,1) + nx(i  ,j  ,1) )
    nx_S = - HALF*( nx(i-1,j  ,1) + nx(i  ,j  ,1) )
    nx_O =   HALF*( nx(i-1,j+1,2) + nx(i  ,j+1,2) )
    nx_E = - HALF*( nx(i-1,j  ,2) + nx(i  ,j  ,2) )
    !
    ny_N =   HALF*( ny(i+1,j  ,1) + ny(i  ,j  ,1) )
    ny_S = - HALF*( ny(i-1,j  ,1) + ny(i  ,j  ,1) )
    ny_O =   HALF*( ny(i-1,j+1,2) + ny(i  ,j+1,2) )
    ny_E = - HALF*( ny(i-1,j  ,2) + ny(i  ,j  ,2) )
    
    val_N = TWENTYFOURTH * (- velx(i+1,j  ) + 26.d0 * velx(i  ,j  ) - velx(i-1,j  ))
    val_S = TWENTYFOURTH * (- velx(i  ,j  ) + 26.d0 * velx(i-1,j  ) - velx(i-2,j  ))
    
    val_E = ccross * (  -        (-velx(i-2,j-2)+ 9.d0 * velx(i-1,j-2) + 9.d0 * velx(i  ,j-2) - velx(i+1,j-2)) &
                        + 7.d0 * (-velx(i-2,j-1)+ 9.d0 * velx(i-1,j-1) + 9.d0 * velx(i  ,j-1) - velx(i+1,j-1)) &
                        + 7.d0 * (-velx(i-2,j  )+ 9.d0 * velx(i-1,j  ) + 9.d0 * velx(i  ,j  ) - velx(i+1,j  )) &
                        -        (-velx(i-2,j+1)+ 9.d0 * velx(i-1,j+1) + 9.d0 * velx(i  ,j+1) - velx(i+1,j+1)) )
                     
    val_O = ccross * (  -        (-velx(i-2,j-1)+ 9.d0 * velx(i-1,j-1) + 9.d0 * velx(i  ,j-1) - velx(i+1,j-1)) &
                        + 7.d0 * (-velx(i-2,j  )+ 9.d0 * velx(i-1,j  ) + 9.d0 * velx(i  ,j  ) - velx(i+1,j  )) &
                        + 7.d0 * (-velx(i-2,j+1)+ 9.d0 * velx(i-1,j+1) + 9.d0 * velx(i  ,j+1) - velx(i+1,j+1)) &
                        -        (-velx(i-2,j+2)+ 9.d0 * velx(i-1,j+2) + 9.d0 * velx(i  ,j+2) - velx(i+1,j+2)) )
            
    
    ux = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    uy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = TWENTYFOURTH * (- vely(i+1,j  ) + 26.d0 * vely(i  ,j  ) - vely(i-1,j  ))
    val_S = TWENTYFOURTH * (- vely(i  ,j  ) + 26.d0 * vely(i-1,j  ) - vely(i-2,j  ))
    
    val_E = ccross * (  -        (-vely(i-2,j-2)+ 9.d0 * vely(i-1,j-2) + 9.d0 * vely(i  ,j-2) - vely(i+1,j-2)) &
                        + 7.d0 * (-vely(i-2,j-1)+ 9.d0 * vely(i-1,j-1) + 9.d0 * vely(i  ,j-1) - vely(i+1,j-1)) &
                        + 7.d0 * (-vely(i-2,j  )+ 9.d0 * vely(i-1,j  ) + 9.d0 * vely(i  ,j  ) - vely(i+1,j  )) &
                        -        (-vely(i-2,j+1)+ 9.d0 * vely(i-1,j+1) + 9.d0 * vely(i  ,j+1) - vely(i+1,j+1)) )
                     
    val_O = ccross * (  -        (-vely(i-2,j-1)+ 9.d0 * vely(i-1,j-1) + 9.d0 * vely(i  ,j-1) - vely(i+1,j-1)) &
                        + 7.d0 * (-vely(i-2,j  )+ 9.d0 * vely(i-1,j  ) + 9.d0 * vely(i  ,j  ) - vely(i+1,j  )) &
                        + 7.d0 * (-vely(i-2,j+1)+ 9.d0 * vely(i-1,j+1) + 9.d0 * vely(i  ,j+1) - vely(i+1,j+1)) &
                        -        (-vely(i-2,j+2)+ 9.d0 * vely(i-1,j+2) + 9.d0 * vely(i  ,j+2) - vely(i+1,j+2)) )
    
    vx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    vy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = TWENTYFOURTH * (- velz(i+1,j  ) + 26.d0 * velz(i  ,j  ) - velz(i-1,j  ))
    val_S = TWENTYFOURTH * (- velz(i  ,j  ) + 26.d0 * velz(i-1,j  ) - velz(i-2,j  ))
    
    val_E = ccross * (  -        (-velz(i-2,j-2)+ 9.d0 * velz(i-1,j-2) + 9.d0 * velz(i  ,j-2) - velz(i+1,j-2)) &
                        + 7.d0 * (-velz(i-2,j-1)+ 9.d0 * velz(i-1,j-1) + 9.d0 * velz(i  ,j-1) - velz(i+1,j-1)) &
                        + 7.d0 * (-velz(i-2,j  )+ 9.d0 * velz(i-1,j  ) + 9.d0 * velz(i  ,j  ) - velz(i+1,j  )) &
                        -        (-velz(i-2,j+1)+ 9.d0 * velz(i-1,j+1) + 9.d0 * velz(i  ,j+1) - velz(i+1,j+1)) )
                     
    val_O = ccross * (  -        (-velz(i-2,j-1)+ 9.d0 * velz(i-1,j-1) + 9.d0 * velz(i  ,j-1) - velz(i+1,j-1)) &
                        + 7.d0 * (-velz(i-2,j  )+ 9.d0 * velz(i-1,j  ) + 9.d0 * velz(i  ,j  ) - velz(i+1,j  )) &
                        + 7.d0 * (-velz(i-2,j+1)+ 9.d0 * velz(i-1,j+1) + 9.d0 * velz(i  ,j+1) - velz(i+1,j+1)) &
                        -        (-velz(i-2,j+2)+ 9.d0 * velz(i-1,j+2) + 9.d0 * velz(i  ,j+2) - velz(i+1,j+2)) )
    
    wx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    wy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1 
                 
    val_N = TWENTYFOURTH * (- tloc(i+1,j  ) + 26.d0 * tloc(i  ,j  ) - tloc(i-1,j  ))
    val_S = TWENTYFOURTH * (- tloc(i  ,j  ) + 26.d0 * tloc(i-1,j  ) - tloc(i-2,j  ))
    
    val_E = ccross * (  -        (-tloc(i-2,j-2)+ 9.d0 * tloc(i-1,j-2) + 9.d0 * tloc(i  ,j-2) - tloc(i+1,j-2)) &
                        + 7.d0 * (-tloc(i-2,j-1)+ 9.d0 * tloc(i-1,j-1) + 9.d0 * tloc(i  ,j-1) - tloc(i+1,j-1)) &
                        + 7.d0 * (-tloc(i-2,j  )+ 9.d0 * tloc(i-1,j  ) + 9.d0 * tloc(i  ,j  ) - tloc(i+1,j  )) &
                        -        (-tloc(i-2,j+1)+ 9.d0 * tloc(i-1,j+1) + 9.d0 * tloc(i  ,j+1) - tloc(i+1,j+1)) )
                     
    val_O = ccross * (  -        (-tloc(i-2,j-1)+ 9.d0 * tloc(i-1,j-1) + 9.d0 * tloc(i  ,j-1) - tloc(i+1,j-1)) &
                        + 7.d0 * (-tloc(i-2,j  )+ 9.d0 * tloc(i-1,j  ) + 9.d0 * tloc(i  ,j  ) - tloc(i+1,j  )) &
                        + 7.d0 * (-tloc(i-2,j+1)+ 9.d0 * tloc(i-1,j+1) + 9.d0 * tloc(i  ,j+1) - tloc(i+1,j+1)) &
                        -        (-tloc(i-2,j+2)+ 9.d0 * tloc(i-1,j+2) + 9.d0 * tloc(i  ,j+2) - tloc(i+1,j+2)) )
    
    tx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    ty = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    !
    ! Computation of viscous fluxes
    !
    uu     = 0.0625d0* (-velx(i-2,j  ) + 9.d0 * velx(i-1,j  ) + 9.d0 * velx(i  ,j  ) - velx(i+1,j  ))
    vv     = 0.0625d0* (-vely(i-2,j  ) + 9.d0 * vely(i-1,j  ) + 9.d0 * vely(i  ,j  ) - vely(i+1,j  ))
    ww     = 0.0625d0* (-velz(i-2,j  ) + 9.d0 * velz(i-1,j  ) + 9.d0 * velz(i  ,j  ) - velz(i+1,j  ))
    !
    mmu    = 0.0625d0* (  -mu(i-2,j  ) + 9.d0 *   mu(i-1,j  ) + 9.d0 *   mu(i  ,j  ) -   mu(i+1,j  ))
 
    lambda = mmu *cpprandtl
    !
    fvrou1 = TWOTHIRD*mmu*( TWO*ux -     vy        )    
    fvrov1 =          mmu*(     uy +     vx        )    
    fvrow1 =          mmu*(                    wx  )    
    fvroe1 = (lambda*tx + uu*fvrou1 + vv*fvrov1 + ww * fvrow1) 
    !
    gvrou1 =          mmu*(     uy +     vx        )
    gvrov1 = TWOTHIRD*mmu*(    -ux + TWO*vy        )
    gvrow1 =          mmu*(                    wy  )    
    gvroe1 = lambda*ty + uu*gvrou1 + vv*gvrov1 + ww * gvrow1
    
    
    !!  fvrou1 = TWOTHIRD*mmu*( TWO*ux -     vy  -   wz  )    
    !!  fvrov1 =          mmu*(     uy +     vx          )    
    !!  fvrow1 =          mmu*(     uz           +   wx  )    
    !!  fvroe1 = lambda*tx + uu*fvrou1 + vv*fvrov1 + ww * fvrow1 
    !!  !
    !!  gvrou1 =          mmu*(     uy +     vx          )
    !!  gvrov1 = TWOTHIRD*mmu*(    -ux + TWO*vy   -  wz  )
    !!  gvrow1 =          mmu*(              vz   +  wy  )    
    !!  gvroe1 = lambda*ty + uu*gvrou1 + vv*gvrov1 + ww * gvrow1
