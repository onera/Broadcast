                 
    ! grad for viscous fluxes o2 - 3p
    !
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
    
    val_N = FOURTH*( velx(i,j  ) + velx(i+1,j) + velx(i,j-1) + velx(i+1,j-1 ) )
    val_S = FOURTH*( velx(i,j  ) + velx(i-1,j) + velx(i,j-1) + velx(i-1,j-1 ) )
    val_O =          velx(i,j  )
    val_E =          velx(i,j-1)
    
    ux = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    uy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = FOURTH*( vely(i,j  ) + vely(i+1,j) + vely(i,j-1) + vely(i+1,j-1 ) )
    val_S = FOURTH*( vely(i,j  ) + vely(i-1,j) + vely(i,j-1) + vely(i-1,j-1 ) )
    val_O =          vely(i,j  )
    val_E =          vely(i,j-1)
    
    vx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    vy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = FOURTH*( velz(i,j  ) + velz(i+1,j) + velz(i,j-1) + velz(i+1,j-1 ) )
    val_S = FOURTH*( velz(i,j  ) + velz(i-1,j) + velz(i,j-1) + velz(i-1,j-1 ) )
    val_O =          velz(i,j  )
    val_E =          velz(i,j-1)
    
    wx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    wy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N = FOURTH*( tloc(i,j  ) + tloc(i+1,j) + tloc(i,j-1) + tloc(i+1,j-1 ) )
    val_S = FOURTH*( tloc(i,j  ) + tloc(i-1,j) + tloc(i,j-1) + tloc(i-1,j-1 ) )
    val_O =          tloc(i,j  )
    val_E =          tloc(i,j-1)
    
    tx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    ty = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    
    uu  = HALF*( velx(i,j) + velx(i,j-1) )
    vv  = HALF*( vely(i,j) + vely(i,j-1) )
    ww  = HALF*( velz(i,j) + velz(i,j-1) )
    !
    mmu    = HALF*( mu(i  ,j)   + mu(i  ,j-1) )
    lambda = HALF*( mu(i  ,j)   + mu(i  ,j-1) )*cpprandtl
    !
    fvrou2 = TWOTHIRD*mmu*( TWO*ux -     vy        )    
    fvrov2 =          mmu*(     uy +     vx        )    
    fvrow2 =          mmu*(                    wx  )    
    fvroe2 = (lambda*tx + uu*fvrou2 + vv*fvrov2 + ww * fvrow2)
    
    gvrou2 =          mmu*(     uy +     vx      )  
    gvrov2 = TWOTHIRD*mmu*(    -ux + TWO*vy      )  
    gvrow2 =          mmu*(                   wy )  
    gvroe2 = (lambda*ty + uu*gvrou2 + vv*gvrov2 + ww*gvrow2) 
    
    !!  fvrou2 = TWOTHIRD*mmu*( TWO*ux -     vy  -  wz  )    
    !!  fvrov2 =          mmu*(     uy +     vx         )    
    !!  fvrow2 =          mmu*(     uz           +  wx  )    
    !!  fvroe2 = (lambda*tx + uu*fvrou2 + vv*fvrov2 + ww * fvrow2)
    !!  
    !!  gvrou2 =          mmu*(     uy +     vx       )  
    !!  gvrov2 = TWOTHIRD*mmu*(    -ux + TWO*vy -  wz )  
    !!  gvrow2 =          mmu*(              vz +  wy )  
    !!  gvroe2 = (lambda*ty + uu*gvrou2 + vv*gvrov2 + ww*gvrow2) 
