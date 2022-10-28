                 
    ! grad for viscous fluxes o2 - 3p
    !

    volm1   = volf(i,j,1)
    
    distm1  = HALF*(y0(i,j) + y0(i,j+1))-ym
    if (distm1.gt.1.d-15) then
        distm1 = ONE/distm1
    endif

    ! Normals
    
    nx_N =   HALF*( nx(i+1,j  ,1) + nx(i  ,j  ,1) )
    nx_S = - HALF*( nx(i-1,j  ,1) + nx(i  ,j  ,1) )
    nx_O =   HALF*( nx(i-1,j+1,2) + nx(i  ,j+1,2) )
    nx_E = - HALF*( nx(i-1,j  ,2) + nx(i  ,j  ,2) )
    !
    ny_N =   HALF*( ny(i+1,j  ,1) + ny(i  ,j  ,1) )
    ny_S = - HALF*( ny(i-1,j  ,1) + ny(i  ,j  ,1) )
    ny_O =   HALF*( ny(i-1,j+1,2) + ny(i  ,j+1,2) )
    ny_E = - HALF*( ny(i-1,j  ,2) + ny(i  ,j  ,2) )
    
    
    
    val_N =          velx(i  ,j)
    val_S =          velx(i-1,j)
    val_E = FOURTH*( velx(i  ,j) + velx(i,j-1) + velx(i-1,j) + velx(i-1,j-1) )
    val_O = FOURTH*( velx(i  ,j) + velx(i,j+1) + velx(i-1,j) + velx(i-1,j+1) )
    
    ux = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    uy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    
    val_N =          vely(i  ,j)
    val_S =          vely(i-1,j)
    val_E = FOURTH*( vely(i  ,j) + vely(i,j-1) + vely(i-1,j) + vely(i-1,j-1) )
    val_O = FOURTH*( vely(i  ,j) + vely(i,j+1) + vely(i-1,j) + vely(i-1,j+1) )
    
    vx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    vy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
    
    val_N =          velz(i  ,j)
    val_S =          velz(i-1,j)
    val_E = FOURTH*( velz(i  ,j) + velz(i,j-1) + velz(i-1,j) + velz(i-1,j-1) )
    val_O = FOURTH*( velz(i  ,j) + velz(i,j+1) + velz(i-1,j) + velz(i-1,j+1) )
    
    wx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    wy = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1
                 
    val_N =          tloc(i  ,j)
    val_S =          tloc(i-1,j)
    val_E = FOURTH*( tloc(i  ,j) + tloc(i,j-1) + tloc(i-1,j) + tloc(i-1,j-1) )
    val_O = FOURTH*( tloc(i  ,j) + tloc(i,j+1) + tloc(i-1,j) + tloc(i-1,j+1) )
    !
    tx = ( val_N*nx_N + val_S*nx_S + val_O*nx_O + val_E*nx_E )*volm1
    ty = ( val_N*ny_N + val_S*ny_S + val_O*ny_O + val_E*ny_E )*volm1  
    !
    ! Computation of viscous fluxes
    !
    uu   = HALF*( velx(i,j) + velx(i-1,j) )
    vv   = HALF*( vely(i,j) + vely(i-1,j) )
    ww   = HALF*( velz(i,j) + velz(i-1,j) )
    !
    mmu    = HALF*(mu(i  ,j  ) + mu(i-1,j  ))
    lambda = HALF*(mu(i  ,j  ) + mu(i-1,j  ))*cpprandtl
    divu   = ux + vy + vv*distm1
    !
    fvrou1 = mmu * (TWO     * ux - TWOTHIRD * divu)    
    fvrov1 = mmu * (          uy +     vx         )  
    fvrow1 = mmu * (                   wx         ) 
    fvroe1 = (lambda*tx + uu*fvrou1 + vv*fvrov1 + ww * fvrow1) 
    !
    gvrou1 = mmu * (          uy +     vx         )
    gvrov1 = mmu * (TWO     * vy - TWOTHIRD * divu)
    gvrow1 = mmu * (          wy - ww*distm1      )    
    gvroe1 = lambda*ty + uu*gvrou1 + vv*gvrov1 + ww * gvrow1
    
    
