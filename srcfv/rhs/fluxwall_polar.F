    !jdir
    !Euler part
    pw = ct0*p(i,j  ) +&
         ct1*p(i,j+1) 
    !!  ! order 3
    !!  pw = ct0*p(i,j  ) +&
    !!       ct1*p(i,j+1) +&
    !!       ct2*p(i,j+2)
    ! fxro2  = ZERO
    ! fxrou2 = ZERO
    fxrov2 = pw 
    ! fxrow2 = ZERO
    ! fxroe2 = ZERO
    
    !ns part
    ! viscosity is taken at the first center
    mmu = mu(i,j)
    !
    ! compute gradient-------------------------------------
    ! vel_wall = 0 ==> dveli/dxi = 2 * velcenter_i *nxi/vol
    ux = TWO * velx(i,j) * nx(i,j,2) * volf(i,j,2)
    vx = TWO * vely(i,j) * nx(i,j,2) * volf(i,j,2)
    wx = TWO * velz(i,j) * nx(i,j,2) * volf(i,j,2)
    
    uy = TWO * velx(i,j) * ny(i,j,2) * volf(i,j,2)
    vy = TWO * vely(i,j) * ny(i,j,2) * volf(i,j,2)
    wy = TWO * velz(i,j) * ny(i,j,2) * volf(i,j,2)
    ! ux = ZERO
    ! vx = ZERO
    ! wx = ZERO
    divu   = ux + vy !+ vv*distm1 vv= O on wall
    !
    ! set viscous fluxes at wall----------------------------
    ! dir normal to the wall
    
    
    fvrou2 = mmu * (TWO     * ux - TWOTHIRD * divu)    
    fvrov2 = mmu * (          uy +     vx         )  
    fvrow2 = mmu * (                   wx         ) 
    !
    gvrou2 = mmu * (          uy +     vx         )
    gvrov2 = mmu * (TWO     * vy - TWOTHIRD * divu)
    gvrow2 = mmu * (                    wy        )    
    
    !
    fvroe2 = ZERO !adiab
    gvroe2 = ZERO !adiab
    !
    hn(i,j,1,2) =   ZERO
    hn(i,j,2,2) = fxrov2*nx(i,j,2) - (fvrou2 * nx(i,j,2) + gvrou2 * ny(i,j,2))
    hn(i,j,3,2) = fxrov2*ny(i,j,2) - (fvrov2 * nx(i,j,2) + gvrov2 * ny(i,j,2))
    hn(i,j,4,2) =                  - (fvrow2 * nx(i,j,2) + gvrow2 * ny(i,j,2))
    hn(i,j,5,2) =   ZERO

