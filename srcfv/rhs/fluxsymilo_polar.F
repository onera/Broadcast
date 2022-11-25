    !jdir
    !Euler part
    pw = ct0*p(i,j  ) +&
         ct1*p(i+i0,j+j0) 
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
    !mmu = mu(i,j)
    !impose d(var)/dr = 0
    mmu = ct0*mu(i   ,j   ) +&
          ct1*mu(i+i0,j+i0)
    !
    ! compute gradient-------------------------------------
    ! velysym = 0 ==> dvely/dxi = 2 * velcenter_i *nxi/vol
    ! dval.n = 0 
    ux = gradu(i,j,1)
    vx = ZERO
    wx = ZERO
    
    uy = ZERO
    vy = TWO * vely(i,j) * ny(i,j,1) * volf(i,j,1)
    wy = ZERO
    !ux = ZERO
    ! vx = ZERO
    ! wx = ZERO
    divu   = ux + vy !+ vv*distm1 vv= O on sym
    !
    ! set viscous fluxes at axisym ----------------------------
    ! dir normal to axisym
    
    
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
    hn(i,j,1,1) =   ZERO
    hn(i,j,2,1) = fxrov2*nx(i,j,1) - (fvrou2 * nx(i,j,1) + gvrou2 * ny(i,j,1))
    hn(i,j,3,1) = fxrov2*ny(i,j,1) - (fvrov2 * nx(i,j,1) + gvrov2 * ny(i,j,1))
    hn(i,j,4,1) =                  - (fvrow2 * nx(i,j,1) + gvrow2 * ny(i,j,1))
    hn(i,j,5,1) =   ZERO

