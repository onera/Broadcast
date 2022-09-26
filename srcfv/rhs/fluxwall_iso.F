    !jdir
    !Euler part
    pw = ct0*p(i,j  ) +&
         ct1*p(i,j+1) 
    !pw = p(i,j)
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

    tx = (tloc(i,j) - twall) * nx(i,j,2) * volf(i,j,2)
    ty = (tloc(i,j) - twall) * ny(i,j,2) * volf(i,j,2)
    ! ux = ZERO
    ! vx = ZERO
    ! wx = ZERO
    !
    ! set viscous fluxes at wall----------------------------
    ! dir normal to the wall
    fvrou2 = TWOTHIRD*mmu*( TWO*ux -     vy        )    
    fvrov2 =          mmu*(     uy +     vx        )    
    fvrow2 =          mmu*(                    wx  )    
    fvroe2 = lambda*tx
    !
    gvrou2 =          mmu*(     uy +     vx       )
    gvrov2 = TWOTHIRD*mmu*(    -ux + TWO*vy       )
    gvrow2 =          mmu*(                  wy   )
    !
    gvroe2 = lambda*ty
    !
    hn(i,j,1,2) =   ZERO
    hn(i,j,2,2) = fxrov2*nx(i,j,2) - (fvrou2 * nx(i,j,2) + gvrou2 * ny(i,j,2))
    hn(i,j,3,2) = fxrov2*ny(i,j,2) - (fvrov2 * nx(i,j,2) + gvrov2 * ny(i,j,2))
    hn(i,j,4,2) =                  - (fvrow2 * nx(i,j,2) + gvrow2 * ny(i,j,2))
    hn(i,j,5,2) = - (fvroe2 * nx(i,j,2) + gvroe2 * ny(i,j,2))

