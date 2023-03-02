! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/

! =============================================================================
!          consistent fluxes for DNC3 2D
! =============================================================================
!
subroutine flux_num_dnc3_nowall_2d(residu,w,x0,y0,nx,ny,xc,yc,vol,volf,gh,cp,cv,prandtl,&
                            gam,rgaz,cs,muref,tref,s_suth,k2,k4,im,jm)
!
  implicit none
! variables for dimension -----------------------------------------
  integer :: im,jm,gh
! required arguments ----------------------------------------------
  real(8),intent(in) :: cp,cv,prandtl,gam,rgaz ! thermo
  real(8),intent(in) :: k2,k4 ! dissipation
  real(8),intent(in) :: cs,muref,tref,s_suth ! viscosity
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1      ),intent(in) :: x0
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1      ),intent(in) :: y0
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1, 2   ),intent(in) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1, 2   ),intent(in) :: ny
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: xc
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: yc
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 2   ),intent(in) :: volf
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: w
! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: residu
! Non-required arguments -------------------------------------------
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1, 5, 2) :: hn
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ) :: f,g
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: velx
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: vely
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: velz
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: tloc
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: p
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: mu
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     , 2) :: gradu,gradv,gradw,gradmu
  real(8),dimension(5)    :: deltar,deltal,t
  integer :: i,j,h
  real(8) :: nx_N, nx_S, nx_E, nx_O, volm1, ux,vx,wx,tx
  real(8) :: ny_N, ny_S, ny_E, ny_O, uy,vy,wy,ty
  real(8) :: val_N, val_S, val_E, val_O
  real(8) :: fxro1,fxro2, fxrou1,fxrou2,fxrov1,fxrov2,fxrow1,fxrow2,fxroe1,fxroe2
  real(8) :: dissro1,dissro2, dissrou1,dissrou2,dissrov1,dissrov2,dissrow1,dissrow2,dissroe1,dissroe2
  real(8) :: fvro1,fvro2, fvrou1,fvrou2,fvrov1,fvrov2,fvrow1,fvrow2,fvroe1,fvroe2
  real(8) :: gvro1,gvro2, gvrou1,gvrou2,gvrov1,gvrov2,gvrow1,gvrow2,gvroe1,gvroe2
  real(8) :: cpprandtl, mmu, lambda, uu, vv, ww, ro,rom1,htot,eloc,ec
  real(8) :: predro1,predrou1,predrov1,predrow1,predroe1, eps2,eps4
  real(8) :: predro2,predrou2,predrov2,predrow2,predroe2, rspec
  real(8) :: divu,divu2,vort2,dx,dy,dxm1,dym1,dxm2,dym2
  real(8) :: gui,gvi,gwi,gmui
  real(8) :: guj,gvj,gwj,gmuj
  real(8) :: rhom,rhomr,rhoml,rhom1l,c2l,c2r,rr,r,u,ur,ul,vr,vl,wr,wl,c2x,nx2,ny2
  real(8) :: ab, sq,ducros1,ducros2,k_sensor1,k_sensor2
  real(8) :: b1,b2,b3,b4,b5,c1,c2,c3,c4,c5,wiggle, denom,betas
  real(8) :: nxloc, nyloc, sn, invsn, sc1, sc2
  real(8) :: d1, d2, cvm1
  real(8) :: coef,omrr,test,diffro,diffrou,diffrov,diffrow,diffroe,v
  real(8) :: HALF,ONE,ZERO,TWO,TWOTHIRD,FOURTH,TWELFTH 
  real(8) :: TWENTYFOURTH,ccross
  real(8) :: mach2,alpha,vn2,cprim 
! -----------------------------------------------------------------
!
  HALF     = 0.5d0
  ONE      = 1.d0
  ZERO     = 0.d0
  TWO      = 2.d0
  TWOTHIRD = 2.d0/3.d0
  FOURTH   = 0.25d0
  TWELFTH  = 0.25d0/3.d0  
  TWENTYFOURTH = ONE/24.d0
  ccross = TWELFTH*0.0625d0
  
  cpprandtl = cp/prandtl
  cvm1      = ONE/cv 
  
! k2 = 1.d0
! k4 = 1.d0/12d0 already in coef for predictor
  
  eps2 = ZERO
  eps4 = k4
  
  diffro   = ZERO
  diffrou  = ZERO
  diffrov  = ZERO
  diffrow  = ZERO
  diffroe  = ZERO
  
  dissro1   = ZERO
  dissrou1  = ZERO
  dissrov1  = ZERO
  dissrow1  = ZERO
  dissroe1  = ZERO
  
  dissro2   = ZERO
  dissrou2  = ZERO
  dissrov2  = ZERO
  dissrow2  = ZERO
  dissroe2  = ZERO
  
  
! Coef for grad o4
  denom = 1.d0/12.d0
  
! b1 =  8.d0 * denom
! b2 = -       denom
  b1 = 0.5d0
  
!expression for FV

! c2 = b2
! c1 = c2 + b1
  
  c1 =  7.d0 * denom
  c2 = -       denom
  
!predictor
  
  d1 = 3.d0 * denom
  d2 =        denom
  
! Primitives
  betas = muref*(tref + cs)/(sqrt(tref)*tref) ! for Sutherland
  
!$AD II-LOOP
      do j=1-gh,jm+gh
!$AD II-LOOP
!DIR$ IVDEP
      do i=1-gh,im+gh

    ro = w(i,j,1)
    rom1 = ONE/ro
    velx(i,j) = w(i,j,2) * rom1
    vely(i,j) = w(i,j,3) * rom1
    velz(i,j) = w(i,j,4) * rom1
!
    ec  = HALF*( velx(i,j)*velx(i,j) &
               + vely(i,j)*vely(i,j) &
               + velz(i,j)*velz(i,j))

! ec  = HALF*( velx(i,j)*velx(i,j) &
!            + vely(i,j)*vely(i,j))
!
    eloc = (w(i,j,5) - ec*ro)*rom1
!
    tloc(i,j) = eloc*cvm1
!
    p(i,j)  = (gam-ONE)*ro*eloc
!  p(i,j)  = ro*rgaz*tloc(i,j)
    
    htot= (w(i,j,5) + p(i,j))*rom1
    
    f(i,j,1)  = w(i,j,2) 
    f(i,j,2)  = w(i,j,2) * velx(i,j) + p(i,j)
    f(i,j,3)  = w(i,j,2) * vely(i,j)
    f(i,j,4)  = w(i,j,2) * velz(i,j)
    f(i,j,5)  = w(i,j,2) * htot
    
    g(i,j,1)  = w(i,j,3) 
    g(i,j,2)  = w(i,j,3) * velx(i,j) 
    g(i,j,3)  = w(i,j,3) * vely(i,j) + p(i,j)
    g(i,j,4)  = w(i,j,3) * velz(i,j)
    g(i,j,5)  = w(i,j,3) * htot
    
!h(i,j,1)  = w(i,j,4)
!h(i,j,2)  = w(i,j,4) * velx(i,j)
!h(i,j,3)  = w(i,j,4) * vely(i,j)
!h(i,j,4)  = w(i,j,4) * velz(i,j) + p(i,j)
!h(i,j,5)  = w(i,j,4) * htot
    
        mu(i,j) = betas/(tloc(i,j) + s_suth) * sqrt(tloc(i,j)) * tloc(i,j)
    
    
  
      enddo
      enddo
  
  
! Work on interior domain minus one cell
!$AD II-LOOP
  do j = 1, jm
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1, im
!idir
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) )
                               
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) )
       
                           
       
!jdir
       guj   = ( b1*(velx(i,j+1   ) - velx(i,j-1   )) )
       
       gvj   = ( b1*(vely(i,j+1   ) - vely(i,j-1   )) )
       
       gwj   = ( b1*(velz(i,j+1   ) - velz(i,j-1   )) )

       gmuj  = ( b1*(  mu(i,j+1   ) -   mu(i,j-1   )) )
       
                           
       
      volm1 = ONE/vol(i,j)
      dxm1 = HALF*(nx(i,j,1)+nx(i+1,j  ,1))*volm1
      dxm2 = HALF*(nx(i,j,2)+nx(i  ,j+1,2))*volm1
      
      dym1 = HALF*(ny(i,j,1)+ny(i+1,j  ,1))*volm1
      dym2 = HALF*(ny(i,j,2)+ny(i  ,j+1,2))*volm1

       gradu( i,j,1)  = dxm1*gui   + dxm2*guj
       gradv( i,j,1)  = dxm1*gvi   + dxm2*gvj
       gradw( i,j,1)  = dxm1*gwi   + dxm2*gwj 
       gradmu(i,j,1)  = dxm1*gmui  + dxm2*gmuj
                         
       gradu( i,j,2)  = dym1*gui   + dym2*guj
       gradv( i,j,2)  = dym1*gvi   + dym2*gvj
       gradw( i,j,2)  = dym1*gwi   + dym2*gwj 
       gradmu(i,j,2)  = dym1*gmui  + dym2*gmuj
                         
                         
  enddo
  enddo
  
  
! trick for ducros in dissipation
  do i=1,2
  do h=1,gh
    gradu(:,1-h,:) =  TWO*gradu(:,2-h,:) - gradu(:,3-h,:)
    gradv(:,1-h,:) =  TWO*gradv(:,2-h,:) - gradv(:,3-h,:)
    gradw(:,1-h,:) =  TWO*gradw(:,2-h,:) - gradw(:,3-h,:)
  
    gradu(:,jm+h,:) = TWO*gradu(:,jm-1+h,:) - gradu(:,jm-2+h,:)
    gradv(:,jm+h,:) = TWO*gradv(:,jm-1+h,:) - gradv(:,jm-2+h,:)
    gradw(:,jm+h,:) = TWO*gradw(:,jm-1+h,:) - gradw(:,jm-2+h,:)
  
    gradu(1-h,:,:) =  TWO*gradu(2-h,:,:) - gradu(3-h,:,:)
    gradv(1-h,:,:) =  TWO*gradv(2-h,:,:) - gradv(3-h,:,:)
    gradw(1-h,:,:) =  TWO*gradw(2-h,:,:) - gradw(3-h,:,:)
  
    gradu(im+h,:,:) = TWO*gradu(im-1+h,:,:) - gradu(im-2+h,:,:)
    gradv(im+h,:,:) = TWO*gradv(im-1+h,:,:) - gradv(im-2+h,:,:)
    gradw(im+h,:,:) = TWO*gradw(im-1+h,:,:) - gradw(im-2+h,:,:)
  enddo
  enddo
  
  
!$AD II-LOOP
  do j = 1, jm+1
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1 , im + 1
!
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) ) * ny(i,j,1)
      
  
  fxro2    =  ( c1* (f(i  ,j    , 1) + f(i  ,j-1  , 1) ) +               &
                c2* (f(i  ,j+1  , 1) + f(i  ,j-2  , 1) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 1) + g(i  ,j-1  , 1) ) +               &
                c2* (g(i  ,j+1  , 1) + g(i  ,j-2  , 1) ) ) * ny(i,j,2)
      
  fxrou2   =  ( c1* (f(i  ,j    , 2) + f(i  ,j-1  , 2) ) +               &
                c2* (f(i  ,j+1  , 2) + f(i  ,j-2  , 2) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 2) + g(i  ,j-1  , 2) ) +               &
                c2* (g(i  ,j+1  , 2) + g(i  ,j-2  , 2) ) ) * ny(i,j,2)
      
  fxrov2   =  ( c1* (f(i  ,j    , 3) + f(i  ,j-1  , 3) ) +               &
                c2* (f(i  ,j+1  , 3) + f(i  ,j-2  , 3) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 3) + g(i  ,j-1  , 3) ) +               &
                c2* (g(i  ,j+1  , 3) + g(i  ,j-2  , 3) ) ) * ny(i,j,2)
      
  fxrow2   =  ( c1* (f(i  ,j    , 4) + f(i  ,j-1  , 4) ) +               &
                c2* (f(i  ,j+1  , 4) + f(i  ,j-2  , 4) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 4) + g(i  ,j-1  , 4) ) +               &
                c2* (g(i  ,j+1  , 4) + g(i  ,j-2  , 4) ) ) * ny(i,j,2)
      
  fxroe2   =  ( c1* (f(i  ,j    , 5) + f(i  ,j-1  , 5) ) +               &
                c2* (f(i  ,j+1  , 5) + f(i  ,j-2  , 5) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 5) + g(i  ,j-1  , 5) ) +               &
                c2* (g(i  ,j+1  , 5) + g(i  ,j-2  , 5) ) ) * ny(i,j,2)
      
  
    predro1   =      d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) 
                   
    predrou1  =      d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2)
                  
    predrov1  =      d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3)
                   
    predrow1  =      d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4)
                   
    predroe1  =      d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5)
                                  
                   
                                                     
    predro2   =      d2*w(i, j- 2,1) &
                   - d1*w(i, j- 1,1) &
                   + d1*w(i, j   ,1) &
                   - d2*w(i, j+ 1,1)
                   
    predrou2  =      d2*w(i, j- 2,2) &
                   - d1*w(i, j- 1,2) &
                   + d1*w(i, j   ,2) &
                   - d2*w(i, j+ 1,2)
                  
    predrov2  =      d2*w(i, j- 2,3) &
                   - d1*w(i, j- 1,3) &
                   + d1*w(i, j   ,3) &
                   - d2*w(i, j+ 1,3)
                   
    predrow2  =      d2*w(i, j- 2,4) &
                   - d1*w(i, j- 1,4) &
                   + d1*w(i, j   ,4) &
                   - d2*w(i, j+ 1,4)
                   
    predroe2  =      d2*w(i, j- 2,5) &
                   - d1*w(i, j- 1,5) &
                   + d1*w(i, j   ,5) &
                   - d2*w(i, j+ 1,5)
    
                 
! grad for viscous fluxes o2 - 3p
!

    volm1 = volf(i,j,1)
    
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
! #include "rhs/spectralradius_i.F"
! correction Rossow JCP 2000
!#include "rhs/spectralradiusRossow_i.F"
!1st direction
    
                 
    rhomr     = w(i,j,1)
    ur        = w(i,j,2)/rhomr
    vr        = w(i,j,3)/rhomr
    c2r       = gam*rgaz*tloc(i,j)
!
    rhoml     = w(i-1,j,1)
    ul        = w(i-1,j,2)/rhoml
    vl        = w(i-1,j,3)/rhoml
    c2l       = gam*rgaz*tloc(i-1,j)
!
    r         = sqrt( rhomr/rhoml)
    rr        = ONE/(ONE+r)
    omrr      = ONE-rr
!
    u         =  ul*rr + ur*omrr
    v         =  vl*rr + vr*omrr
!
    c2x       = c2l*rr + c2r*omrr
    nx2       = nx(i,j,1)*nx(i,j,1)+ny(i,j,1)*ny(i,j,1)
!
    ab        = abs(nx(i,j,1)*u+ny(i,j,1)*v)
    sq        = sqrt(c2x*nx2)
!
    rspec     = ab + sq            
                            
    
!
!1st direction
    
    k_sensor1 = ABS(p(i-1,j) - TWO*p(i,j) + p(i+1,j)) / &
                ABS(p(i-1,j) + TWO*p(i,j) + p(i+1,j))
                
    k_sensor2 = ABS(p(i-2,j) - TWO*p(i-1,j) + p(i,j)) / &
                ABS(p(i-2,j) + TWO*p(i-1,j) + p(i,j))
    
                
    divu      = (gradu(i,j,1)+gradv(i,j,2))            
    divu2     = divu * divu
    vort2     = (gradv(i,j,1)-gradu(i,j,2)) * (gradv(i,j,1)-gradu(i,j,2))
    ducros1   = divu2/(divu2+vort2+1d-15)
! dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sq+1.d-15)*divu) )
    dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sqrt(c2r*nx2)+1.d-15)*divu) )
 
     
    divu      = (gradu(i-1,j,1)+gradv(i-1,j,2))
    divu2     = divu * divu
    vort2     = (gradv(i-1,j,1)-gradu(i-1,j,2)) * (gradv(i-1,j,1)-gradu(i-1,j,2))
    ducros2   = divu2/(divu2+vort2+1d-15)
! dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i-1,j)/(sq+1.d-15)*divu) )
    dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i-1,j)/(sqrt(c2l*nx2)+1.d-15)*divu) )                            
    
!coef      = MAX( k_sensor1*ducros1*dxm1, &
!                 k_sensor2*ducros2*dxm2)
    
    coef      = max(k_sensor1, k_sensor2) * max(ducros1, ducros2) * max(dxm1, dxm2)
   
!     ! extension of shock capturing

    
!
!k_sensor1 = ABS(p(i,j) - TWO*p(i+1,j) + p(i+2,j)) / &
!            ABS(p(i,j) + TWO*p(i+1,j) + p(i+2,j))
!
!k_sensor2 = ABS(p(i-3,j) - TWO*p(i-2,j) + p(i-1,j)) / &
!            ABS(p(i-3,j) + TWO*p(i-2,j) + p(i-1,j))
!
!divu      = gradu(i+1,j,1)+gradv(i+1,j,2)
!divu2     = divu*divu
!vort2     = (gradv(i+1,j,1)-gradu(i+1,j,2)) * (gradv(i+1,j,1)-gradu(i+1,j,2))
!ducros1   = divu2/(divu2+vort2+1d-15)
!dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i+1,j)/(sq+1.d-15)*divu) )
!
!
!divu      = gradu(i-2,j,1)+gradv(i-2,j,2)
!divu2     = divu*divu
!vort2     = (gradv(i-2,j,1)-gradu(i-2,j,2)) * (gradv(i-2,j,1)-gradu(i-2,j,2))
!ducros2   = divu2/(divu2+vort2+1d-15)
!dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i-2,j)/(sq+1.d-15)*divu) )
    
!coef      = MAX( k_sensor1*ducros1*dxm1, &
!                 k_sensor2*ducros2*dxm2, coef)
    
!coef      = MAX(sensor(i-1,j  ,1), sensor(i,j  ,1), &
!                sensor(i-1,j-1,1), sensor(i,j-1,1), &
!                sensor(i-1,j+1,1), sensor(i,j+1,1)  )
!
!coef      = MAX(sensor(i-3,j,1), sensor(i-2,j,1), &
!                sensor(i-1,j,1), sensor(i  ,j,1), &
!                sensor(i+1,j,1), sensor(i+2,j,1)  )
!
!coef      = MAX(sensor(i-2,j,1), &
!                sensor(i-1,j,1), sensor(i  ,j,1), &
!                sensor(i+1,j,1))
!
!rspec      = rconv(i,j,1)
    
!
    eps2      = k2*coef
    eps4      = MAX(ZERO,k4-eps2 * 12.d0) ! to follow Sciacovelli CF 2021

    diffro   = HALF * (w(i,j,1) - w(i-1,j,1))
    diffrou  = HALF * (w(i,j,2) - w(i-1,j,2))
    diffrov  = HALF * (w(i,j,3) - w(i-1,j,3))
    diffrow  = HALF * (w(i,j,4) - w(i-1,j,4))
    diffroe  = HALF * (w(i,j,5) - w(i-1,j,5))
!
!
!     if (eps4.gt.1.d-12) then
! #include "rhs/wiggle_diri.F"
!     endif
    
    dissro1  = rspec * (eps2*diffro  + eps4*predro1 )
    dissrou1 = rspec * (eps2*diffrou + eps4*predrou1)
    dissrov1 = rspec * (eps2*diffrov + eps4*predrov1)
    dissrow1 = rspec * (eps2*diffrow + eps4*predrow1)
    dissroe1 = rspec * (eps2*diffroe + eps4*predroe1)
    

! #include "rhs/spectralradius_j.F"
! correction Rossow JCP 2000
!#include "rhs/spectralradiusRossow_j.F"
!1st direction
    
                 
    rhomr     = w(i,j,1)
    ur        = w(i,j,2)/rhomr
    vr        = w(i,j,3)/rhomr
    c2r       = gam*rgaz*tloc(i,j)
!
    rhoml     = w(i,j-1,1)
    ul        = w(i,j-1,2)/rhoml
    vl        = w(i,j-1,3)/rhoml
    c2l       = gam*rgaz*tloc(i,j-1)
!
    r         = sqrt( rhomr/rhoml)
    rr        = ONE/(ONE+r)
    omrr      = ONE-rr
!
    u         =  ul*rr + ur*omrr
    v         =  vl*rr + vr*omrr
!
    c2x       = c2l*rr + c2r*omrr
    nx2       = nx(i,j,2)*nx(i,j,2)+ny(i,j,2)*ny(i,j,2)
!
    ab        = abs(nx(i,j,2)*u+ny(i,j,2)*v)
    sq        = sqrt(c2x*nx2)
!
    rspec     = ab + sq            
                            
    
!
    
! 2nd direction

    k_sensor1 = ABS(p(i,j-1) - TWO*p(i,j) + p(i,j+1)) / &
                ABS(p(i,j-1) + TWO*p(i,j) + p(i,j+1))
                
    k_sensor2 = ABS(p(i,j-2) - TWO*p(i,j-1) + p(i,j)) / &
                ABS(p(i,j-2) + TWO*p(i,j-1) + p(i,j))
                
!     ducros1 is done in dissipation_ducros_x

    divu      = (gradu(i,j,1)+gradv(i,j,2))            
    divu2     = divu * divu
    vort2     = (gradv(i,j,1)-gradu(i,j,2)) * (gradv(i,j,1)-gradu(i,j,2))
    ducros1   = divu2/(divu2+vort2+1d-15)
! dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sq+1.d-15)*divu) )
    dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j)/(sqrt(c2r*nx2)+1.d-15)*divu) )

    divu      = gradu(i,j-1,1)+gradv(i,j-1,2)                           
    divu2     = divu * divu
    vort2     = (gradv(i,j-1,1)-gradu(i,j-1,2)) * (gradv(i,j-1,1)-gradu(i,j-1,2))
    ducros2   = divu2/(divu2+vort2+1.d-15)
    
!sq = sound/dx
! dxm2       = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j-1)/(sq+1.d-15)*divu) )
    dxm2       = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j-1)/(sqrt(c2l*nx2)+1.d-15)*divu) )              
    
!coef      = MAX( k_sensor1*ducros1*dxm1, &
!                 k_sensor2*ducros2*dxm2)
    
    coef      = max(k_sensor1, k_sensor2) * max(ducros1, ducros2) * max(dxm1, dxm2)

    
! extension of shock capturing
    
!k_sensor1 = ABS(p(i,j) - TWO*p(i,j+1) + p(i,j+2)) / &
!            ABS(p(i,j) + TWO*p(i,j+1) + p(i,j+2))
!
!k_sensor2 = ABS(p(i,j-3) - TWO*p(i,j-2) + p(i,j-1)) / &
!            ABS(p(i,j-3) + TWO*p(i,j-2) + p(i,j-1))
!
!divu      = gradu(i,j+1,1)+gradv(i,j+1,2)
!divu2     = divu*divu
!vort2     = (gradv(i,j+1,1)-gradu(i,j+1,2)) * (gradv(i,j+1,1)-gradu(i,j+1,2))
!ducros1   = divu2/(divu2+vort2+1d-15)
!dxm1      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j+1)/(sq+1.d-15)*divu) )
!
!
!divu      = gradu(i,j-2,1)+gradv(i,j-2,2)
!divu2     = divu*divu
!vort2     = (gradv(i,j-2,1)-gradu(i,j-2,2)) * (gradv(i,j-2,1)-gradu(i,j-2,2))
!ducros2   = divu2/(divu2+vort2+1d-15)
!dxm2      = HALF*(ONE-tanh(2.5d0+10.d0*vol(i,j-2)/(sq+1.d-15)*divu) )
!
!coef      = MAX( k_sensor1*ducros1*dxm1, &
!                 k_sensor2*ducros2*dxm2, coef)

!coef      = MAX(sensor(i  ,j-1,2), sensor(i  ,j,2), &
!                sensor(i-1,j-1,2), sensor(i-1,j,2), &
!                sensor(i+1,j-1,2), sensor(i+1,j,2)  )
!
!coef      = MAX(sensor(i,j-3,2), sensor(i,j-2,2), &
!                sensor(i,j-1,2), sensor(i,j  ,2), &
!                sensor(i,j+1,2), sensor(i,j+2,2)  )
!
!coef      = MAX(sensor(i,j-2,2), &
!                sensor(i,j-1,2), sensor(i,j,2), &
!                sensor(i,j+1,2))
!
!rspec      = rconv(i,j,2)

    
    eps2      = k2*coef
    eps4      = MAX(ZERO, k4 - 12.d0 * eps2)  ! 12.d0 to follow Sciacovelli CF 2021
    

    diffro   = HALF * (w(i,j,1) - w(i,j-1,1))
    diffrou  = HALF * (w(i,j,2) - w(i,j-1,2))
    diffrov  = HALF * (w(i,j,3) - w(i,j-1,3))
    diffrow  = HALF * (w(i,j,4) - w(i,j-1,4))
    diffroe  = HALF * (w(i,j,5) - w(i,j-1,5))
!
!
!     if (eps4.gt.1.d-12) then
! #include "rhs/wiggle_dirj.F"
!     endif
    
!
    dissro2  = rspec * (eps2*diffro  + eps4*predro2 )
    dissrou2 = rspec * (eps2*diffrou + eps4*predrou2)
    dissrov2 = rspec * (eps2*diffrov + eps4*predrov2)
    dissrow2 = rspec * (eps2*diffrow + eps4*predrow2)
    dissroe2 = rspec * (eps2*diffroe + eps4*predroe2)
       sc1 = nx(i,j,1)
       sc2 = ny(i,j,1)
       sn  = sqrt(sc1*sc1 + sc2*sc2)
       invsn = ONE/sn
       nxloc = sc1*invsn
       nyloc = sc2*invsn
       
       hn(i,j,1,1) = fxro1  -  dissro1  
!
       hn(i,j,2,1) = fxrou1 - dissrou1 - (fvrou1 * nxloc + gvrou1 * nyloc)*sn 
!
       hn(i,j,3,1) = fxrov1 - dissrov1 - (fvrov1 * nxloc + gvrov1 * nyloc)*sn
!
       hn(i,j,4,1) = fxrow1 - dissrow1 - (fvrow1 * nxloc + gvrow1 * nyloc)*sn
!
       hn(i,j,5,1) = fxroe1 - dissroe1 - (fvroe1 * nxloc + gvroe1 * nyloc)*sn

       sc1 = nx(i,j,2)
       sc2 = ny(i,j,2)
       sn  = sqrt(sc1*sc1 + sc2*sc2)
       invsn = ONE/sn
       nxloc = sc1*invsn
       nyloc = sc2*invsn
       
       hn(i,j,1,2) = fxro2  -  dissro2  
!
       hn(i,j,2,2) = fxrou2 - dissrou2 - (fvrou2 * nxloc + gvrou2 * nyloc)*sn 
!
       hn(i,j,3,2) = fxrov2 - dissrov2 - (fvrov2 * nxloc + gvrov2 * nyloc)*sn
!
       hn(i,j,4,2) = fxrow2 - dissrow2 - (fvrow2 * nxloc + gvrow2 * nyloc)*sn
!
       hn(i,j,5,2) = fxroe2 - dissroe2 - (fvroe2 * nxloc + gvroe2 * nyloc)*sn

!
  enddo
  enddo
  
  
! fluxes balance

!$AD II-LOOP
  do j = 1 , jm
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1 , im

residu(i  ,j,1) = - ( hn(i+1,j  ,1,1) - hn(i,j,1,1) ) &
                  - ( hn(i  ,j+1,1,2) - hn(i,j,1,2) ) 
                  
residu(i  ,j,2) = - ( hn(i+1,j  ,2,1) - hn(i,j,2,1) ) &
                  - ( hn(i  ,j+1,2,2) - hn(i,j,2,2) ) 

residu(i  ,j,3) = - ( hn(i+1,j  ,3,1) - hn(i,j,3,1) ) &
                  - ( hn(i  ,j+1,3,2) - hn(i,j,3,2) ) 

residu(i  ,j,4) = - ( hn(i+1,j  ,4,1) - hn(i,j,4,1) ) &
                  - ( hn(i  ,j+1,4,2) - hn(i,j,4,2) ) 

residu(i  ,j,5) = - ( hn(i+1,j  ,5,1) - hn(i,j,5,1) ) &
                  - ( hn(i  ,j+1,5,2) - hn(i,j,5,2) ) 
                  
      
  enddo
  enddo
      


end subroutine flux_num_dnc3_nowall_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!