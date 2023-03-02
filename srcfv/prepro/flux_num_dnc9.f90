! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/

! =============================================================================
!          consistent fluxes for DNC9 2D
! =============================================================================
!
subroutine flux_num_dnc9_2d(residu,w,x0,y0,nx,ny,xc,yc,vol,volf,gh,cp,cv,prandtl,&
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
  real(8) :: fvro1,fvro2, fvrou1,fvrou2,fvrov1,fvrov2,fvrow1,fvrow2,fvroe1,fvroe2
  real(8) :: gvro1,gvro2, gvrou1,gvrou2,gvrov1,gvrov2,gvrow1,gvrow2,gvroe1,gvroe2
  real(8) :: dissro1,dissro2, dissrou1,dissrou2,dissrov1,dissrov2,dissrow1,dissrow2,dissroe1,dissroe2
  real(8) :: cpprandtl, mmu, lambda, uu, vv, ww, ro,rom1,htot,eloc,ec
  real(8) :: predro1,predrou1,predrov1,predrow1,predroe1, eps2,eps4
  real(8) :: predro2,predrou2,predrov2,predrow2,predroe2, rspec
  real(8) :: nxloc, nyloc, sn, invsn, sc1, sc2
  real(8) :: divu,divu2,vort2,dxm1,dym1,dxm2,dym2
  real(8) :: gui,gvi,gwi,gmui
  real(8) :: guj,gvj,gwj,gmuj
  real(8) :: rhom,rhomr,rhoml,rhom1l,c2l,c2r,rr,r,u,ur,ul,vr,vl,wr,wl,c2x,nx2,ny2
  real(8) :: ab, sq,ducros1,ducros2,k_sensor1,k_sensor2
  real(8) :: b1,b2,b3,b4,b5,c1,c2,c3,c4,c5,d1,d2,d3,d4,d5,wiggle,denom,betas
  real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 ! for bnd off-centered schemes
  real(8) :: c9d0,c9d1,c9d2,c9d3,c9d4,c9d5,c9d6,c9d7,c9d8,c9d9,c9d10
  real(8) :: c7d0,c7d1,c7d2,c7d3,c7d4,c7d5,c7d6,c7d7,c7d8,c7d9,c7d10
  real(8) :: c5d0,c5d1,c5d2,c5d3,c5d4,c5d5,c5d6,c5d7,c5d8,c5d9,c5d10
  real(8) :: c3d0,c3d1,c3d2,c3d3,c3d4,c3d5,c3d6,c3d7,c3d8,c3d9,c3d10
  real(8) :: pw,ct0,ct1,ct2,coef,omrr,test,diffro,diffrou,diffrov,diffrow,diffroe,v
  real(8) :: HALF,ONE,ZERO,TWO,TWOTHIRD,FOURTH,TWELFTH,cvm1
  real(8) :: TWENTYFOURTH,ccross
  real(8) :: mach2,alpha,vn2,cprim
! -----------------------------------------------------------------
!
  HALF     = 0.5d0
  ONE      = 1.d0
  ZERO     = 0.d0
  TWO      = 2.d0
  TWOTHIRD = TWO/3.d0
  FOURTH   = 0.25d0
  TWELFTH  = 0.25d0/3.d0 
  TWENTYFOURTH = ONE/24.d0
  ccross = TWELFTH*0.0625d0 
 
  cpprandtl = cp/prandtl
  cvm1      = ONE/cv
  
! k2 = 1.d0
! k4 = 1.d0/1260.d0 already in coef for predictor
  
  eps2 = ZERO
  eps4 = k4
  
  diffro   = ZERO
  diffrou  = ZERO
  diffrov  = ZERO
  diffrow  = ZERO
  diffroe  = ZERO
!
  
! Coef for grad o10
! b1 =  5.d0/6.d0
! b2 = -5.d0/21.d0
! b3 =  5.d0/84.d0
! b4 = -5.d0/504.d0
! b5 =  1.d0/1260.d0
  denom = 1.d0/840.d0
!Coef for grad o8
  b1 =  672.d0 * denom! 4.d0/5.d0
  b2 = -168.d0 * denom!-1.d0/5.d0
  b3 =  32.d0  * denom! 4.d0/105.d0
  b4 = - 3.d0  * denom!-1.d0/280.d0
  
!expression for FV
! c5 = a5
! c4 = c5 + a4
! c3 = c4 + a3
! c2 = c3 + a2
! c1 = c2 + a1
  denom = 1.d0/2520.d0
!
  c1 =  1627.d0 * denom
  c2 = - 473.d0 * denom
  c3 =   127.d0 * denom
  c4 = -  23.d0 * denom
  c5 =     2.d0 * denom
  
! coef for predictor
  denom = 1.d0/1260.d0
!
  d1 = 126.d0 * denom
  d2 =  84.d0 * denom
  d3 =  36.d0 * denom
  d4 =   9.d0 * denom
  d5 =          denom
  
  
  betas = muref*(tref + cs)/(sqrt(tref)*tref)
    
! Primitives
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
  
  
!$AD II-LOOP
  do j = 1, jm
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1, im
!idir
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) &
               + b2*(velx(i+2,j     ) - velx(i-2,j     )) &
               + b3*(velx(i+3,j     ) - velx(i-3,j     )) &
               + b4*(velx(i+4,j     ) - velx(i-4,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) &
               + b2*(vely(i+2,j     ) - vely(i-2,j     )) &
               + b3*(vely(i+3,j     ) - vely(i-3,j     )) &
               + b4*(vely(i+4,j     ) - vely(i-4,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) &
               + b2*(velz(i+2,j     ) - velz(i-2,j     )) &
               + b3*(velz(i+3,j     ) - velz(i-3,j     )) &
               + b4*(velz(i+4,j     ) - velz(i-4,j     )) )
                               
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) &
               + b2*(  mu(i+2,j     ) -   mu(i-2,j     )) &
               + b3*(  mu(i+3,j     ) -   mu(i-3,j     )) &
               + b4*(  mu(i+4,j     ) -   mu(i-4,j     )) )
       
                           
       
!jdir
       guj   = ( b1*(velx(i,j+1   ) - velx(i,j-1   )) &
               + b2*(velx(i,j+2   ) - velx(i,j-2   )) &
               + b3*(velx(i,j+3   ) - velx(i,j-3   )) &
               + b4*(velx(i,j+4   ) - velx(i,j-4   )) )
       
       gvj   = ( b1*(vely(i,j+1   ) - vely(i,j-1   )) &
               + b2*(vely(i,j+2   ) - vely(i,j-2   )) &
               + b3*(vely(i,j+3   ) - vely(i,j-3   )) &
               + b4*(vely(i,j+4   ) - vely(i,j-4   )) )
       
       gwj   = ( b1*(velz(i,j+1   ) - velz(i,j-1   )) &
               + b2*(velz(i,j+2   ) - velz(i,j-2   )) &
               + b3*(velz(i,j+3   ) - velz(i,j-3   )) &
               + b4*(velz(i,j+4   ) - velz(i,j-4   )) )

       gmuj  = ( b1*(  mu(i,j+1   ) -   mu(i,j-1   )) &
               + b2*(  mu(i,j+2   ) -   mu(i,j-2   )) &
               + b3*(  mu(i,j+3   ) -   mu(i,j-3   )) &
               + b4*(  mu(i,j+4   ) -   mu(i,j-4   )) )
       
                           
       
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
  do j= 6,jm + 1
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1 , im + 1
!
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) +               & 
                c3* (f(i+2,j  , 1) + f(i-3,j  , 1) ) +               &
                c4* (f(i+3,j  , 1) + f(i-4,j  , 1) ) +               &
                c5* (f(i+4,j  , 1) + f(i-5,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) +               & 
                c3* (g(i+2,j  , 1) + g(i-3,j  , 1) ) +               &
                c4* (g(i+3,j  , 1) + g(i-4,j  , 1) ) +               &
                c5* (g(i+4,j  , 1) + g(i-5,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) +               & 
                c3* (f(i+2,j  , 2) + f(i-3,j  , 2) ) +               &
                c4* (f(i+3,j  , 2) + f(i-4,j  , 2) ) +               &
                c5* (f(i+4,j  , 2) + f(i-5,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) +               & 
                c3* (g(i+2,j  , 2) + g(i-3,j  , 2) ) +               &
                c4* (g(i+3,j  , 2) + g(i-4,j  , 2) ) +               &
                c5* (g(i+4,j  , 2) + g(i-5,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) +               & 
                c3* (f(i+2,j  , 3) + f(i-3,j  , 3) ) +               &
                c4* (f(i+3,j  , 3) + f(i-4,j  , 3) ) +               &
                c5* (f(i+4,j  , 3) + f(i-5,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) +               & 
                c3* (g(i+2,j  , 3) + g(i-3,j  , 3) ) +               &
                c4* (g(i+3,j  , 3) + g(i-4,j  , 3) ) +               &
                c5* (g(i+4,j  , 3) + g(i-5,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) +               & 
                c3* (f(i+2,j  , 4) + f(i-3,j  , 4) ) +               &
                c4* (f(i+3,j  , 4) + f(i-4,j  , 4) ) +               &
                c5* (f(i+4,j  , 4) + f(i-5,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) +               & 
                c3* (g(i+2,j  , 4) + g(i-3,j  , 4) ) +               &
                c4* (g(i+3,j  , 4) + g(i-4,j  , 4) ) +               &
                c5* (g(i+4,j  , 4) + g(i-5,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) +               & 
                c3* (f(i+2,j  , 5) + f(i-3,j  , 5) ) +               &
                c4* (f(i+3,j  , 5) + f(i-4,j  , 5) ) +               &
                c5* (f(i+4,j  , 5) + f(i-5,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) +               & 
                c3* (g(i+2,j  , 5) + g(i-3,j  , 5) ) +               &
                c4* (g(i+3,j  , 5) + g(i-4,j  , 5) ) +               &
                c5* (g(i+4,j  , 5) + g(i-5,j  , 5) ) ) * ny(i,j,1)
      
  
  fxro2    =  ( c1* (f(i  ,j    , 1) + f(i  ,j-1  , 1) ) +               &
                c2* (f(i  ,j+1  , 1) + f(i  ,j-2  , 1) ) +               & 
                c3* (f(i  ,j+2  , 1) + f(i  ,j-3  , 1) ) +               &
                c4* (f(i  ,j+3  , 1) + f(i  ,j-4  , 1) ) +               &
                c5* (f(i  ,j+4  , 1) + f(i  ,j-5  , 1) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 1) + g(i  ,j-1  , 1) ) +               &
                c2* (g(i  ,j+1  , 1) + g(i  ,j-2  , 1) ) +               & 
                c3* (g(i  ,j+2  , 1) + g(i  ,j-3  , 1) ) +               &
                c4* (g(i  ,j+3  , 1) + g(i  ,j-4  , 1) ) +               &
                c5* (g(i  ,j+4  , 1) + g(i  ,j-5  , 1) ) ) * ny(i,j,2)
      
  fxrou2   =  ( c1* (f(i  ,j    , 2) + f(i  ,j-1  , 2) ) +               &
                c2* (f(i  ,j+1  , 2) + f(i  ,j-2  , 2) ) +               & 
                c3* (f(i  ,j+2  , 2) + f(i  ,j-3  , 2) ) +               &
                c4* (f(i  ,j+3  , 2) + f(i  ,j-4  , 2) ) +               &
                c5* (f(i  ,j+4  , 2) + f(i  ,j-5  , 2) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 2) + g(i  ,j-1  , 2) ) +               &
                c2* (g(i  ,j+1  , 2) + g(i  ,j-2  , 2) ) +               & 
                c3* (g(i  ,j+2  , 2) + g(i  ,j-3  , 2) ) +               &
                c4* (g(i  ,j+3  , 2) + g(i  ,j-4  , 2) ) +               &
                c5* (g(i  ,j+4  , 2) + g(i  ,j-5  , 2) ) ) * ny(i,j,2)
      
  fxrov2   =  ( c1* (f(i  ,j    , 3) + f(i  ,j-1  , 3) ) +               &
                c2* (f(i  ,j+1  , 3) + f(i  ,j-2  , 3) ) +               & 
                c3* (f(i  ,j+2  , 3) + f(i  ,j-3  , 3) ) +               &
                c4* (f(i  ,j+3  , 3) + f(i  ,j-4  , 3) ) +               &
                c5* (f(i  ,j+4  , 3) + f(i  ,j-5  , 3) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 3) + g(i  ,j-1  , 3) ) +               &
                c2* (g(i  ,j+1  , 3) + g(i  ,j-2  , 3) ) +               & 
                c3* (g(i  ,j+2  , 3) + g(i  ,j-3  , 3) ) +               &
                c4* (g(i  ,j+3  , 3) + g(i  ,j-4  , 3) ) +               &
                c5* (g(i  ,j+4  , 3) + g(i  ,j-5  , 3) ) ) * ny(i,j,2)
      
  fxrow2   =  ( c1* (f(i  ,j    , 4) + f(i  ,j-1  , 4) ) +               &
                c2* (f(i  ,j+1  , 4) + f(i  ,j-2  , 4) ) +               & 
                c3* (f(i  ,j+2  , 4) + f(i  ,j-3  , 4) ) +               &
                c4* (f(i  ,j+3  , 4) + f(i  ,j-4  , 4) ) +               &
                c5* (f(i  ,j+4  , 4) + f(i  ,j-5  , 4) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 4) + g(i  ,j-1  , 4) ) +               &
                c2* (g(i  ,j+1  , 4) + g(i  ,j-2  , 4) ) +               & 
                c3* (g(i  ,j+2  , 4) + g(i  ,j-3  , 4) ) +               &
                c4* (g(i  ,j+3  , 4) + g(i  ,j-4  , 4) ) +               &
                c5* (g(i  ,j+4  , 4) + g(i  ,j-5  , 4) ) ) * ny(i,j,2)
      
  fxroe2   =  ( c1* (f(i  ,j    , 5) + f(i  ,j-1  , 5) ) +               &
                c2* (f(i  ,j+1  , 5) + f(i  ,j-2  , 5) ) +               & 
                c3* (f(i  ,j+2  , 5) + f(i  ,j-3  , 5) ) +               &
                c4* (f(i  ,j+3  , 5) + f(i  ,j-4  , 5) ) +               &
                c5* (f(i  ,j+4  , 5) + f(i  ,j-5  , 5) ) ) * nx(i,j,2) + &
              ( c1* (g(i  ,j    , 5) + g(i  ,j-1  , 5) ) +               &
                c2* (g(i  ,j+1  , 5) + g(i  ,j-2  , 5) ) +               & 
                c3* (g(i  ,j+2  , 5) + g(i  ,j-3  , 5) ) +               &
                c4* (g(i  ,j+3  , 5) + g(i  ,j-4  , 5) ) +               &
                c5* (g(i  ,j+4  , 5) + g(i  ,j-5  , 5) ) ) * ny(i,j,2)
      
  
    predro1   =    - d5*w(i - 5,j,1) &
                   + d4*w(i - 4,j,1) &
                   - d3*w(i - 3,j,1) &
                   + d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) &
                   + d3*w(i + 2,j,1) &
                   - d4*w(i + 3,j,1) &
                   + d5*w(i + 4,j,1)
                   
    predrou1  =    - d5*w(i - 5,j,2) &
                   + d4*w(i - 4,j,2) &
                   - d3*w(i - 3,j,2) &
                   + d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2) &
                   + d3*w(i + 2,j,2) &
                   - d4*w(i + 3,j,2) &
                   + d5*w(i + 4,j,2)
                  
    predrov1  =    - d5*w(i - 5,j,3) &
                   + d4*w(i - 4,j,3) &
                   - d3*w(i - 3,j,3) &
                   + d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3) &
                   + d3*w(i + 2,j,3) &
                   - d4*w(i + 3,j,3) &
                   + d5*w(i + 4,j,3)
                   
    predrow1  =    - d5*w(i - 5,j,4) &
                   + d4*w(i - 4,j,4) &
                   - d3*w(i - 3,j,4) &
                   + d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4) &
                   + d3*w(i + 2,j,4) &
                   - d4*w(i + 3,j,4) &
                   + d5*w(i + 4,j,4)
                   
    predroe1  =    - d5*w(i - 5,j,5) &
                   + d4*w(i - 4,j,5) &
                   - d3*w(i - 3,j,5) &
                   + d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5) &
                   + d3*w(i + 2,j,5) &
                   - d4*w(i + 3,j,5) &
                   + d5*w(i + 4,j,5)
                                  
                   
    
    
                                                     
    predro2   =    - d5*w(i, j- 5,1) &
                   + d4*w(i, j- 4,1) &
                   - d3*w(i, j- 3,1) &
                   + d2*w(i, j- 2,1) &
                   - d1*w(i, j- 1,1) &
                   + d1*w(i, j   ,1) &
                   - d2*w(i, j+ 1,1) &
                   + d3*w(i, j+ 2,1) &
                   - d4*w(i, j+ 3,1) &
                   + d5*w(i, j+ 4,1)
                   
    predrou2  =    - d5*w(i, j- 5,2) &
                   + d4*w(i, j- 4,2) &
                   - d3*w(i, j- 3,2) &
                   + d2*w(i, j- 2,2) &
                   - d1*w(i, j- 1,2) &
                   + d1*w(i, j   ,2) &
                   - d2*w(i, j+ 1,2) &
                   + d3*w(i, j+ 2,2) &
                   - d4*w(i, j+ 3,2) &
                   + d5*w(i, j+ 4,2)
                  
    predrov2  =    - d5*w(i, j- 5,3) &
                   + d4*w(i, j- 4,3) &
                   - d3*w(i, j- 3,3) &
                   + d2*w(i, j- 2,3) &
                   - d1*w(i, j- 1,3) &
                   + d1*w(i, j   ,3) &
                   - d2*w(i, j+ 1,3) &
                   + d3*w(i, j+ 2,3) &
                   - d4*w(i, j+ 3,3) &
                   + d5*w(i, j+ 4,3)
                   
    predrow2  =    - d5*w(i, j- 5,4) &
                   + d4*w(i, j- 4,4) &
                   - d3*w(i, j- 3,4) &
                   + d2*w(i, j- 2,4) &
                   - d1*w(i, j- 1,4) &
                   + d1*w(i, j   ,4) &
                   - d2*w(i, j+ 1,4) &
                   + d3*w(i, j+ 2,4) &
                   - d4*w(i, j+ 3,4) &
                   + d5*w(i, j+ 4,4)
                   
    predroe2  =    - d5*w(i, j- 5,5) &
                   + d4*w(i, j- 4,5) &
                   - d3*w(i, j- 3,5) &
                   + d2*w(i, j- 2,5) &
                   - d1*w(i, j- 1,5) &
                   + d1*w(i, j   ,5) &
                   - d2*w(i, j+ 1,5) &
                   + d3*w(i, j+ 2,5) &
                   - d4*w(i, j+ 3,5) &
                   + d5*w(i, j+ 4,5)
    
                 
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
  
! Manage Wall Fluxes
   
! off centered schemes for fluxes near wall
! coef for near wall off-centered fluxes
!  coef are obtained with sympy for FD scheme
!face 9/2:
!{a8: -3/56, a1: -1/63, a5: 6/5, a7: 4/21, a2: 3/28, a9: 1/105, a10: -1/1260, a6: -1/2, a3: -4/7, a0: 1/840, a4: -11/30}
  a0  =   1.d0/840.d0
  a1  = - 1.d0/63.d0
  a2  =   3.d0/28.d0
  a3  = - 4.d0/7.d0
  a4  = -11.d0/30.d0
  a5  =   6.d0/5.d0
  a6  = - 0.5d0
  a7  =   4.d0/21.d0
  a8  = - 3.d0/56.d0
  a9  =   1.d0/105.d0
  a10 = - 1.d0/1260.d0
! Phi(9/2) = Phi(11/2) - Sum( ai*fi)
  c9d0  = c5 - a0
  c9d1  = c4 - a1
  c9d2  = c3 - a2
  c9d3  = c2 - a3
  c9d4  = c1 - a4
  c9d5  = c1 - a5
  c9d6  = c2 - a6
  c9d7  = c3 - a7
  c9d8  = c4 - a8
  c9d9  = c5 - a9
  c9d10 = -a10
! face 7/2
! {a8: 3/40, a1: 1/24, a5: -21/20, a7: -1/4, a2: -3/8, a9: -1/72, a10: 1/840, a6: 7/12, a3: -319/420, a0: -1/360, a4: 7/4}
! FD37
  a0  = -  1.d0/360.d0
  a1  =    1.d0/24.d0
  a2  = -  3.d0/8.d0
  a3  = -319.d0/420.d0
  a4  =    1.75d0
  a5  = - 21.d0/20.d0
  a6  =    7.d0/12.d0
  a7  = -  0.25d0
  a8  =    0.075d0
  a9  = -  1.d0/72.d0
  a10 =    1.d0/840.d0
! Phi(7/2) = Phi(9/2) - Sum( ai*fi)
  c7d0  = c9d0  - a0
  c7d1  = c9d1  - a1
  c7d2  = c9d2  - a2
  c7d3  = c9d3  - a3
  c7d4  = c9d4  - a4
  c7d5  = c9d5  - a5
  c7d6  = c9d6  - a6
  c7d7  = c9d7  - a7
  c7d8  = c9d8  - a8
  c7d9  = c9d9  - a9
  c7d10 = c9d10 - a10
! face 5/2
! {a8: -1/6, a1: -2/9, a5: 28/15, a7: 8/15, a2: -341/280, a9: 2/63, a10: -1/360, a6: -7/6, a3: 8/3, a0: 1/90, a4: -7/3}
! FD28
  a0=    1.d0/90.d0
  a1= -  2.d0/9.d0
  a2= -341.d0/280.d0
  a3=    8.d0/3.d0
  a4= -  7.d0/3.d0
  a5=   28.d0/15.d0
  a6= -  7.d0/6.d0
  a7=    8.d0/15.d0
  a8= -  1.d0/6.d0
  a9=    2.d0/63.d0
  a10= - 1.d0/360.d0
! Phi(5/2) = Phi(7/2) - Sum( ai*fi) given in Berland
  c5d0  = c7d0  - a0
  c5d1  = c7d1  - a1
  c5d2  = c7d2  - a2
  c5d3  = c7d3  - a3
  c5d4  = c7d4  - a4
  c5d5  = c7d5  - a5
  c5d6  = c7d6  - a6
  c5d7  = c7d7  - a7
  c5d8  = c7d8  - a8
  c5d9  = c7d9  - a9
  c5d10 = c7d10 - a10
! face 3/2
! {a8: 9/14, a1: -4609/2520, a5: -63/10, a7: -2, a2: 9/2, a9: -1/8, a10: 1/90, a6: 21/5, a3: -6, a0: -1/10, a4: 7}
! FD19
  a0= -   0.1d0
  a1= -4609.d0/2520.d0
  a2=     4.5d0
  a3= -   6.d0
  a4=     7.d0
  a5= -  6.3d0
  a6=    21.d0*0.2d0
  a7= -   2.d0
  a8=     9.d0/14.d0
  a9= -   0.125d0
  a10=    1.d0/90.d0
! Phi(3/2) = Phi(5/2) - Sum( ai*fi) given in Berland
  c3d0  = c5d0  - a0
  c3d1  = c5d1  - a1
  c3d2  = c5d2  - a2
  c3d3  = c5d3  - a3
  c3d4  = c5d4  - a4
  c3d5  = c5d5  - a5
  c3d6  = c5d6  - a6
  c3d7  = c5d7  - a7
  c3d8  = c5d8  - a8
  c3d9  = c5d9  - a9
  c3d10 = c5d10 - a10
  
  j = 5
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1,im+1
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) +               & 
                c3* (f(i+2,j  , 1) + f(i-3,j  , 1) ) +               &
                c4* (f(i+3,j  , 1) + f(i-4,j  , 1) ) +               &
                c5* (f(i+4,j  , 1) + f(i-5,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) +               & 
                c3* (g(i+2,j  , 1) + g(i-3,j  , 1) ) +               &
                c4* (g(i+3,j  , 1) + g(i-4,j  , 1) ) +               &
                c5* (g(i+4,j  , 1) + g(i-5,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) +               & 
                c3* (f(i+2,j  , 2) + f(i-3,j  , 2) ) +               &
                c4* (f(i+3,j  , 2) + f(i-4,j  , 2) ) +               &
                c5* (f(i+4,j  , 2) + f(i-5,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) +               & 
                c3* (g(i+2,j  , 2) + g(i-3,j  , 2) ) +               &
                c4* (g(i+3,j  , 2) + g(i-4,j  , 2) ) +               &
                c5* (g(i+4,j  , 2) + g(i-5,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) +               & 
                c3* (f(i+2,j  , 3) + f(i-3,j  , 3) ) +               &
                c4* (f(i+3,j  , 3) + f(i-4,j  , 3) ) +               &
                c5* (f(i+4,j  , 3) + f(i-5,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) +               & 
                c3* (g(i+2,j  , 3) + g(i-3,j  , 3) ) +               &
                c4* (g(i+3,j  , 3) + g(i-4,j  , 3) ) +               &
                c5* (g(i+4,j  , 3) + g(i-5,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) +               & 
                c3* (f(i+2,j  , 4) + f(i-3,j  , 4) ) +               &
                c4* (f(i+3,j  , 4) + f(i-4,j  , 4) ) +               &
                c5* (f(i+4,j  , 4) + f(i-5,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) +               & 
                c3* (g(i+2,j  , 4) + g(i-3,j  , 4) ) +               &
                c4* (g(i+3,j  , 4) + g(i-4,j  , 4) ) +               &
                c5* (g(i+4,j  , 4) + g(i-5,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) +               & 
                c3* (f(i+2,j  , 5) + f(i-3,j  , 5) ) +               &
                c4* (f(i+3,j  , 5) + f(i-4,j  , 5) ) +               &
                c5* (f(i+4,j  , 5) + f(i-5,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) +               & 
                c3* (g(i+2,j  , 5) + g(i-3,j  , 5) ) +               &
                c4* (g(i+3,j  , 5) + g(i-4,j  , 5) ) +               &
                c5* (g(i+4,j  , 5) + g(i-5,j  , 5) ) ) * ny(i,j,1)
      
  
  fxro2  = ( c9d0 * f(i, 1,1) + c9d1  * f(i, 2,1) + c9d2 * f(i,3,1) &
           + c9d3 * f(i, 4,1) + c9d4  * f(i, 5,1) + c9d5 * f(i,6,1) &
           + c9d6 * f(i, 7,1) + c9d7  * f(i, 8,1) + c9d8 * f(i,9,1) &
           + c9d9 * f(i,10,1) + c9d10 * f(i,11,1) ) * nx(i,j,2)   + &
           ( c9d0 * g(i, 1,1) + c9d1  * g(i, 2,1) + c9d2 * g(i,3,1) &
           + c9d3 * g(i, 4,1) + c9d4  * g(i, 5,1) + c9d5 * g(i,6,1) &
           + c9d6 * g(i, 7,1) + c9d7  * g(i, 8,1) + c9d8 * g(i,9,1) &
           + c9d9 * g(i,10,1) + c9d10 * g(i,11,1) ) * ny(i,j,2)
           
  fxrou2 = ( c9d0 * f(i, 1,2) + c9d1  * f(i, 2,2) + c9d2 * f(i,3,2) &
           + c9d3 * f(i, 4,2) + c9d4  * f(i, 5,2) + c9d5 * f(i,6,2) &
           + c9d6 * f(i, 7,2) + c9d7  * f(i, 8,2) + c9d8 * f(i,9,2) &
           + c9d9 * f(i,10,2) + c9d10 * f(i,11,2) ) * nx(i,j,2)   + &
           ( c9d0 * g(i, 1,2) + c9d1  * g(i, 2,2) + c9d2 * g(i,3,2) &
           + c9d3 * g(i, 4,2) + c9d4  * g(i, 5,2) + c9d5 * g(i,6,2) &
           + c9d6 * g(i, 7,2) + c9d7  * g(i, 8,2) + c9d8 * g(i,9,2) &
           + c9d9 * g(i,10,2) + c9d10 * g(i,11,2) ) * ny(i,j,2)
          
  fxrov2 = ( c9d0 * f(i, 1,3) + c9d1  * f(i, 2,3) + c9d2 * f(i,3,3) &
           + c9d3 * f(i, 4,3) + c9d4  * f(i, 5,3) + c9d5 * f(i,6,3) &
           + c9d6 * f(i, 7,3) + c9d7  * f(i, 8,3) + c9d8 * f(i,9,3) &
           + c9d9 * f(i,10,3) + c9d10 * f(i,11,3) ) * nx(i,j,2)   + &
           ( c9d0 * g(i, 1,3) + c9d1  * g(i, 2,3) + c9d2 * g(i,3,3) &
           + c9d3 * g(i, 4,3) + c9d4  * g(i, 5,3) + c9d5 * g(i,6,3) &
           + c9d6 * g(i, 7,3) + c9d7  * g(i, 8,3) + c9d8 * g(i,9,3) &
           + c9d9 * g(i,10,3) + c9d10 * g(i,11,3) ) * ny(i,j,2)
           
  fxrow2 = ( c9d0 * f(i, 1,4) + c9d1  * f(i, 2,4) + c9d2 * f(i,3,4) &
           + c9d3 * f(i, 4,4) + c9d4  * f(i, 5,4) + c9d5 * f(i,6,4) &
           + c9d6 * f(i, 7,4) + c9d7  * f(i, 8,4) + c9d8 * f(i,9,4) &
           + c9d9 * f(i,10,4) + c9d10 * f(i,11,4) ) * nx(i,j,2)   + &
           ( c9d0 * g(i, 1,4) + c9d1  * g(i, 2,4) + c9d2 * g(i,3,4) &
           + c9d3 * g(i, 4,4) + c9d4  * g(i, 5,4) + c9d5 * g(i,6,4) &
           + c9d6 * g(i, 7,4) + c9d7  * g(i, 8,4) + c9d8 * g(i,9,4) &
           + c9d9 * g(i,10,4) + c9d10 * g(i,11,4) ) * ny(i,j,2)
           
  fxroe2 = ( c9d0 * f(i, 1,5) + c9d1  * f(i, 2,5) + c9d2 * f(i,3,5) &
           + c9d3 * f(i, 4,5) + c9d4  * f(i, 5,5) + c9d5 * f(i,6,5) &
           + c9d6 * f(i, 7,5) + c9d7  * f(i, 8,5) + c9d8 * f(i,9,5) &
           + c9d9 * f(i,10,5) + c9d10 * f(i,11,5) ) * nx(i,j,2)   + &
           ( c9d0 * g(i, 1,5) + c9d1  * g(i, 2,5) + c9d2 * g(i,3,5) &
           + c9d3 * g(i, 4,5) + c9d4  * g(i, 5,5) + c9d5 * g(i,6,5) &
           + c9d6 * g(i, 7,5) + c9d7  * g(i, 8,5) + c9d8 * g(i,9,5) &
           + c9d9 * g(i,10,5) + c9d10 * g(i,11,5) ) * ny(i,j,2)
          
    predro1   =    - d5*w(i - 5,j,1) &
                   + d4*w(i - 4,j,1) &
                   - d3*w(i - 3,j,1) &
                   + d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) &
                   + d3*w(i + 2,j,1) &
                   - d4*w(i + 3,j,1) &
                   + d5*w(i + 4,j,1)
                   
    predrou1  =    - d5*w(i - 5,j,2) &
                   + d4*w(i - 4,j,2) &
                   - d3*w(i - 3,j,2) &
                   + d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2) &
                   + d3*w(i + 2,j,2) &
                   - d4*w(i + 3,j,2) &
                   + d5*w(i + 4,j,2)
                  
    predrov1  =    - d5*w(i - 5,j,3) &
                   + d4*w(i - 4,j,3) &
                   - d3*w(i - 3,j,3) &
                   + d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3) &
                   + d3*w(i + 2,j,3) &
                   - d4*w(i + 3,j,3) &
                   + d5*w(i + 4,j,3)
                   
    predrow1  =    - d5*w(i - 5,j,4) &
                   + d4*w(i - 4,j,4) &
                   - d3*w(i - 3,j,4) &
                   + d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4) &
                   + d3*w(i + 2,j,4) &
                   - d4*w(i + 3,j,4) &
                   + d5*w(i + 4,j,4)
                   
    predroe1  =    - d5*w(i - 5,j,5) &
                   + d4*w(i - 4,j,5) &
                   - d3*w(i - 3,j,5) &
                   + d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5) &
                   + d3*w(i + 2,j,5) &
                   - d4*w(i + 3,j,5) &
                   + d5*w(i + 4,j,5)
                                  
                   
    
    
                                                     
    predro2   =    - d5*w(i, j- 5,1) &
                   + d4*w(i, j- 4,1) &
                   - d3*w(i, j- 3,1) &
                   + d2*w(i, j- 2,1) &
                   - d1*w(i, j- 1,1) &
                   + d1*w(i, j   ,1) &
                   - d2*w(i, j+ 1,1) &
                   + d3*w(i, j+ 2,1) &
                   - d4*w(i, j+ 3,1) &
                   + d5*w(i, j+ 4,1)
                   
    predrou2  =    - d5*w(i, j- 5,2) &
                   + d4*w(i, j- 4,2) &
                   - d3*w(i, j- 3,2) &
                   + d2*w(i, j- 2,2) &
                   - d1*w(i, j- 1,2) &
                   + d1*w(i, j   ,2) &
                   - d2*w(i, j+ 1,2) &
                   + d3*w(i, j+ 2,2) &
                   - d4*w(i, j+ 3,2) &
                   + d5*w(i, j+ 4,2)
                  
    predrov2  =    - d5*w(i, j- 5,3) &
                   + d4*w(i, j- 4,3) &
                   - d3*w(i, j- 3,3) &
                   + d2*w(i, j- 2,3) &
                   - d1*w(i, j- 1,3) &
                   + d1*w(i, j   ,3) &
                   - d2*w(i, j+ 1,3) &
                   + d3*w(i, j+ 2,3) &
                   - d4*w(i, j+ 3,3) &
                   + d5*w(i, j+ 4,3)
                   
    predrow2  =    - d5*w(i, j- 5,4) &
                   + d4*w(i, j- 4,4) &
                   - d3*w(i, j- 3,4) &
                   + d2*w(i, j- 2,4) &
                   - d1*w(i, j- 1,4) &
                   + d1*w(i, j   ,4) &
                   - d2*w(i, j+ 1,4) &
                   + d3*w(i, j+ 2,4) &
                   - d4*w(i, j+ 3,4) &
                   + d5*w(i, j+ 4,4)
                   
    predroe2  =    - d5*w(i, j- 5,5) &
                   + d4*w(i, j- 4,5) &
                   - d3*w(i, j- 3,5) &
                   + d2*w(i, j- 2,5) &
                   - d1*w(i, j- 1,5) &
                   + d1*w(i, j   ,5) &
                   - d2*w(i, j+ 1,5) &
                   + d3*w(i, j+ 2,5) &
                   - d4*w(i, j+ 3,5) &
                   + d5*w(i, j+ 4,5)
    
                 
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

  enddo
  j=4                           
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1,im+1
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) +               & 
                c3* (f(i+2,j  , 1) + f(i-3,j  , 1) ) +               &
                c4* (f(i+3,j  , 1) + f(i-4,j  , 1) ) +               &
                c5* (f(i+4,j  , 1) + f(i-5,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) +               & 
                c3* (g(i+2,j  , 1) + g(i-3,j  , 1) ) +               &
                c4* (g(i+3,j  , 1) + g(i-4,j  , 1) ) +               &
                c5* (g(i+4,j  , 1) + g(i-5,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) +               & 
                c3* (f(i+2,j  , 2) + f(i-3,j  , 2) ) +               &
                c4* (f(i+3,j  , 2) + f(i-4,j  , 2) ) +               &
                c5* (f(i+4,j  , 2) + f(i-5,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) +               & 
                c3* (g(i+2,j  , 2) + g(i-3,j  , 2) ) +               &
                c4* (g(i+3,j  , 2) + g(i-4,j  , 2) ) +               &
                c5* (g(i+4,j  , 2) + g(i-5,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) +               & 
                c3* (f(i+2,j  , 3) + f(i-3,j  , 3) ) +               &
                c4* (f(i+3,j  , 3) + f(i-4,j  , 3) ) +               &
                c5* (f(i+4,j  , 3) + f(i-5,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) +               & 
                c3* (g(i+2,j  , 3) + g(i-3,j  , 3) ) +               &
                c4* (g(i+3,j  , 3) + g(i-4,j  , 3) ) +               &
                c5* (g(i+4,j  , 3) + g(i-5,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) +               & 
                c3* (f(i+2,j  , 4) + f(i-3,j  , 4) ) +               &
                c4* (f(i+3,j  , 4) + f(i-4,j  , 4) ) +               &
                c5* (f(i+4,j  , 4) + f(i-5,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) +               & 
                c3* (g(i+2,j  , 4) + g(i-3,j  , 4) ) +               &
                c4* (g(i+3,j  , 4) + g(i-4,j  , 4) ) +               &
                c5* (g(i+4,j  , 4) + g(i-5,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) +               & 
                c3* (f(i+2,j  , 5) + f(i-3,j  , 5) ) +               &
                c4* (f(i+3,j  , 5) + f(i-4,j  , 5) ) +               &
                c5* (f(i+4,j  , 5) + f(i-5,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) +               & 
                c3* (g(i+2,j  , 5) + g(i-3,j  , 5) ) +               &
                c4* (g(i+3,j  , 5) + g(i-4,j  , 5) ) +               &
                c5* (g(i+4,j  , 5) + g(i-5,j  , 5) ) ) * ny(i,j,1)
      
  
  fxro2  = ( c7d0 * f(i, 1,1) + c7d1  * f(i, 2,1) + c7d2 * f(i,3,1) &
           + c7d3 * f(i, 4,1) + c7d4  * f(i, 5,1) + c7d5 * f(i,6,1) &
           + c7d6 * f(i, 7,1) + c7d7  * f(i, 8,1) + c7d8 * f(i,9,1) &
           + c7d9 * f(i,10,1) + c7d10 * f(i,11,1) ) * nx(i,j,2)   + &
           ( c7d0 * g(i, 1,1) + c7d1  * g(i, 2,1) + c7d2 * g(i,3,1) &
           + c7d3 * g(i, 4,1) + c7d4  * g(i, 5,1) + c7d5 * g(i,6,1) &
           + c7d6 * g(i, 7,1) + c7d7  * g(i, 8,1) + c7d8 * g(i,9,1) &
           + c7d9 * g(i,10,1) + c7d10 * g(i,11,1) ) * ny(i,j,2)
           
  fxrou2 = ( c7d0 * f(i, 1,2) + c7d1  * f(i, 2,2) + c7d2 * f(i,3,2) &
           + c7d3 * f(i, 4,2) + c7d4  * f(i, 5,2) + c7d5 * f(i,6,2) &
           + c7d6 * f(i, 7,2) + c7d7  * f(i, 8,2) + c7d8 * f(i,9,2) &
           + c7d9 * f(i,10,2) + c7d10 * f(i,11,2) ) * nx(i,j,2)   + &
           ( c7d0 * g(i, 1,2) + c7d1  * g(i, 2,2) + c7d2 * g(i,3,2) &
           + c7d3 * g(i, 4,2) + c7d4  * g(i, 5,2) + c7d5 * g(i,6,2) &
           + c7d6 * g(i, 7,2) + c7d7  * g(i, 8,2) + c7d8 * g(i,9,2) &
           + c7d9 * g(i,10,2) + c7d10 * g(i,11,2) ) * ny(i,j,2)
          
  fxrov2 = ( c7d0 * f(i, 1,3) + c7d1  * f(i, 2,3) + c7d2 * f(i,3,3) &
           + c7d3 * f(i, 4,3) + c7d4  * f(i, 5,3) + c7d5 * f(i,6,3) &
           + c7d6 * f(i, 7,3) + c7d7  * f(i, 8,3) + c7d8 * f(i,9,3) &
           + c7d9 * f(i,10,3) + c7d10 * f(i,11,3) ) * nx(i,j,2)   + &
           ( c7d0 * g(i, 1,3) + c7d1  * g(i, 2,3) + c7d2 * g(i,3,3) &
           + c7d3 * g(i, 4,3) + c7d4  * g(i, 5,3) + c7d5 * g(i,6,3) &
           + c7d6 * g(i, 7,3) + c7d7  * g(i, 8,3) + c7d8 * g(i,9,3) &
           + c7d9 * g(i,10,3) + c7d10 * g(i,11,3) ) * ny(i,j,2)
           
  fxrow2 = ( c7d0 * f(i, 1,4) + c7d1  * f(i, 2,4) + c7d2 * f(i,3,4) &
           + c7d3 * f(i, 4,4) + c7d4  * f(i, 5,4) + c7d5 * f(i,6,4) &
           + c7d6 * f(i, 7,4) + c7d7  * f(i, 8,4) + c7d8 * f(i,9,4) &
           + c7d9 * f(i,10,4) + c7d10 * f(i,11,4) ) * nx(i,j,2)   + &
           ( c7d0 * g(i, 1,4) + c7d1  * g(i, 2,4) + c7d2 * g(i,3,4) &
           + c7d3 * g(i, 4,4) + c7d4  * g(i, 5,4) + c7d5 * g(i,6,4) &
           + c7d6 * g(i, 7,4) + c7d7  * g(i, 8,4) + c7d8 * g(i,9,4) &
           + c7d9 * g(i,10,4) + c7d10 * g(i,11,4) ) * ny(i,j,2)
           
  fxroe2 = ( c7d0 * f(i, 1,5) + c7d1  * f(i, 2,5) + c7d2 * f(i,3,5) &
           + c7d3 * f(i, 4,5) + c7d4  * f(i, 5,5) + c7d5 * f(i,6,5) &
           + c7d6 * f(i, 7,5) + c7d7  * f(i, 8,5) + c7d8 * f(i,9,5) &
           + c7d9 * f(i,10,5) + c7d10 * f(i,11,5) ) * nx(i,j,2)   + &
           ( c7d0 * g(i, 1,5) + c7d1  * g(i, 2,5) + c7d2 * g(i,3,5) &
           + c7d3 * g(i, 4,5) + c7d4  * g(i, 5,5) + c7d5 * g(i,6,5) &
           + c7d6 * g(i, 7,5) + c7d7  * g(i, 8,5) + c7d8 * g(i,9,5) &
           + c7d9 * g(i,10,5) + c7d10 * g(i,11,5) ) * ny(i,j,2)
  
  
    predro1   =    - d5*w(i - 5,j,1) &
                   + d4*w(i - 4,j,1) &
                   - d3*w(i - 3,j,1) &
                   + d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) &
                   + d3*w(i + 2,j,1) &
                   - d4*w(i + 3,j,1) &
                   + d5*w(i + 4,j,1)
                   
    predrou1  =    - d5*w(i - 5,j,2) &
                   + d4*w(i - 4,j,2) &
                   - d3*w(i - 3,j,2) &
                   + d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2) &
                   + d3*w(i + 2,j,2) &
                   - d4*w(i + 3,j,2) &
                   + d5*w(i + 4,j,2)
                  
    predrov1  =    - d5*w(i - 5,j,3) &
                   + d4*w(i - 4,j,3) &
                   - d3*w(i - 3,j,3) &
                   + d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3) &
                   + d3*w(i + 2,j,3) &
                   - d4*w(i + 3,j,3) &
                   + d5*w(i + 4,j,3)
                   
    predrow1  =    - d5*w(i - 5,j,4) &
                   + d4*w(i - 4,j,4) &
                   - d3*w(i - 3,j,4) &
                   + d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4) &
                   + d3*w(i + 2,j,4) &
                   - d4*w(i + 3,j,4) &
                   + d5*w(i + 4,j,4)
                   
    predroe1  =    - d5*w(i - 5,j,5) &
                   + d4*w(i - 4,j,5) &
                   - d3*w(i - 3,j,5) &
                   + d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5) &
                   + d3*w(i + 2,j,5) &
                   - d4*w(i + 3,j,5) &
                   + d5*w(i + 4,j,5)
                                  
                   
    
    
                                                     
    predro2   =    - d5*w(i, j- 5,1) &
                   + d4*w(i, j- 4,1) &
                   - d3*w(i, j- 3,1) &
                   + d2*w(i, j- 2,1) &
                   - d1*w(i, j- 1,1) &
                   + d1*w(i, j   ,1) &
                   - d2*w(i, j+ 1,1) &
                   + d3*w(i, j+ 2,1) &
                   - d4*w(i, j+ 3,1) &
                   + d5*w(i, j+ 4,1)
                   
    predrou2  =    - d5*w(i, j- 5,2) &
                   + d4*w(i, j- 4,2) &
                   - d3*w(i, j- 3,2) &
                   + d2*w(i, j- 2,2) &
                   - d1*w(i, j- 1,2) &
                   + d1*w(i, j   ,2) &
                   - d2*w(i, j+ 1,2) &
                   + d3*w(i, j+ 2,2) &
                   - d4*w(i, j+ 3,2) &
                   + d5*w(i, j+ 4,2)
                  
    predrov2  =    - d5*w(i, j- 5,3) &
                   + d4*w(i, j- 4,3) &
                   - d3*w(i, j- 3,3) &
                   + d2*w(i, j- 2,3) &
                   - d1*w(i, j- 1,3) &
                   + d1*w(i, j   ,3) &
                   - d2*w(i, j+ 1,3) &
                   + d3*w(i, j+ 2,3) &
                   - d4*w(i, j+ 3,3) &
                   + d5*w(i, j+ 4,3)
                   
    predrow2  =    - d5*w(i, j- 5,4) &
                   + d4*w(i, j- 4,4) &
                   - d3*w(i, j- 3,4) &
                   + d2*w(i, j- 2,4) &
                   - d1*w(i, j- 1,4) &
                   + d1*w(i, j   ,4) &
                   - d2*w(i, j+ 1,4) &
                   + d3*w(i, j+ 2,4) &
                   - d4*w(i, j+ 3,4) &
                   + d5*w(i, j+ 4,4)
                   
    predroe2  =    - d5*w(i, j- 5,5) &
                   + d4*w(i, j- 4,5) &
                   - d3*w(i, j- 3,5) &
                   + d2*w(i, j- 2,5) &
                   - d1*w(i, j- 1,5) &
                   + d1*w(i, j   ,5) &
                   - d2*w(i, j+ 1,5) &
                   + d3*w(i, j+ 2,5) &
                   - d4*w(i, j+ 3,5) &
                   + d5*w(i, j+ 4,5)
    
                 
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

  enddo
      
  j = 3                           
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1,im+1
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) +               & 
                c3* (f(i+2,j  , 1) + f(i-3,j  , 1) ) +               &
                c4* (f(i+3,j  , 1) + f(i-4,j  , 1) ) +               &
                c5* (f(i+4,j  , 1) + f(i-5,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) +               & 
                c3* (g(i+2,j  , 1) + g(i-3,j  , 1) ) +               &
                c4* (g(i+3,j  , 1) + g(i-4,j  , 1) ) +               &
                c5* (g(i+4,j  , 1) + g(i-5,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) +               & 
                c3* (f(i+2,j  , 2) + f(i-3,j  , 2) ) +               &
                c4* (f(i+3,j  , 2) + f(i-4,j  , 2) ) +               &
                c5* (f(i+4,j  , 2) + f(i-5,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) +               & 
                c3* (g(i+2,j  , 2) + g(i-3,j  , 2) ) +               &
                c4* (g(i+3,j  , 2) + g(i-4,j  , 2) ) +               &
                c5* (g(i+4,j  , 2) + g(i-5,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) +               & 
                c3* (f(i+2,j  , 3) + f(i-3,j  , 3) ) +               &
                c4* (f(i+3,j  , 3) + f(i-4,j  , 3) ) +               &
                c5* (f(i+4,j  , 3) + f(i-5,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) +               & 
                c3* (g(i+2,j  , 3) + g(i-3,j  , 3) ) +               &
                c4* (g(i+3,j  , 3) + g(i-4,j  , 3) ) +               &
                c5* (g(i+4,j  , 3) + g(i-5,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) +               & 
                c3* (f(i+2,j  , 4) + f(i-3,j  , 4) ) +               &
                c4* (f(i+3,j  , 4) + f(i-4,j  , 4) ) +               &
                c5* (f(i+4,j  , 4) + f(i-5,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) +               & 
                c3* (g(i+2,j  , 4) + g(i-3,j  , 4) ) +               &
                c4* (g(i+3,j  , 4) + g(i-4,j  , 4) ) +               &
                c5* (g(i+4,j  , 4) + g(i-5,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) +               & 
                c3* (f(i+2,j  , 5) + f(i-3,j  , 5) ) +               &
                c4* (f(i+3,j  , 5) + f(i-4,j  , 5) ) +               &
                c5* (f(i+4,j  , 5) + f(i-5,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) +               & 
                c3* (g(i+2,j  , 5) + g(i-3,j  , 5) ) +               &
                c4* (g(i+3,j  , 5) + g(i-4,j  , 5) ) +               &
                c5* (g(i+4,j  , 5) + g(i-5,j  , 5) ) ) * ny(i,j,1)
      
  
  fxro2  = ( c5d0 * f(i, 1,1) + c5d1  * f(i, 2,1) + c5d2 * f(i,3,1) &
           + c5d3 * f(i, 4,1) + c5d4  * f(i, 5,1) + c5d5 * f(i,6,1) &
           + c5d6 * f(i, 7,1) + c5d7  * f(i, 8,1) + c5d8 * f(i,9,1) &
           + c5d9 * f(i,10,1) + c5d10 * f(i,11,1) ) * nx(i,j,2)   + &
           ( c5d0 * g(i, 1,1) + c5d1  * g(i, 2,1) + c5d2 * g(i,3,1) &
           + c5d3 * g(i, 4,1) + c5d4  * g(i, 5,1) + c5d5 * g(i,6,1) &
           + c5d6 * g(i, 7,1) + c5d7  * g(i, 8,1) + c5d8 * g(i,9,1) &
           + c5d9 * g(i,10,1) + c5d10 * g(i,11,1) ) * ny(i,j,2)
           
  fxrou2 = ( c5d0 * f(i, 1,2) + c5d1  * f(i, 2,2) + c5d2 * f(i,3,2) &
           + c5d3 * f(i, 4,2) + c5d4  * f(i, 5,2) + c5d5 * f(i,6,2) &
           + c5d6 * f(i, 7,2) + c5d7  * f(i, 8,2) + c5d8 * f(i,9,2) &
           + c5d9 * f(i,10,2) + c5d10 * f(i,11,2) ) * nx(i,j,2)   + &
           ( c5d0 * g(i, 1,2) + c5d1  * g(i, 2,2) + c5d2 * g(i,3,2) &
           + c5d3 * g(i, 4,2) + c5d4  * g(i, 5,2) + c5d5 * g(i,6,2) &
           + c5d6 * g(i, 7,2) + c5d7  * g(i, 8,2) + c5d8 * g(i,9,2) &
           + c5d9 * g(i,10,2) + c5d10 * g(i,11,2) ) * ny(i,j,2)
          
  fxrov2 = ( c5d0 * f(i, 1,3) + c5d1  * f(i, 2,3) + c5d2 * f(i,3,3) &
           + c5d3 * f(i, 4,3) + c5d4  * f(i, 5,3) + c5d5 * f(i,6,3) &
           + c5d6 * f(i, 7,3) + c5d7  * f(i, 8,3) + c5d8 * f(i,9,3) &
           + c5d9 * f(i,10,3) + c5d10 * f(i,11,3) ) * nx(i,j,2)   + &
           ( c5d0 * g(i, 1,3) + c5d1  * g(i, 2,3) + c5d2 * g(i,3,3) &
           + c5d3 * g(i, 4,3) + c5d4  * g(i, 5,3) + c5d5 * g(i,6,3) &
           + c5d6 * g(i, 7,3) + c5d7  * g(i, 8,3) + c5d8 * g(i,9,3) &
           + c5d9 * g(i,10,3) + c5d10 * g(i,11,3) ) * ny(i,j,2)
           
  fxrow2 = ( c5d0 * f(i, 1,4) + c5d1  * f(i, 2,4) + c5d2 * f(i,3,4) &
           + c5d3 * f(i, 4,4) + c5d4  * f(i, 5,4) + c5d5 * f(i,6,4) &
           + c5d6 * f(i, 7,4) + c5d7  * f(i, 8,4) + c5d8 * f(i,9,4) &
           + c5d9 * f(i,10,4) + c5d10 * f(i,11,4) ) * nx(i,j,2)   + &
           ( c5d0 * g(i, 1,4) + c5d1  * g(i, 2,4) + c5d2 * g(i,3,4) &
           + c5d3 * g(i, 4,4) + c5d4  * g(i, 5,4) + c5d5 * g(i,6,4) &
           + c5d6 * g(i, 7,4) + c5d7  * g(i, 8,4) + c5d8 * g(i,9,4) &
           + c5d9 * g(i,10,4) + c5d10 * g(i,11,4) ) * ny(i,j,2)
           
  fxroe2 = ( c5d0 * f(i, 1,5) + c5d1  * f(i, 2,5) + c5d2 * f(i,3,5) &
           + c5d3 * f(i, 4,5) + c5d4  * f(i, 5,5) + c5d5 * f(i,6,5) &
           + c5d6 * f(i, 7,5) + c5d7  * f(i, 8,5) + c5d8 * f(i,9,5) &
           + c5d9 * f(i,10,5) + c5d10 * f(i,11,5) ) * nx(i,j,2)   + &
           ( c5d0 * g(i, 1,5) + c5d1  * g(i, 2,5) + c5d2 * g(i,3,5) &
           + c5d3 * g(i, 4,5) + c5d4  * g(i, 5,5) + c5d5 * g(i,6,5) &
           + c5d6 * g(i, 7,5) + c5d7  * g(i, 8,5) + c5d8 * g(i,9,5) &
           + c5d9 * g(i,10,5) + c5d10 * g(i,11,5) ) * ny(i,j,2)
    predro1   =    - d5*w(i - 5,j,1) &
                   + d4*w(i - 4,j,1) &
                   - d3*w(i - 3,j,1) &
                   + d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) &
                   + d3*w(i + 2,j,1) &
                   - d4*w(i + 3,j,1) &
                   + d5*w(i + 4,j,1)
                   
    predrou1  =    - d5*w(i - 5,j,2) &
                   + d4*w(i - 4,j,2) &
                   - d3*w(i - 3,j,2) &
                   + d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2) &
                   + d3*w(i + 2,j,2) &
                   - d4*w(i + 3,j,2) &
                   + d5*w(i + 4,j,2)
                  
    predrov1  =    - d5*w(i - 5,j,3) &
                   + d4*w(i - 4,j,3) &
                   - d3*w(i - 3,j,3) &
                   + d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3) &
                   + d3*w(i + 2,j,3) &
                   - d4*w(i + 3,j,3) &
                   + d5*w(i + 4,j,3)
                   
    predrow1  =    - d5*w(i - 5,j,4) &
                   + d4*w(i - 4,j,4) &
                   - d3*w(i - 3,j,4) &
                   + d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4) &
                   + d3*w(i + 2,j,4) &
                   - d4*w(i + 3,j,4) &
                   + d5*w(i + 4,j,4)
                   
    predroe1  =    - d5*w(i - 5,j,5) &
                   + d4*w(i - 4,j,5) &
                   - d3*w(i - 3,j,5) &
                   + d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5) &
                   + d3*w(i + 2,j,5) &
                   - d4*w(i + 3,j,5) &
                   + d5*w(i + 4,j,5)
                                  
                   
    
    
                                                     
    predro2   =    - d5*w(i, j- 5,1) &
                   + d4*w(i, j- 4,1) &
                   - d3*w(i, j- 3,1) &
                   + d2*w(i, j- 2,1) &
                   - d1*w(i, j- 1,1) &
                   + d1*w(i, j   ,1) &
                   - d2*w(i, j+ 1,1) &
                   + d3*w(i, j+ 2,1) &
                   - d4*w(i, j+ 3,1) &
                   + d5*w(i, j+ 4,1)
                   
    predrou2  =    - d5*w(i, j- 5,2) &
                   + d4*w(i, j- 4,2) &
                   - d3*w(i, j- 3,2) &
                   + d2*w(i, j- 2,2) &
                   - d1*w(i, j- 1,2) &
                   + d1*w(i, j   ,2) &
                   - d2*w(i, j+ 1,2) &
                   + d3*w(i, j+ 2,2) &
                   - d4*w(i, j+ 3,2) &
                   + d5*w(i, j+ 4,2)
                  
    predrov2  =    - d5*w(i, j- 5,3) &
                   + d4*w(i, j- 4,3) &
                   - d3*w(i, j- 3,3) &
                   + d2*w(i, j- 2,3) &
                   - d1*w(i, j- 1,3) &
                   + d1*w(i, j   ,3) &
                   - d2*w(i, j+ 1,3) &
                   + d3*w(i, j+ 2,3) &
                   - d4*w(i, j+ 3,3) &
                   + d5*w(i, j+ 4,3)
                   
    predrow2  =    - d5*w(i, j- 5,4) &
                   + d4*w(i, j- 4,4) &
                   - d3*w(i, j- 3,4) &
                   + d2*w(i, j- 2,4) &
                   - d1*w(i, j- 1,4) &
                   + d1*w(i, j   ,4) &
                   - d2*w(i, j+ 1,4) &
                   + d3*w(i, j+ 2,4) &
                   - d4*w(i, j+ 3,4) &
                   + d5*w(i, j+ 4,4)
                   
    predroe2  =    - d5*w(i, j- 5,5) &
                   + d4*w(i, j- 4,5) &
                   - d3*w(i, j- 3,5) &
                   + d2*w(i, j- 2,5) &
                   - d1*w(i, j- 1,5) &
                   + d1*w(i, j   ,5) &
                   - d2*w(i, j+ 1,5) &
                   + d3*w(i, j+ 2,5) &
                   - d4*w(i, j+ 3,5) &
                   + d5*w(i, j+ 4,5)
    
                 
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


  enddo
     
  j = 2                          
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1,im+1
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) +               & 
                c3* (f(i+2,j  , 1) + f(i-3,j  , 1) ) +               &
                c4* (f(i+3,j  , 1) + f(i-4,j  , 1) ) +               &
                c5* (f(i+4,j  , 1) + f(i-5,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) +               & 
                c3* (g(i+2,j  , 1) + g(i-3,j  , 1) ) +               &
                c4* (g(i+3,j  , 1) + g(i-4,j  , 1) ) +               &
                c5* (g(i+4,j  , 1) + g(i-5,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) +               & 
                c3* (f(i+2,j  , 2) + f(i-3,j  , 2) ) +               &
                c4* (f(i+3,j  , 2) + f(i-4,j  , 2) ) +               &
                c5* (f(i+4,j  , 2) + f(i-5,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) +               & 
                c3* (g(i+2,j  , 2) + g(i-3,j  , 2) ) +               &
                c4* (g(i+3,j  , 2) + g(i-4,j  , 2) ) +               &
                c5* (g(i+4,j  , 2) + g(i-5,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) +               & 
                c3* (f(i+2,j  , 3) + f(i-3,j  , 3) ) +               &
                c4* (f(i+3,j  , 3) + f(i-4,j  , 3) ) +               &
                c5* (f(i+4,j  , 3) + f(i-5,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) +               & 
                c3* (g(i+2,j  , 3) + g(i-3,j  , 3) ) +               &
                c4* (g(i+3,j  , 3) + g(i-4,j  , 3) ) +               &
                c5* (g(i+4,j  , 3) + g(i-5,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) +               & 
                c3* (f(i+2,j  , 4) + f(i-3,j  , 4) ) +               &
                c4* (f(i+3,j  , 4) + f(i-4,j  , 4) ) +               &
                c5* (f(i+4,j  , 4) + f(i-5,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) +               & 
                c3* (g(i+2,j  , 4) + g(i-3,j  , 4) ) +               &
                c4* (g(i+3,j  , 4) + g(i-4,j  , 4) ) +               &
                c5* (g(i+4,j  , 4) + g(i-5,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) +               & 
                c3* (f(i+2,j  , 5) + f(i-3,j  , 5) ) +               &
                c4* (f(i+3,j  , 5) + f(i-4,j  , 5) ) +               &
                c5* (f(i+4,j  , 5) + f(i-5,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) +               & 
                c3* (g(i+2,j  , 5) + g(i-3,j  , 5) ) +               &
                c4* (g(i+3,j  , 5) + g(i-4,j  , 5) ) +               &
                c5* (g(i+4,j  , 5) + g(i-5,j  , 5) ) ) * ny(i,j,1)
      
  
  fxro2  = ( c3d0 * f(i, 1,1) + c3d1  * f(i, 2,1) + c3d2 * f(i,3,1) &
           + c3d3 * f(i, 4,1) + c3d4  * f(i, 5,1) + c3d5 * f(i,6,1) &
           + c3d6 * f(i, 7,1) + c3d7  * f(i, 8,1) + c3d8 * f(i,9,1) &
           + c3d9 * f(i,10,1) + c3d10 * f(i,11,1) ) * nx(i,j,2)   + &
           ( c3d0 * g(i, 1,1) + c3d1  * g(i, 2,1) + c3d2 * g(i,3,1) &
           + c3d3 * g(i, 4,1) + c3d4  * g(i, 5,1) + c3d5 * g(i,6,1) &
           + c3d6 * g(i, 7,1) + c3d7  * g(i, 8,1) + c3d8 * g(i,9,1) &
           + c3d9 * g(i,10,1) + c3d10 * g(i,11,1) ) * ny(i,j,2)
           
  fxrou2 = ( c3d0 * f(i, 1,2) + c3d1  * f(i, 2,2) + c3d2 * f(i,3,2) &
           + c3d3 * f(i, 4,2) + c3d4  * f(i, 5,2) + c3d5 * f(i,6,2) &
           + c3d6 * f(i, 7,2) + c3d7  * f(i, 8,2) + c3d8 * f(i,9,2) &
           + c3d9 * f(i,10,2) + c3d10 * f(i,11,2) ) * nx(i,j,2)   + &
           ( c3d0 * g(i, 1,2) + c3d1  * g(i, 2,2) + c3d2 * g(i,3,2) &
           + c3d3 * g(i, 4,2) + c3d4  * g(i, 5,2) + c3d5 * g(i,6,2) &
           + c3d6 * g(i, 7,2) + c3d7  * g(i, 8,2) + c3d8 * g(i,9,2) &
           + c3d9 * g(i,10,2) + c3d10 * g(i,11,2) ) * ny(i,j,2)
          
  fxrov2 = ( c3d0 * f(i, 1,3) + c3d1  * f(i, 2,3) + c3d2 * f(i,3,3) &
           + c3d3 * f(i, 4,3) + c3d4  * f(i, 5,3) + c3d5 * f(i,6,3) &
           + c3d6 * f(i, 7,3) + c3d7  * f(i, 8,3) + c3d8 * f(i,9,3) &
           + c3d9 * f(i,10,3) + c3d10 * f(i,11,3) ) * nx(i,j,2)   + &
           ( c3d0 * g(i, 1,3) + c3d1  * g(i, 2,3) + c3d2 * g(i,3,3) &
           + c3d3 * g(i, 4,3) + c3d4  * g(i, 5,3) + c3d5 * g(i,6,3) &
           + c3d6 * g(i, 7,3) + c3d7  * g(i, 8,3) + c3d8 * g(i,9,3) &
           + c3d9 * g(i,10,3) + c3d10 * g(i,11,3) ) * ny(i,j,2)
           
  fxrow2 = ( c3d0 * f(i, 1,4) + c3d1  * f(i, 2,4) + c3d2 * f(i,3,4) &
           + c3d3 * f(i, 4,4) + c3d4  * f(i, 5,4) + c3d5 * f(i,6,4) &
           + c3d6 * f(i, 7,4) + c3d7  * f(i, 8,4) + c3d8 * f(i,9,4) &
           + c3d9 * f(i,10,4) + c3d10 * f(i,11,4) ) * nx(i,j,2)   + &
           ( c3d0 * g(i, 1,4) + c3d1  * g(i, 2,4) + c3d2 * g(i,3,4) &
           + c3d3 * g(i, 4,4) + c3d4  * g(i, 5,4) + c3d5 * g(i,6,4) &
           + c3d6 * g(i, 7,4) + c3d7  * g(i, 8,4) + c3d8 * g(i,9,4) &
           + c3d9 * g(i,10,4) + c3d10 * g(i,11,4) ) * ny(i,j,2)
           
  fxroe2 = ( c3d0 * f(i, 1,5) + c3d1  * f(i, 2,5) + c3d2 * f(i,3,5) &
           + c3d3 * f(i, 4,5) + c3d4  * f(i, 5,5) + c3d5 * f(i,6,5) &
           + c3d6 * f(i, 7,5) + c3d7  * f(i, 8,5) + c3d8 * f(i,9,5) &
           + c3d9 * f(i,10,5) + c3d10 * f(i,11,5) ) * nx(i,j,2)   + &
           ( c3d0 * g(i, 1,5) + c3d1  * g(i, 2,5) + c3d2 * g(i,3,5) &
           + c3d3 * g(i, 4,5) + c3d4  * g(i, 5,5) + c3d5 * g(i,6,5) &
           + c3d6 * g(i, 7,5) + c3d7  * g(i, 8,5) + c3d8 * g(i,9,5) &
           + c3d9 * g(i,10,5) + c3d10 * g(i,11,5) ) * ny(i,j,2)
          
    predro1   =    - d5*w(i - 5,j,1) &
                   + d4*w(i - 4,j,1) &
                   - d3*w(i - 3,j,1) &
                   + d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) &
                   + d3*w(i + 2,j,1) &
                   - d4*w(i + 3,j,1) &
                   + d5*w(i + 4,j,1)
                   
    predrou1  =    - d5*w(i - 5,j,2) &
                   + d4*w(i - 4,j,2) &
                   - d3*w(i - 3,j,2) &
                   + d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2) &
                   + d3*w(i + 2,j,2) &
                   - d4*w(i + 3,j,2) &
                   + d5*w(i + 4,j,2)
                  
    predrov1  =    - d5*w(i - 5,j,3) &
                   + d4*w(i - 4,j,3) &
                   - d3*w(i - 3,j,3) &
                   + d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3) &
                   + d3*w(i + 2,j,3) &
                   - d4*w(i + 3,j,3) &
                   + d5*w(i + 4,j,3)
                   
    predrow1  =    - d5*w(i - 5,j,4) &
                   + d4*w(i - 4,j,4) &
                   - d3*w(i - 3,j,4) &
                   + d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4) &
                   + d3*w(i + 2,j,4) &
                   - d4*w(i + 3,j,4) &
                   + d5*w(i + 4,j,4)
                   
    predroe1  =    - d5*w(i - 5,j,5) &
                   + d4*w(i - 4,j,5) &
                   - d3*w(i - 3,j,5) &
                   + d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5) &
                   + d3*w(i + 2,j,5) &
                   - d4*w(i + 3,j,5) &
                   + d5*w(i + 4,j,5)
                                  
                   
    
    
                                                     
    predro2   =    - d5*w(i, j- 5,1) &
                   + d4*w(i, j- 4,1) &
                   - d3*w(i, j- 3,1) &
                   + d2*w(i, j- 2,1) &
                   - d1*w(i, j- 1,1) &
                   + d1*w(i, j   ,1) &
                   - d2*w(i, j+ 1,1) &
                   + d3*w(i, j+ 2,1) &
                   - d4*w(i, j+ 3,1) &
                   + d5*w(i, j+ 4,1)
                   
    predrou2  =    - d5*w(i, j- 5,2) &
                   + d4*w(i, j- 4,2) &
                   - d3*w(i, j- 3,2) &
                   + d2*w(i, j- 2,2) &
                   - d1*w(i, j- 1,2) &
                   + d1*w(i, j   ,2) &
                   - d2*w(i, j+ 1,2) &
                   + d3*w(i, j+ 2,2) &
                   - d4*w(i, j+ 3,2) &
                   + d5*w(i, j+ 4,2)
                  
    predrov2  =    - d5*w(i, j- 5,3) &
                   + d4*w(i, j- 4,3) &
                   - d3*w(i, j- 3,3) &
                   + d2*w(i, j- 2,3) &
                   - d1*w(i, j- 1,3) &
                   + d1*w(i, j   ,3) &
                   - d2*w(i, j+ 1,3) &
                   + d3*w(i, j+ 2,3) &
                   - d4*w(i, j+ 3,3) &
                   + d5*w(i, j+ 4,3)
                   
    predrow2  =    - d5*w(i, j- 5,4) &
                   + d4*w(i, j- 4,4) &
                   - d3*w(i, j- 3,4) &
                   + d2*w(i, j- 2,4) &
                   - d1*w(i, j- 1,4) &
                   + d1*w(i, j   ,4) &
                   - d2*w(i, j+ 1,4) &
                   + d3*w(i, j+ 2,4) &
                   - d4*w(i, j+ 3,4) &
                   + d5*w(i, j+ 4,4)
                   
    predroe2  =    - d5*w(i, j- 5,5) &
                   + d4*w(i, j- 4,5) &
                   - d3*w(i, j- 3,5) &
                   + d2*w(i, j- 2,5) &
                   - d1*w(i, j- 1,5) &
                   + d1*w(i, j   ,5) &
                   - d2*w(i, j+ 1,5) &
                   + d3*w(i, j+ 2,5) &
                   - d4*w(i, j+ 3,5) &
                   + d5*w(i, j+ 4,5)
    
                 
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

                                                           
  enddo
! force dp/dn = 0. at 2nd oder
  ct0 = 1.125d0! 9/8
  ct1 =-0.125d0!-1/8
  
! force dp/dn = 0. at 3rd oder
!ct0 = 225.d0/184.d0
!ct1 = -25.d0/92.d0
!ct2 = 9.d0/184.d0
  j = 1
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1,im+1
! for idir (only works for infinite Wall)
  fxro1    =  ( c1* (f(i  ,j  , 1) + f(i-1,j  , 1) ) +               &
                c2* (f(i+1,j  , 1) + f(i-2,j  , 1) ) +               & 
                c3* (f(i+2,j  , 1) + f(i-3,j  , 1) ) +               &
                c4* (f(i+3,j  , 1) + f(i-4,j  , 1) ) +               &
                c5* (f(i+4,j  , 1) + f(i-5,j  , 1) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 1) + g(i-1,j  , 1) ) +               &
                c2* (g(i+1,j  , 1) + g(i-2,j  , 1) ) +               & 
                c3* (g(i+2,j  , 1) + g(i-3,j  , 1) ) +               &
                c4* (g(i+3,j  , 1) + g(i-4,j  , 1) ) +               &
                c5* (g(i+4,j  , 1) + g(i-5,j  , 1) ) ) * ny(i,j,1)
      
  fxrou1   =  ( c1* (f(i  ,j  , 2) + f(i-1,j  , 2) ) +               &
                c2* (f(i+1,j  , 2) + f(i-2,j  , 2) ) +               & 
                c3* (f(i+2,j  , 2) + f(i-3,j  , 2) ) +               &
                c4* (f(i+3,j  , 2) + f(i-4,j  , 2) ) +               &
                c5* (f(i+4,j  , 2) + f(i-5,j  , 2) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 2) + g(i-1,j  , 2) ) +               &
                c2* (g(i+1,j  , 2) + g(i-2,j  , 2) ) +               & 
                c3* (g(i+2,j  , 2) + g(i-3,j  , 2) ) +               &
                c4* (g(i+3,j  , 2) + g(i-4,j  , 2) ) +               &
                c5* (g(i+4,j  , 2) + g(i-5,j  , 2) ) ) * ny(i,j,1)
      
  fxrov1   =  ( c1* (f(i  ,j  , 3) + f(i-1,j  , 3) ) +               &
                c2* (f(i+1,j  , 3) + f(i-2,j  , 3) ) +               & 
                c3* (f(i+2,j  , 3) + f(i-3,j  , 3) ) +               &
                c4* (f(i+3,j  , 3) + f(i-4,j  , 3) ) +               &
                c5* (f(i+4,j  , 3) + f(i-5,j  , 3) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 3) + g(i-1,j  , 3) ) +               &
                c2* (g(i+1,j  , 3) + g(i-2,j  , 3) ) +               & 
                c3* (g(i+2,j  , 3) + g(i-3,j  , 3) ) +               &
                c4* (g(i+3,j  , 3) + g(i-4,j  , 3) ) +               &
                c5* (g(i+4,j  , 3) + g(i-5,j  , 3) ) ) * ny(i,j,1)
      
  fxrow1   =  ( c1* (f(i  ,j  , 4) + f(i-1,j  , 4) ) +               &
                c2* (f(i+1,j  , 4) + f(i-2,j  , 4) ) +               & 
                c3* (f(i+2,j  , 4) + f(i-3,j  , 4) ) +               &
                c4* (f(i+3,j  , 4) + f(i-4,j  , 4) ) +               &
                c5* (f(i+4,j  , 4) + f(i-5,j  , 4) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 4) + g(i-1,j  , 4) ) +               &
                c2* (g(i+1,j  , 4) + g(i-2,j  , 4) ) +               & 
                c3* (g(i+2,j  , 4) + g(i-3,j  , 4) ) +               &
                c4* (g(i+3,j  , 4) + g(i-4,j  , 4) ) +               &
                c5* (g(i+4,j  , 4) + g(i-5,j  , 4) ) ) * ny(i,j,1)
      
  fxroe1   =  ( c1* (f(i  ,j  , 5) + f(i-1,j  , 5) ) +               &
                c2* (f(i+1,j  , 5) + f(i-2,j  , 5) ) +               & 
                c3* (f(i+2,j  , 5) + f(i-3,j  , 5) ) +               &
                c4* (f(i+3,j  , 5) + f(i-4,j  , 5) ) +               &
                c5* (f(i+4,j  , 5) + f(i-5,j  , 5) ) ) * nx(i,j,1) + &
              ( c1* (g(i  ,j  , 5) + g(i-1,j  , 5) ) +               &
                c2* (g(i+1,j  , 5) + g(i-2,j  , 5) ) +               & 
                c3* (g(i+2,j  , 5) + g(i-3,j  , 5) ) +               &
                c4* (g(i+3,j  , 5) + g(i-4,j  , 5) ) +               &
                c5* (g(i+4,j  , 5) + g(i-5,j  , 5) ) ) * ny(i,j,1)
      
  
    predro1   =    - d5*w(i - 5,j,1) &
                   + d4*w(i - 4,j,1) &
                   - d3*w(i - 3,j,1) &
                   + d2*w(i - 2,j,1) &
                   - d1*w(i - 1,j,1) &
                   + d1*w(i    ,j,1) &
                   - d2*w(i + 1,j,1) &
                   + d3*w(i + 2,j,1) &
                   - d4*w(i + 3,j,1) &
                   + d5*w(i + 4,j,1)
                   
    predrou1  =    - d5*w(i - 5,j,2) &
                   + d4*w(i - 4,j,2) &
                   - d3*w(i - 3,j,2) &
                   + d2*w(i - 2,j,2) &
                   - d1*w(i - 1,j,2) &
                   + d1*w(i    ,j,2) &
                   - d2*w(i + 1,j,2) &
                   + d3*w(i + 2,j,2) &
                   - d4*w(i + 3,j,2) &
                   + d5*w(i + 4,j,2)
                  
    predrov1  =    - d5*w(i - 5,j,3) &
                   + d4*w(i - 4,j,3) &
                   - d3*w(i - 3,j,3) &
                   + d2*w(i - 2,j,3) &
                   - d1*w(i - 1,j,3) &
                   + d1*w(i    ,j,3) &
                   - d2*w(i + 1,j,3) &
                   + d3*w(i + 2,j,3) &
                   - d4*w(i + 3,j,3) &
                   + d5*w(i + 4,j,3)
                   
    predrow1  =    - d5*w(i - 5,j,4) &
                   + d4*w(i - 4,j,4) &
                   - d3*w(i - 3,j,4) &
                   + d2*w(i - 2,j,4) &
                   - d1*w(i - 1,j,4) &
                   + d1*w(i    ,j,4) &
                   - d2*w(i + 1,j,4) &
                   + d3*w(i + 2,j,4) &
                   - d4*w(i + 3,j,4) &
                   + d5*w(i + 4,j,4)
                   
    predroe1  =    - d5*w(i - 5,j,5) &
                   + d4*w(i - 4,j,5) &
                   - d3*w(i - 3,j,5) &
                   + d2*w(i - 2,j,5) &
                   - d1*w(i - 1,j,5) &
                   + d1*w(i    ,j,5) &
                   - d2*w(i + 1,j,5) &
                   + d3*w(i + 2,j,5) &
                   - d4*w(i + 3,j,5) &
                   + d5*w(i + 4,j,5)
                                  
                   
    
    
                 
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
! ux = ZERO
! vx = ZERO
! wx = ZERO
!
! set viscous fluxes at wall----------------------------
! dir normal to the wall
    fvrou2 = TWOTHIRD*mmu*( TWO*ux -     vy        )    
    fvrov2 =          mmu*(     uy +     vx        )    
    fvrow2 =          mmu*(                    wx  )    
    fvroe2 = ZERO ! adiab
!
    gvrou2 =          mmu*(     uy +     vx       )
    gvrov2 = TWOTHIRD*mmu*(    -ux + TWO*vy       )
    gvrow2 =          mmu*(                  wy   )
!
    fvroe2 = ZERO !adiab
!
    hn(i,j,1,2) =   ZERO
    hn(i,j,2,2) = fxrov2*nx(i,j,2) - (fvrou2 * nx(i,j,2) + gvrou2 * ny(i,j,2))
    hn(i,j,3,2) = fxrov2*ny(i,j,2) - (fvrov2 * nx(i,j,2) + gvrov2 * ny(i,j,2))
    hn(i,j,4,2) =                  - (fvrow2 * nx(i,j,2) + gvrow2 * ny(i,j,2))
    hn(i,j,5,2) =   ZERO


        
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
      
  

end subroutine flux_num_dnc9_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
