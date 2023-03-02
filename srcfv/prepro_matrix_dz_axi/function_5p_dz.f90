! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

! =============================================================================
!          consistent fluxes for DNC3 2D
! =============================================================================
!
subroutine function_5p_dz(func0,func1,func2,func3,func4,func5,func6,func7,func8,&
                            func9,func10,func11,func12,func13,func14,func15,&
                            w,ym,x0,y0,nx,ny,xc,yc,vol,volf,gh,cp,cv,prandtl,&
                            gam,rgaz,cs,muref,tref,s_suth,im,jm)
!
  implicit none
! variables for dimension -----------------------------------------
  integer :: im,jm,gh
! required arguments ----------------------------------------------
  real(8),intent(in) :: cp,cv,prandtl,gam,rgaz ! thermo
  real(8),intent(in) :: cs,muref,tref,s_suth ! viscosity
  real(8),intent(in) :: ym ! axis position
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
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func0
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func1
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func2
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func3
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func4
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func5
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func6
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func7
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func8
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func9
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func10
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func11
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func12
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func13
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func14
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: func15
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
  real(8) :: divu2,vort2,dx,dy,dxm1,dym1,dxm2,dym2,distm1
  real(8) :: gui,gvi,gwi,gmui
  real(8) :: guj,gvj,gwj,gmuj
  real(8) :: rhom,rhomr,rhoml,rhom1l,c2l,c2r,rr,r,u,ur,ul,vr,vl,wr,wl,c2x,nx2,ny2
  real(8) :: ab, sq,ducros1,ducros2,k_sensor1,k_sensor2
  real(8) :: b1,b2,b3,b4,b5,c1,c2,c3,c4,c5,wiggle, denom,betas
  real(8) :: nxloc, nyloc, sn, invsn, sc1, sc2
  real(8) :: d1, d2
  real(8) :: c3d0,c3d1
  real(8) :: pw,ct0,ct1,ct2,cvm1
  real(8) :: coef,omrr,test,diffro,diffrou,diffrov,diffrow,diffroe,v
  real(8) :: HALF,ONE,ZERO,TWO,TWOTHIRD,FOURTH,TWELFTH 
  real(8) :: TWENTYFOURTH,ccross
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
   
! Coef for grad o4
  denom = 1.d0/12.d0
  
  b1 =  8.d0 * denom
  b2 = -       denom 
  
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
!DIR$ IVDEP
!$AD II-LOOP
  do i = 1, im
!idir
       gui   = ( b1*(velx(i+1,j     ) - velx(i-1,j     )) &
               + b2*(velx(i+2,j     ) - velx(i-2,j     )) )
       
       gvi   = ( b1*(vely(i+1,j     ) - vely(i-1,j     )) &
               + b2*(vely(i+2,j     ) - vely(i-2,j     )) )
       
       gwi   = ( b1*(velz(i+1,j     ) - velz(i-1,j     )) &
               + b2*(velz(i+2,j     ) - velz(i-2,j     )) )
                               
       gmui  = ( b1*(  mu(i+1,j     ) -   mu(i-1,j     )) &
               + b2*(  mu(i+2,j     ) -   mu(i-2,j     )) )
       
                           
       
!jdir
       guj   = ( b1*(velx(i,j+1   ) - velx(i,j-1   )) &
               + b2*(velx(i,j+2   ) - velx(i,j-2   )) )
       
       gvj   = ( b1*(vely(i,j+1   ) - vely(i,j-1   )) &
               + b2*(vely(i,j+2   ) - vely(i,j-2   )) )
       
       gwj   = ( b1*(velz(i,j+1   ) - velz(i,j-1   )) &
               + b2*(velz(i,j+2   ) - velz(i,j-2   )) )

       gmuj  = ( b1*(  mu(i,j+1   ) -   mu(i,j-1   )) &
               + b2*(  mu(i,j+2   ) -   mu(i,j-2   )) )
       
                           
       
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
  do j = 1 , jm 
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1 , im 
!

    distm1    = ONE/abs( yc(i,j) - ym )

    func0(i,j,1) = w(i,j,4)
    func0(i,j,2) = w(i,j,2) * velz(i,j) - mu(i,j) * gradw(i,j,1)
    func0(i,j,3) = w(i,j,3) * velz(i,j) - mu(i,j) * gradw(i,j,2) + mu(i,j) * velz(i,j) * distm1
    func0(i,j,4) = w(i,j,4) * velz(i,j) + p(i,j) + TWOTHIRD * mu(i,j) * ( gradu(i,j,1) + gradv(i,j,2) - TWO * vely(i,j) * distm1 )
    func0(i,j,5) = ( w(i,j,5) + p(i,j) ) * velz(i,j)

    func1(i,j,1) = ZERO
    func1(i,j,2) = velz(i,j)
    func1(i,j,3) = velz(i,j)
    func1(i,j,4) = velx(i,j)
    func1(i,j,5) = velz(i,j)

    func2(i,j,1) = ZERO
    func2(i,j,2) = gradw(i,j,1)
    func2(i,j,3) = gradw(i,j,2)
    func2(i,j,4) = vely(i,j)
    func2(i,j,5) = gradw(i,j,1)

    func3(i,j,1) = ZERO
    func3(i,j,2) = ZERO
    func3(i,j,3) = ZERO
    func3(i,j,4) = gradu(i,j,1) + gradv(i,j,2)
    func3(i,j,5) = velx(i,j)

    func4(i,j,1) = ZERO
    func4(i,j,2) = ZERO
    func4(i,j,3) = ZERO
    func4(i,j,4) = ZERO
    func4(i,j,5) = gradu(i,j,1)

    func5(i,j,1) = ZERO
    func5(i,j,2) = ZERO
    func5(i,j,3) = ZERO
    func5(i,j,4) = ZERO
    func5(i,j,5) = velz(i,j)

    func6(i,j,1) = ZERO
    func6(i,j,2) = ZERO
    func6(i,j,3) = ZERO
    func6(i,j,4) = ZERO
    func6(i,j,5) = gradw(i,j,2)

    func7(i,j,1) = ZERO
    func7(i,j,2) = ZERO
    func7(i,j,3) = ZERO
    func7(i,j,4) = ZERO
    func7(i,j,5) = vely(i,j)

    func8(i,j,1) = ZERO
    func8(i,j,2) = ZERO
    func8(i,j,3) = ZERO
    func8(i,j,4) = ZERO
    func8(i,j,5) = gradv(i,j,2)

    func9(i,j,1) = ZERO
    func9(i,j,2) = ZERO
    func9(i,j,3) = ZERO
    func9(i,j,4) = ZERO
    func9(i,j,5) = mu(i,j) * velx(i,j)

    func10(i,j,1) = ZERO
    func10(i,j,2) = ZERO
    func10(i,j,3) = ZERO
    func10(i,j,4) = ZERO
    func10(i,j,5) = mu(i,j) * vely(i,j)

    func11(i,j,1) = ZERO
    func11(i,j,2) = ZERO
    func11(i,j,3) = ZERO
    func11(i,j,4) = ZERO
    func11(i,j,5) = gradw(i,j,2) - velz(i,j) * distm1

    func12(i,j,1) = ZERO
    func12(i,j,2) = ZERO
    func12(i,j,3) = ZERO
    func12(i,j,4) = ZERO
    func12(i,j,5) = mu(i,j) * velz(i,j)

    func13(i,j,1) = ZERO
    func13(i,j,2) = ZERO
    func13(i,j,3) = ZERO
    func13(i,j,4) = ZERO
    func13(i,j,5) = gradu(i,j,1) + gradv(i,j,2) - TWO * vely(i,j) * distm1

    func14(i,j,1) = ZERO
    func14(i,j,2) = ZERO
    func14(i,j,3) = ZERO
    func14(i,j,4) = ZERO
    func14(i,j,5) = ZERO

    func15(i,j,1) = ZERO
    func15(i,j,2) = ZERO
    func15(i,j,3) = ZERO
    func15(i,j,4) = ZERO
    func15(i,j,5) = ZERO




!
  enddo
  enddo


end subroutine function_5p_dz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
