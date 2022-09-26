! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

! =============================================================================
!          consistent fluxes for DNC7 2D
! =============================================================================
!
subroutine function_9p_dz(func0,func1,func2,func3,func4,func5,func6,func7,func8,&
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
  real(8) :: fvro1,fvro2, fvrou1,fvrou2,fvrov1,fvrov2,fvrow1,fvrow2,fvroe1,fvroe2
  real(8) :: gvro1,gvro2, gvrou1,gvrou2,gvrov1,gvrov2,gvrow1,gvrow2,gvroe1,gvroe2
  real(8) :: dissro1,dissro2, dissrou1,dissrou2,dissrov1,dissrov2,dissrow1,dissrow2,dissroe1,dissroe2
  real(8) :: cpprandtl, mmu, lambda, uu, vv, ww, ro,rom1,htot,eloc,ec
  real(8) :: predro1,predrou1,predrov1,predrow1,predroe1, eps2,eps4
  real(8) :: predro2,predrou2,predrov2,predrow2,predroe2, rspec
  real(8) :: divu2,vort2,dxm1,dym1,dxm2,dym2,distm1
  real(8) :: gui,gvi,gwi,gmui
  real(8) :: guj,gvj,gwj,gmuj
  real(8) :: nxloc, nyloc, sn, invsn, sc1, sc2
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
 
  denom = 1.d0/840.d0
  ! Coef for grad o8
  b1 =  672.d0 * denom! 4.d0/5.d0
  b2 = -168.d0 * denom!-1.d0/5.d0
  b3 =  32.d0  * denom! 4.d0/105.d0
  b4 = - 3.d0  * denom!-1.d0/280.d0
   
  betas = muref*(tref + cs)/(sqrt(tref)*tref)
    
  ! Primitives
#include "rhs/primvisc.F" 
  
  ! Work on interior domain minus one cell
  
!$AD II-LOOP
  do j = 1, jm
!$AD II-LOOP
!DIR$ IVDEP      
  do i = 1, im
#include "rhs/gradop_9pi.F"
#include "rhs/gradop_9pj.F"
#include "rhs/gradient.F"
  enddo
  enddo
  
  ! trick for ducros in dissipation
#include "rhs/gradveloingh.F"
  
  
!$AD II-LOOP
  do j = 1 , jm 
!$AD II-LOOP
!DIR$ IVDEP      
  do i = 1 , im 
    !
#include "matrix_dz_axi/function_dz.F"
    !
  enddo
  enddo
  

end subroutine function_9p_dz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
