! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

! =============================================================================
!          consistent fluxes for DNC3 2D
! =============================================================================
!
subroutine flux_num_dnc3_polar_2d(residu,w,ym,x0,y0,nx,ny,xc,yc,vol,volf,gh,cp,cv,prandtl,&
                            gam,rgaz,cs,muref,tref,s_suth,k2,k4,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  real(8),intent(in) :: cp,cv,prandtl,gam,rgaz ! thermo
  real(8),intent(in) :: k2,k4,ym ! dissipation, axis position 
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
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     , 2) :: grad_rho,grad_temp
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
  real(8) :: divu2,vort2,dx,dy,dxm1,dym1,dxm2,dym2
  real(8) :: grhoi,gui,gvi,gwi,gti,gmui
  real(8) :: grhoj,guj,gvj,gwj,gtj,gmuj
  real(8) :: rhom,rhomr,rhoml,rhom1l,c2l,c2r,rr,r,u,ur,ul,vr,vl,wr,wl,c2x,nx2,ny2
  real(8) :: ab, sq,ducros1,ducros2,k_sensor1,k_sensor2
  real(8) :: b1,b2,b3,b4,b5,c1,c2,c3,c4,c5,wiggle, denom,betas
  real(8) :: nxloc, nyloc, sn, invsn, sc1, sc2
  real(8) :: d1, d2
  real(8) :: c3d0,c3d1
  real(8) :: pw,ct0,ct1,ct2,cvm1
  real(8) :: coef,omrr,test,diffro,diffrou,diffrov,diffrow,diffroe,v
  real(8) :: HALF,ONE,ZERO,TWO,TWOTHIRD,FOURTH,TWELFTH 
  real(8) :: TWENTYFOURTH,ccross,divu,utheta,voldistm1,distm1
  real(8) :: src1,src2,src3,src4,src5
  real(8) :: srce1,srce2,srce3,srce4,srce5
  real(8) :: srcv1,srcv2,srcv3,srcv4,srcv5
  real(8) :: t_xy, t_yy, t_yz  
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
  
  b1 =  8.d0 * denom
  b2 = -       denom 
  
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
  
#include "rhs/primvisc.F"  
  
  ! Work on interior domain minus one cell
  
!$AD II-LOOP
  do j = 1, jm
!$AD II-LOOP
!DIR$ IVDEP            
  do i = 1, im
#include "rhs/gradop_prim_5pi.F"
#include "rhs/gradop_prim_5pj.F"
#include "rhs/gradient_prim.F"
  enddo
  enddo
  
  
  ! trick for ducros in dissipation
#include "rhs/gradveloingh.F"
  
  ! do j = 3 , jm + 1
  ! do i = 1 , im + 1
  ! nbjb = 11 
!$AD II-LOOP
  do j = 3, jm+1
!DIR$ IVDEP            
!$AD II-LOOP
  do i = 1 , im + 1
    !
#include "rhs/euler_o4.F"  
#include "rhs/predictor_5p.F"
#include "rhs/flux_visqueux_polar_o2_i.F"
#include "rhs/flux_visqueux_polar_o2_j.F"
#include "rhs/dissipation_ducros.F"
#include "rhs/fluxnumassembly_i.F"
#include "rhs/fluxnumassembly_j.F"
    !
  enddo
  enddo
  
  ! Manage Wall Fluxes
   
  ! coef for near wall off-centered fluxes
! #include rhs/coefnearbnd_5p
  c3d0 = HALF
  c3d1 = HALF
   ! off centered schemes for fluxes near wall
      
  j = 2       
!$AD II-LOOP
!DIR$ IVDEP            
  do i = 1,im+1
    
#include "rhs/euler_o4_i.F"  
#include "rhs/nearbndfluxes3demi_5p.F"
#include "rhs/predictor_5p.F"
#include "rhs/flux_visqueux_polar_o2_i.F"
#include "rhs/flux_visqueux_polar_o2_j.F"
#include "rhs/dissipation_ducros.F"    
#include "rhs/fluxnumassembly_i.F"
#include "rhs/fluxnumassembly_j.F"
                                                           
  enddo
  ! force dp/dn = 0. at 2nd oder
  ct0 = 1.125d0! 9/8
  ct1 =-0.125d0!-1/8
  
  ! force dp/dn = 0. at 3rd oder
  ! ct0 = 225.d0/184.d0
  ! ct1 = -25.d0/92.d0
  ! ct2 = 9.d0/184.d0
  
  j = 1
!$AD II-LOOP
!DIR$ IVDEP            
  do i = 1,im+1
    ! for idir (only works for infinite Wall)   
#include "rhs/euler_o4_i.F"
#include "rhs/predictor_5p_i.F"
#include "rhs/flux_visqueux_polar_o2_i.F"
#include "rhs/dissipation_ducros_i.F" 
#include "rhs/fluxnumassembly_i.F"
#include "rhs/fluxwall_polar.F"   

    
  enddo     
   
  
  ! fluxes balance
#include "rhs/balance_polar.F"


end subroutine flux_num_dnc3_polar_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
