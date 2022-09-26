! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

!
!==============================================================================
!     BC No-Ref_2D:  Euler & NS
!==============================================================================
!
subroutine bc_no_reflexion_forhessian_2d(wout,w,wbd,loc,interf,nx,ny,gam,gh,im,jm,lm)
  !
  implicit none
  ! Variables for dimension ---------------------------------------------------
  integer,intent(in) :: im,jm,lm
  integer,intent(in) :: gh
  ! Input variables -----------------------------------------------------------
  real(8),intent(in) :: gam
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(lm,5),intent(in) :: wbd
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(in) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(in) :: ny
  ! Returned objects ----------------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,5),intent(inout) :: w
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,5),intent(inout) :: wout
  ! Local variables -----------------------------------------------------------
  real(8) :: ro0,uu,vv,ww,roc0,rovn0,epsm,epsp,eps0,p0
  real(8) :: ros,us,vs,ws,ps,uts,vts,vns
  real(8) :: rod,ud,vd,wdi,pd,utd,vtd,vnd
  real(8) :: ro,p,e,ut,vt,wt,vn,am,ap,bs,b0
  real(8) :: ro0m1,rodm1,roc0m1,rosm1,rom1,roe
  real(8) :: nxloc,nyloc,nsum,nsumi,sens
  real(8) :: nxnorm,nynorm,gam1,ONE,HALF,TWO
  real(8) :: rou,rov,row,roe1,roe2
  integer :: da,i1,j1
  ! ---------------------------------------------------------------------------
#include "init_2d.F"
  
  ONE  =  1.d0
  HALF = 0.5d0
  TWO  =  2.d0

  wout = w

  gam1 = gam - ONE
  !
  sens = float(i0+j0)
  !
  nxloc = 0.d0
  nyloc = 0.d0
  i1 = i0 * i0
  j1 = j0 * j0
!$AD II-LOOP
  do l = lmin , lmax
    i =  imin + (l-lmin)*j1
    j =  jmin + (l-lmin)*i1
    ! normals------------------------------------------------------------------
    nxloc = nx(i+high*i1,j+high*j1,kdir)
    nyloc = ny(i+high*i1,j+high*j1,kdir)
    ! highis used to select the last normal when the face is of type *hi
    nsum    = sqrt(nxloc*nxloc + nyloc*nyloc)
    nsumi   = ONE/nsum
    nxnorm  = nxloc * nsumi * sens
    nynorm  = nyloc * nsumi * sens

    ! 0-State --------------------------------------------------------------
    ro0   = w(i,j,1)
    ro0m1 = ONE/ro0
    uu    = w(i,j,2)*ro0m1
    vv    = w(i,j,3)*ro0m1
    ww    = w(i,j,4)*ro0m1
    roe1  = w(i,j,5) - HALF*ro0*(uu*uu + vv*vv + ww*ww)
    p0    = gam1*(roe1)
    roc0  = ro0*sqrt(gam*p0*ro0m1)
    roc0m1= ONE/roc0
    rovn0 = (w(i,j,2)*nxnorm + w(i,j,3)*nynorm )
    epsm  = HALF + sign(HALF,roc0-rovn0)
    eps0  = HALF + sign(HALF,-rovn0)
    epsp  = HALF + sign(HALF,-roc0-rovn0)
    
    
    ! d-State ------------for computational domain border: dState = state_inf ---
    rod   = wbd(l,1)
    rodm1 = ONE/rod
    ud    = wbd(l,2)*rodm1
    vd    = wbd(l,3)*rodm1
    wdi    = wbd(l,4)*rodm1
    pd    = gam1*(wbd(l,5) - rod*HALF*(ud*ud + vd*vd + wdi*wdi))
    vnd   = (ud*nxnorm + vd*nynorm)
    utd   = ud - vnd * nxnorm
    vtd   = vd - vnd * nynorm
    

    ! sch-State ---------------------------------------------------------------
    ros   = w(i,j,1)
    rosm1 = ONE/ros
    us    = w(i,j,2)*rosm1
    vs    = w(i,j,3)*rosm1
    ws    = w(i,j,4)*rosm1
    ps    = gam1*(w(i,j,5) - ros*HALF*(us*us + vs*vs + ws*ws))
    vns   = (us*nxnorm + vs*nynorm)
    uts   = us - vns * nxnorm
    vts   = vs - vns * nynorm
    ! wts   = ws - vns * nznorm
    ! updated State -----------------------------------------------------------
    ut = eps0*uts + (ONE-eps0)*utd
    vt = eps0*vts + (ONE-eps0)*vtd
    wt = eps0*ws  + (ONE-eps0)*wdi 
    ! wt = TWO *wd  - ws !extrap o2

    am = epsm*(ps-roc0*vns) + (ONE-epsm)*(pd-roc0*vnd)
    ap = epsp*(ps+roc0*vns) + (ONE-epsp)*(pd+roc0*vnd)
    vn = (ap-am) * HALF *roc0m1
    p  = (ap+am) * HALF

    bs = (p-ps)*ro0*ro0*roc0m1*roc0m1 + ros
    b0 = (p-pd)*ro0*ro0*roc0m1*roc0m1 + rod
    ro = eps0*bs + (ONE-eps0)*b0
    rom1 = ONE/ro
    roe  = p/gam1
    
    ! rou = ro  * (TWO*(ut + vn*nxnorm) - us)
    ! rov = ro  * (TWO*(vt + vn*nynorm) - vs)
    rou = ro  * (ut + vn*nxnorm)
    rov = ro  * (vt + vn*nynorm)

    row = ro  * wt

    roe2 = roe + HALF * rom1 * (rou*rou+rov*rov+row*row)
    
    do de = 1, gh
      da = de - 1

      wout(i-de*i0,j-de*j0,1) = ro
      wout(i-de*i0,j-de*j0,2) = rou
      wout(i-de*i0,j-de*j0,3) = rov
      wout(i-de*i0,j-de*j0,4) = row
      ! w(i-de*i0,j-de*j0,5) = roe + HALF * rom1 * (rou*rou+rov*rov+row*row)

      wout(i-de*i0,j-de*j0,5) = roe2
             
      !extrap o2 
      ! roe2 =  roe
      ro   =  TWO * ro  - wout(i-da*i0,j-da*j0,1)  
      rou  =  TWO * rou - wout(i-da*i0,j-da*j0,2)  
      rov  =  TWO * rov - wout(i-da*i0,j-da*j0,3)  
      row  =  TWO * row - wout(i-da*i0,j-da*j0,4)  
      ! roe  =  TWO * roe - roe1
      ! roe1 =  roe2
      ! rom1 =  ONE/ro

      roe2 =  TWO * roe2- wout(i-da*i0,j-da*j0,5)  

    enddo
  enddo
  
end subroutine bc_no_reflexion_forhessian_2d
!
