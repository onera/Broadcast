! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!
!==============================================================================
!     BC Pressure No-Ref_2D:  Euler & NS
!==============================================================================
!
subroutine bc_pressure_2d(w,loc,interf,pext,noref,gam,nx,ny,im,jm,gh,em)
!
  implicit none
! Variables for dimension ---------------------------------------------------
  integer,intent(in) :: im,jm
  integer,intent(in) :: gh,em
! Input variables -----------------------------------------------------------
  real(8),intent(in) :: pext
  logical,intent(in) :: noref
  real(8),intent(in) :: gam
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(in) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(in) :: ny
! Returned objects ----------------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(inout) :: w
! Local variables -----------------------------------------------------------
  real(8) :: ro0,ro0m1,uu,vv,ww,rovn0,epsm,epsp,eps0,p0
  real(8) :: ros,rosm1,us,vs,ws,ps,uts,vts,vns,c2s,c2sm1,rocs,rocsm1
  real(8) :: rod,rodm1,ud,vd,wd,pd,utd,vtd,vnd
  real(8) :: ro,rom1,p,e,ut,vt,wt,vn,am,ap,bs,b0
  real(8) :: rou,rov,row,roe,roe1,roe2
  real(8) :: c20,c20m1,roc0,roc0m1,alpha
  real(8) :: nxloc,nyloc,nsum,nsumi,sens
  real(8) :: nxnorm,nynorm,gam1,ONE,HALF,TWO
  integer :: da,i1,j1,m
! ---------------------------------------------------------------------------
! Local variables -----------------------------------------------------------
  integer :: kdir,de,i,j,imin,imax,jmin,jmax,lmin,lmax,l,i0,j0,high
! ---------------------------------------------------------------------------
!
  imin=interf(1,1)
  jmin=interf(1,2)
  imax=interf(2,1)
  jmax=interf(2,2)
  
!   write(200,*) loc, 'imin = ', interf(1,1)
!   write(200,*) loc, 'jmin = ', interf(1,2)
!   write(200,*) loc, 'imax = ', interf(2,1)
!   write(200,*) loc, 'jmax = ', interf(2,2)
  
  i0 = 0
  j0 = 0
  high = 0
  
  i = imin
  j = jmin
  
  lmin = 1
  
  if (loc.eq.'Ilo') then
    kdir = 1
    i0 = 1
    lmax = jmax - jmin + 1
  else if (loc.eq.'Ihi')then
    kdir = 1
    i0 = -1
    lmax = jmax - jmin + 1
    high = 1
  else if (loc.eq.'Jlo')then
    kdir=2
    j0 = 1
    lmax = imax - imin + 1
  else if (loc.eq.'Jhi')then
    kdir=2
    j0 = -1
    lmax = imax - imin + 1
    high = 1
  endif
  
  
  
  ONE  =  1.d0
  HALF = 0.5d0
  TWO  =  2.d0
  
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
! NB: normals into the domain
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
    c20   = gam*p0*ro0m1
    c20m1 = ONE/(c20)
    roc0  = ro0*sqrt(c20)
    roc0m1= ONE/roc0
    

! sch-State ---------------------------------------------------------------
    ros   = w(i,j,1)
    rosm1 = ONE/ros
    us    = w(i,j,2)*rosm1
    vs    = w(i,j,3)*rosm1
    ws    = w(i,j,4)*rosm1
    ps    = gam1*(w(i,j,5) - ros*HALF*(us*us + vs*vs + ws*ws))
    vns   = (us*nxnorm + vs*nynorm)
    uts   = us - (vns * nxnorm)
    vts   = vs - (vns * nynorm)
    
    
! d-State --------------------------------------------------------------
    pd  = pext
    rod = ros + (pd - ps)*c20m1
    vnd = vns + (pd - ps)*roc0m1 ! normal into the domain
    utd = uts
    vtd = vts
    ud  = utd + (vnd * nxnorm)
    vd  = vtd + (vnd * nynorm)
    wd  = ws


! updated State -----------------------------------------------------------
    IF (noref) THEN
      rovn0 = (w(i,j,2)*nxnorm + w(i,j,3)*nynorm )
      epsm  = HALF + sign(HALF,roc0-rovn0)
      eps0  = HALF + sign(HALF,-rovn0)
      epsp  = HALF + sign(HALF,-roc0-rovn0)

      ut = eps0*uts + (ONE-eps0)*utd
      vt = eps0*vts + (ONE-eps0)*vtd
      wt = eps0*ws  + (ONE-eps0)*wd 
      
      am = epsm*(ps-roc0*vns) + (ONE-epsm)*(pd-roc0*vnd)
      ap = epsp*(ps+roc0*vns) + (ONE-epsp)*(pd+roc0*vnd)
      vn = (ap-am) * HALF *roc0m1
      p  = (ap+am) * HALF

      bs = (p-ps)*ro0*ro0*roc0m1*roc0m1 + ros
      b0 = (p-pd)*ro0*ro0*roc0m1*roc0m1 + rod
      ro = eps0*bs + (ONE-eps0)*b0
      rom1 = ONE/ro
      roe  = p/gam1
      
      rou = ro * (ut + vn*nxnorm)
      rov = ro * (vt + vn*nynorm)
      row = ro * wt
    ELSE
      p = pd
      ro   = rod
      rom1 = ONE/ro
      rou  = ro * ud
      rov  = ro * vd
      row  = ro * wd
      roe  = p/gam1
    END IF
    
    do de = 1, gh
      da = de - 1

      w(i-de*i0,j-de*j0,1) = ro
      w(i-de*i0,j-de*j0,2) = rou
      w(i-de*i0,j-de*j0,3) = rov
      w(i-de*i0,j-de*j0,4) = row
      w(i-de*i0,j-de*j0,5) = roe + HALF * rom1 * (rou*rou+rov*rov+row*row)
      do m = 6, em
        w(i-de*i0,j-de*j0,m) = w(i,j,m) * ro0m1 * ro
      enddo

    enddo
  enddo
end subroutine bc_pressure_2d
!
