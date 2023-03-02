! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
subroutine bc_supandsubinlet_2d(w,loc,interf,field,nx,ny,gam,im,jm,lm,gh)
! works for inlet boundary condition on imin interface
  implicit none
! variables for dimension -----------------------------------------
  integer,intent(in) :: im,jm,gh,lm
! required arguments ----------------------------------------------
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),intent(in) :: gam
  real(8),dimension(     lm     ,     gh     , 5 ),intent(in   ) :: field
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1, 2 ),intent(in   ) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1, 2 ),intent(in   ) :: ny
! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5 ),intent(inout) :: w
! Non-required arguments -------------------------------------------
  real(8) :: ro,uu,vv,ww,velo2,p0,sound2,mach,rom1
  real(8) :: ros,us,vs,ws,ps,uts,vts,vns
  real(8) :: rod,ud,vd,wd,pd,utd,vtd,vnd
  real(8) :: p,e,ut,vt,vn,am,ap,bs,b0
  real(8) :: ro0m1,rodm1,roc0m1,rosm1
  real(8) :: nsum,nxloc,nyloc,nxnorm,nynorm
  real(8) :: gam1,epsm,eps0,epsp,roc0,rovn0
  real(8) :: r, rr, oneonrplusone
  real(8) :: wr1,wr2,wr3,wr4,wr5,pr,hr
  real(8) :: wl1,wl2,wl3,wl4,wl5,pl,hl,ee,hh
  real(8) :: rou,rov,row,roe
  real(8) :: HALF,ONE
  integer :: da
  
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
  
  
  
  
  HALF  = 0.5d0
  ONE   = 1.d0  
  gam1  = gam - ONE
  
  
!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
! face state is taken as Roe State
    wr1    = w(i, j,1)
    rom1   = ONE/wr1
    wr2    = w(i, j,2)*rom1
    wr3    = w(i, j,3)*rom1
    wr4    = w(i, j,4)*rom1
    pr     = gam1 * (w(i, j,5)-wr1*(wr2*wr2+wr3*wr3+wr4*wr4))
    hr     = (w(i, j, 5) + pr)*rom1
  
    wl1    = field(j,i,1)
    rosm1  = ONE/wl1
    wl2    = field(j,i,2)*rosm1  
    wl3    = field(j,i,3)*rosm1  
    wl4    = field(j,i,4)*rosm1  
    pl     = gam1 * (field(j,i,5)-wl1*(wl2*wl2+wl3*wl3+wl4*wl4))
    hl     = (field(j,i,5) + pl)*rosm1
  
    r  = sqrt(wr1/wl1)
    rr = sqrt(wr1*wl1)

    oneonrplusone = ONE/(r+ONE)

    uu     = (wr2*r+wl2)*oneonrplusone
    vv     = (wr3*r+wl3)*oneonrplusone
    ww     = (wr4*r+wl4)*oneonrplusone
    hh     = (hr*r +hl )*oneonrplusone
    velo2  = uu*uu + vv*vv + ww
    ee     = HALF*velo2
  
    sound2 = gam1 * (hh-ee)
    mach   = sqrt(velo2/sound2)
    if (mach.lt.ONE) then
!
      nxloc  = nx(i,j,kdir)
      nyloc  = ny(i,j,kdir)
      nsum   = ONE/ sqrt(nxloc*nxloc + nyloc*nyloc)
      nxnorm = nxloc*nsum
      nynorm = nyloc*nsum
      roc0   = rr*sqrt(sound2)
      rovn0  = rr*(uu*nxnorm+vv*nynorm)
      epsm   = HALF + sign(HALF,roc0-rovn0)
      eps0   = HALF + sign(HALF,-rovn0)
      epsp   = HALF + sign(HALF,-roc0-rovn0)
    endif
        
  
    do de = 1,gh
! dirichlet
      w(i-de*i0,j-de*j0,1) = field(l,de,1)
      w(i-de*i0,j-de*j0,2) = field(l,de,2)
      w(i-de*i0,j-de*j0,3) = field(l,de,3)
      w(i-de*i0,j-de*j0,4) = field(l,de,4)
      w(i-de*i0,j-de*j0,5) = field(l,de,5)
      
      da = de - 1
      if (mach.lt.ONE) then
! l-state
        rod   = w(i-de*i0,j-de*j0,1)
        rom1  = ONE/rod
        ud    = w(i-de*i0,j-de*j0,2)*rom1
        vd    = w(i-de*i0,j-de*j0,3)*rom1
        wd    = w(i-de*i0,j-de*j0,4)*rom1
        pd    = gam1*(w(i-de*i0,j-de*j0,5) - rod*HALF*(ud*ud + vd*vd + wd*wd))
        vnd   = (ud*nxnorm + vd*nynorm)
        utd   = ud - vnd * nxnorm
        vtd   = vd - vnd * nynorm
        hl    = (w(i-de*i0,j-de*j0,5)+pd)*rom1
        
! r-state
        
        ros   = w(i-da*i0,j-da*j0,1)
        rosm1 = ONE/ros
        us    = w(i-da*i0,j-da*j0,2)*rosm1
        vs    = w(i-da*i0,j-da*j0,3)*rosm1
        ws    = w(i-da*i0,j-da*j0,4)*rosm1
        ps    = gam1*(w(i-da*i0,j-da*j0,5) - ros*HALF*(us*us + vs*vs + ws*ws))
        vns   = (us*nxnorm + vs*nynorm)
        uts   = us - vns * nxnorm
        vts   = vs - vns * nynorm
        hr    = (w(i-da*i0,j-da*j0,5)+ps)*rosm1
        
        r  = SQRT(ros*rom1)
        rr = SQRT(ros*rod)

        oneonrplusone = ONE/(r+ONE)

        uu     = (us*r+ud)*oneonrplusone
        vv     = (vs*r+vd)*oneonrplusone
        ww     = (ws*r+wd)*oneonrplusone
        ee     = HALF*(uu*uu + vv*vv + ww*ww)
        hh     = (hr*r+hl)*oneonrplusone
        sound2 = gam1*(hh-ee)
        
        roc0   = rr*sqrt(sound2)
        roc0m1 = ONE/roc0
        
! updated State -----------------------------------------------------------
! Wave directions (eps) from Bnd are used in the ghost cells
        
        ut = eps0*uts + (ONE-eps0)*utd
        vt = eps0*vts + (ONE-eps0)*vtd
        

        am = epsm*(ps-roc0*vns) + (ONE-epsm)*(pd-roc0*vnd)
        ap = epsp*(ps+roc0*vns) + (ONE-epsp)*(pd+roc0*vnd)
        vn = (ap-am) * HALF *roc0m1
        p  = (ap+am) * HALF

        bs = (p-ps)*rr*rr*roc0m1*roc0m1 + ros
        b0 = (p-pd)*rr*rr*roc0m1*roc0m1 + rod
        ro = eps0*bs + (ONE-eps0)*b0
        rom1 = ONE/ro
        roe  = p/gam1
      
        rou = ro  * (ut + vn*nxnorm)
        rov = ro  * (vt + vn*nynorm)
        
        row =field(j, de ,4)
        
        w(i-de*i0,j-de*j0,1) = ro
        w(i-de*i0,j-de*j0,2) = rou
        w(i-de*i0,j-de*j0,3) = rov
        w(i-de*i0,j-de*j0,4) = row
        w(i-de*i0,j-de*j0,5) = roe+HALF*(rou*rou+rov*rov+row*row)/ro  
      
      endif 
    enddo
  enddo
end subroutine bc_supandsubinlet_2d

subroutine bc_filteringilo_2d(w,coef,coefbnd,im,jm,gh,cl,em)
  implicit none
! variables for dimension -----------------------------------------
  integer,intent(in) :: im,jm,gh,cl,em
  real(8),dimension(cl), intent(in) :: coef
  real(8),dimension(gh,2*gh+1), intent(in) :: coefbnd
!-------------------------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , em ),intent(inout) :: w
!-------------------------------------------------------------------
  real(8),dimension(1-gh:gh+1  ,1-gh:jm+gh  , em ) :: wt
  real(8),dimension(5 ) :: phi
  integer :: i,j,g,l,sten
  real(8) :: sig
!-------------------------------------------------------------------
  sten = 2*gh
  sig = 0.4d0
  do i=0,sten
    wt(1-gh+i,:,:) = w(1-gh+i,:,:)
  enddo 
  
  do j=1+gh,jm-gh
! dir j
    do i= 1-gh,0
      phi(:) = coef(1)*wt(i,j,:)
        do l = 1,cl-1
          phi(:) = phi(:) + coef(l+1)*(wt(i,j+l,:)+wt(i,j-l,:))
        enddo
      w(i,j,:) = w(i,j,:) - sig*phi(:) 
    enddo
! dir i
    do i= 1-gh,0
        do l = 0,sten
          phi(:) = phi(:) + coefbnd(i+gh,l+1)*wt(1-gh+l,j,:)
        enddo
      w(i,j,:) = w(i,j,:) - sig*phi(:) 
    enddo
  enddo
  
end subroutine bc_filteringilo_2d 


