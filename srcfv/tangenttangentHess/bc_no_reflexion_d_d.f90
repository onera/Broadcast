! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_no_reflexion_2d_d in forward (tangent) mode (with options with!SliceDeadControl with!SliceDeadInstrs wit
!h!StaticTaping):
!   variations   of useful results: w wd
!   with respect to varying inputs: w
!   RW status of diff variables: w:in-out wd:out
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_no_reflexion_2d in forward (tangent) mode (with options with!SliceDeadControl with!SliceDeadInstrs with!
!StaticTaping):
!   variations   of useful results: w
!   with respect to varying inputs: w
!   RW status of diff variables: w:in-out
!
!==============================================================================
!     BC No-Ref_2D:  Euler & NS
!==============================================================================
!
SUBROUTINE BC_NO_REFLEXION_2D_D_D(w, wd0, wd, wdd, wbd, loc, interf, nx&
& , ny, gam, gh, im, jm, lm)
  IMPLICIT NONE
! Variables for dimension ---------------------------------------------------
  INTEGER, INTENT(IN) :: im, jm, lm
  INTEGER, INTENT(IN) :: gh
! Input variables -----------------------------------------------------------
  REAL*8, INTENT(IN) :: gam
  CHARACTER(len=3), INTENT(IN) :: loc
  INTEGER, DIMENSION(2, 2), INTENT(IN) :: interf
  REAL*8, DIMENSION(lm, 5), INTENT(IN) :: wbd
  REAL*8, DIMENSION(1-gh:im+gh+1, 1-gh:jm+gh+1, 2), INTENT(IN) :: nx
  REAL*8, DIMENSION(1-gh:im+gh+1, 1-gh:jm+gh+1, 2), INTENT(IN) :: ny
! Returned objects ----------------------------------------------------------
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: w
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wd0
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wd
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wdd
! Local variables -----------------------------------------------------------
  REAL*8 :: ro0, uu, vv, ww, roc0, rovn0, epsm, epsp, eps0, p0
  REAL*8 :: ro0d0, uud0, vvd0, wwd0, roc0d0, p0d0
  REAL*8 :: ro0d, uud, vvd, wwd, roc0d, p0d
  REAL*8 :: ro0dd, uudd, vvdd, wwdd, roc0dd, p0dd
  REAL*8 :: ros, us, vs, ws, ps, uts, vts, vns
  REAL*8 :: rosd0, usd0, vsd0, wsd0, psd0, utsd0, vtsd0, vnsd0
  REAL*8 :: rosd, usd, vsd, wsd, psd, utsd, vtsd, vnsd
  REAL*8 :: rosdd, usdd, vsdd, wsdd, psdd, utsdd, vtsdd, vnsdd
  REAL*8 :: rod, ud, vd, wdi, pd, utd, vtd, vnd
  REAL*8 :: ro, p, e, ut, vt, wt, vn, am, ap, bs, b0
  REAL*8 :: rod1, pd1, utd1, vtd1, wtd0, vnd1, amd0, apd0, bsd0, b0d0
  REAL*8 :: rod0, pd0, utd0, vtd0, wtd, vnd0, amd, apd, bsd, b0d
  REAL*8 :: rod0d, pd0d, utd0d, vtd0d, wtdd, vnd0d, amdd, apdd, bsdd, &
& b0dd
  REAL*8 :: ro0m1, rodm1, roc0m1, rosm1, rom1, roe
  REAL*8 :: ro0m1d0, roc0m1d0, rosm1d0, rom1d0, roed0
  REAL*8 :: ro0m1d, roc0m1d, rosm1d, rom1d, roed
  REAL*8 :: ro0m1dd, roc0m1dd, rosm1dd, rom1dd, roedd
  REAL*8 :: nxloc, nyloc, nsum, nsumi, sens
  REAL*8 :: nxnorm, nynorm, gam1, one, half, two
  REAL*8 :: rou, rov, row, roe1, roe2
  REAL*8 :: roud0, rovd0, rowd0, roe1d0, roe2d0
  REAL*8 :: roud, rovd, rowd, roe1d, roe2d
  REAL*8 :: roudd, rovdd, rowdd, roe1dd, roe2dd
  INTEGER :: da, i1, j1
! ---------------------------------------------------------------------------
! Local variables -----------------------------------------------------------
  INTEGER :: kdir, de, i, j, imin, imax, jmin, jmax, lmin, lmax, l, i0, &
& j0, high
  INTRINSIC FLOAT
  INTRINSIC SQRT
  INTRINSIC SIGN
  REAL*8 :: arg1
  REAL*8 :: arg1d0
  REAL*8 :: arg1d
  REAL*8 :: arg1dd
  REAL*8 :: result1
  REAL*8 :: result1d0
  REAL*8 :: result1d
  REAL*8 :: result1dd
  REAL*8 :: temp
  REAL*8 :: tempd
  REAL*8 :: temp0
  REAL*8 :: temp1
! ---------------------------------------------------------------------------
!
  imin = interf(1, 1)
  jmin = interf(1, 2)
  imax = interf(2, 1)
  jmax = interf(2, 2)
!   write(200,*) loc, 'imin = ', interf(1,1)
!   write(200,*) loc, 'jmin = ', interf(1,2)
!   write(200,*) loc, 'imax = ', interf(2,1)
!   write(200,*) loc, 'jmax = ', interf(2,2)
  i0 = 0
  j0 = 0
  high = 0
  lmin = 1
  IF (loc .EQ. 'Ilo') THEN
    kdir = 1
    i0 = 1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Ihi') THEN
    kdir = 1
    i0 = -1
    lmax = jmax - jmin + 1
    high = 1
  ELSE IF (loc .EQ. 'Jlo') THEN
    kdir = 2
    j0 = 1
    lmax = imax - imin + 1
  ELSE IF (loc .EQ. 'Jhi') THEN
    kdir = 2
    j0 = -1
    lmax = imax - imin + 1
    high = 1
  END IF
  one = 1.d0
  half = 0.5d0
  two = 2.d0
  gam1 = gam - one
!
  sens = FLOAT(i0 + j0)
!
  i1 = i0*i0
  j1 = j0*j0
  DO l=lmin,lmax
    i = imin + (l-lmin)*j1
    j = jmin + (l-lmin)*i1
! normals------------------------------------------------------------------
    nxloc = nx(i+high*i1, j+high*j1, kdir)
    nyloc = ny(i+high*i1, j+high*j1, kdir)
! highis used to select the last normal when the face is of type *hi
    arg1 = nxloc*nxloc + nyloc*nyloc
    nsum = SQRT(arg1)
    nsumi = one/nsum
    nxnorm = nxloc*nsumi*sens
    nynorm = nyloc*nsumi*sens
! 0-State --------------------------------------------------------------
    ro0dd = wdd(i, j, 1)
    ro0d = wd(i, j, 1)
    ro0d0 = wd0(i, j, 1)
    ro0 = w(i, j, 1)
    temp0 = ro0d/(ro0*ro0)
    ro0m1dd = -(one*(ro0dd-temp0*2*ro0*ro0d0)/ro0**2)
    ro0m1d = -(one*temp0)
    ro0m1d0 = -(one*ro0d0/ro0**2)
    ro0m1 = one/ro0
    uudd = wd(i, j, 2)*ro0m1d0 + ro0m1*wdd(i, j, 2) + ro0m1d*wd0(i, j, 2&
&     ) + w(i, j, 2)*ro0m1dd
    uud = ro0m1*wd(i, j, 2) + w(i, j, 2)*ro0m1d
    uud0 = ro0m1*wd0(i, j, 2) + w(i, j, 2)*ro0m1d0
    uu = w(i, j, 2)*ro0m1
    vvdd = wd(i, j, 3)*ro0m1d0 + ro0m1*wdd(i, j, 3) + ro0m1d*wd0(i, j, 3&
&     ) + w(i, j, 3)*ro0m1dd
    vvd = ro0m1*wd(i, j, 3) + w(i, j, 3)*ro0m1d
    vvd0 = ro0m1*wd0(i, j, 3) + w(i, j, 3)*ro0m1d0
    vv = w(i, j, 3)*ro0m1
    wwdd = wd(i, j, 4)*ro0m1d0 + ro0m1*wdd(i, j, 4) + ro0m1d*wd0(i, j, 4&
&     ) + w(i, j, 4)*ro0m1dd
    wwd = ro0m1*wd(i, j, 4) + w(i, j, 4)*ro0m1d
    wwd0 = ro0m1*wd0(i, j, 4) + w(i, j, 4)*ro0m1d0
    ww = w(i, j, 4)*ro0m1
    tempd = 2*uu*uud0 + 2*vv*vvd0 + 2*ww*wwd0
    temp = uu*uu + vv*vv + ww*ww
    temp0 = 2*uu*uud + 2*vv*vvd + 2*ww*wwd
    roe1dd = wdd(i, j, 5) - half*(ro0d*tempd+temp*ro0dd+temp0*ro0d0+ro0*&
&     (uud*2*uud0+2*uu*uudd+vvd*2*vvd0+2*vv*vvdd+wwd*2*wwd0+2*ww*wwdd))
    roe1d = wd(i, j, 5) - half*(temp*ro0d+ro0*temp0)
    roe1d0 = wd0(i, j, 5) - half*(temp*ro0d0+ro0*tempd)
    roe1 = w(i, j, 5) - half*(ro0*temp)
    p0dd = gam1*roe1dd
    p0d = gam1*roe1d
    p0d0 = gam1*roe1d0
    p0 = gam1*roe1
    arg1dd = gam*(p0d*ro0m1d0+ro0m1*p0dd+ro0m1d*p0d0+p0*ro0m1dd)
    arg1d = gam*(ro0m1*p0d+p0*ro0m1d)
    arg1d0 = gam*(ro0m1*p0d0+p0*ro0m1d0)
    arg1 = gam*p0*ro0m1
    temp0 = SQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      tempd = 0.0_8
    ELSE
      tempd = arg1d0/(2.0*temp0)
    END IF
    temp = temp0
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.0_8
      result1dd = 0.0_8
    ELSE
      temp0 = arg1d/(2.0*temp)
      result1dd = (arg1dd-temp0*2.0*tempd)/(2.0*temp)
      result1d = temp0
    END IF
    result1d0 = tempd
    result1 = temp
    roc0dd = ro0d*result1d0 + result1*ro0dd + result1d*ro0d0 + ro0*&
&     result1dd
    roc0d = result1*ro0d + ro0*result1d
    roc0d0 = result1*ro0d0 + ro0*result1d0
    roc0 = ro0*result1
    temp0 = roc0d/(roc0*roc0)
    roc0m1dd = -(one*(roc0dd-temp0*2*roc0*roc0d0)/roc0**2)
    roc0m1d = -(one*temp0)
    roc0m1d0 = -(one*roc0d0/roc0**2)
    roc0m1 = one/roc0
    rovn0 = w(i, j, 2)*nxnorm + w(i, j, 3)*nynorm
    epsm = half + SIGN(half, roc0 - rovn0)
    eps0 = half + SIGN(half, -rovn0)
    epsp = half + SIGN(half, -roc0 - rovn0)
! d-State ------------for computational domain border: dState = state_inf ---
    rod = wbd(l, 1)
    rodm1 = one/rod
    ud = wbd(l, 2)*rodm1
    vd = wbd(l, 3)*rodm1
    wdi = wbd(l, 4)*rodm1
    pd = gam1*(wbd(l, 5)-rod*half*(ud*ud+vd*vd+wdi*wdi))
    vnd = ud*nxnorm + vd*nynorm
    utd = ud - vnd*nxnorm
    vtd = vd - vnd*nynorm
! sch-State ---------------------------------------------------------------
    rosdd = wdd(i, j, 1)
    rosd = wd(i, j, 1)
    rosd0 = wd0(i, j, 1)
    ros = w(i, j, 1)
    temp0 = rosd/(ros*ros)
    rosm1dd = -(one*(rosdd-temp0*2*ros*rosd0)/ros**2)
    rosm1d = -(one*temp0)
    rosm1d0 = -(one*rosd0/ros**2)
    rosm1 = one/ros
    usdd = wd(i, j, 2)*rosm1d0 + rosm1*wdd(i, j, 2) + rosm1d*wd0(i, j, 2&
&     ) + w(i, j, 2)*rosm1dd
    usd = rosm1*wd(i, j, 2) + w(i, j, 2)*rosm1d
    usd0 = rosm1*wd0(i, j, 2) + w(i, j, 2)*rosm1d0
    us = w(i, j, 2)*rosm1
    vsdd = wd(i, j, 3)*rosm1d0 + rosm1*wdd(i, j, 3) + rosm1d*wd0(i, j, 3&
&     ) + w(i, j, 3)*rosm1dd
    vsd = rosm1*wd(i, j, 3) + w(i, j, 3)*rosm1d
    vsd0 = rosm1*wd0(i, j, 3) + w(i, j, 3)*rosm1d0
    vs = w(i, j, 3)*rosm1
    wsdd = wd(i, j, 4)*rosm1d0 + rosm1*wdd(i, j, 4) + rosm1d*wd0(i, j, 4&
&     ) + w(i, j, 4)*rosm1dd
    wsd = rosm1*wd(i, j, 4) + w(i, j, 4)*rosm1d
    wsd0 = rosm1*wd0(i, j, 4) + w(i, j, 4)*rosm1d0
    ws = w(i, j, 4)*rosm1
    tempd = 2*us*usd0 + 2*vs*vsd0 + 2*ws*wsd0
    temp = us*us + vs*vs + ws*ws
    temp0 = 2*us*usd + 2*vs*vsd + 2*ws*wsd
    psdd = gam1*(wdd(i, j, 5)-half*(rosd*tempd+temp*rosdd+temp0*rosd0+&
&     ros*(usd*2*usd0+2*us*usdd+vsd*2*vsd0+2*vs*vsdd+wsd*2*wsd0+2*ws*&
&     wsdd)))
    psd = gam1*(wd(i, j, 5)-half*(temp*rosd+ros*temp0))
    psd0 = gam1*(wd0(i, j, 5)-half*(temp*rosd0+ros*tempd))
    ps = gam1*(w(i, j, 5)-half*(ros*temp))
    vnsdd = nxnorm*usdd + nynorm*vsdd
    vnsd = nxnorm*usd + nynorm*vsd
    vnsd0 = nxnorm*usd0 + nynorm*vsd0
    vns = us*nxnorm + vs*nynorm
    utsdd = usdd - nxnorm*vnsdd
    utsd = usd - nxnorm*vnsd
    utsd0 = usd0 - nxnorm*vnsd0
    uts = us - vns*nxnorm
    vtsdd = vsdd - nynorm*vnsdd
    vtsd = vsd - nynorm*vnsd
    vtsd0 = vsd0 - nynorm*vnsd0
    vts = vs - vns*nynorm
! wts   = ws - vns * nznorm
! updated State -----------------------------------------------------------
    utd0d = eps0*utsdd
    utd0 = eps0*utsd
    utd1 = eps0*utsd0
    ut = eps0*uts + (one-eps0)*utd
    vtd0d = eps0*vtsdd
    vtd0 = eps0*vtsd
    vtd1 = eps0*vtsd0
    vt = eps0*vts + (one-eps0)*vtd
    wtdd = eps0*wsdd
    wtd = eps0*wsd
    wtd0 = eps0*wsd0
    wt = eps0*ws + (one-eps0)*wdi
! wt = TWO *wd  - ws !extrap o2
    amdd = epsm*(psdd-roc0d*vnsd0-vns*roc0dd-vnsd*roc0d0-roc0*vnsdd) - (&
&     one-epsm)*vnd*roc0dd
    amd = epsm*(psd-vns*roc0d-roc0*vnsd) - (one-epsm)*vnd*roc0d
    amd0 = epsm*(psd0-vns*roc0d0-roc0*vnsd0) - (one-epsm)*vnd*roc0d0
    am = epsm*(ps-roc0*vns) + (one-epsm)*(pd-roc0*vnd)
    apdd = epsp*(psdd+roc0d*vnsd0+vns*roc0dd+vnsd*roc0d0+roc0*vnsdd) + (&
&     one-epsp)*vnd*roc0dd
    apd = epsp*(psd+vns*roc0d+roc0*vnsd) + (one-epsp)*vnd*roc0d
    apd0 = epsp*(psd0+vns*roc0d0+roc0*vnsd0) + (one-epsp)*vnd*roc0d0
    ap = epsp*(ps+roc0*vns) + (one-epsp)*(pd+roc0*vnd)
    vnd0d = half*((apd-amd)*roc0m1d0+roc0m1*(apdd-amdd)+roc0m1d*(apd0-&
&     amd0)+(ap-am)*roc0m1dd)
    vnd0 = half*(roc0m1*(apd-amd)+(ap-am)*roc0m1d)
    vnd1 = half*(roc0m1*(apd0-amd0)+(ap-am)*roc0m1d0)
    vn = (ap-am)*half*roc0m1
    pd0d = half*(apdd+amdd)
    pd0 = half*(apd+amd)
    pd1 = half*(apd0+amd0)
    p = (ap+am)*half
    tempd = ro0**2*(pd1-psd0) + (p-ps)*2*ro0*ro0d0
    temp = (p-ps)*(ro0*ro0)
    temp0 = ro0*ro0*(pd0-psd) + 2*(p-ps)*ro0*ro0d
    temp1 = 2*temp*roc0m1
    bsdd = temp0*2*roc0m1*roc0m1d0 + roc0m1**2*((pd0-psd)*2*ro0*ro0d0+&
&     ro0**2*(pd0d-psdd)+ro0*ro0d*2*(pd1-psd0)+2*(p-ps)*(ro0d*ro0d0+ro0*&
&     ro0dd)) + roc0m1d*(roc0m1*2*tempd+2*temp*roc0m1d0) + temp1*&
&     roc0m1dd + rosdd
    bsd = roc0m1*roc0m1*temp0 + temp1*roc0m1d + rosd
    bsd0 = roc0m1**2*tempd + temp*2*roc0m1*roc0m1d0 + rosd0
    bs = temp*(roc0m1*roc0m1) + ros
    tempd = ro0**2*pd1 + (p-pd)*2*ro0*ro0d0
    temp = (p-pd)*(ro0*ro0)
    temp1 = ro0*ro0*pd0 + 2*(p-pd)*ro0*ro0d
    temp0 = 2*temp*roc0m1
    b0dd = temp1*2*roc0m1*roc0m1d0 + roc0m1**2*(pd0*2*ro0*ro0d0+ro0**2*&
&     pd0d+ro0*ro0d*2*pd1+2*(p-pd)*(ro0d*ro0d0+ro0*ro0dd)) + roc0m1d*(&
&     roc0m1*2*tempd+2*temp*roc0m1d0) + temp0*roc0m1dd
    b0d = roc0m1*roc0m1*temp1 + temp0*roc0m1d
    b0d0 = roc0m1**2*tempd + temp*2*roc0m1*roc0m1d0
    b0 = rod + temp*(roc0m1*roc0m1)
    rod0d = eps0*bsdd + (one-eps0)*b0dd
    rod0 = eps0*bsd + (one-eps0)*b0d
    rod1 = eps0*bsd0 + (one-eps0)*b0d0
    ro = eps0*bs + (one-eps0)*b0
    temp1 = rod0/(ro*ro)
    rom1dd = -(one*(rod0d-temp1*2*ro*rod1)/ro**2)
    rom1d = -(one*temp1)
    rom1d0 = -(one*rod1/ro**2)
    rom1 = one/ro
    roedd = pd0d/gam1
    roed = pd0/gam1
    roed0 = pd1/gam1
    roe = p/gam1
! rou = ro  * (TWO*(ut + vn*nxnorm) - us)
! rov = ro  * (TWO*(vt + vn*nynorm) - vs)
    roudd = rod0*(utd1+nxnorm*vnd1) + (ut+nxnorm*vn)*rod0d + (utd0+&
&     nxnorm*vnd0)*rod1 + ro*(utd0d+nxnorm*vnd0d)
    roud = (ut+nxnorm*vn)*rod0 + ro*(utd0+nxnorm*vnd0)
    roud0 = (ut+nxnorm*vn)*rod1 + ro*(utd1+nxnorm*vnd1)
    rou = ro*(ut+vn*nxnorm)
    rovdd = rod0*(vtd1+nynorm*vnd1) + (vt+nynorm*vn)*rod0d + (vtd0+&
&     nynorm*vnd0)*rod1 + ro*(vtd0d+nynorm*vnd0d)
    rovd = (vt+nynorm*vn)*rod0 + ro*(vtd0+nynorm*vnd0)
    rovd0 = (vt+nynorm*vn)*rod1 + ro*(vtd1+nynorm*vnd1)
    rov = ro*(vt+vn*nynorm)
    rowdd = rod0*wtd0 + wt*rod0d + wtd*rod1 + ro*wtdd
    rowd = wt*rod0 + ro*wtd
    rowd0 = wt*rod1 + ro*wtd0
    row = ro*wt
    tempd = 2*rou*roud0 + 2*rov*rovd0 + 2*row*rowd0
    temp = rou*rou + rov*rov + row*row
    temp1 = 2*rou*roud + 2*rov*rovd + 2*row*rowd
    roe2dd = roedd + half*(rom1d*tempd+temp*rom1dd+temp1*rom1d0+rom1*(&
&     roud*2*roud0+2*rou*roudd+rovd*2*rovd0+2*rov*rovdd+rowd*2*rowd0+2*&
&     row*rowdd))
    roe2d = roed + half*(temp*rom1d+rom1*temp1)
    roe2d0 = roed0 + half*(temp*rom1d0+rom1*tempd)
    roe2 = roe + half*(rom1*temp)
    DO de=1,gh
      da = de - 1
      wdd(i-de*i0, j-de*j0, 1) = rod0d
      wd(i-de*i0, j-de*j0, 1) = rod0
      wd0(i-de*i0, j-de*j0, 1) = rod1
      w(i-de*i0, j-de*j0, 1) = ro
      wdd(i-de*i0, j-de*j0, 2) = roudd
      wd(i-de*i0, j-de*j0, 2) = roud
      wd0(i-de*i0, j-de*j0, 2) = roud0
      w(i-de*i0, j-de*j0, 2) = rou
      wdd(i-de*i0, j-de*j0, 3) = rovdd
      wd(i-de*i0, j-de*j0, 3) = rovd
      wd0(i-de*i0, j-de*j0, 3) = rovd0
      w(i-de*i0, j-de*j0, 3) = rov
      wdd(i-de*i0, j-de*j0, 4) = rowdd
      wd(i-de*i0, j-de*j0, 4) = rowd
      wd0(i-de*i0, j-de*j0, 4) = rowd0
      w(i-de*i0, j-de*j0, 4) = row
! w(i-de*i0,j-de*j0,5) = roe + HALF * rom1 * (rou*rou+rov*rov+row*row)
      wdd(i-de*i0, j-de*j0, 5) = roe2dd
      wd(i-de*i0, j-de*j0, 5) = roe2d
      wd0(i-de*i0, j-de*j0, 5) = roe2d0
      w(i-de*i0, j-de*j0, 5) = roe2
!extrap o2
! roe2 =  roe
      rod0d = two*rod0d - wdd(i-da*i0, j-da*j0, 1)
      rod0 = two*rod0 - wd(i-da*i0, j-da*j0, 1)
      rod1 = two*rod1 - wd0(i-da*i0, j-da*j0, 1)
      ro = two*ro - w(i-da*i0, j-da*j0, 1)
      roudd = two*roudd - wdd(i-da*i0, j-da*j0, 2)
      roud = two*roud - wd(i-da*i0, j-da*j0, 2)
      roud0 = two*roud0 - wd0(i-da*i0, j-da*j0, 2)
      rou = two*rou - w(i-da*i0, j-da*j0, 2)
      rovdd = two*rovdd - wdd(i-da*i0, j-da*j0, 3)
      rovd = two*rovd - wd(i-da*i0, j-da*j0, 3)
      rovd0 = two*rovd0 - wd0(i-da*i0, j-da*j0, 3)
      rov = two*rov - w(i-da*i0, j-da*j0, 3)
      rowdd = two*rowdd - wdd(i-da*i0, j-da*j0, 4)
      rowd = two*rowd - wd(i-da*i0, j-da*j0, 4)
      rowd0 = two*rowd0 - wd0(i-da*i0, j-da*j0, 4)
      row = two*row - w(i-da*i0, j-da*j0, 4)
! roe  =  TWO * roe - roe1
! roe1 =  roe2
! rom1 =  ONE/ro
      roe2dd = two*roe2dd - wdd(i-da*i0, j-da*j0, 5)
      roe2d = two*roe2d - wd(i-da*i0, j-da*j0, 5)
      roe2d0 = two*roe2d0 - wd0(i-da*i0, j-da*j0, 5)
      roe2 = two*roe2 - w(i-da*i0, j-da*j0, 5)
    END DO
  END DO
END SUBROUTINE BC_NO_REFLEXION_2D_D_D
!

