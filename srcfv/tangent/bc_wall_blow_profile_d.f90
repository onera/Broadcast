! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_wall_blow_profile_2d in forward (tangent) mode (with options with!SliceDeadControl with!SliceDeadInstrs 
!with!StaticTaping):
!   variations   of useful results: w
!   with respect to varying inputs: w velprof gam
!   RW status of diff variables: w:in-out velprof:in gam:in
SUBROUTINE BC_WALL_BLOW_PROFILE_2D_D(w, wd, velprof, velprofd, loc, gam&
& , gamd, interf, gh, im, jm, lm)
  IMPLICIT NONE
!
!
! Variable for dimension ------------------------------------------
  INTEGER, INTENT(IN) :: im, jm, lm
  INTEGER, INTENT(IN) :: gh
! Input variables -------------------------------------------------
  CHARACTER(len=3), INTENT(IN) :: loc
  INTEGER, DIMENSION(2, 2), INTENT(IN) :: interf
  REAL*8, INTENT(IN) :: gam
  REAL*8, INTENT(IN) :: gamd
  REAL*8, DIMENSION(lm), INTENT(IN) :: velprof
  REAL*8, DIMENSION(lm), INTENT(IN) :: velprofd
! Returned variables ----------------------------------------------
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: w
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wd
! Local variables -------------------------------------------------
  INTEGER :: da
  REAL*8 :: pe, roe, ue, ve, we
  REAL*8 :: ped, roed, ued, ved, wed
  REAL*8 :: pi, roi, ui, vi, wi, ei, roiei, pw
  REAL*8 :: pid, roid, uid, vid, wid, roieid, pwd
  REAL*8 :: roem1, roe1m1
  REAL*8 :: roem1d, roe1m1d
  REAL*8 :: pe1, roe1, ue1, ve1, ve21, we1
  REAL*8 :: pe1d, roe1d, ue1d, ve1d, ve21d, we1d
  REAL*8 :: ve2, gami, gam1, ct0, ct1, third, four, half, one, two
  REAL*8 :: ve2d, gamid, gam1d
! -----------------------------------------------------------------
! Local variables -----------------------------------------------------------
  INTEGER :: kdir, de, i, j, imin, imax, jmin, jmax, lmin, lmax, l, i0, &
& j0, high
  REAL*8 :: pwr1
  REAL*8 :: pwr1d
  REAL*8 :: pwx2
  REAL*8 :: pwx2d
  REAL*8 :: temp
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
  lmin = 1
  IF (loc .EQ. 'Ilo') THEN
    i0 = 1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Ihi') THEN
    i0 = -1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Jlo') THEN
    j0 = 1
    lmax = imax - imin + 1
  ELSE IF (loc .EQ. 'Jhi') THEN
    j0 = -1
    lmax = imax - imin + 1
  END IF
!bid
  gam1d = gamd
  gam1 = gam - 1.d0
  gamid = -(gamd/gam**2)
  gami = 1.d0/gam
! 9/8
  ct0 = 1.125d0
!-1/8
  ct1 = -0.125d0
  third = 1.d0/3.d0
  four = 4.d0
  half = 0.5d0
  one = 1.d0
  two = 2.d0
!bid
!
  DO l=lmin,lmax
    i = imin + (l-lmin)*j0*j0
    j = jmin + (l-lmin)*i0*i0
!
    roed = wd(i, j, 1)
    roe = w(i, j, 1)
    roem1d = -(one*roed/roe**2)
    roem1 = one/roe
    ued = roem1*wd(i, j, 2) + w(i, j, 2)*roem1d
    ue = w(i, j, 2)*roem1
    ved = roem1*wd(i, j, 3) + w(i, j, 3)*roem1d
    ve = w(i, j, 3)*roem1
    wed = roem1*wd(i, j, 4) + w(i, j, 4)*roem1d
    we = w(i, j, 4)*roem1
    ve2d = 2*ue*ued + 2*ve*ved + 2*we*wed
    ve2 = ue*ue + ve*ve + we*we
! pe    = gam1 * (w(i+i0,j+j0,5) - HALF*roe*ve2)
    temp = w(i, j, 5) - half*roe*ve2
    ped = temp*gam1d + gam1*(wd(i, j, 5)-half*(ve2*roed+roe*ve2d))
    pe = gam1*temp
!
    roe1d = wd(i+i0, j+j0, 1)
    roe1 = w(i+i0, j+j0, 1)
    roe1m1d = -(one*roe1d/roe1**2)
    roe1m1 = one/roe1
    temp = w(i+i0, j+j0, 2)
    ue1d = roe1m1*wd(i+i0, j+j0, 2) + temp*roe1m1d
    ue1 = temp*roe1m1
    temp = w(i+i0, j+j0, 3)
    ve1d = roe1m1*wd(i+i0, j+j0, 3) + temp*roe1m1d
    ve1 = temp*roe1m1
    temp = w(i+i0, j+j0, 4)
    we1d = roe1m1*wd(i+i0, j+j0, 4) + temp*roe1m1d
    we1 = temp*roe1m1
    ve21d = 2*ue1*ue1d + 2*ve1*ve1d + 2*we1*we1d
    ve21 = ue1*ue1 + ve1*ve1 + we1*we1
    temp = w(i+i0, j+j0, 5) - half*roe1*ve21
    pe1d = temp*gam1d + gam1*(wd(i+i0, j+j0, 5)-half*(ve21*roe1d+roe1*&
&     ve21d))
    pe1 = gam1*temp
!
! extrap dt/dn = 0 and dp/dn = 0 o2
    pwd = ct0*ped + ct1*pe1d
    pw = ct0*pe + ct1*pe1
!
! ???
    pid = third*(four*pwd-ped)
    pi = third*(four*pw-pe)
! pi  = pe
! roi = (pe*roe**gam / pi)**gami
    temp = roe**gam
    IF (roe .LE. 0.0 .AND. (gam .EQ. 0.0 .OR. gam .NE. INT(gam))) THEN
      pwr1d = 0.0_8
    ELSE IF (roe .LE. 0.0) THEN
      pwr1d = gam*roe**(gam-1)*roed
    ELSE
      pwr1d = gam*roe**(gam-1)*roed + temp*LOG(roe)*gamd
    END IF
    pwr1 = temp
    temp = pi*pwr1/pe
    pwx2d = (pwr1*pid+pi*pwr1d-temp*ped)/pe
    pwx2 = temp
    temp = pwx2**gami
    IF (pwx2 .LE. 0.0 .AND. (gami .EQ. 0.0 .OR. gami .NE. INT(gami))) &
&   THEN
      roid = 0.0_8
    ELSE IF (pwx2 .LE. 0.0) THEN
      roid = gami*pwx2**(gami-1)*pwx2d
    ELSE
      roid = gami*pwx2**(gami-1)*pwx2d + temp*LOG(pwx2)*gamid
    END IF
    roi = temp
!
    roieid = (pid-pi*gam1d/gam1)/gam1
    roiei = pi/gam1
!
    uid = -ued
    ui = -ue
    vid = two*velprofd(l) - ved
    vi = two*velprof(l) - ve
    wid = -wed
    wi = -we
!
    DO de=1,gh
      wd(i-de*i0, j-de*j0, 1) = roid
      w(i-de*i0, j-de*j0, 1) = roi
      wd(i-de*i0, j-de*j0, 2) = ui*roid + roi*uid
      w(i-de*i0, j-de*j0, 2) = roi*ui
      wd(i-de*i0, j-de*j0, 3) = vi*roid + roi*vid
      w(i-de*i0, j-de*j0, 3) = roi*vi
      wd(i-de*i0, j-de*j0, 4) = wi*roid + roi*wid
      w(i-de*i0, j-de*j0, 4) = roi*wi
      temp = ui*ui + vi*vi + wi*wi
      wd(i-de*i0, j-de*j0, 5) = roieid + half*(temp*roid+roi*(2*ui*uid+2&
&       *vi*vid+2*wi*wid))
      w(i-de*i0, j-de*j0, 5) = roiei + half*(roi*temp)
!
      temp = one/w(i+de*i0, j+de*j0, 1)
      roem1d = -(temp*wd(i+de*i0, j+de*j0, 1)/w(i+de*i0, j+de*j0, 1))
      roem1 = temp
      temp = w(i+de*i0, j+de*j0, 2)
      ue1d = roem1*wd(i+de*i0, j+de*j0, 2) + temp*roem1d
      ue1 = temp*roem1
      temp = w(i+de*i0, j+de*j0, 3)
      ve1d = roem1*wd(i+de*i0, j+de*j0, 3) + temp*roem1d
      ve1 = temp*roem1
      temp = w(i+de*i0, j+de*j0, 4)
      we1d = roem1*wd(i+de*i0, j+de*j0, 4) + temp*roem1d
      we1 = temp*roem1
      uid = -ue1d
      ui = -ue1
      vid = two*velprofd(l) - ve1d
      vi = two*velprof(l) - ve1
      wid = -we1d
      wi = -we1
    END DO
  END DO
END SUBROUTINE BC_WALL_BLOW_PROFILE_2D_D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

