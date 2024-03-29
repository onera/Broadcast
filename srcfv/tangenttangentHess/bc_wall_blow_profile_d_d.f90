!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_wall_blow_profile_2d_d in forward (tangent) mode (with options with!SliceDeadControl with!SliceDeadInstr
!s with!StaticTaping):
!   variations   of useful results: w wd
!   with respect to varying inputs: w velprof
!   RW status of diff variables: w:in-out velprof:in wd:out
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_wall_blow_profile_2d in forward (tangent) mode (with options with!SliceDeadControl with!SliceDeadInstrs 
!with!StaticTaping):
!   variations   of useful results: w
!   with respect to varying inputs: w velprof gam
!   RW status of diff variables: w:in-out velprof:in gam:in
SUBROUTINE BC_WALL_BLOW_PROFILE_2D_D_D(w, wd0, wd, wdd, velprof, &
& velprofd0, velprofd, loc, gam, gamd, interf, gh, im, jm, lm)
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
  REAL*8, DIMENSION(lm), INTENT(IN) :: velprofd0
  REAL*8, DIMENSION(lm), INTENT(IN) :: velprofd
! Returned variables ----------------------------------------------
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: w
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wd0
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wd
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wdd
! Local variables -------------------------------------------------
  INTEGER :: da
  REAL*8 :: pe, roe, ue, ve, we
  REAL*8 :: ped0, roed0, ued0, ved0, wed0
  REAL*8 :: ped, roed, ued, ved, wed
  REAL*8 :: pedd, roedd, uedd, vedd, wedd
  REAL*8 :: pi, roi, ui, vi, wi, ei, roiei, pw
  REAL*8 :: pid0, roid0, uid0, vid0, wid0, roieid0, pwd0
  REAL*8 :: pid, roid, uid, vid, wid, roieid, pwd
  REAL*8 :: pidd, roidd, uidd, vidd, widd, roieidd, pwdd
  REAL*8 :: roem1, roe1m1
  REAL*8 :: roem1d0, roe1m1d0
  REAL*8 :: roem1d, roe1m1d
  REAL*8 :: roem1dd, roe1m1dd
  REAL*8 :: pe1, roe1, ue1, ve1, ve21, we1
  REAL*8 :: pe1d0, roe1d0, ue1d0, ve1d0, ve21d0, we1d0
  REAL*8 :: pe1d, roe1d, ue1d, ve1d, ve21d, we1d
  REAL*8 :: pe1dd, roe1dd, ue1dd, ve1dd, ve21dd, we1dd
  REAL*8 :: ve2, gami, gam1, ct0, ct1, third, four, half, one, two
  REAL*8 :: ve2d0
  REAL*8 :: ve2d, gamid, gam1d
  REAL*8 :: ve2dd
! -----------------------------------------------------------------
! Local variables -----------------------------------------------------------
  INTEGER :: kdir, de, i, j, imin, imax, jmin, jmax, lmin, lmax, l, i0, &
& j0, high
  REAL*8 :: pwr1
  REAL*8 :: pwr1d0
  REAL*8 :: pwr1d
  REAL*8 :: pwr1dd
  REAL*8 :: pwx2
  REAL*8 :: pwx2d0
  REAL*8 :: pwx2d
  REAL*8 :: pwx2dd
  REAL*8 :: temp
  REAL*8 :: tempd
  INTRINSIC INT
  INTRINSIC LOG
  REAL*8 :: pwy1
  REAL*8 :: pwr10
  REAL*8 :: pwr10d
  REAL*8 :: temp0
  REAL*8 :: temp1
  REAL*8 :: temp2
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
  wdd = 0.0_8
!bid
!
  DO l=lmin,lmax
    i = imin + (l-lmin)*j0*j0
    j = jmin + (l-lmin)*i0*i0
!
    roedd = wdd(i, j, 1)
    roed = wd(i, j, 1)
    roed0 = wd0(i, j, 1)
    roe = w(i, j, 1)
    temp0 = roed/(roe*roe)
    roem1dd = -(one*(roedd-temp0*2*roe*roed0)/roe**2)
    roem1d = -(one*temp0)
    roem1d0 = -(one*roed0/roe**2)
    roem1 = one/roe
    uedd = wd(i, j, 2)*roem1d0 + roem1*wdd(i, j, 2) + roem1d*wd0(i, j, 2&
&     ) + w(i, j, 2)*roem1dd
    ued = roem1*wd(i, j, 2) + w(i, j, 2)*roem1d
    ued0 = roem1*wd0(i, j, 2) + w(i, j, 2)*roem1d0
    ue = w(i, j, 2)*roem1
    vedd = wd(i, j, 3)*roem1d0 + roem1*wdd(i, j, 3) + roem1d*wd0(i, j, 3&
&     ) + w(i, j, 3)*roem1dd
    ved = roem1*wd(i, j, 3) + w(i, j, 3)*roem1d
    ved0 = roem1*wd0(i, j, 3) + w(i, j, 3)*roem1d0
    ve = w(i, j, 3)*roem1
    wedd = wd(i, j, 4)*roem1d0 + roem1*wdd(i, j, 4) + roem1d*wd0(i, j, 4&
&     ) + w(i, j, 4)*roem1dd
    wed = roem1*wd(i, j, 4) + w(i, j, 4)*roem1d
    wed0 = roem1*wd0(i, j, 4) + w(i, j, 4)*roem1d0
    we = w(i, j, 4)*roem1
    ve2dd = ued*2*ued0 + 2*ue*uedd + ved*2*ved0 + 2*ve*vedd + wed*2*wed0&
&     + 2*we*wedd
    ve2d = 2*ue*ued + 2*ve*ved + 2*we*wed
    ve2d0 = 2*ue*ued0 + 2*ve*ved0 + 2*we*wed0
    ve2 = ue*ue + ve*ve + we*we
! pe    = gam1 * (w(i+i0,j+j0,5) - HALF*roe*ve2)
    tempd = wd0(i, j, 5) - half*(ve2*roed0+roe*ve2d0)
    temp = w(i, j, 5) - half*roe*ve2
    pedd = gam1d*tempd + gam1*(wdd(i, j, 5)-half*(roed*ve2d0+ve2*roedd+&
&     ve2d*roed0+roe*ve2dd))
    ped = temp*gam1d + gam1*(wd(i, j, 5)-half*(ve2*roed+roe*ve2d))
    ped0 = gam1*tempd
    pe = gam1*temp
!
    roe1dd = wdd(i+i0, j+j0, 1)
    roe1d = wd(i+i0, j+j0, 1)
    roe1d0 = wd0(i+i0, j+j0, 1)
    roe1 = w(i+i0, j+j0, 1)
    temp0 = roe1d/(roe1*roe1)
    roe1m1dd = -(one*(roe1dd-temp0*2*roe1*roe1d0)/roe1**2)
    roe1m1d = -(one*temp0)
    roe1m1d0 = -(one*roe1d0/roe1**2)
    roe1m1 = one/roe1
    tempd = wd0(i+i0, j+j0, 2)
    temp = w(i+i0, j+j0, 2)
    temp0 = wd(i+i0, j+j0, 2)
    ue1dd = temp0*roe1m1d0 + roe1m1*wdd(i+i0, j+j0, 2) + roe1m1d*tempd +&
&     temp*roe1m1dd
    ue1d = roe1m1*temp0 + temp*roe1m1d
    ue1d0 = roe1m1*tempd + temp*roe1m1d0
    ue1 = temp*roe1m1
    tempd = wd0(i+i0, j+j0, 3)
    temp = w(i+i0, j+j0, 3)
    temp0 = wd(i+i0, j+j0, 3)
    ve1dd = temp0*roe1m1d0 + roe1m1*wdd(i+i0, j+j0, 3) + roe1m1d*tempd +&
&     temp*roe1m1dd
    ve1d = roe1m1*temp0 + temp*roe1m1d
    ve1d0 = roe1m1*tempd + temp*roe1m1d0
    ve1 = temp*roe1m1
    tempd = wd0(i+i0, j+j0, 4)
    temp = w(i+i0, j+j0, 4)
    temp0 = wd(i+i0, j+j0, 4)
    we1dd = temp0*roe1m1d0 + roe1m1*wdd(i+i0, j+j0, 4) + roe1m1d*tempd +&
&     temp*roe1m1dd
    we1d = roe1m1*temp0 + temp*roe1m1d
    we1d0 = roe1m1*tempd + temp*roe1m1d0
    we1 = temp*roe1m1
    ve21dd = ue1d*2*ue1d0 + 2*ue1*ue1dd + ve1d*2*ve1d0 + 2*ve1*ve1dd + &
&     we1d*2*we1d0 + 2*we1*we1dd
    ve21d = 2*ue1*ue1d + 2*ve1*ve1d + 2*we1*we1d
    ve21d0 = 2*ue1*ue1d0 + 2*ve1*ve1d0 + 2*we1*we1d0
    ve21 = ue1*ue1 + ve1*ve1 + we1*we1
    tempd = wd0(i+i0, j+j0, 5) - half*(ve21*roe1d0+roe1*ve21d0)
    temp = w(i+i0, j+j0, 5) - half*roe1*ve21
    pe1dd = gam1d*tempd + gam1*(wdd(i+i0, j+j0, 5)-half*(roe1d*ve21d0+&
&     ve21*roe1dd+ve21d*roe1d0+roe1*ve21dd))
    pe1d = temp*gam1d + gam1*(wd(i+i0, j+j0, 5)-half*(ve21*roe1d+roe1*&
&     ve21d))
    pe1d0 = gam1*tempd
    pe1 = gam1*temp
!
! extrap dt/dn = 0 and dp/dn = 0 o2
    pwdd = ct0*pedd + ct1*pe1dd
    pwd = ct0*ped + ct1*pe1d
    pwd0 = ct0*ped0 + ct1*pe1d0
    pw = ct0*pe + ct1*pe1
!
! ???
    pidd = third*(four*pwdd-pedd)
    pid = third*(four*pwd-ped)
    pid0 = third*(four*pwd0-ped0)
    pi = third*(four*pw-pe)
! pi  = pe
! roi = (pe*roe**gam / pi)**gami
    IF (roe .LE. 0.0 .AND. (gam .EQ. 0.0 .OR. gam .NE. INT(gam))) THEN
      tempd = 0.0_8
    ELSE
      tempd = gam*roe**(gam-1)*roed0
    END IF
    temp = roe**gam
    IF (roe .LE. 0.0 .AND. (gam .EQ. 0.0 .OR. gam .NE. INT(gam))) THEN
      pwr1d = 0.0_8
      pwr1dd = 0.0_8
    ELSE IF (roe .LE. 0.0) THEN
      pwy1 = gam - 1
      IF (roe .LE. 0.0 .AND. (pwy1 .EQ. 0.0 .OR. pwy1 .NE. INT(pwy1))) &
&     THEN
        pwr10d = 0.0_8
      ELSE
        pwr10d = pwy1*roe**(pwy1-1)*roed0
      END IF
      pwr10 = roe**pwy1
      pwr1dd = gam*(roed*pwr10d+pwr10*roedd)
      pwr1d = gam*pwr10*roed
    ELSE
      pwy1 = gam - 1
      IF (roe .LE. 0.0 .AND. (pwy1 .EQ. 0.0 .OR. pwy1 .NE. INT(pwy1))) &
&     THEN
        pwr10d = 0.0_8
      ELSE
        pwr10d = pwy1*roe**(pwy1-1)*roed0
      END IF
      pwr10 = roe**pwy1
      temp0 = LOG(roe)
      pwr1dd = gam*(roed*pwr10d+pwr10*roedd) + gamd*(temp0*tempd+temp*&
&       roed0/roe)
      pwr1d = gam*(pwr10*roed) + gamd*(temp*temp0)
    END IF
    pwr1d0 = tempd
    pwr1 = temp
    temp0 = pi*pwr1/pe
    tempd = (pwr1*pid0+pi*pwr1d0-temp0*ped0)/pe
    temp = temp0
    temp0 = (pwr1*pid+pi*pwr1d-temp*ped)/pe
    pwx2dd = (pid*pwr1d0+pwr1*pidd+pwr1d*pid0+pi*pwr1dd-ped*tempd-temp*&
&     pedd-temp0*ped0)/pe
    pwx2d = temp0
    pwx2d0 = tempd
    pwx2 = temp
    IF (pwx2 .LE. 0.0 .AND. (gami .EQ. 0.0 .OR. gami .NE. INT(gami))) &
&   THEN
      tempd = 0.0_8
    ELSE
      tempd = gami*pwx2**(gami-1)*pwx2d0
    END IF
    temp = pwx2**gami
    IF (pwx2 .LE. 0.0 .AND. (gami .EQ. 0.0 .OR. gami .NE. INT(gami))) &
&   THEN
      roid = 0.0_8
      roidd = 0.0_8
    ELSE IF (pwx2 .LE. 0.0) THEN
      pwy1 = gami - 1
      IF (pwx2 .LE. 0.0 .AND. (pwy1 .EQ. 0.0 .OR. pwy1 .NE. INT(pwy1))) &
&     THEN
        pwr10d = 0.0_8
      ELSE
        pwr10d = pwy1*pwx2**(pwy1-1)*pwx2d0
      END IF
      pwr10 = pwx2**pwy1
      roidd = gami*(pwx2d*pwr10d+pwr10*pwx2dd)
      roid = gami*pwr10*pwx2d
    ELSE
      pwy1 = gami - 1
      IF (pwx2 .LE. 0.0 .AND. (pwy1 .EQ. 0.0 .OR. pwy1 .NE. INT(pwy1))) &
&     THEN
        pwr10d = 0.0_8
      ELSE
        pwr10d = pwy1*pwx2**(pwy1-1)*pwx2d0
      END IF
      pwr10 = pwx2**pwy1
      temp0 = LOG(pwx2)
      roidd = gami*(pwx2d*pwr10d+pwr10*pwx2dd) + gamid*(temp0*tempd+temp&
&       *pwx2d0/pwx2)
      roid = gami*(pwr10*pwx2d) + gamid*(temp*temp0)
    END IF
    roid0 = tempd
    roi = temp
!
    roieidd = (pidd-gam1d*pid0/gam1)/gam1
    roieid = (pid-pi*gam1d/gam1)/gam1
    roieid0 = pid0/gam1
    roiei = pi/gam1
!
    uidd = -uedd
    uid = -ued
    uid0 = -ued0
    ui = -ue
    vidd = -vedd
    vid = two*velprofd(l) - ved
    vid0 = two*velprofd0(l) - ved0
    vi = two*velprof(l) - ve
    widd = -wedd
    wid = -wed
    wid0 = -wed0
    wi = -we
!
    DO de=1,gh
      wdd(i-de*i0, j-de*j0, 1) = roidd
      wd(i-de*i0, j-de*j0, 1) = roid
      wd0(i-de*i0, j-de*j0, 1) = roid0
      w(i-de*i0, j-de*j0, 1) = roi
      wdd(i-de*i0, j-de*j0, 2) = roid*uid0 + ui*roidd + uid*roid0 + roi*&
&       uidd
      wd(i-de*i0, j-de*j0, 2) = ui*roid + roi*uid
      wd0(i-de*i0, j-de*j0, 2) = ui*roid0 + roi*uid0
      w(i-de*i0, j-de*j0, 2) = roi*ui
      wdd(i-de*i0, j-de*j0, 3) = roid*vid0 + vi*roidd + vid*roid0 + roi*&
&       vidd
      wd(i-de*i0, j-de*j0, 3) = vi*roid + roi*vid
      wd0(i-de*i0, j-de*j0, 3) = vi*roid0 + roi*vid0
      w(i-de*i0, j-de*j0, 3) = roi*vi
      wdd(i-de*i0, j-de*j0, 4) = roid*wid0 + wi*roidd + wid*roid0 + roi*&
&       widd
      wd(i-de*i0, j-de*j0, 4) = wi*roid + roi*wid
      wd0(i-de*i0, j-de*j0, 4) = wi*roid0 + roi*wid0
      w(i-de*i0, j-de*j0, 4) = roi*wi
      tempd = 2*ui*uid0 + 2*vi*vid0 + 2*wi*wid0
      temp = ui*ui + vi*vi + wi*wi
      temp0 = 2*ui*uid + 2*vi*vid + 2*wi*wid
      wdd(i-de*i0, j-de*j0, 5) = roieidd + half*(roid*tempd+temp*roidd+&
&       temp0*roid0+roi*(uid*2*uid0+2*ui*uidd+vid*2*vid0+2*vi*vidd+wid*2&
&       *wid0+2*wi*widd))
      wd(i-de*i0, j-de*j0, 5) = roieid + half*(temp*roid+roi*temp0)
      wd0(i-de*i0, j-de*j0, 5) = roieid0 + half*(temp*roid0+roi*tempd)
      w(i-de*i0, j-de*j0, 5) = roiei + half*(roi*temp)
!
      temp0 = one/w(i+de*i0, j+de*j0, 1)
      tempd = -(temp0*wd0(i+de*i0, j+de*j0, 1)/w(i+de*i0, j+de*j0, 1))
      temp = temp0
      temp0 = w(i+de*i0, j+de*j0, 1)
      temp1 = wd(i+de*i0, j+de*j0, 1)
      temp2 = temp*temp1/temp0
      roem1dd = -((temp1*tempd+temp*wdd(i+de*i0, j+de*j0, 1)-temp2*wd0(i&
&       +de*i0, j+de*j0, 1))/temp0)
      roem1d = -temp2
      roem1d0 = tempd
      roem1 = temp
      tempd = wd0(i+de*i0, j+de*j0, 2)
      temp = w(i+de*i0, j+de*j0, 2)
      temp2 = wd(i+de*i0, j+de*j0, 2)
      ue1dd = temp2*roem1d0 + roem1*wdd(i+de*i0, j+de*j0, 2) + roem1d*&
&       tempd + temp*roem1dd
      ue1d = roem1*temp2 + temp*roem1d
      ue1d0 = roem1*tempd + temp*roem1d0
      ue1 = temp*roem1
      tempd = wd0(i+de*i0, j+de*j0, 3)
      temp = w(i+de*i0, j+de*j0, 3)
      temp2 = wd(i+de*i0, j+de*j0, 3)
      ve1dd = temp2*roem1d0 + roem1*wdd(i+de*i0, j+de*j0, 3) + roem1d*&
&       tempd + temp*roem1dd
      ve1d = roem1*temp2 + temp*roem1d
      ve1d0 = roem1*tempd + temp*roem1d0
      ve1 = temp*roem1
      tempd = wd0(i+de*i0, j+de*j0, 4)
      temp = w(i+de*i0, j+de*j0, 4)
      temp2 = wd(i+de*i0, j+de*j0, 4)
      we1dd = temp2*roem1d0 + roem1*wdd(i+de*i0, j+de*j0, 4) + roem1d*&
&       tempd + temp*roem1dd
      we1d = roem1*temp2 + temp*roem1d
      we1d0 = roem1*tempd + temp*roem1d0
      we1 = temp*roem1
      uidd = -ue1dd
      uid = -ue1d
      uid0 = -ue1d0
      ui = -ue1
      vidd = -ve1dd
      vid = two*velprofd(l) - ve1d
      vid0 = two*velprofd0(l) - ve1d0
      vi = two*velprof(l) - ve1
      widd = -we1dd
      wid = -we1d
      wid0 = -we1d0
      wi = -we1
    END DO
  END DO
END SUBROUTINE BC_WALL_BLOW_PROFILE_2D_D_D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

