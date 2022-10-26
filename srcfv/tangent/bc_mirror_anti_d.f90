! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_mirror_anti_2d in forward (tangent) mode (with options with!SliceDeadControl with!SliceDeadInstrs with!S
!taticTaping):
!   variations   of useful results: w
!   with respect to varying inputs: w
!   RW status of diff variables: w:in-out
! =============================================================================
!                             BC mirror antisym 2D
! =============================================================================
SUBROUTINE BC_MIRROR_ANTI_2D_D(w, wd, loc, interf, nx, ny, gh, im, jm)
  IMPLICIT NONE
! write(1200,*) "i, j, w, field1", i-de*i0, j-de*j0, field(l,de,1)
! write(1200,*) "i, j, w, field2", i-de*i0, j-de*j0, field(l,de,2)
! write(1200,*) "i, j, w, field3", i-de*i0, j-de*j0, field(l,de,3)
! write(1200,*) "i, j, w, field4", i-de*i0, j-de*j0, field(l,de,4)
! write(1200,*) "i, j, w, field5", i-de*i0, j-de*j0, field(l,de,5)
!
! Variables for dimension ---------------------------------------------------
  INTEGER, INTENT(IN) :: im, jm
  INTEGER, INTENT(IN) :: gh
! Input variables -----------------------------------------------------------
  CHARACTER(len=3), INTENT(IN) :: loc
  INTEGER, DIMENSION(2, 2), INTENT(IN) :: interf
  REAL*8, DIMENSION(1-gh:im+gh+1, 1-gh:jm+gh+1, 2), INTENT(IN) :: nx
  REAL*8, DIMENSION(1-gh:im+gh+1, 1-gh:jm+gh+1, 2), INTENT(IN) :: ny
! Output variables ----------------------------------------------------------
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: w
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, 5), INTENT(INOUT) :: wd
  INTEGER :: da, i1, j1
  REAL*8 :: zero, one, vi, vj, nsum, nsumi, nxnorm, nynorm, nxloc, nyloc&
& , sens
  REAL*8 :: vjd
! Local variables -----------------------------------------------------------
  INTEGER :: kdir, de, i, j, imin, imax, jmin, jmax, lmin, lmax, l, i0, &
& j0, high
  INTRINSIC FLOAT
  INTRINSIC SQRT
  REAL*8 :: arg1
  REAL*8 :: arg1d
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
  zero = 0.d0
  one = 1.d0
  sens = FLOAT(i0 + j0)
!
  i1 = i0*i0
  j1 = j0*j0
  DO l=lmin,lmax
    i = imin + (l-lmin)*j0*j0
    j = jmin + (l-lmin)*i0*i0
! write(1200,*) ' '
! write(1200,*) ' =====================', j
    nxloc = nx(i+high*i1, j+high*j1, kdir)
    nyloc = ny(i+high*i1, j+high*j1, kdir)
! highis used to select the last normal when the face is of type *hi
    arg1 = nxloc*nxloc + nyloc*nyloc
    nsum = SQRT(arg1)
    nsumi = one/nsum
    nxnorm = nxloc*nsumi*sens
    nynorm = nyloc*nsumi*sens
    DO de=1,gh
      da = de - 1
      vi = zero
      arg1d = nynorm**2*2*w(i+da*i0, j+da*j0, 2)*wd(i+da*i0, j+da*j0, 2)&
&       + nxnorm**2*2*w(i+da*i0, j+da*j0, 3)*wd(i+da*i0, j+da*j0, 3)
      arg1 = w(i+da*i0, j+da*j0, 2)*w(i+da*i0, j+da*j0, 2)*nynorm*nynorm&
&       + w(i+da*i0, j+da*j0, 3)*w(i+da*i0, j+da*j0, 3)*nxnorm*nxnorm
      temp = SQRT(arg1)
      IF (arg1 .EQ. 0.0) THEN
        vjd = 0.0_8
      ELSE
        vjd = arg1d/(2.0*temp)
      END IF
      vj = temp
      wd(i-de*i0, j-de*j0, 1) = -wd(i+da*i0, j+da*j0, 1)
      w(i-de*i0, j-de*j0, 1) = -w(i+da*i0, j+da*j0, 1)
! w(i-de*i0,j-de*j0,2) = w(i+da*i0,j+da*j0,2)
      wd(i-de*i0, j-de*j0, 2) = nynorm*vjd
      w(i-de*i0, j-de*j0, 2) = -(nxnorm*vi-nynorm*vj)
! w(i-de*i0,j-de*j0,3) = w(i+da*i0,j+da*j0,3)
! w(i-de*i0,j-de*j0,3) = ZERO
      wd(i-de*i0, j-de*j0, 3) = nxnorm*vjd
      w(i-de*i0, j-de*j0, 3) = nynorm*vi + nxnorm*vj
      wd(i-de*i0, j-de*j0, 4) = -wd(i+da*i0, j+da*j0, 4)
      w(i-de*i0, j-de*j0, 4) = -w(i+da*i0, j+da*j0, 4)
      wd(i-de*i0, j-de*j0, 5) = -wd(i+da*i0, j+da*j0, 5)
      w(i-de*i0, j-de*j0, 5) = -w(i+da*i0, j+da*j0, 5)
    END DO
  END DO
END SUBROUTINE BC_MIRROR_ANTI_2D_D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

