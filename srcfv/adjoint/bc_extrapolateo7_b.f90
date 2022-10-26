! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_extrapolate_o7_2d in reverse (adjoint) mode (with options with!SliceDeadControl with!SliceDeadInstrs wit
!h!StaticTaping):
!   gradient     of useful results: w
!   with respect to varying inputs: w
!   RW status of diff variables: w:in-out
SUBROUTINE BC_EXTRAPOLATE_O7_2D_B(w, wb, loc, interf, im, jm, gh, em)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: im, jm, gh, em
  CHARACTER(len=3), INTENT(IN) :: loc
  INTEGER, DIMENSION(2, 2), INTENT(IN) :: interf
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, em), INTENT(INOUT) :: w
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, em), INTENT(INOUT) :: wb
! Local variables -----------------------------------------------------------
  INTEGER :: d1, d2, d3, d4, d5, d6, d7
! Local variables -----------------------------------------------------------
  INTEGER :: kdir, de, i, j, imin, imax, jmin, jmax, lmin, lmax, l, i0, &
& j0, high
  REAL*8, DIMENSION(em) :: tmp
  REAL*8, DIMENSION(em) :: tmpb
  INTEGER :: branch
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
    CALL PUSHCONTROL2B(3)
    i0 = 1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Ihi') THEN
    CALL PUSHCONTROL2B(2)
    i0 = -1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Jlo') THEN
    CALL PUSHCONTROL2B(1)
    j0 = 1
    lmax = imax - imin + 1
  ELSE IF (loc .EQ. 'Jhi') THEN
    CALL PUSHCONTROL2B(0)
    j0 = -1
    lmax = imax - imin + 1
  ELSE
    CALL PUSHCONTROL2B(0)
  END IF
  DO l=lmin,lmax
    CALL PUSHINTEGER4(i)
    i = imin + (l-lmin)*j0*j0
    CALL PUSHINTEGER4(j)
    j = jmin + (l-lmin)*i0*i0
  END DO
  DO l=lmax,lmin,-1
    DO de=gh,1,-1
      d1 = de - 1
      d2 = de - 2
      d3 = de - 3
      d4 = de - 4
      d5 = de - 5
      d6 = de - 6
      d7 = de - 7
      tmpb(:) = wb(i-i0*de, j-j0*de, :)
      wb(i-i0*de, j-j0*de, :) = 0.0_8
      wb(i-i0*d2, j-j0*d2, :) = wb(i-i0*d2, j-j0*d2, :) + 21.d0*tmpb
      wb(i-i0*d1, j-j0*d1, :) = wb(i-i0*d1, j-j0*d1, :) - 7.d0*tmpb
      wb(i-i0*d4, j-j0*d4, :) = wb(i-i0*d4, j-j0*d4, :) + 35.d0*tmpb
      wb(i-i0*d3, j-j0*d3, :) = wb(i-i0*d3, j-j0*d3, :) - 35.d0*tmpb
      wb(i-i0*d6, j-j0*d6, :) = wb(i-i0*d6, j-j0*d6, :) + 7.d0*tmpb
      wb(i-i0*d5, j-j0*d5, :) = wb(i-i0*d5, j-j0*d5, :) - 21.d0*tmpb
      wb(i-i0*d7, j-j0*d7, :) = wb(i-i0*d7, j-j0*d7, :) - tmpb
    END DO
    CALL POPINTEGER4(j)
    CALL POPINTEGER4(i)
  END DO
  CALL POPCONTROL2B(branch)
END SUBROUTINE BC_EXTRAPOLATE_O7_2D_B

