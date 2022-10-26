! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bc_general_2d in reverse (adjoint) mode (with options with!SliceDeadControl with!SliceDeadInstrs with!Stati
!cTaping):
!   gradient     of useful results: w
!   with respect to varying inputs: w
!   RW status of diff variables: w:in-out
! =============================================================================
!                             BC general 2D
! =============================================================================
SUBROUTINE BC_GENERAL_2D_B(w, wb, loc, interf, field, gh, im, jm, lm, em&
&)
  IMPLICIT NONE
!
! Variables for dimension ---------------------------------------------------
  INTEGER, INTENT(IN) :: lm, im, jm, em
  INTEGER, INTENT(IN) :: gh
! Input variables -----------------------------------------------------------
  CHARACTER(len=3), INTENT(IN) :: loc
  INTEGER, DIMENSION(2, 2), INTENT(IN) :: interf
  REAL*8, DIMENSION(lm, gh, em), INTENT(IN) :: field
! Output variables ----------------------------------------------------------
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, em), INTENT(INOUT) :: w
  REAL*8, DIMENSION(1-gh:im+gh, 1-gh:jm+gh, em), INTENT(INOUT) :: wb
! Local variables -----------------------------------------------------------
  INTEGER :: kdir, de, i, j, imin, imax, jmin, jmax, lmin, lmax, l, i0, &
& j0, high
  INTEGER :: ad_branch
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
    ad_branch = 0
    i0 = 1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Ihi') THEN
    ad_branch = 1
    i0 = -1
    lmax = jmax - jmin + 1
  ELSE IF (loc .EQ. 'Jlo') THEN
    ad_branch = 2
    j0 = 1
    lmax = imax - imin + 1
  ELSE IF (loc .EQ. 'Jhi') THEN
    ad_branch = 3
    j0 = -1
    lmax = imax - imin + 1
  ELSE
    ad_branch = 3
  END IF
!$BWD-OF II-LOOP 
  DO l=lmin,lmax
    i = imin + (l-lmin)*j0*j0
    j = jmin + (l-lmin)*i0*i0
    DO de=gh,1,-1
      wb(i-de*i0, j-de*j0, :) = 0.0_8
    END DO
  END DO
END SUBROUTINE BC_GENERAL_2D_B
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

