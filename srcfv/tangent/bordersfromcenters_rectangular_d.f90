! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bordersfromcenters_rectangular_2d in forward (tangent) mode (with options with!SliceDeadControl with!SliceD
!eadInstrs with!StaticTaping):
!   variations   of useful results: x0 y0
!   with respect to varying inputs: xc ymin xmin yc x0 y0
!   RW status of diff variables: xc:in ymin:in xmin:in yc:in x0:in-out
!                y0:in-out
SUBROUTINE BORDERSFROMCENTERS_RECTANGULAR_2D_D(x0, x0d, y0, y0d, xc, xcd&
& , yc, ycd, im, jm, xmin, xmind, ymin, ymind)
  IMPLICIT NONE
! variables for dimension -----------------------------------------
  INTEGER :: im, jm, gh
! required arguments ----------------------------------------------
! Returned objects ------------------------------------------------
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: x0
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: x0d
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: y0
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: y0d
  REAL*8, DIMENSION(im, jm), INTENT(IN) :: xc
  REAL*8, DIMENSION(im, jm), INTENT(IN) :: xcd
  REAL*8, DIMENSION(im, jm), INTENT(IN) :: yc
  REAL*8, DIMENSION(im, jm), INTENT(IN) :: ycd
  REAL*8, INTENT(IN) :: xmin, ymin
  REAL*8, INTENT(IN) :: xmind, ymind
! Non-required arguments -------------------------------------------
  INTEGER :: i, j, g, dummy
! ------------------------------------------------------------------
  REAL*8 :: two, half, fourth
  two = 2.d0
  DO j=1,jm+1
    x0d(1, j) = xmind
  END DO
  DO i=1,im+1
    y0d(i, 1) = ymind
  END DO
  DO j=1,jm
    DO i=2,im+1
      x0d(i, j) = two*xcd(i-1, j) - x0d(i-1, j)
    END DO
  END DO
  DO j=2,jm+1
    DO i=1,im
      y0d(i, j) = two*ycd(i, j-1) - y0d(i, j-1)
    END DO
  END DO
  DO i=1,im+1
    x0d(i, jm+1) = x0d(i, jm)
  END DO
  DO j=1,jm+1
    y0d(im+1, j) = y0d(im, j)
  END DO
END SUBROUTINE BORDERSFROMCENTERS_RECTANGULAR_2D_D

