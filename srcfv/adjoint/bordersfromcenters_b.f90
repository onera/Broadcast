! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
!  Differentiation of bordersfromcenters_2d in reverse (adjoint) mode (with options with!SliceDeadControl with!SliceDeadInstrs wi
!th!StaticTaping):
!   gradient     of useful results: xc yc x0 y0
!   with respect to varying inputs: xc yc x0 y0
!   RW status of diff variables: xc:incr yc:incr x0:in-out y0:in-out
SUBROUTINE BORDERSFROMCENTERS_2D_B(x0, x0b, y0, y0b, xc, xcb, yc, ycb, &
& im, jm)
  IMPLICIT NONE
! variables for dimension -----------------------------------------
  INTEGER :: im, jm, gh
! required arguments ----------------------------------------------
! Returned objects ------------------------------------------------
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: x0
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: x0b
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: y0
  REAL*8, DIMENSION(im+1, jm+1), INTENT(INOUT) :: y0b
  REAL*8, DIMENSION(im, jm), INTENT(IN) :: xc
  REAL*8, DIMENSION(im, jm) :: xcb
  REAL*8, DIMENSION(im, jm), INTENT(IN) :: yc
  REAL*8, DIMENSION(im, jm) :: ycb
! Non-required arguments -------------------------------------------
  INTEGER :: i, j, g, dummy
! ------------------------------------------------------------------
  REAL*8 :: two, half, fourth
  REAL*8 :: tempb
  fourth = 0.25d0
  y0b(1, 2) = y0b(1, 2) + y0b(1, 1)
  y0b(2, 1) = y0b(2, 1) + y0b(1, 1)
  y0b(2, 2) = y0b(2, 2) - y0b(1, 1)
  y0b(1, 1) = 0.0_8
  x0b(1, 2) = x0b(1, 2) + x0b(1, 1)
  x0b(2, 1) = x0b(2, 1) + x0b(1, 1)
  x0b(2, 2) = x0b(2, 2) - x0b(1, 1)
  x0b(1, 1) = 0.0_8
  y0b(im, 1) = y0b(im, 1) + y0b(im+1, 1)
  y0b(im+1, 2) = y0b(im+1, 2) + y0b(im+1, 1)
  y0b(im, 2) = y0b(im, 2) - y0b(im+1, 1)
  y0b(im+1, 1) = 0.0_8
  x0b(im, 1) = x0b(im, 1) + x0b(im+1, 1)
  x0b(im+1, 2) = x0b(im+1, 2) + x0b(im+1, 1)
  x0b(im, 2) = x0b(im, 2) - x0b(im+1, 1)
  x0b(im+1, 1) = 0.0_8
  y0b(1, jm) = y0b(1, jm) + y0b(1, jm+1)
  y0b(2, jm+1) = y0b(2, jm+1) + y0b(1, jm+1)
  y0b(2, jm) = y0b(2, jm) - y0b(1, jm+1)
  y0b(1, jm+1) = 0.0_8
  x0b(1, jm) = x0b(1, jm) + x0b(1, jm+1)
  x0b(2, jm+1) = x0b(2, jm+1) + x0b(1, jm+1)
  x0b(2, jm) = x0b(2, jm) - x0b(1, jm+1)
  x0b(1, jm+1) = 0.0_8
  y0b(im, jm+1) = y0b(im, jm+1) + y0b(im+1, jm+1)
  y0b(im+1, jm) = y0b(im+1, jm) + y0b(im+1, jm+1)
  y0b(im, jm) = y0b(im, jm) - y0b(im+1, jm+1)
  y0b(im+1, jm+1) = 0.0_8
  x0b(im, jm+1) = x0b(im, jm+1) + x0b(im+1, jm+1)
  x0b(im+1, jm) = x0b(im+1, jm) + x0b(im+1, jm+1)
  x0b(im, jm) = x0b(im, jm) - x0b(im+1, jm+1)
  x0b(im+1, jm+1) = 0.0_8
!$BWD-OF II-LOOP 
  DO i=2,im
    ycb(i, 1) = ycb(i, 1) + y0b(i, 1)
    ycb(i-1, 1) = ycb(i-1, 1) + y0b(i, 1)
    y0b(i, 2) = y0b(i, 2) - y0b(i, 1)
    y0b(i, 1) = 0.0_8
    xcb(i, 1) = xcb(i, 1) + x0b(i, 1)
    xcb(i-1, 1) = xcb(i-1, 1) + x0b(i, 1)
    x0b(i, 2) = x0b(i, 2) - x0b(i, 1)
    x0b(i, 1) = 0.0_8
    ycb(i, jm) = ycb(i, jm) + y0b(i, jm+1)
    ycb(i-1, jm) = ycb(i-1, jm) + y0b(i, jm+1)
    y0b(i, jm) = y0b(i, jm) - y0b(i, jm+1)
    y0b(i, jm+1) = 0.0_8
    xcb(i, jm) = xcb(i, jm) + x0b(i, jm+1)
    xcb(i-1, jm) = xcb(i-1, jm) + x0b(i, jm+1)
    x0b(i, jm) = x0b(i, jm) - x0b(i, jm+1)
    x0b(i, jm+1) = 0.0_8
  END DO
!$BWD-OF II-LOOP 
  DO j=2,jm
    ycb(1, j) = ycb(1, j) + y0b(1, j)
    ycb(1, j-1) = ycb(1, j-1) + y0b(1, j)
    y0b(2, j) = y0b(2, j) - y0b(1, j)
    y0b(1, j) = 0.0_8
    xcb(1, j) = xcb(1, j) + x0b(1, j)
    xcb(1, j-1) = xcb(1, j-1) + x0b(1, j)
    x0b(2, j) = x0b(2, j) - x0b(1, j)
    x0b(1, j) = 0.0_8
    ycb(im, j) = ycb(im, j) + y0b(im+1, j)
    ycb(im, j-1) = ycb(im, j-1) + y0b(im+1, j)
    y0b(im, j) = y0b(im, j) - y0b(im+1, j)
    y0b(im+1, j) = 0.0_8
    xcb(im, j) = xcb(im, j) + x0b(im+1, j)
    xcb(im, j-1) = xcb(im, j-1) + x0b(im+1, j)
    x0b(im, j) = x0b(im, j) - x0b(im+1, j)
    x0b(im+1, j) = 0.0_8
  END DO
!$BWD-OF II-LOOP 
  DO j=2,jm
!$BWD-OF II-LOOP 
!DIR$ IVDEP
    DO i=2,im
      tempb = fourth*y0b(i, j)
      y0b(i, j) = 0.0_8
      ycb(i-1, j-1) = ycb(i-1, j-1) + tempb
      ycb(i, j-1) = ycb(i, j-1) + tempb
      ycb(i-1, j) = ycb(i-1, j) + tempb
      ycb(i, j) = ycb(i, j) + tempb
      tempb = fourth*x0b(i, j)
      x0b(i, j) = 0.0_8
      xcb(i-1, j-1) = xcb(i-1, j-1) + tempb
      xcb(i, j-1) = xcb(i, j-1) + tempb
      xcb(i-1, j) = xcb(i-1, j) + tempb
      xcb(i, j) = xcb(i, j) + tempb
    END DO
  END DO
END SUBROUTINE BORDERSFROMCENTERS_2D_B

