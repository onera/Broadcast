/* Copyright (C) 1991-2018 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */
! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
subroutine bordersfromcenters_2d(x0,y0,xc,yc,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1:im+1,1:jm+1 ),intent(inout) :: x0
  real(8),dimension(1:im+1,1:jm+1 ),intent(inout) :: y0
  real(8),dimension(1:im ,1:jm ),intent(in) :: xc
  real(8),dimension(1:im ,1:jm ),intent(in) :: yc
  ! Non-required arguments -------------------------------------------
  integer :: i,j,g,dummy
  real(8) :: TWO,HALF,FOURTH! ------------------------------------------------------------------
  TWO = 2.d0
  HALF = 0.5d0
  FOURTH = 0.25d0
!$AD II-LOOP
  do j=2, jm
!$AD II-LOOP
!DIR$ IVDEP
  do i=2, im
    x0(i,j) = FOURTH*(xc(i-1,j-1) + xc(i,j-1) + xc(i-1,j) + xc(i,j))
    y0(i,j) = FOURTH*(yc(i-1,j-1) + yc(i,j-1) + yc(i-1,j) + yc(i,j))
  enddo
  enddo
!$AD II-LOOP
  do j=2, jm
    x0(im+1,j) = xc(im,j) + xc(im,j-1) - x0(im,j)
    y0(im+1,j) = yc(im,j) + yc(im,j-1) - y0(im,j)
    x0(1 ,j) = xc(1 ,j) + xc(1 ,j-1) - x0(2 ,j)
    y0(1 ,j) = yc(1 ,j) + yc(1 ,j-1) - y0(2 ,j)
  enddo
!$AD II-LOOP
  do i=2, im
    x0(i,jm+1) = xc(i,jm) + xc(i-1,jm) - x0(i,jm)
    y0(i,jm+1) = yc(i,jm) + yc(i-1,jm) - y0(i,jm)
    x0(i,1 ) = xc(i,1 ) + xc(i-1,1 ) - x0(i,2 )
    y0(i,1 ) = yc(i,1 ) + yc(i-1,1 ) - y0(i,2 )
  enddo
  x0(im+1,jm+1) = x0(im,jm+1) + x0(im+1,jm ) - x0(im,jm)
  y0(im+1,jm+1) = y0(im,jm+1) + y0(im+1,jm ) - y0(im,jm)
  x0(1 ,jm+1) = x0(1 ,jm ) + x0(2 ,jm+1) - x0(2 ,jm)
  y0(1 ,jm+1) = y0(1 ,jm ) + y0(2 ,jm+1) - y0(2 ,jm)
  x0(im+1,1 ) = x0(im,1 ) + x0(im+1,2 ) - x0(im,2 )
  y0(im+1,1 ) = y0(im,1 ) + y0(im+1,2 ) - y0(im,2 )
  x0(1 ,1 ) = x0(1 ,2 ) + x0(2 ,1 ) - x0(2 ,2 )
  y0(1 ,1 ) = y0(1 ,2 ) + y0(2 ,1 ) - y0(2 ,2 )
end subroutine bordersfromcenters_2d
