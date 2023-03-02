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
subroutine computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1 ),intent(inout) :: x0
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1 ),intent(inout) :: y0
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(inout) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(inout) :: ny
  real(8),dimension(1-gh:im+gh ,1-gh:jm+gh ),intent(inout) :: xc
  real(8),dimension(1-gh:im+gh ,1-gh:jm+gh ),intent(inout) :: yc
  real(8),dimension(1-gh:im+gh ,1-gh:jm+gh ),intent(inout) :: vol
  real(8),dimension(1-gh:im+gh ,1-gh:jm+gh ,2),intent(inout) :: volf
  ! Non-required arguments -------------------------------------------
  integer :: i,j,g,dummy
  real(8) :: TWO,HALF,FOURTH
  real(8) :: xa,xb,xc1,xd,ya,yb,yc1,yd
  real(8) :: abx,aby,acx,acy,q1
  real(8) :: dbx,dby,dcx,dcy,q2
  ! ------------------------------------------------------------------
  TWO = 2.d0
  HALF = 0.5d0
  FOURTH = 0.25d0
  !extend in ghost
  do g=1,gh
    do dummy = 1,2
      x0( 1-g,:) = TWO * x0( 2-g,:) - x0( 3-g,:)
      x0(im+1+g,:) = TWO * x0(im +g,:) - x0(im-1+g,:)
      x0(: , 1-g) = TWO * x0(: , 2-g) - x0(: , 3-g)
      x0(: ,jm+1+g) = TWO * x0(: ,jm+g) - x0(: ,jm-1+g)
      y0(:, 1-g) = TWO * y0(:, 2-g) - y0(:, 3-g)
      y0(:,jm+1+g) = TWO * y0(:,jm +g) - y0(:,jm-1+g)
      y0( 1-g,:) = TWO * y0( 2-g, :) - y0( 3-g, :)
      y0(im+1+g,:) = TWO * y0(im+g, :) - y0(im-1+g, :)
    enddo
  enddo
!$AD II-LOOP
  do j=1, jm+1
!$AD II-LOOP
!DIR$ IVDEP
  do i=1, im+1
      xc(i,j) = FOURTH*(x0(i,j) + x0(i+1,j) + x0(i,j+1) + x0(i+1,j+1))
      yc(i,j) = FOURTH*(y0(i,j) + y0(i+1,j) + y0(i,j+1) + y0(i+1,j+1))
      xa = x0(i ,j ) ; ya = y0(i ,j )
      xb = x0(i+1,j ) ; yb = y0(i+1,j )
      xc1 = x0(i ,j+1) ; yc1 = y0(i ,j+1)
      xd = x0(i+1,j+1) ; yd = y0(i+1,j+1)
      !
      abx = xb -xa ; aby = yb -ya
      acx = xc1-xa ; acy = yc1-ya
      q1 = HALF*abs( abx*acy - acx*aby )
      !
      dcx = xc1-xd ; dcy = yc1-yd
      dbx = xb -xd ; dby = yb -yd
      q2 = HALF*abs( dcx*dby - dbx*dcy )
      !
      vol(i,j) = q1 + q2
      nx(i,j,1) = y0(i,j+1) - y0(i,j)
      ny(i,j,1) = x0(i,j ) - x0(i,j+1)
      nx(i,j,2) = y0(i ,j) - y0(i+1,j)
      ny(i,j,2) = x0(i+1,j) - x0(i ,j)
  enddo
  enddo
  !extend in ghost
  do g=1,gh
    do dummy = 1,2
      xc( 1-g,:) = TWO * xc( 2-g, :) - xc( 3-g, :)
      xc(im+g,:) = TWO * xc(im-1+g, :) - xc(im-2+g, :)
      xc(: , 1-g) = TWO * xc(: , 2-g) - xc(: , 3-g)
      xc(: ,jm+g) = TWO * xc(: ,jm-1+g) - xc(: ,jm-2+g)
      yc( 1-g,:) = TWO * yc( 2-g, :) - yc( 3-g, :)
      yc(im+g,:) = TWO * yc(im-1+g, :) - yc(im-2+g, :)
      yc(: , 1-g) = TWO * yc(: , 2-g) - yc(: , 3-g)
      yc(: ,jm+g) = TWO * yc(: ,jm-1+g) - yc(: ,jm-2+g)
      vol( 1-g,1:jm) = TWO * vol( 2-g, 1:jm) - vol( 3-g, 1:jm)
      vol(im+g,1:jm) = TWO * vol(im-1+g, 1:jm) - vol(im-2+g, 1:jm)
      vol(: , 1-g) = TWO * vol(: , 2-g) - vol(: , 3-g)
      vol(: ,jm+g) = TWO * vol(: ,jm-1+g) - vol(: ,jm-2+g)
      nx( 1-g,:,:) = TWO * nx( 2-g,:,:) - nx( 3-g,:,:)
      ny( 1-g,:,:) = TWO * ny( 2-g,:,:) - ny( 3-g,:,:)
      nx(im+1+g,:,:) = TWO * nx(im +g,:,:) - nx(im-1+g,:,:)
      ny(im+1+g,:,:) = TWO * ny(im +g,:,:) - ny(im-1+g,:,:)
      nx(: , 1-g,:) = TWO * nx(:, 2-g,:) - nx(: , 3-g,:)
      ny(: , 1-g,:) = TWO * ny(:, 2-g,:) - ny(: , 3-g,:)
      nx(: ,1+jm+g,:) = TWO * nx(:, jm+g,:) - nx(: ,g+jm-1,:)
      ny(: ,1+jm+g,:) = TWO * ny(:, jm+g,:) - ny(: ,g+jm-1,:)
    enddo
  enddo
!$AD II-LOOP
  do j=1, jm+1
!$AD II-LOOP
!DIR$ IVDEP
  do i=1, im+1
    volf(i,j,1) = TWO/(vol(i,j)+vol(i-1,j))
    volf(i,j,2) = TWO/(vol(i,j)+vol(i,j-1))
  enddo
  enddo
end subroutine computegeom_2d
