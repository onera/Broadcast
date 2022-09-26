! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
! =============================================================================
!                             BC symmetry 2D
! =============================================================================


subroutine bc_symmetry_2d(w,loc,interf,nx,ny,gh,im,jm)
  !
  implicit none
  ! Variables for dimension ---------------------------------------------------
  integer,intent(in) :: im,jm
  integer,intent(in) :: gh
  ! Input variables -----------------------------------------------------------
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(in) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(in) :: ny
  ! Output variables ----------------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,5),intent(inout) :: w
  ! Local variables -------------------------------------------------
  integer :: da, i1, j1
  real(8) :: ONE, TWO, HALF, nsum, nsumi, nxnorm, nynorm, nxloc, nyloc, sens
  real(8) :: rho,rhoinv,velx,vely,veln
  ! -----------------------------------------------------------------
#include "init_2d.F"
ONE = 1.d0
TWO = 2.d0
HALF= 0.5d0
sens = float(i0+j0)
!
nxloc = 0.d0
nyloc = 0.d0
i1 = i0 * i0
j1 = j0 * j0

!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    ! write(1200,*) ' '
    ! write(1200,*) ' =====================', j
    nxloc = nx(i+high*i1,j+high*j1,kdir)
    nyloc = ny(i+high*i1,j+high*j1,kdir)
    ! highis used to select the last normal when the face is of type *hi
    nsum    = sqrt(nxloc*nxloc + nyloc*nyloc)
    nsumi   = ONE/nsum
    nxnorm  = nxloc * nsumi * sens
    nynorm  = nyloc * nsumi * sens
    !
    do de = 1,gh
      da = de-1
      !
      rho = w(i+da*i0,j+da*j0,1)
      rhoinv = ONE/rho
      velx = w(i+da*i0,j+da*j0,2)*rhoinv
      vely = w(i+da*i0,j+da*j0,3)*rhoinv
      ! veln = ( velx*nxloc + vely*nyloc ) * nsumi
      veln = velx*nxnorm + vely*nynorm
      !
      ! velx = velx - TWO * veln * nxloc*nsumi
      velx = velx - TWO * veln * nxnorm
      ! vely = vely - TWO * veln * nyloc*nsumi
      vely = vely - TWO * veln * nynorm
      !
      w(i-de*i0,j-de*j0,1) = rho
      w(i-de*i0,j-de*j0,2) = rho * velx
      w(i-de*i0,j-de*j0,3) = rho * vely
      w(i-de*i0,j-de*j0,5) = w(i+da*i0,j+da*j0,5) &
                    -HALF* (  w(i+da*i0,j+da*j0,2)**2  &
                            + w(i+da*i0,j+da*j0,3)**2) &
                            / w(i+da*i0,j+da*j0,1) &
                    +HALF* (  w(i-de*i0,j-de*j0,2)**2  &
                            + w(i-de*i0,j-de*j0,3)**2) &
                            / w(i-de*i0,j-de*j0,1)
    enddo
  end do
end subroutine bc_symmetry_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
