! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
subroutine initblfv(muinf,xc,yc,w,field,im,gh,jm)
  !
  implicit none
   ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  real(8),intent(in) :: muinf
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: xc 
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: yc
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: w
  real(8),dimension(    jm      ,     gh     , 5   ),intent(inout) :: field
  ! Non-required arguments -------------------------------------------
  real(8) :: eta,ubl,ue ! BL parameter
  real(8) :: ro, roE, rou, xdeltam1,rou2
  integer :: i,j

  !
  ! w contains state inf adim
  ro     = w(3,3,1)
  rou    = w(3,3,2)
  roE    = w(3,3,5)
  ue     = rou/ro
  rou2   = 0.5d0*rou*ue
  xdeltam1  = sqrt(rou/muinf)
  !init w
  do j=1,jm+gh
  do i=1-gh, im+gh
    eta = yc(i,j)*xdeltam1/sqrt(xc(i,j))
    if (eta.lt.1.d0) then
      ubl = 2.d0*eta -2.d0*eta*eta*eta + eta*eta*eta*eta
    else
      ubl = 1.d0
    endif
    ubl = ubl * ue
    w(i,j,2) = ro*ubl
    w(i,j,5) = roE - rou2 + 0.5d0*ro*ubl*ubl
  enddo
  enddo
  w(:,0,2) = - w(:,1,2) 
  !set field
  do j=1,jm
  do i=1-gh,0
    eta = yc(i,j)*xdeltam1*sqrt(xc(i,j))
    if (eta.lt.1.d0) then
      ubl = 2.d0*eta -2.d0*eta*eta*eta + eta*eta*eta*eta
    else
      ubl = 1.d0
    endif
    ubl = ubl * ue
    field(j,1-i,2) = ro*ubl
    field(j,1-i,5) = roE - rou2 + 0.5d0*ro*ubl*ubl
    ! write(209876,*) "i,j, field(j,i)", i,j, 
  enddo
  enddo

end subroutine initblfv
subroutine initbl(muinf,xc,yc,w,field,im,gh,jm)
  !
  implicit none
   ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  real(8),intent(in) :: muinf
  real(8),dimension(1-gh:im+gh                     ),intent(in) :: xc
  real(8),dimension(               1-gh:jm+gh      ),intent(in) :: yc 
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: w
  real(8),dimension(    jm      ,     gh     , 5   ),intent(inout) :: field
  ! Non-required arguments -------------------------------------------
  real(8) :: eta,ubl,ue ! BL parameter
  real(8) :: ro, roE, rou, xdeltam1,rou2
  integer :: i,j
  !
  ! w contains state inf adim
  ro     = w(3,3,1)
  rou    = w(3,3,2)
  roE    = w(3,3,5)
  ue     = rou/ro
  rou2   = 0.5d0*rou*ue
  xdeltam1  = sqrt(rou/muinf)
  !init w
  do j=1,jm+gh
  do i=1-gh, im+gh
    eta = yc(j)*xdeltam1/sqrt(xc(i))
    if (eta.lt.1.d0) then
      ubl = 2.d0*eta -2.d0*eta*eta*eta + eta*eta*eta*eta
    else
      ubl = 1.d0
    endif
    ubl = ubl * ue
    w(i,j,2) = ro*ubl
    w(i,j,5) = roE - rou2 + 0.5d0*ro*ubl*ubl
  enddo
  enddo
  w(:,0,2) = - w(:,1,2) 
  !set field
  do j=1,jm
  do i=1-gh,0
    eta = yc(j)*xdeltam1*sqrt(xc(i))
    if (eta.lt.1.d0) then
      ubl = 2.d0*eta -2.d0*eta*eta*eta + eta*eta*eta*eta
    else
      ubl = 1.d0
    endif
    ubl = ubl * ue
    field(j,1-i,2) = ro*ubl
    field(j,1-i,5) = roE - rou2 + 0.5d0*ro*ubl*ubl
    ! write(209876,*) "i,j, field(j,i)", i,j, 
  enddo
  enddo

end subroutine initbl
