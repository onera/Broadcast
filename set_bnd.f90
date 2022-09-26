! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
subroutine set_bndbl_2d(w,field,wbd,im,jm,gh)
  ! we suppose initfield correct in gh
  implicit none
   ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: w
  ! Returned objects ------------------------------------------------
  real(8),dimension(    jm      ,     gh     , 5   ),intent(inout) :: field
  real(8),dimension(    im+gh                   , 5   ),intent(inout) :: wbd
  ! Non-required arguments -------------------------------------------
  integer :: i,j,depth
  
  do j = 1, jm
  do depth = 1,gh
      field(j,depth,:) = w(1-depth, j,:)  
  enddo
  enddo
  do i = 1,im+gh
    ! wbd(i,:) = 0.5d0 * (w(i-gh,jm,:) + w(i-gh,jm+1,:))
    wbd(i,:) = w(i-gh,jm,:)
  enddo
end subroutine set_bndbl_2d
