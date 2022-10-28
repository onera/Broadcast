! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
! =============================================================================
!                             BC general 2D
! =============================================================================

subroutine bc_general_2d(w,loc,interf,field,gh,im,jm,lm,em,gh1)
  !
  implicit none
  ! Variables for dimension ---------------------------------------------------
  integer,intent(in) :: lm,im,jm,em
  integer,intent(in) :: gh,gh1
  ! Input variables -----------------------------------------------------------
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(lm,gh1,em),intent(in) :: field
  ! Output variables ----------------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
#include "init_2d.F"
 
!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    ! write(1200,*) ' '
    ! write(1200,*) ' =====================', j
    do de = 1,gh
      w(i-de*i0,j-de*j0,:) = field(l,de,:)
    end do
  end do
  !
end subroutine bc_general_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
