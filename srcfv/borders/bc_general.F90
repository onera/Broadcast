! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
! =============================================================================
!                             BC general 2D
! =============================================================================

subroutine bc_general_2d(w,loc,interf,field,gh,im,jm,lm)
  !
  implicit none
  ! Variables for dimension ---------------------------------------------------
  integer,intent(in) :: lm,im,jm
  integer,intent(in) :: gh
  ! Input variables -----------------------------------------------------------
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(lm,gh,5),intent(in) :: field
  ! Output variables ----------------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,5),intent(inout) :: w
#include "init_2d.F"
 
!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    ! write(1200,*) ' '
    ! write(1200,*) ' =====================', j
    do de = 1,gh
      w(i-de*i0,j-de*j0,1) = field(l,de,1)
      w(i-de*i0,j-de*j0,2) = field(l,de,2)
      w(i-de*i0,j-de*j0,3) = field(l,de,3)
      w(i-de*i0,j-de*j0,4) = field(l,de,4)
      w(i-de*i0,j-de*j0,5) = field(l,de,5)
      ! write(1200,*) "i, j, w, field1", i-de*i0, j-de*j0, field(l,de,1)
      ! write(1200,*) "i, j, w, field2", i-de*i0, j-de*j0, field(l,de,2)
      ! write(1200,*) "i, j, w, field3", i-de*i0, j-de*j0, field(l,de,3)
      ! write(1200,*) "i, j, w, field4", i-de*i0, j-de*j0, field(l,de,4)
      ! write(1200,*) "i, j, w, field5", i-de*i0, j-de*j0, field(l,de,5)
    end do
  end do
  !
end subroutine bc_general_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
