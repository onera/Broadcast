! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/

subroutine filteringlhs_2d(dwi,w,gh,im,jm,em)
  implicit none
! variables for dimension -----------------------------------------
  integer,intent(in) :: em,im,jm
  integer,intent(in) :: gh
!requiered
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(in) :: w
! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: dwi
! Local variables -------------------------------------------------
  integer :: i,j,eq,iro,iroe
  real(8) :: alpha,ZERO,ONE
  iro  = 1
  iroe = 5
  ZERO = 0.d0
  ONE  = 1.d0
!===========================================================
! filtering of conservative variables
!===========================================================
  do j=1,jm
  do i=1,im
    alpha = max(ZERO, -dwi(i,j,iro)/w(i,j,iro), -dwi(i,j,iroe)/w(i,j,iroe))
    alpha = ONE / (ONE + alpha)
    do eq=1,5
      dwi(i,j,eq) = alpha*dwi(i,j,eq)
    enddo
    if(em.gt.5) then
      do eq=6,em
        alpha       = max(ZERO, -dwi(i,j,eq)/max(w(i,j,eq), 1.d-15))
        alpha       = ONE / (ONE + alpha)
        dwi(i,j,eq) = alpha*dwi(i,j,eq)
      enddo
    endif
  enddo
  enddo

end subroutine filteringlhs_2d
!


