! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/
! add source term due to time integration at 2nd order (for Gear, BDF2 algo)

subroutine source_t_2d(res,w,w_save,dwr,Phi,dt,dtloc,vol,gh,im,jm,em)
  !
  implicit none
  ! variables for dimension ---------------------------------------------------
  integer,intent(in) :: im,jm,em,gh
  ! required arguments --------------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh, em),intent(in)  :: w,w_save,dwr
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh    ),intent(in)  :: vol
  real(8),intent(in) :: phi,dt,dtloc
  ! returned objects ----------------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh, em),intent(inout)  :: res
  ! local variables -----------------------------------------------------------
  integer :: i,j
  real(8) :: dtm1, dtlocm1, ONE
  ! ---------------------------------------------------------------------------
  ONE     = 1.d0
  dtm1    = (ONE+phi)/dt
  dtlocm1 = phi/dtloc
  do j = 1-gh , jm+gh
  do i = 1-gh , im+gh
    !
    !ts(:,i,j) = (ONE/dt)*( (ONE+phi)*( w(:,i,j)-w_save(:,i,j) ) &
    !                       -      phi *  dwr(:,i,j) &
    !                     )
    !
    res(i,j,:) = res(i,j,:) - vol(i,j) * & 
                 (dtm1 * ( w(i,j,:)-w_save(i,j,:))- dtlocm1 * dwr(i,j,:))
                  
                    
    !
  enddo
  enddo
  !
end subroutine source_t_2d
