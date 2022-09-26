! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
subroutine compute_norml2(norm,nmoy,rhs,im,jm,gh)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh, 5),intent( in) :: rhs
  real(8),dimension(                       5),intent(out) :: norm, nmoy
  ! Non-required arguments -------------------------------------------
  integer :: i,j
  real(8),dimension(5) :: norml2
  ! ------------------------------------------------------------------
  norm   = 0.d0
  nmoy   = 0.d0
  norml2 = 0.d0
  
  do j=1, jm
  do i=1, im
      norml2(1) = norml2(1) + rhs(i,j,1) * rhs(i,j,1)
      norml2(2) = norml2(2) + rhs(i,j,2) * rhs(i,j,2)
      norml2(3) = norml2(3) + rhs(i,j,3) * rhs(i,j,3)
      norml2(4) = norml2(4) + rhs(i,j,4) * rhs(i,j,4)
      norml2(5) = norml2(5) + rhs(i,j,5) * rhs(i,j,5)
  enddo
  enddo  

  norm = sqrt(norml2)
  nmoy = norm/(im*jm)

end subroutine compute_norml2

subroutine compute_norml2inf(norm,ninf,rhs,im,jm,gh)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh, 5),intent( in) :: rhs
  real(8),dimension(                       5),intent(out) :: norm, ninf
  ! Non-required arguments -------------------------------------------
  integer :: i,j
  real(8),dimension(5) :: norml2
  real(8),dimension(5) :: norminf
  real(8):: expo,expom1
  ! ------------------------------------------------------------------
  norm    = 0.d0
  ninf    = 0.d0
  norml2  = 0.d0
  norminf = 0.d0
  expo    = 10.d0
  expom1  = 0.1d0
  
  do j=1, jm
  do i=1, im
      norml2(1) = norml2(1) + rhs(i,j,1) * rhs(i,j,1)
      norml2(2) = norml2(2) + rhs(i,j,2) * rhs(i,j,2)
      norml2(3) = norml2(3) + rhs(i,j,3) * rhs(i,j,3)
      norml2(4) = norml2(4) + rhs(i,j,4) * rhs(i,j,4)
      norml2(5) = norml2(5) + rhs(i,j,5) * rhs(i,j,5)
      
      norminf(1) = norminf(1) + rhs(i,j,1) ** expo
      norminf(2) = norminf(2) + rhs(i,j,2) ** expo
      norminf(3) = norminf(3) + rhs(i,j,3) ** expo
      norminf(4) = norminf(4) + rhs(i,j,4) ** expo
      norminf(5) = norminf(5) + rhs(i,j,5) ** expo
       
      
  enddo
  enddo  

  norm = sqrt(norml2)
  ninf = (norminf)**expom1

end subroutine compute_norml2inf

subroutine compute_maxnorm(norm,imax,jmax,rhs,im,jm,gh)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh, 5),intent( in) :: rhs
  real(8),dimension(                       5),intent(out) :: norm
  integer,dimension(                       5),intent(out) :: imax,jmax
  ! Non-required arguments -------------------------------------------
  integer :: i,j
  real(8),dimension(5) :: rhsloc
  ! ------------------------------------------------------------------
  norm    = 0.d0
  
  do j=1, jm
  do i=1, im
      rhsloc(:)  = abs(rhs(i,j,:))
      norm(:)    = abs(norm(:))
      if (rhsloc(1).gt.norm(1)) then
          norm(1) =  rhs(i,j,1)
          imax(1) = i 
          jmax(1) = j 
      endif
      if (rhsloc(2).gt.norm(2)) then
          norm(2) =  rhs(i,j,2)
          imax(2) = i 
          jmax(2) = j
      endif
      if (rhsloc(3).gt.norm(3)) then
          norm(3) =  rhs(i,j,3)
          imax(3) = i 
          jmax(3) = j
      endif
      if (rhsloc(4).gt.norm(4)) then
          norm(4) =  rhs(i,j,4)
          imax(4) = i 
          jmax(4) = j
      endif
      if (rhsloc(5).gt.norm(5)) then
          norm(5) =  rhs(i,j,5)
          imax(5) = i 
          jmax(5) = j
      endif
  enddo
  enddo  

end subroutine compute_maxnorm
