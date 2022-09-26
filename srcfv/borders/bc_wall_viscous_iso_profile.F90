! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
subroutine bc_wall_viscous_iso_profile_2d(w,twallprof,loc,gam,rgaz,interf,gh,im,jm,lm)
  !
  implicit none
  ! Variable for dimension ------------------------------------------
  integer,intent(in) :: im,jm,lm
  integer,intent(in) :: gh
  ! Input variables -------------------------------------------------
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8), intent(in) :: gam
  real(8), intent(in) :: rgaz
  real(8),dimension(lm),intent(in) :: twallprof
  ! Returned variables ----------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,5),intent(inout) :: w
  ! Local variables -------------------------------------------------
  integer :: da
  real(8) :: pe,roe,ue,ve,we
  real(8) :: pi,roi,ui,vi,wi,ei,roiei,pw
  real(8) :: roem1,roe1m1
  real(8) :: pe1,roe1,ue1,ve1,ve21,we1
  real(8) :: ve2,gami,gam1,ct0,ct1,THIRD,FOUR,HALF,ONE,TWO
  ! -----------------------------------------------------------------
#include "init_2d.F"
  
  !bid
  gam1  = gam - 1.d0
  gami  = 1.d0/gam
  ct0   = 1.125d0! 9/8
  ct1   =-0.125d0!-1/8
  THIRD = 1.d0/3.d0
  FOUR  = 4.d0
  HALF  = 0.5d0
  ONE   = 1.d0
  TWO   = 2.d0
  !bid
  !
!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    !
    roe   = w(i,j,1)
    roem1 = ONE/roe
    ue    = w(i,j,2)*roem1
    ve    = w(i,j,3)*roem1
    we    = w(i,j,4)*roem1
    ve2   = ue*ue + ve*ve + we*we
    pe    = gam1 * (w(i,j,5) - HALF*roe*ve2)
    !
    roe1   = w(i+i0,j+j0,1)
    roe1m1 = ONE/roe1
    ue1    = w(i+i0,j+j0,2)*roe1m1
    ve1    = w(i+i0,j+j0,3)*roe1m1
    we1    = w(i+i0,j+j0,4)*roe1m1
    ve21   = ue1*ue1 + ve1*ve1 + we1*we1
    pe1    = gam1 * (w(i+i0,j+j0,5) - HALF*roe1*ve21)
    !
    ! extrap dp/dn = 0 o2
    pw = ct0 * pe + ct1 * pe1
    !! pw = pe
    !!
    !! ???
    roi = TWO * pw / (rgaz * twallprof(l)) - roe
    ! pi  = pw
    pi  = THIRD * (FOUR * pw - pe)
    ! roi = (pi/pe*roe**gam )**gami
    !
    roiei  = pi/gam1
    !
    ui = - ue
    vi = - ve
    wi = - we
    !
    ! roi = roe
    ! pi  = pe
    ! roiei  = pi/gam1 + HALF*roi*(ui*ui+vi*vi+wi*wi)
    
    do de = 1, gh
      
      w(i-de*i0,j-de*j0,1) = roi
      w(i-de*i0,j-de*j0,2) = roi*ui
      w(i-de*i0,j-de*j0,3) = roi*vi
      w(i-de*i0,j-de*j0,4) = roi*wi
      w(i-de*i0,j-de*j0,5) = roiei + HALF*roi*(ui*ui+vi*vi+wi*wi)
      !
      ! roi   = w(i+de*i0,j+de*j0,1)
      ! roem1 = ONE/roi
      da = de - 1
      roi   = TWO * w(i-de*i0,j-de*j0,1) - w(i-da*i0,j-da*j0,1)
      roem1 = ONE/w(i+de*i0,j+de*j0,1)
      ue1   = w(i+de*i0,j+de*j0,2)*roem1
      ve1   = w(i+de*i0,j+de*j0,3)*roem1
      we1   = w(i+de*i0,j+de*j0,4)*roem1
      ! ve21  = ue1*ue1 + ve1*ve1 + we1*we1
      ! roiei = w(i+de*i0,j+de*j0,5)
      
      ui    = - ue1
      vi    = - ve1
      wi    = - we1
      
    enddo
    !
  end do
  !
end subroutine bc_wall_viscous_iso_profile_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
