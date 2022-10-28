! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! @summary:       These subroutines compute extrapolation of a function
!                 at various orders of accuracy and sometimes using a
!                 mesh geometry to ponderate the extrapolation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bc_extrapolate_o0_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1
#include "init_2d.F"
  
             
!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
    
      !ordre 0
      w(i-i0*de,j-j0*de,:) = w(i-i0*d1,j-j0*d1,:) 
      
      
                                                             
  enddo
  enddo
  
end subroutine bc_extrapolate_o0_2d

subroutine bc_extrapolate_o2_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1,d2
  real(8) :: TWO
#include "init_2d.F"
  
  TWO  = 2.d0
             
!$AD II-LOOP
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
      d2 = de - 2
    
      !ordre 2
      w(i-i0*de,j-j0*de,:) = TWO * w(i-i0*d1,j-j0*d1,:) &
                           -       w(i-i0*d2,j-j0*d2,:)
      
      
                                                             
  enddo
  enddo
  
end subroutine bc_extrapolate_o2_2d

subroutine bc_extrapolate_o3_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1,d2,d3
  real(8) :: THREE
  
#include "init_2d.F"
  
  THREE = 3.d0
              
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
      d2 = de - 2
      d3 = de - 3
      !ordre 3
      w(i-i0*de,j-j0*de,:) = THREE * w(i-i0*d1,j-j0*d1,:) &
                           - THREE * w(i-i0*d2,j-j0*d2,:) &
                           +         w(i-i0*d3,j-j0*d3,:) 
      
                                                       
  enddo
  enddo
end subroutine bc_extrapolate_o3_2d

subroutine bc_extrapolate_o4_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1,d2,d3,d4
  real(8) :: FOUR, SIX
  
#include "init_2d.F"
  
  FOUR = 4.d0
  SIX  = 6.d0
              
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
      d2 = de - 2
      d3 = de - 3
      d4 = de - 4
      !ordre 4
      w(i-i0*de,j-j0*de,:) = FOUR*w(i-i0*d1,j-j0*d1,:) &
                           - SIX *w(i-i0*d2,j-j0*d2,:) &
                           + FOUR*w(i-i0*d3,j-j0*d3,:) &
                           -      w(i-i0*d4,j-j0*d4,:) 
      
                                                       
  enddo
  enddo
end subroutine bc_extrapolate_o4_2d


subroutine bc_extrapolate_o5_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1,d2,d3,d4,d5
  real(8) :: FIVE,TEN
  
#include "init_2d.F"
  
  FIVE = 5.d0
  TEN  = 10.d0
              
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
      d2 = de - 2
      d3 = de - 3
      d4 = de - 4
      d5 = de - 5
      !ordre 5
      w(i-i0*de,j-j0*de,:) = FIVE * w(i-i0*d1,j-j0*d1,:) &
                           - TEN  * w(i-i0*d2,j-j0*d2,:) &
                           + TEN  * w(i-i0*d3,j-j0*d3,:) &
                           - FIVE * w(i-i0*d4,j-j0*d4,:) &
                           +        w(i-i0*d5,j-j0*d5,:)
      
                                                       
  enddo
  enddo
end subroutine bc_extrapolate_o5_2d

subroutine bc_extrapolate_o7_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1,d2,d3,d4,d5,d6,d7  
  
#include "init_2d.F"
              
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
      d2 = de - 2
      d3 = de - 3
      d4 = de - 4
      d5 = de - 5
      d6 = de - 6
      d7 = de - 7
      !ordre 7
      w(i-i0*de,j-j0*de,:) = -  7.d0 * w(i-i0*d1,j-j0*d1,:) &
                             + 21.d0 * w(i-i0*d2,j-j0*d2,:) &
                             - 35.d0 * w(i-i0*d3,j-j0*d3,:) &
                             + 35.d0 * w(i-i0*d4,j-j0*d4,:) &
                             - 21.d0 * w(i-i0*d5,j-j0*d5,:) &
                             +  7.d0 * w(i-i0*d6,j-j0*d6,:) &
                             -         w(i-i0*d7,j-j0*d7,:)
      
      
      
      
                                                       
  enddo
  enddo
end subroutine bc_extrapolate_o7_2d
subroutine bc_extrapolate_o9_2d(w,loc,interf,im,jm,gh,em)
  !
  !
  implicit none
  integer,intent(in) :: im,jm,gh,em
  character(len=3),intent(in) :: loc
  integer,dimension(2,2),intent(in) :: interf
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em),intent(inout) :: w
  ! Local variables -----------------------------------------------------------
  integer :: d1,d2,d3,d4,d5,d6,d7,d8,d9
    
#include "init_2d.F"
              
  do l = lmin,lmax
    i =  imin + (l-lmin)*j0*j0
    j =  jmin + (l-lmin)*i0*i0
    do de = 1,gh
      d1 = de - 1
      d2 = de - 2
      d3 = de - 3
      d4 = de - 4
      d5 = de - 5
      d6 = de - 6
      d7 = de - 7
      d8 = de - 8
      d9 = de - 9
      !ordre 9
      w(i-i0*de,j-j0*de,:) = +   9.d0 * w(i-i0*d1,j-j0*d1,:) &
                             -  36.d0 * w(i-i0*d2,j-j0*d2,:) &
                             +  84.d0 * w(i-i0*d3,j-j0*d3,:) &
                             - 126.d0 * w(i-i0*d4,j-j0*d4,:) &
                             + 126.d0 * w(i-i0*d5,j-j0*d5,:) &
                             -  84.d0 * w(i-i0*d6,j-j0*d6,:) &
                             +  36.d0 * w(i-i0*d7,j-j0*d7,:) &
                             -   9.d0 * w(i-i0*d8,j-j0*d8,:) &
                             +          w(i-i0*d9,j-j0*d9,:)
        
      
      
      
                                                       
  enddo
  enddo
end subroutine bc_extrapolate_o9_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of file
