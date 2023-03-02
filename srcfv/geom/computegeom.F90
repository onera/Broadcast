! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

subroutine computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1  ),intent(inout) :: x0
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1  ),intent(inout) :: y0
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(inout) :: nx
  real(8),dimension(1-gh:im+gh+1,1-gh:jm+gh+1,2),intent(inout) :: ny
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ),intent(inout) :: xc
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ),intent(inout) :: yc 
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ),intent(inout) :: vol
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,2),intent(inout) :: volf
  ! Non-required arguments -------------------------------------------
  integer :: i,j,g,dummy
  real(8) :: TWO,HALF,FOURTH
  real(8) :: xa,xb,xc1,xd,ya,yb,yc1,yd
  real(8) :: abx,aby,acx,acy,q1
  real(8) :: dbx,dby,dcx,dcy,q2
  ! ------------------------------------------------------------------
  TWO    = 2.d0
  HALF   = 0.5d0
  FOURTH = 0.25d0
  !extend in ghost 
  do g=1,gh
    do dummy = 1,2

    
      x0(   1-g,:) = TWO * x0(   2-g,:) - x0(   3-g,:)
      x0(im+1+g,:) = TWO * x0(im  +g,:) - x0(im-1+g,:)

      x0(:   , 1-g) = TWO * x0(:     ,   2-g) - x0(:     ,   3-g)
      x0(:   ,jm+1+g) = TWO * x0(:     ,jm+g) - x0(:     ,jm-1+g)
       
      y0(:,   1-g) = TWO * y0(:,   2-g) - y0(:,   3-g)
      y0(:,jm+1+g) = TWO * y0(:,jm  +g) - y0(:,jm-1+g)

      y0( 1-g,:) = TWO * y0(   2-g,  :) - y0(   3-g,  :)
      y0(im+1+g,:) = TWO * y0(im+g,  :) - y0(im-1+g,  :)
    enddo
  
  enddo
  
!$AD II-LOOP
  do j=1, jm+1
!$AD II-LOOP
!DIR$ IVDEP      
  do i=1, im+1
#include "geom/centers.F"
#include "geom/volumes.F"
#include "geom/normals_idir.F"
#include "geom/normals_jdir.F"
  enddo
  enddo  
  
  !extend in ghost 
  do g=1,gh
    do dummy = 1,2
      xc( 1-g,:) = TWO * xc(   2-g,  :) - xc(   3-g,  :)
      xc(im+g,:) = TWO * xc(im-1+g,  :) - xc(im-2+g,  :)
      
      xc(:   , 1-g) = TWO * xc(:     ,   2-g) - xc(:     ,   3-g)
      xc(:   ,jm+g) = TWO * xc(:     ,jm-1+g) - xc(:     ,jm-2+g)
     
      yc( 1-g,:) = TWO * yc(   2-g,  :) - yc(   3-g,  :)
      yc(im+g,:) = TWO * yc(im-1+g,  :) - yc(im-2+g,  :)
      
      yc(:   , 1-g) = TWO * yc(:     ,   2-g) - yc(:     ,   3-g)
      yc(:   ,jm+g) = TWO * yc(:     ,jm-1+g) - yc(:     ,jm-2+g)
   
       
      vol( 1-g,1:jm) = TWO * vol(   2-g,  1:jm) - vol(   3-g,  1:jm)
      vol(im+g,1:jm) = TWO * vol(im-1+g,  1:jm) - vol(im-2+g,  1:jm)
      
      vol(:   , 1-g) = TWO * vol(:     ,   2-g) - vol(:     ,   3-g)
      vol(:   ,jm+g) = TWO * vol(:     ,jm-1+g) - vol(:     ,jm-2+g)
   
      nx(   1-g,:,:) = TWO * nx(   2-g,:,:) - nx(   3-g,:,:)
      ny(   1-g,:,:) = TWO * ny(   2-g,:,:) - ny(   3-g,:,:)
      nx(im+1+g,:,:) = TWO * nx(im  +g,:,:) - nx(im-1+g,:,:)
      ny(im+1+g,:,:) = TWO * ny(im  +g,:,:) - ny(im-1+g,:,:)
      
      nx(: ,   1-g,:) = TWO * nx(:,   2-g,:) - nx(: ,   3-g,:)
      ny(: ,   1-g,:) = TWO * ny(:,   2-g,:) - ny(: ,   3-g,:)
      nx(: ,1+jm+g,:) = TWO * nx(:,  jm+g,:) - nx(: ,g+jm-1,:)
      ny(: ,1+jm+g,:) = TWO * ny(:,  jm+g,:) - ny(: ,g+jm-1,:)
    
    enddo

  enddo
!$AD II-LOOP
  do j=1, jm+1
!$AD II-LOOP
!DIR$ IVDEP      
  do i=1, im+1
    volf(i,j,1) = TWO/(vol(i,j)+vol(i-1,j))
    volf(i,j,2) = TWO/(vol(i,j)+vol(i,j-1))
  enddo
  enddo 
end subroutine computegeom_2d
