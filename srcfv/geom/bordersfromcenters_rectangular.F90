
subroutine bordersfromcenters_rectangular_2d(x0,y0,xc,yc,im,jm,xmin,ymin)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  ! Returned objects ------------------------------------------------
  real(8),dimension(1:im+1,1:jm+1  ),intent(inout) :: x0
  real(8),dimension(1:im+1,1:jm+1  ),intent(inout) :: y0
  real(8),dimension(1:im  ,1:jm    ),intent(in) :: xc
  real(8),dimension(1:im  ,1:jm    ),intent(in) :: yc 
  real(8),intent(in) :: xmin, ymin
  ! Non-required arguments -------------------------------------------
  integer :: i,j,g,dummy
  real(8) :: TWO,HALF,FOURTH! ---------------------------------------
  TWO    = 2.d0
  HALF   = 0.5d0
  FOURTH = 0.25d0   


  do j=1, jm+1
    x0(1,j) = xmin
  enddo 


  do i=1, im +1 
    y0(i,1) = ymin
  enddo


  do j=1, jm
   
  do i=2, im+1

    x0(i,j) = TWO*xc(i-1,j) - x0(i-1,j)
    
  enddo
  enddo


  do j=2, jm+1
     
  do i=1, im

    y0(i,j) = TWO*yc(i,j-1) - y0(i,j-1)

  enddo
  enddo  

   
  do i=1, im+1

    x0(i,jm+1) = x0(i,jm)

  enddo

 
  do j=1, jm+1

    y0(im+1,j) = y0(im,j)

  enddo  
  


end subroutine bordersfromcenters_rectangular_2d
