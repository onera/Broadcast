subroutine bordersfromcenters_2d(x0,y0,xc,yc,im,jm)
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
  ! Non-required arguments -------------------------------------------
  integer :: i,j,g,dummy
  real(8) :: TWO,HALF,FOURTH! ----------------------------------------
  TWO    = 2.d0
  HALF   = 0.5d0
  FOURTH = 0.25d0    

!$AD II-LOOP
  do j=2, jm
!$AD II-LOOP
!DIR$ IVDEP      
  do i=2, im

    x0(i,j) = FOURTH*(xc(i-1,j-1) + xc(i,j-1) + xc(i-1,j) + xc(i,j))
    y0(i,j) = FOURTH*(yc(i-1,j-1) + yc(i,j-1) + yc(i-1,j) + yc(i,j))

  enddo
  enddo

!$AD II-LOOP
  do j=2, jm

    x0(im+1,j) = xc(im,j) + xc(im,j-1) - x0(im,j)
    y0(im+1,j) = yc(im,j) + yc(im,j-1) - y0(im,j)
    x0(1   ,j) = xc(1 ,j) + xc(1 ,j-1) - x0(2 ,j)
    y0(1   ,j) = yc(1 ,j) + yc(1 ,j-1) - y0(2 ,j)

  enddo

!$AD II-LOOP  
  do i=2, im

    x0(i,jm+1) = xc(i,jm) + xc(i-1,jm) - x0(i,jm)
    y0(i,jm+1) = yc(i,jm) + yc(i-1,jm) - y0(i,jm)
    x0(i,1   ) = xc(i,1 ) + xc(i-1,1 ) - x0(i,2 )
    y0(i,1   ) = yc(i,1 ) + yc(i-1,1 ) - y0(i,2 )

  enddo  

  x0(im+1,jm+1) = x0(im,jm+1) + x0(im+1,jm  ) - x0(im,jm)
  y0(im+1,jm+1) = y0(im,jm+1) + y0(im+1,jm  ) - y0(im,jm)

  x0(1   ,jm+1) = x0(1 ,jm  ) + x0(2   ,jm+1) - x0(2 ,jm)
  y0(1   ,jm+1) = y0(1 ,jm  ) + y0(2   ,jm+1) - y0(2 ,jm)

  x0(im+1,1   ) = x0(im,1   ) + x0(im+1,2   ) - x0(im,2 )
  y0(im+1,1   ) = y0(im,1   ) + y0(im+1,2   ) - y0(im,2 )

  x0(1   ,1   ) = x0(1 ,2   ) + x0(2   ,1   ) - x0(2 ,2 )
  y0(1   ,1   ) = y0(1 ,2   ) + y0(2   ,1   ) - y0(2 ,2 )


end subroutine bordersfromcenters_2d
