
! =============================================================================
!          Chu matrix-product 2D
! =============================================================================
!
subroutine chuproduct_2d(product,w,wf,vol,gh,mach,cv,&
                            gam,rgaz,cs,muref,tref,s_suth,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer :: im,jm,gh
  ! required arguments ----------------------------------------------
  real(8),intent(in) :: cv,gam,rgaz ! thermo
  real(8),intent(in) :: cs,muref,tref,s_suth ! viscosity
  real(8),intent(in) :: mach
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: w
  real(8),dimension(1:im  ,1:jm  , 5   ),intent(in) :: wf
  ! Returned objects ------------------------------------------------
  real(8),dimension(1:im  ,1:jm  , 5   ),intent(inout) :: product
  ! Non-required arguments ------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ) :: f,g
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: velx
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: vely
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: velz
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: tloc
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: p
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ) :: mu
  integer :: i,j,h
  real(8) :: ro,rom1,htot,eloc,ec
  real(8) :: betas
  real(8) :: cvm1
  real(8) :: a1,a2,gam1
  real(8) :: HALF,ONE,ZERO,TWO
  ! -----------------------------------------------------------------
  !
  HALF     = 0.5d0
  ONE      = 1.d0
  ZERO     = 0.d0
  TWO      = 2.d0
  
  cvm1     = ONE/cv 
  gam1     = gam - ONE
  
  ! Primitives
  betas = muref*(tref + cs)/(sqrt(tref)*tref) ! for Sutherland
  
#include "rhs/primvisc.F"  
  
  ! Work on interior domain minus one cell
!$AD II-LOOP  
  do j = 1, jm
!$AD II-LOOP
!DIR$ IVDEP                
  do i = 1, im
#include "chu/prod.F"
  enddo
  enddo
  


end subroutine chuproduct_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
