
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
  
!$AD II-LOOP
      do j=1-gh,jm+gh
!$AD II-LOOP
!DIR$ IVDEP
      do i=1-gh,im+gh

    ro = w(i,j,1)
    rom1 = ONE/ro
    velx(i,j) = w(i,j,2) * rom1
    vely(i,j) = w(i,j,3) * rom1
    velz(i,j) = w(i,j,4) * rom1
!
    ec  = HALF*( velx(i,j)*velx(i,j) &
               + vely(i,j)*vely(i,j) &
               + velz(i,j)*velz(i,j))

! ec  = HALF*( velx(i,j)*velx(i,j) &
!            + vely(i,j)*vely(i,j))
!
    eloc = (w(i,j,5) - ec*ro)*rom1
!
    tloc(i,j) = eloc*cvm1
!
    p(i,j)  = (gam-ONE)*ro*eloc
!  p(i,j)  = ro*rgaz*tloc(i,j)
    
    htot= (w(i,j,5) + p(i,j))*rom1
    
    f(i,j,1)  = w(i,j,2) 
    f(i,j,2)  = w(i,j,2) * velx(i,j) + p(i,j)
    f(i,j,3)  = w(i,j,2) * vely(i,j)
    f(i,j,4)  = w(i,j,2) * velz(i,j)
    f(i,j,5)  = w(i,j,2) * htot
    
    g(i,j,1)  = w(i,j,3) 
    g(i,j,2)  = w(i,j,3) * velx(i,j) 
    g(i,j,3)  = w(i,j,3) * vely(i,j) + p(i,j)
    g(i,j,4)  = w(i,j,3) * velz(i,j)
    g(i,j,5)  = w(i,j,3) * htot
    
!h(i,j,1)  = w(i,j,4)
!h(i,j,2)  = w(i,j,4) * velx(i,j)
!h(i,j,3)  = w(i,j,4) * vely(i,j)
!h(i,j,4)  = w(i,j,4) * velz(i,j) + p(i,j)
!h(i,j,5)  = w(i,j,4) * htot
    
        mu(i,j) = betas/(tloc(i,j) + s_suth) * sqrt(tloc(i,j)) * tloc(i,j)
    
    
  
      enddo
      enddo
  
  
! Work on interior domain minus one cell
!$AD II-LOOP
  do j = 1, jm
!$AD II-LOOP
!DIR$ IVDEP
  do i = 1, im

a1 = gam1*gam*mach*mach*w(i,j,1)/tloc(i,j)
a2 = (HALF*(velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j))-cv*tloc(i,j))/w(i,j,1)

product(i,j,1) = HALF * vol(i,j) * (((velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j))/w(i,j,1)+tloc(i,j)/(w(i,j,1)*gam*mach*mach) &
                                    +a1*a2*a2)*wf(i,j,1) &
                                    - velx(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,2)   &
                                    - vely(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,3)    &
                                    + ZERO*wf(i,j,4)   &
                                    + a1*a2/w(i,j,1)*wf(i,j,5) )

product(i,j,2) = HALF * vol(i,j) * (- velx(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,1) &
                                    + (ONE/w(i,j,1)+velx(i,j)*velx(i,j)*a1/(w(i,j,1)*w(i,j,1)))*wf(i,j,2)   &
                                    + velx(i,j)*vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,3)    &
                                    + ZERO*wf(i,j,4)   &
                                    - velx(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,5) )

product(i,j,3) = HALF * vol(i,j) * (- vely(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,1) &
                                    + velx(i,j)*vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,2)    &
                                    + (ONE/w(i,j,1)+vely(i,j)*vely(i,j)*a1/(w(i,j,1)*w(i,j,1)))*wf(i,j,3)   & 
                                    + ZERO*wf(i,j,4)   &
                                    - vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,5) )

product(i,j,4) = HALF * vol(i,j) * (ZERO*wf(i,j,1) &
                                    + ZERO*wf(i,j,2)   &
                                    + ZERO*wf(i,j,3)    &
                                    + ONE/w(i,j,1)*wf(i,j,4)   &
                                    + ZERO*wf(i,j,5)  )

product(i,j,5) = HALF * vol(i,j) * (a1*a2/w(i,j,1)*wf(i,j,1) &
                                    - velx(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,2)   &
                                    - vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,3)   &
                                    + ZERO*wf(i,j,4)   &
                                    + a1/(w(i,j,1)*w(i,j,1))*wf(i,j,5) )
  enddo
  enddo
  


end subroutine chuproduct_2d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
