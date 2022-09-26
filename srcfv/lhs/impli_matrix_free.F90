! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.


! =============================================================================
!                Implicit matrix free phase 2D NS
! =============================================================================
subroutine impli_matrix_free_2d(dwi,nx,ny,w,dw,vol,volf,dtcoef,cfl, &
                                gam,rgaz,prandtl,lmax,gh,cv,cs, &
                                muref, tref, s_suth,im,jm,em)

  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer,intent(in) :: em,im,jm,lmax
  integer,intent(in) :: gh
  ! required arguments ----------------------------------------------
  real(8),intent(in) :: dtcoef,gam,rgaz,prandtl,cfl,cs, muref, tref, s_suth, cv
  real(8),dimension(1-gh:im+1+gh,1-gh:jm+1+gh, 2),intent(in) :: nx,ny
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 2),intent(in) :: volf
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     ),intent(in) :: vol
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(in) :: w
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(in) :: dw
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em    ),intent(inout) :: dwi
  ! Local variables -------------------------------------------------
  integer :: i,j,l,equa,kdir,i0,j0,ipt,le,eq2
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ) :: tloc,mu
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh,em ) :: d2w,dfi,dgi,hn,f,g
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ) :: coefdiag
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh ,2 ) :: coef !specrad
  real(8),dimension(em) :: wi,fi,gi
  real(8) :: velxi,velyi,velzi,pint,tint,htot
  real(8) :: rhom1,roki,norm,specdiff,iro,rol,ror,rom1
  real(8) :: ONE,HALF, ZERO,TWO,rom1l,rom1r,gampr,gamprt
  real(8) :: uu,vv,cc,gam1,mutot
  real(8) :: dtm1, dist,vitc,dt_euler,dt_ns,cvm1
  real(8) :: ec,eloc,betas,diagm1,cflm1,dtm1save
  ! -----------------------------------------------------------------
  HALF = 0.5d0
  ZERO = 0.d0
  TWO  = 2.d0
  ONE  = 1.d0
  !
  dwi = ZERO
  dfi = ZERO
  dgi = ZERO
  !
  gampr  = 2.d0*gam/prandtl
  gam1   = gam - ONE
  !
  cflm1  = ONE/cfl
  cvm1   = ONE/cv
  betas = muref*(tref + cs)/(sqrt(tref)*tref)
  
  do j = 1-gh , jm+gh
  do i = 1-gh , im+gh
#include "lhs/PrimitivesLhs.F"
#include "phys/viscosity.F"
#include "lhs/time_step.F"
  ! compute diagonal coefficient------------------------------------------------------
  coefdiag(i,j) = dtm1 * dtcoef * vol(i,j) ! dtcoef/tau0   
  enddo      
  enddo 
  
!   !for global timeStep
!   dtm1save = 1.d-15
!   do j = 1 , jm
!   do i = 1 , im
! #include "lhs/PrimitivesLhs.F"
! #include "phys/viscosity.F"
! #include "lhs/time_step.F"
!     IF (dtm1.gt.dtm1save) THEN
!       dtm1save = dtm1
!     ENDIF
!     dtm1 = dtm1save
!   enddo      
!   enddo 
!     !
!   coefdiag(1:im,1:jm) = dtm1 * dtcoef * vol(1:im,1:jm) ! dtcoef/tau0 
  

  !
  loop_kdir: do kdir=1,2
    i0 = -kdir+2
    j0 = kdir-1
    !
    do j = 1 , jm+j0
    do i = 1 , im+i0
#include "lhs/specrad_ns.F"
    enddo
    enddo
    do j = 1 , jm
    do i = 1 , im
      coefdiag(i,j)=coefdiag(i,j) + coef(i,j,kdir) + coef(i+i0,j+j0,kdir)
    enddo
    enddo
  enddo loop_kdir
  
  loop_subite: do l=1,lmax

    d2w = dw
    !
    ! Computation of the left hand side
    !
    loop_kdir_inner: do kdir=1,2
      i0 = -kdir+2
      j0 = kdir-1
      !
      do j = 1 , jm + j0
      do i = 1 , im + i0
        !norm = dsqrt( nx(i,j,kdir)**2 + ny(i,j,kdir)**2)
        hn(i,j,:)= &
             + HALF *(dfi(i-i0,j-j0,:)+dfi(i,j,:))*nx(i,j,kdir) &
             + HALF *(dgi(i-i0,j-j0,:)+dgi(i,j,:))*ny(i,j,kdir) 
      enddo
      enddo
      !
      do j = 1 , jm
      do i = 1 , im
        d2w(i,j,:) =  d2w(i,j,:)  + hn(i,j,:) - hn(i+i0,j+j0,:) &
                     + coef(i   ,j   ,kdir) * dwi(i-i0,j-j0,:) &
                     + coef(i+i0,j+j0,kdir) * dwi(i+i0,j+j0,:)
        
      enddo
      enddo
      
    enddo loop_kdir_inner
    !
    ! Computation of the intermediate increment
    !
    do j=1,jm
    do i=1,im
        diagm1 = ONE/coefdiag(i,j)
        dwi(i,j,1) = d2w(i,j,1) * diagm1
        dwi(i,j,2) = d2w(i,j,2) * diagm1
        dwi(i,j,3) = d2w(i,j,3) * diagm1
        dwi(i,j,4) = d2w(i,j,4) * diagm1
        dwi(i,j,5) = d2w(i,j,5) * diagm1
    enddo
    enddo
    !
    ! Actualisation de wi
    !

    do j = 1 , jm
    do i = 1 , im
      wi(:) = w(i,j,:) + dwi(i,j,:)
      rhom1 = ONE/wi(1)
      velxi = wi(2) * rhom1
      velyi = wi(3) * rhom1
      velzi = wi(4) * rhom1
      pint = gam1*( wi(5) - HALF*wi(1)*( velxi*velxi + velyi*velyi +velzi*velzi) )
      !
      ! Actualisation des flux intermediaires
      !
      fi(1) = wi(2)
      fi(2) = wi(2)*velxi + pint
      fi(3) = wi(2)*velyi
      fi(4) = wi(2)*velzi
      fi(5) = velxi*( wi(5) + pint )
      gi(1) = wi(3)
      gi(2) = wi(3)*velxi
      gi(3) = wi(3)*velyi + pint
      gi(4) = wi(3)*velzi
      gi(5) = velyi*( wi(5) + pint )
      !
      ! Calcul des increments des flux
      !
      dfi(i,j,:) = fi(:) - f(i,j,:)
      dgi(i,j,:) = gi(:) - g(i,j,:)
    enddo
    enddo
    !
  enddo loop_subite
  !
end subroutine impli_matrix_free_2d
!

