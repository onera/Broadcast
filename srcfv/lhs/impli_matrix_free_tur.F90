! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.


! =============================================================================
!                Implicit matrix free phase 2D RANS coupled
! =============================================================================
subroutine impli_matrix_free_rans_2d(dwi,nx,ny,src_i,w,mut,dw,vol,volf,dtcoef,cfl, &
                                    sigma,gam,rgaz,prandtl,prandtlturb,lmax,gh, &
                                    cv,cs, muref, tref, s_suth,im,jm,em)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer,intent(in) :: em,im,jm,lmax
  integer,intent(in) :: gh
  ! required arguments ----------------------------------------------
  real(8),intent(in) :: dtcoef,gam,rgaz,prandtl,prandtlturb,sigma
  real(8),intent(in) :: cv,cfl,cs, muref, tref, s_suth
  real(8),dimension(1-gh:im+1+gh,1-gh:jm+1+gh, 2),intent(in) :: nx,ny
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 2),intent(in) :: volf
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     ),intent(in) :: vol
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(in) :: w
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(in) :: dw
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     ),intent(in) :: mut
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     ),intent(in) :: src_i
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em    ),intent(inout) :: dwi
  ! Local variables -------------------------------------------------
  integer :: i,j,l,equa,kdir,i0,j0,ipt,le,eq2
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    )  :: mu,tloc
  real(8),dimension(   0:im+1   ,   0:jm+1,em  )  :: d2w,dfi,dgi,hn,f,g
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,2)  :: coefdiag
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh,2,2)  :: coef !specrad
  real(8),dimension(em) :: wi,fi,gi
  real(8) :: velxi,velyi,velzi,pint,tint,htot
  real(8) :: rhom1,roki,norm,specdiff,iro,rol,ror,sigmainv
  real(8) :: ONE,HALF, ZERO,TWO,rom1l,rom1r,gampr,gamprt,mutot,difftur
  real(8) :: uu,vv,cc,gam1,itur
  real(8) :: dist,vitc,dt_euler,dt_ns,cflm1,distm1,cvm1 
  real(8) :: ec,eloc,betas,rom1,dtm1,dt
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
  gamprt = 2.d0*gam/prandtlturb
  gam1   = gam - ONE
  !
  cflm1 = ONE/cfl
  cvm1  = ONE/cv
  betas = muref*(tref + cs)/(sqrt(tref)*tref)  
  
  sigmainv = ONE/sigma
  
  do j = 1-gh , jm+gh
  do i = 1-gh , im+gh
#include "lhs/PrimitivesLhs.F"
#include "phys/viscosity.F"
   do itur = 6,em
#include "lhs/physfluxtur.F"
   enddo
#include "lhs/time_step.F"
  ! compute diagonal coefficient------------------------------------------------------
  coefdiag(1:im,1:jm,1) = dtcoef*dtm1*vol(1:im,1:jm) ! 1/dt ou 1.5/dt
  coefdiag(1:im,1:jm,2) = dtcoef*dtm1*vol(1:im,1:jm) ! 1/dt ou 1.5/dt
  enddo      
  enddo      
  !

  loop_kdir: do kdir=1,2
    i0 = -kdir+2
    j0 = kdir-1
    !
    do j = 1 , jm+1
    do i = 1 , im+1
#include "lhs/specradtur.F"
    enddo
    enddo
    do j = 1 , jm
    do i = 1 , im
      coefdiag(i,j,1)=coefdiag(i,j,1) + coef(i,j,1,kdir) + coef(i+i0,j+j0,1,kdir)
      coefdiag(i,j,2)=coefdiag(i,j,2) + coef(i,j,2,kdir) + coef(i+i0,j+j0,2,kdir)
    enddo
    enddo
  enddo loop_kdir
  !
  loop_subite: do l=1,lmax
    !
    do equa = 1,em
        d2w(0:im+1,0:jm+1,equa) = dw(0:im+1,0:jm+1,equa)*vol(0:im+1,0:jm+1)
    enddo
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
        d2w(i,j,1:5) = d2w(i,j,1:5) + hn(i,j,1:5) - hn(i+i0,j+j0,1:5) &
                     + coef(i,j,1,kdir)       * dwi(i-i0,j-j0,1:5) &
                     + coef(i+i0,j+j0,1,kdir) * dwi(i+i0,j+j0,1:5)
        
      enddo
      enddo
      do equa = 6,em
      do j = 1 , jm
      do i = 1 , im
         d2w(i,j,equa) = d2w(i,j,equa) + hn(i,j,equa) - hn(i+i0,j+j0,equa) &
                     + coef(i,j,2,kdir)       * dwi(i-i0,j-j0,equa) &
                     + coef(i+i0,j+j0,2,kdir) * dwi(i+i0,j+j0,equa)
      enddo
      enddo
      enddo
    enddo loop_kdir_inner
    !
    ! Computation of the intermediate increment
    !
    dwi(1:im,1:jm,1   ) = d2w(1:im,1:jm,1   ) /  coefdiag(1:im,1:jm,1)
    dwi(1:im,1:jm,2   ) = d2w(1:im,1:jm,2   ) /  coefdiag(1:im,1:jm,1)
    dwi(1:im,1:jm,3   ) = d2w(1:im,1:jm,3   ) /  coefdiag(1:im,1:jm,1)
    dwi(1:im,1:jm,4   ) = d2w(1:im,1:jm,4   ) /  coefdiag(1:im,1:jm,1)
    dwi(1:im,1:jm,5   ) = d2w(1:im,1:jm,5   ) /  coefdiag(1:im,1:jm,1)
    do equa = 6,em
      dwi(1:im,1:jm,equa) = d2w(1:im,1:jm,equa) / (coefdiag(1:im,1:jm,2)+src_i(1:im,1:jm))
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
      ! turbulent quantities
      fi(6:em) = velxi * wi(6:em)
      gi(6:em) = velyi * wi(6:em)
      ! Calcul des increments des flux
      dfi(i,j,:) = fi(:) - f(i,j,:)
      dgi(i,j,:) = gi(:) - g(i,j,:)
    enddo
    enddo
    !
  enddo loop_subite
  !
end subroutine impli_matrix_free_rans_2d
!

