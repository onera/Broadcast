! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.


! =============================================================================
!                Implicit matrix free phase 2D RANS coupled
! =============================================================================
subroutine impli_matrix_free_rans_2d(dwi,nx,ny,src_i,w,dw,vol,volf,dtcoef,cfl, &
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
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     ),intent(in) :: src_i
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(inout) :: dwi
  ! Local variables -------------------------------------------------
  integer :: i,j,l,equa,kdir,i0,j0,ipt,le,eq2,ip,jp
  integer :: ipas,imin,jmin,imax,jmax,irspec,jrspec
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    )  :: tau0,mu,tloc,mut
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh,em )  :: d1w,f,g,ssor
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,2)  :: coefdiag
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh,2,2)  :: coef !specrad
  real(8),dimension(em) :: wi,fi,gi,d2w,dfi,dgi
  real(8) :: velxi,velyi,velzi,pint,tint,htot
  real(8) :: rhom1,roki,norm,specdiff,iro,rol,ror,sigmainv
  real(8) :: ONE,HALF, ZERO,TWO,rom1l,rom1r,gampr,gamprt,mutot,difftur
  real(8) :: uu,vv,cc,gam1,itur
  real(8) :: dist,vitc,dt_euler,dt_ns,cflm1,distm1,cvm1 
  real(8) :: ec,eloc,betas,rom1,dtm1,dt, alpha
  real(8) :: xkhi, xkhi3, safv1,sacv13,sacv1,cdiagm1,step
  real(8) :: multi, multi2,nxloc, nyloc, cdiag2m1
  ! -----------------------------------------------------------------
  HALF = 0.5d0
  ZERO = 0.d0
  TWO  = 2.d0
  ONE  = 1.d0
  !
  sacv1  = 7.1d0 !SA NuTilde constante
  sacv13 = sacv1*sacv1*sacv1
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
  itur = 6 
  
  do j = 1-gh , jm+gh
  do i = 1-gh , im+gh
#include "lhs/PrimitivesLhs.F"
#include "phys/viscosity.F"
#include "tur/mutspalart.F"
#include "lhs/physfluxtur.F"
#include "lhs/time_step_tur.F"
  ! compute diagonal coefficient------------------------------------------------------
  ! coefdiag(i,j,1) = dtcoef*dtm1*vol(i,j) ! 1/dt ou 1.5/dt
  ! coefdiag(i,j,2) = dtcoef*dtm1*vol(i,j) ! 1/dt ou 1.5/dt
  
  tau0(i,j) = dt/vol(i,j)
  
  enddo      
  enddo      
  !

  do kdir=1,2
    i0 = -kdir+2
    j0 = kdir-1
    !
    do j = 1 , jm+1
    do i = 1 , im+1
#include "lhs/specradtur.F"
    enddo
    enddo
  enddo 
  
  do j = 1 , jm
  do i = 1 , im
    coefdiag(i,j,1)=tau0(i,j) * (  coef( i,j,1,1) + coef(i+1,j  ,1,1) &
                                 + coef( i,j,1,2) + coef(i  ,j+1,1,2) ) + dtcoef
    coefdiag(i,j,2)=tau0(i,j) * (  coef( i,j,2,1) + coef(i+1,j  ,2,1) &
                                 + coef( i,j,2,2) + coef(i  ,j+1,2,2) &
                                 + src_i(i,j)                         ) + dtcoef
  
    !copy rhs
    d1w(i,j,:) =  dw(i,j,:) * tau0(i,j)
  enddo
  enddo
  !
  ssor = ZERO
  !
  loop_subite: do l=1,lmax
    !
    ! do equa = 1,em
    !     d2w(0:im+1,0:jm+1,equa) = dw(0:im+1,0:jm+1,equa)
    ! enddo
    step = (-1.d0)**(l+1) 
    if(step.gt.0.5d0) then
      imin =  1
      jmin =  1
      imax = im
      jmax = jm
      ipas =  1
    else
      imin = im
      jmin = jm
      imax =  1
      jmax =  1
      ipas = -1
    endif
    !===========================================================
    !  L phase : L =  1/2 ( df + rspec dw(i-1) )
    !  U phase : U = -1/2 ( df - rspec dw(i+1) )
    !===========================================================
    ! step = (-1.d0)**(l-1)
    ! step = 1.d0
    do j=1,jm
    do i=1,im
      dwi(i,j,:) = d1w(i,j,:) + ssor(i,j,:)
    enddo
    enddo
    ssor = ZERO
    
    
    !
    ! Computation of the left hand side
    !
    !
    do j = jmin , jmax, ipas
    do i = imin , imax, ipas
  
      d2w(:) = ZERO

      do kdir=1,2
          !
          i0 = ipas*(-kdir+2)
          j0 = ipas*( kdir-1)
          !
          ip = i - i0
          jp = j - j0
          !
          irspec = i - i0 * (1-ipas)/2
          jrspec = j - j0 * (1-ipas)/2
          !
          nxloc = nx(irspec,jrspec,kdir)
          nyloc = ny(irspec,jrspec,kdir)
          
          
          !fluxes
          wi(:) = w(ip,jp,:) + dwi(ip,jp,:)
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
          fi(6) = velxi * wi(6)
          gi(6) = velyi * wi(6)
          
          dfi(:) = fi(:) - f(ip,jp,:)
          dgi(:) = gi(:) - g(ip,jp,:)
          

          multi = HALF * step
          multi2 = coef(irspec,jrspec,1,kdir)

          d2w(1:5) = d2w(1:5) + multi  *( dfi(1:5)*nxloc +  &
                                          dgi(1:5)*nyloc  ) &
                              + multi2 *  dwi(ip,jp,1:5)
          
          multi2 = coef(irspec,jrspec,2,kdir)
          d2w(6) = d2w(6) + multi  *(dfi(6)*nxloc +  &
                                     dgi(6)*nyloc  ) &
                          + multi2 *  dwi(ip,jp,6)
          !
          ! end of kdir loop
        enddo
        !
        !===========================================================
        ! solve the system
        !   D Dw = rhs
        !===========================================================
        !
        cdiagm1  = ONE/coefdiag(i,j,1)
        cdiag2m1 = ONE/coefdiag(i,j,2)
        
        d2w(1)       =  d2w(1) * tau0(i,j)
        ssor(i,j,1)  =  d2w(1) + ssor(i,j,1)
        dwi(i,j,1)   = (d2w(1) +  dwi(i,j,1)) * cdiagm1
        
        d2w(2)       =  d2w(2) * tau0(i,j)
        ssor(i,j,2)  =  d2w(2) + ssor(i,j,2)
        dwi(i,j,2)   = (d2w(2) +  dwi(i,j,2)) * cdiagm1
        
        d2w(3)       =  d2w(3) * tau0(i,j)
        ssor(i,j,3)  =  d2w(3) + ssor(i,j,3)
        dwi(i,j,3)   = (d2w(3) +  dwi(i,j,3)) * cdiagm1
        
        d2w(4)       =  d2w(4) * tau0(i,j)
        ssor(i,j,4)  =  d2w(4) + ssor(i,j,4)
        dwi(i,j,4)   = (d2w(4) +  dwi(i,j,4)) * cdiagm1
        
        d2w(5)       =  d2w(5) * tau0(i,j)
        ssor(i,j,5)  =  d2w(5) + ssor(i,j,5)
        dwi(i,j,5)   = (d2w(5) +  dwi(i,j,5)) * cdiagm1
        
        d2w(6)       =  d2w(6) * tau0(i,j)
        ssor(i,j,6)  =  d2w(6) + ssor(i,j,6)
        dwi(i,j,6)   = (d2w(6) +  dwi(i,j,6)) * cdiag2m1
        !/(coefdiag(i,j,2) + src_i(i,j) * tau0(i,j))
        
    enddo
    enddo
    !
    ! Actualisation de wi
    !
! BCs on dwi (Neumann treatment)
    if (l.lt.lmax) then
      j = jm+1
      do i=1,im
        dwi(i,j,1) = dwi(i,j-1,1)
        dwi(i,j,2) = dwi(i,j-1,2)
        dwi(i,j,3) = dwi(i,j-1,3)
        dwi(i,j,4) = dwi(i,j-1,4)
        dwi(i,j,5) = dwi(i,j-1,5)
        dwi(i,j,6) = dwi(i,j-1,6)
      enddo
      j = 0 
      do i=1,im
        dwi(i,j,1) = dwi(i,j+1,1)
        dwi(i,j,2) = dwi(i,j+1,2)
        dwi(i,j,3) = dwi(i,j+1,3)
        dwi(i,j,4) = dwi(i,j+1,4)
        dwi(i,j,5) = dwi(i,j+1,5)
        dwi(i,j,6) = dwi(i,j+1,6)
      enddo
      i = im+1
      do j=1,jm
        dwi(i,j,1) = dwi(i-1,j,1)
        dwi(i,j,2) = dwi(i-1,j,2)
        dwi(i,j,3) = dwi(i-1,j,3)
        dwi(i,j,4) = dwi(i-1,j,4)
        dwi(i,j,5) = dwi(i-1,j,5)
        dwi(i,j,6) = dwi(i-1,j,6)
      enddo
      i = 0 
      do j=1,jm
        dwi(i,j,1) = dwi(i+1,j,1)
        dwi(i,j,2) = dwi(i+1,j,2)
        dwi(i,j,3) = dwi(i+1,j,3)
        dwi(i,j,4) = dwi(i+1,j,4)
        dwi(i,j,5) = dwi(i+1,j,5)
        dwi(i,j,6) = dwi(i+1,j,6)
      enddo
    else
      j = jm+1
      do i=1,im
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
        dwi(i,j,6) = ZERO
      enddo
      j = 0 
      do i=1,im
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
        dwi(i,j,6) = ZERO
      enddo
      i = im+1
      do j=1,jm
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
        dwi(i,j,6) = ZERO
      enddo
      i = 0 
      do j=1,jm
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
        dwi(i,j,6) = ZERO
      enddo
    
    endif
        
    !
  enddo loop_subite
  

  !
end subroutine impli_matrix_free_rans_2d
!

