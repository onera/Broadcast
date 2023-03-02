! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/


! =============================================================================
!                Implicit matrix free phase 2D NS
! =============================================================================
subroutine impli_matrix_free_polar_2d(dwi,ym,yc,nx,ny,w,dw,vol,volf,dtcoef,cfl, &
                                      gam,rgaz,prandtl,lmax,gh,cv,cs, &
                                      muref, tref, s_suth,im,jm,em)

!
  implicit none
! variables for dimension -----------------------------------------
  integer,intent(in) :: em,im,jm,lmax
  integer,intent(in) :: gh
! required arguments ----------------------------------------------
  real(8),intent(in) :: dtcoef,gam,rgaz,prandtl,cfl,cs, muref, tref, s_suth, cv
  real(8),intent(in) :: ym
  real(8),dimension(1-gh:im+1+gh,1-gh:jm+1+gh, 2),intent(in) :: nx,ny
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 2),intent(in) :: volf
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh     ),intent(in) :: vol,yc
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(in) :: w
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  ,em),intent(in) :: dw
! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh,1-gh:jm+gh,em    ),intent(inout) :: dwi
! Local variables -------------------------------------------------
  integer :: i,j,l,equa,kdir,i0,j0,ipt,le,eq2,ip,jp
  integer :: ipas,imin,jmin,imax,jmax,irspec,jrspec
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ) :: tau0,tloc,mu
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh,em ) :: d1w,hn,f,g,ssor
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ) :: coefdiag
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh ,2 ) :: coef !specrad
  real(8),dimension(em) :: wi,fi,gi,d2w,dfi,dgi
  real(8) :: velxi,velyi,velzi,pint,tint,htot
  real(8) :: rhom1,roki,norm,specdiff,iro,rol,ror,rom1
  real(8) :: ONE,HALF, ZERO,TWO,rom1l,rom1r,gampr,gamprt
  real(8) :: uu,vv,cc,gam1,mutot,cdiagm1,step
  real(8) :: dtm1, dist,vitc,dt_euler,dt_ns,cvm1,dt
  real(8) :: ec,eloc,betas,diagm1,cflm1,dtm1save
  real(8) :: multi, multi2,nxloc, nyloc
  real(8) :: srcipol,distm1, voldistm1,ur
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

    ror = w(i,j,1)
    rom1 = ONE/ror
    velxi = w(i,j,2) * rom1
    velyi = w(i,j,3) * rom1
    velzi = w(i,j,4) * rom1
!
    ec  = HALF*( velxi*velxi &
               + velyi*velyi &
               + velzi*velzi)
!
    eloc = (w(i,j,5) - ec*ror)*rom1
!
    tloc(i,j) = eloc*cvm1
!
    pint  = ror*rgaz*tloc(i,j)
    
    htot= (w(i,j,5) + pint)*rom1
    
    f(i,j,1)  = w(i,j,2) 
    f(i,j,2)  = w(i,j,2) * velxi + pint
    f(i,j,3)  = w(i,j,2) * velyi
    f(i,j,4)  = w(i,j,2) * velzi
    f(i,j,5)  = w(i,j,2) * htot
    
    g(i,j,1)  = w(i,j,3) 
    g(i,j,2)  = w(i,j,3) * velxi 
    g(i,j,3)  = w(i,j,3) * velyi + pint
    g(i,j,4)  = w(i,j,3) * velzi
    g(i,j,5)  = w(i,j,3) * htot
    
!h(i,j,1)  = w(i,j,4)
!h(i,j,2)  = w(i,j,4) * velx(i,j)
!h(i,j,3)  = w(i,j,4) * vely(i,j)
!h(i,j,4)  = w(i,j,4) * velz(i,j) + pint
!h(i,j,5)  = w(i,j,4) * htot
    
        mu(i,j) = betas/(tloc(i,j) + s_suth) * sqrt(tloc(i,j)) * tloc(i,j)
    
    
  

    dist =  vol(i,j) / (sqrt((nx(i,j,1)*nx(i,j,1)+ny(i,j,1)*ny(i,j,1)))+ &
                        sqrt((nx(i,j,2)*nx(i,j,2)+ny(i,j,2)*ny(i,j,2)))  )
    

    cc    = sqrt(gam*pint*rom1)
    
    vitc  = sqrt(velxi*velxi+velyi*velyi) + cc
    
    dt_euler  = dist/vitc
    dt_ns     = 2.d0 * dist*dist*ror/(gampr*mu(i,j))
! 2.d0 is linked to the expression of gampr and gamprt
    
    
    dtm1      = cflm1/max(min(dt_euler,dt_ns),1.d-15) ! 1/dt
    dt        = max(min(dt_euler,dt_ns),1.d-15)*cfl ! 1/dt

!write(30101,*) "i,j    =  ", i,j
!write(30101,*) "dtm1   =  ", dtm1
!write(30101,*) "distm1 =  ", distm1
!write(30101,*) "cflm1  =  ", dtm1
!write(30101,*) "mu     =  ", mu(i,j)

    
    
    
    
    
    
    
! compute diagonal coefficient------------------------------------------------------
! coefdiag(i,j) = dtm1 * dtcoef * vol(i,j) ! dtcoef/tau0
    tau0(i,j) = dt/vol(i,j)
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
  do kdir=1,2
    i0 = -kdir+2
    j0 = kdir-1
!
    do j = 1 , jm+j0
    do i = 1 , im+i0
      rol   = w(i-i0,j-j0,1)
      ror   = w(i   ,j   ,1)
      rom1l = ONE/rol
      rom1r = ONE/ror
      iro   = 2.d0/(rol+ror)
      
      uu    = HALF*(w(i-i0,j-j0,2)*rom1l + w(i,j,2)*rom1r)
      vv    = HALF*(w(i-i0,j-j0,3)*rom1l + w(i,j,3)*rom1r)
        
      cc    = gam*rgaz*HALF*(tloc(i-i0,j-j0)+tloc(i,j))
       
      norm  = nx(i,j,kdir)*nx(i,j,kdir)+ny(i,j,kdir)*ny(i,j,kdir)
      
      mutot =     HALF * gampr  * ( mu (i-i0,j-j0) + mu (i,j) ) 

      specdiff = mutot * norm * volf(i,j,kdir) * iro
      
      coef(i,j,kdir) = ONE  * ( abs(  uu*nx(i,j,kdir)  &
                                    + vv*ny(i,j,kdir)) &
                                    + sqrt(cc*norm)    &
                                    + TWO * specdiff   )

    enddo
    enddo
    
  enddo 
  
  do j = 1 , jm
  do i = 1 , im
    
    distm1    = ONE/(yc(i,j) - ym)
! if (distm1.gt.1.d-15) then
!    distm1    = ONE/abs( distm1 )
! endif
    
    voldistm1 = vol(i,j)* distm1
    
    ur = w(i,j,3)/w(i,j,1)
    
    srcipol = gam * ur *voldistm1 ! rspec of Euler part of source term
    coefdiag(i,j)=tau0(i,j) * (  coef(i,j,1) + coef(i+1,j  ,1) &
                               + coef(i,j,2) + coef(i  ,j+1,2) &
                               + srcipol                       ) + dtcoef
  
!copy rhs
    d1w(i,j,:) =  dw(i,j,:) * tau0(i,j)
    
  enddo
  enddo
  
  ssor = ZERO

  loop_subite: do l=1,lmax
!
    step = (-1.d0)**(l-1) 
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
          
          dfi(:) = fi(:) - f(ip,jp,:)
          dgi(:) = gi(:) - g(ip,jp,:)
          

          multi = HALF * step
          multi2 = coef(irspec,jrspec, kdir)

          d2w(1:5) = d2w(1:5) + multi  *( dfi(1:5)*nxloc +  &
                                          dgi(1:5)*nyloc  ) &
                              + multi2 *  dwi(ip,jp,1:5)
          
!
! end of kdir loop
        enddo
!
!===========================================================
! solve the system
!   D Dw = rhs
!===========================================================
!
        cdiagm1 = ONE/coefdiag(i,j)
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
      enddo
      j = 0 
      do i=1,im
        dwi(i,j,1) = dwi(i,j+1,1)
        dwi(i,j,2) = dwi(i,j+1,2)
        dwi(i,j,3) = dwi(i,j+1,3)
        dwi(i,j,4) = dwi(i,j+1,4)
        dwi(i,j,5) = dwi(i,j+1,5)
      enddo
      i = im+1
      do j=1,jm
        dwi(i,j,1) = dwi(i-1,j,1)
        dwi(i,j,2) = dwi(i-1,j,2)
        dwi(i,j,3) = dwi(i-1,j,3)
        dwi(i,j,4) = dwi(i-1,j,4)
        dwi(i,j,5) = dwi(i-1,j,5)
      enddo
      i = 0 
      do j=1,jm
        dwi(i,j,1) = dwi(i+1,j,1)
        dwi(i,j,2) = dwi(i+1,j,2)
        dwi(i,j,3) = dwi(i+1,j,3)
        dwi(i,j,4) = dwi(i+1,j,4)
        dwi(i,j,5) = dwi(i+1,j,5)
      enddo
    else 
      j = jm+1
      do i=1,im
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
      enddo
      j = 0 
      do i=1,im
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
      enddo
      i = im+1
      do j=1,jm
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
      enddo
      i = 0 
      do j=1,jm
        dwi(i,j,1) = ZERO
        dwi(i,j,2) = ZERO
        dwi(i,j,3) = ZERO
        dwi(i,j,4) = ZERO
        dwi(i,j,5) = ZERO
      enddo
    
    endif
    
!
  enddo loop_subite
!
end subroutine impli_matrix_free_polar_2d
!

