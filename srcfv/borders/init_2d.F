! Local variables -----------------------------------------------------------
  integer :: kdir,de,i,j,imin,imax,jmin,jmax,lmin,lmax,l,i0,j0,high
  ! ---------------------------------------------------------------------------
  !
  imin=interf(1,1)
  jmin=interf(1,2)
  imax=interf(2,1)
  jmax=interf(2,2)
  
!   write(200,*) loc, 'imin = ', interf(1,1)
!   write(200,*) loc, 'jmin = ', interf(1,2)
!   write(200,*) loc, 'imax = ', interf(2,1)
!   write(200,*) loc, 'jmax = ', interf(2,2)
  
  i0 = 0
  j0 = 0
  high = 0
  
  i = imin
  j = jmin
  
  lmin = 1
  
  if (loc.eq.'Ilo') then
    kdir = 1
    i0 = 1
    lmax = jmax - jmin + 1
  else if (loc.eq.'Ihi')then
    kdir = 1
    i0 = -1
    lmax = jmax - jmin + 1
    high = 1
  else if (loc.eq.'Jlo')then
    kdir=2
    j0 = 1
    lmax = imax - imin + 1
  else if (loc.eq.'Jhi')then
    kdir=2
    j0 = -1
    lmax = imax - imin + 1
    high = 1
  endif
  
  
