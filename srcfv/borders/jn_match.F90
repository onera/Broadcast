! This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

subroutine jn_match_2d(wr,prr,gh1r,gh2r,gh3r,gh4r,imr,jmr,&
                       wd,prd,gh1d,gh2d,gh3d,gh4d,imd,jmd,&
                       tr,em)
  !
  implicit none
  ! Variables for dimension -----------------------------------------
  integer,intent(in) :: em,imr,jmr,imd,jmd
  integer,intent(in) :: gh1r,gh2r,gh3r,gh4r
  integer,intent(in) :: gh1d,gh2d,gh3d,gh4d
  ! Input variables -------------------------------------------------
  real(8),dimension(1-gh1d:imd+gh2d,1-gh3d:jmd+gh4d,em),intent(in) :: wd
  integer,dimension(2,2),intent(in) :: prr,prd
  integer,dimension(2),intent(in) :: tr
  ! Output variables ------------------------------------------------
  real(8),dimension(1-gh1r:imr+gh2r,1-gh3r:jmr+gh4r,em),intent(inout) :: wr
  ! Local variables -------------------------------------------------
  integer,pointer :: ind1,ind2
  integer :: ir,jr
  integer,target :: id,jd
  integer :: istep,jstep,idir,jdir,idd,jdd,i,j
  integer :: i1,i2,j1,j2
  ! -----------------------------------------------------------------
  istep = sign(1,tr(1))
  jstep = sign(1,tr(2))
  idir = abs(tr(1))
  jdir = abs(tr(2))
  !
  if(istep == 1)then
    i1 = 1 ; i2 = prd(2,idir) - prd(1,idir) + 1
  else
    i1 = prd(2,idir) - prd(1,idir) + 1 ; i2 = 1
  endif
  !
  if(jstep == 1)then
    j1 = 1 ; j2 = prd(2,jdir) - prd(1,jdir) + 1
  else
    j1 = prd(2,jdir) - prd(1,jdir) + 1 ; j2 = 1
  endif
  !
  if(idir == 1)then
    ind1=>id
    ind2=>jd
  elseif(idir == 2)then
    ind1=>jd
    ind2=>id
  endif
  !
  idd = 1 ; jdd = 1
  !
  do j = j1,j2,jstep
    do i = i1,i2,istep
      ir = prr(1,1)+idd-1
      jr = prr(1,2)+jdd-1
      id = prd(1,idir)+i-1
      jd = prd(1,jdir)+j-1
      ! wr(ir,jr,:) = wd(ind1,ind2,:)
      wr(ir,jr,:) = wd(id,jd,:)
      idd = idd + 1
    enddo
    idd = 1
    jdd = jdd + 1
  enddo
  !
end subroutine jn_match_2d
!===============================================================================
!                          JOIN Match for Gradients at centers
!===============================================================================


! =============================================================================
!                              Join match 2-D
! =============================================================================

subroutine jn_match_grad_2d(wr,prr,gh1r,gh2r,gh3r,gh4r,imr,jmr,&
                            wd,prd,gh1d,gh2d,gh3d,gh4d,imd,jmd,&
                            tr)
  !
  implicit none
  ! Variables for dimension -----------------------------------------
  integer,intent(in) :: imr,jmr,imd,jmd
  integer,intent(in) :: gh1r,gh2r,gh3r,gh4r
  integer,intent(in) :: gh1d,gh2d,gh3d,gh4d
  ! Input variables -------------------------------------------------
  real(8),dimension(1-gh1d:imd+gh2d,1-gh3d:jmd+gh4d, 2),intent(in) :: wd
  integer,dimension(2,2),intent(in) :: prr,prd
  integer,dimension(2),intent(in) :: tr
  ! Output variables ------------------------------------------------
  real(8),dimension(1-gh1r:imr+gh2r,1-gh3r:jmr+gh4r, 2),intent(inout) :: wr
  ! Local variables -------------------------------------------------
  integer,pointer :: ind1,ind2
  integer :: ir,jr
  integer,target :: id,jd
  integer :: istep,jstep,idir,jdir,idd,jdd,i,j
  integer :: i1,i2,j1,j2
  ! -----------------------------------------------------------------
  istep = sign(1,tr(1))
  jstep = sign(1,tr(2))
  idir = abs(tr(1))
  jdir = abs(tr(2))
  !
  if(istep == 1)then
    i1 = 1 ; i2 = prd(2,idir) - prd(1,idir) + 1
  else
    i1 = prd(2,idir) - prd(1,idir) + 1 ; i2 = 1
  endif
  !
  if(jstep == 1)then
    j1 = 1 ; j2 = prd(2,jdir) - prd(1,jdir) + 1
  else
    j1 = prd(2,jdir) - prd(1,jdir) + 1 ; j2 = 1
  endif
  !
  if(idir == 1)then
    ind1=>id
    ind2=>jd
  elseif(idir == 2)then
    ind1=>jd
    ind2=>id
  endif
  !
  idd = 1 ; jdd = 1
  !
  do j = j1,j2,jstep
    do i = i1,i2,istep
      ir = prr(1,1)+idd-1
      jr = prr(1,2)+jdd-1
      id = prd(1,idir)+i-1
      jd = prd(1,jdir)+j-1
      wr(ir,jr,:) = wd(ind1,ind2,:)
      idd = idd + 1
    enddo
    idd = 1
    jdd = jdd + 1
  enddo
  !
end subroutine jn_match_grad_2d

!===============================================================================
!                          JOIN Match for Fluxes at faces
!===============================================================================

subroutine jn_match_fx_2d(fr,prr,gh1r,gh2r,gh3r,gh4r,imr,jmr,&
                                 fd,prd,gh1d,gh2d,gh3d,gh4d,imd,jmd,&
                                 locr,locd,tr,em,dimf)
  !
  implicit none
  ! Variables for dimension -----------------------------------------
  integer,intent(in) :: imr,jmr,imd,jmd,dimf,em
  integer,intent(in) :: gh1r,gh2r,gh3r,gh4r
  integer,intent(in) :: gh1d,gh2d,gh3d,gh4d
  ! Input variables -------------------------------------------------
  character(len=3),intent(in) :: locr,locd
  real(8),dimension(1-gh1d:imd+1+gh2d,1-gh3d:jmd+1+gh4d,em,dimf),intent(in) :: fd
  integer,dimension(2,2),intent(in) :: prr,prd
  integer,dimension(2),intent(in) :: tr
  ! Output variables ------------------------------------------------
  real(8),dimension(1-gh1r:imr+1+gh2r,1-gh3r:jmr+1+gh4r,em,dimf),intent(inout) :: fr
  ! Local variables -------------------------------------------------
  integer,pointer :: ind1,ind2
  integer :: ir,jr
  integer,target :: id,jd
  integer :: i,j,istep,jstep,idir,jdir,idd,jdd
  integer :: direction,dird
  integer,dimension(2) :: sgnNx
  integer,dimension(2) :: highir, highjr
  integer,dimension(2) :: highid, highjd,i0,j0,i1,i2,j1,j2
  ! -----------------------------------------------------------------
  istep = sign(1,tr(1))
  jstep = sign(1,tr(2))
  idir = abs(tr(1))
  jdir = abs(tr(2))
  sgnNx(idir) = istep
  sgnNx(jdir) = jstep
  !
  highir = 0; highjr = 0
  highid = 0; highjd = 0
  !
  i0 = 0 ; j0 = 0
  i0(idir) = 1 ; j0(jdir) =1
  !
  if (locr=='Ihi')then
      highir(1)=1
  elseif(locr=='Jhi')then
      highjr(2)=1
  endif
  !
  if(locd=='Ilo')then
      i0(1) = 0 ; j0(1) = 0
      if (idir == 1) then
          highid(1) = 1
      else
          highjd(1) = 1
      endif
   elseif(locd=='Ihi')then
      i0(1) = 0 ; j0(1) = 0
  elseif(locd=='Jlo')then
      i0(2) = 0 ; j0(2) = 0
      if (idir == 1) then
          highjd(2) = 1
      else
          highid(2) = 1
      endif
  elseif(locd=='Jhi')then
      i0(2) = 0 ; j0(2) = 0
  endif
  !
  if(istep == 1)then
    i1(:) = 1 ; i2(:) = prd(2,idir) - prd(1,idir) + 1 + i0(:)
  else
    i1(:) = prd(2,idir) - prd(1,idir) + 1 + i0(:) ; i2(:) = 1
  endif
  !
  if(jstep == 1)then
    j1(:) = 1 ; j2(:) = prd(2,jdir) - prd(1,jdir) + 1 + j0(:)
  else
    j1(:) = prd(2,jdir) - prd(1,jdir) + 1 + j0(:) ; j2(:) = 1
  endif
  !
  if(idir == 1)then
    ind1=>id
    ind2=>jd
  elseif(idir == 2)then
    ind1=>jd
    ind2=>id
  endif
  !
  do direction = 1,2
    !
    idd = 1 ; jdd = 1
    !
    ! get dir for the donor
    dird = abs(tr(direction))
    do j = j1(dird),j2(dird),jstep
      do i = i1(dird),i2(dird),istep
        ir = prr(1,1)+idd-1  + highir(direction)
        jr = prr(1,2)+jdd-1  + highjr(direction)
        id = prd(1,idir)+i-1 + highid(dird)
        jd = prd(1,jdir)+j-1 + highjd(dird)
        !
        fr(:,ir,jr,direction) = sgnNx(dird) * fd(:,ind1,ind2,dird)
        idd = idd + 1
      enddo
      idd = 1
      jdd = jdd + 1
    enddo
  enddo
  !
end subroutine jn_match_fx_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
