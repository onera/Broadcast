
subroutine testvector_poisson(wd,l,k,gh,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,jm,gh,l,k
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh    ),intent(inout) :: wd
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j  
  stepj = 1+2*gh
  wd = 0.d0
  do j = k+1, jm, stepj
!DIR$ IVDEP  
  do i = l+1, im, stepj
    wd(i, j) = 1.d0
  enddo
  enddo
end subroutine testvector_poisson

subroutine computejacobianfromjv_poisson_norelax(jac,ia,ja,resd,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + k * im*jm + l * im*jm*stepj

    ia(current +1)  = (j-1)  + (i-1) * jm

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = valj + vali * jm
        jac(current +1) = resd(i , j)
        ! jac(current +1) = 1.d0
      endif  
    endif    
    
  enddo
  enddo
end subroutine computejacobianfromjv_poisson_norelax

subroutine testvector_dchu(wd,m,gh,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,jm,gh,m
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: wd
  ! Non-required arguments -------------------------------------------
  integer :: i,j  
  wd = 0.d0
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    wd(i, j, m+1) = 1.d0
  enddo
  enddo
end subroutine testvector_dchu

subroutine computederivativechufromjv(jac,ia,ja,resd,m,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m
  ! required arguments ----------------------------------------------
  real(8),dimension(1:im  ,1:jm  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + m * im*jm*5

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    ja(current +1)  = m + (j-1) * 5 + (i-1) * jm*5
    jac(current +1) = resd(i , j, e)
    ! jac(current +1) = 1.d0

  enddo
  enddo
  enddo
end subroutine computederivativechufromjv

subroutine testvector_twall(twallprofd,m,gh,im)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,gh,m
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh),intent(inout) :: twallprofd
  ! Non-required arguments -------------------------------------------
  twallprofd = 0.d0
  twallprofd(m+1) = 1.d0
end subroutine testvector_twall

subroutine testvector_twall_short(twallprofd,m,istart,iend)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: istart,iend,m
  ! Returned objects ------------------------------------------------
  real(8),dimension(istart:iend),intent(inout) :: twallprofd
  ! Non-required arguments -------------------------------------------
  twallprofd = 0.d0
  twallprofd(m) = 1.d0
end subroutine testvector_twall_short

subroutine computejacobian_twall(jac,ia,ja,resd,m,gh,im,jm,nbentry,vol)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + m * im*jm*5

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    ja(current +1)  = m
    jac(current +1) = resd(i , j, e) / vol(i,j)
    ! jac(current +1) = 1.d0

  enddo
  enddo
  enddo
end subroutine computejacobian_twall

subroutine computejacobian_twall_notdivbyvol(jac,ia,ja,resd,m,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + m * im*jm*5

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    ja(current +1)  = m
    jac(current +1) = resd(i , j, e)
    ! jac(current +1) = 1.d0

  enddo
  enddo
  enddo
end subroutine computejacobian_twall_notdivbyvol

subroutine computejacobian_twall_restrict(jac,ia,ja,resd,m,gh,im,jm,nbentry,vol)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, 2*gh+1
!DIR$ IVDEP
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*(2*gh+1) + m * im*(2*gh+1)*5

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    ja(current +1)  = m
    jac(current +1) = resd(i , j, e) / vol(i,j)
    ! jac(current +1) = 1.d0

  enddo
  enddo
  enddo
end subroutine computejacobian_twall_restrict

subroutine computejacobian_twall_restrict_notdivbyvol(jac,ia,ja,resd,m,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, 2*gh+1
!DIR$ IVDEP
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*(2*gh+1) + m * im*(2*gh+1)*5

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    ja(current +1)  = m
    jac(current +1) = resd(i , j, e) 
    ! jac(current +1) = 1.d0

  enddo
  enddo
  enddo
end subroutine computejacobian_twall_restrict_notdivbyvol

subroutine computejacobianfromjv(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = 5*im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        jac(current +1) = - resd(i , j, e)
        ! jac(current +1) = 1.d0
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv

subroutine testvector(wd,m,l,k,gh,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,jm,gh,m,l,k
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: wd
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j  
  stepj = 1+2*gh
  wd = 0.d0
  do j = k+1, jm, stepj
!DIR$ IVDEP  
  do i = l+1, im, stepj
    wd(i, j, m+1) = 1.d0
  enddo
  enddo
end subroutine testvector

subroutine testvector_unity(wd,m,l,k,gh,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,jm,gh,m,l,k
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: wd
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j  
  stepj = 1+2*gh
  wd = 0.d0
  j = k+1
  i = l+1
  wd(i, j, m+1) = 1.d0
end subroutine testvector_unity

subroutine computejacobianfromjv_unity(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*jm + m * im*jm*5*jm*im

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    vali = l
    valj = k
    ja(current +1)  = m + valj * 5 + vali * jm*5
    jac(current +1) = - resd(i , j, e)
    ! jac(current +1) = 1.d0

  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_unity


!!!!! Subroutines close to the wall, offcentered stencil 

subroutine computejacobianfromjv_wall(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  do j = 1, gh
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    valj = k

    if (l <= gh) then
      if (i <= l+1 + gh) then
        vali = l 
      else
        vali = (i - (l+1) + gh) / stepj * stepj + l
      endif   
    else
      if (i <= l - gh) then
        vali = im+1
      else
        vali = (i - (l+1) + gh) / stepj * stepj + l
      endif
    endif    
    if (vali >= im) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else
      ja(current +1)  = m + valj * 5 + vali * jm*5
      jac(current +1) = - resd(i , j, e)
      ! jac(current +1) = 1.d0
    endif
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_wall

subroutine testvector_wall(wd,m,l,k,gh,im,jm)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,jm,gh,m,l,k
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: wd
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j  
  stepj = 1+2*gh
  wd = 0.d0
  do j = k+1, jm+1, stepj
!DIR$ IVDEP  
  do i = l+1, im+1, stepj
     wd(i, j, m+1) = 1.d0
  enddo
  enddo
end subroutine testvector_wall


subroutine computejacobianfromjv_relaxed(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry, coefdiag)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(im,jm),intent(in) :: coefdiag
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = 5*im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        if (ia(current +1) == ja(current +1)) then
          jac(current +1) = coefdiag(i,j) - resd(i , j, e)
        else
          jac(current +1) = - resd(i , j, e)  
        endif  
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_relaxed


subroutine computejacobianfromjv_relaxed_dbyvol(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry, coefdiag, vol)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(im,jm),intent(in) :: coefdiag
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = 5*im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        if (ia(current +1) == ja(current +1)) then
          jac(current +1) = (coefdiag(i,j) - resd(i , j, e)) / vol(i,j)
        else
          jac(current +1) = - resd(i , j, e) / vol(i,j)
        endif  
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_relaxed_dbyvol

subroutine computejacobianfromjv_dbyvol(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry, vol)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = 5*im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        jac(current +1) = - resd(i , j, e) / vol(i,j)
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_dbyvol

subroutine computejacobianfromdz(jac,ia,ja,dz,m,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: dz
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (k <= gh) then
      if (j <= k+1 + gh) then
        valj = k 
      else
        valj = (j - (k+1) + gh) / stepj * stepj + k
      endif   
    else
      if (j <= k - gh) then
        valj = jm+1
      else
        valj = (j - (k+1) + gh) / stepj * stepj + k
      endif
    endif  

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = 5*im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        jac(current +1) = dz(i , j, e)
        ! jac(current +1) = 1.d0
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromdz

subroutine computejacobianfromjv_complex(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  complex*16,dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  complex*16,dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -1
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im+1
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        ia(current +1)  = 5*im*jm -1
        ja(current +1)  = 0
        jac(current +1) = 0.d0
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        jac(current +1) = - resd(i , j, e)
        ! jac(current +1) = 1.d0
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_complex

!! computejacobianfromjv_relaxed_withjn ONLY join in i-direction for im multiple of stepj
subroutine computejacobianfromjv_relaxed_withjn(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry, coefdiag)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(im,jm),intent(in) :: coefdiag
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -2
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im-1 - 2*gh + l 
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        if (vali - im <= gh-1) then
          vali = l
          ja(current +1)  = m + valj * 5 + vali * jm*5
          if (ia(current +1) == ja(current +1)) then
            jac(current +1) = coefdiag(i,j) - resd(i , j, e)
          else
            jac(current +1) = - resd(i , j, e)  
            ! jac(current +1) = 1.d0
          endif  
        else
          ia(current +1)  = 5*im*jm -2
          ja(current +1)  = 0
          jac(current +1) = 0.d0
        endif  
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        if (ia(current +1) == ja(current +1)) then
          jac(current +1) = coefdiag(i,j) - resd(i , j, e)
        else
          jac(current +1) = - resd(i , j, e)  
          ! jac(current +1) = 1.d0 
        endif
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_relaxed_withjn

!! computejacobianfromjv_withjn ONLY join in i-direction for im multiple of stepj
subroutine computejacobianfromjv_withjn(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -2
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im-1 - 2*gh + l 
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        if (vali - im <= gh-1) then
          vali = l
          ja(current +1)  = m + valj * 5 + vali * jm*5
          jac(current +1) = - resd(i , j, e)  
          ! jac(current +1) = 1.d0
        else
          ia(current +1)  = 5*im*jm -2
          ja(current +1)  = 0
          jac(current +1) = 0.d0
        endif  
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        jac(current +1) = - resd(i , j, e)  
        ! jac(current +1) = 1.d0 
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_withjn

!! computejacobianfromjv_withjn divided by the volume ONLY join in i-direction for im multiple of stepj
subroutine computejacobianfromjv_withjn_dbyvol(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry, vol)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh        ),intent(in) :: vol
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -2
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif   
      else
        if (i <= l - gh) then
          vali = im-1 - 2*gh + l 
          ! vali = (i - (l+1) + gh) / stepj * stepj + l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif    
      if (vali >= im) then
        if (vali - im <= gh-1) then
          vali = l
          ja(current +1)  = m + valj * 5 + vali * jm*5
          jac(current +1) = - resd(i , j, e) / vol(i,j) 
          ! jac(current +1) = 1.d0  
        else
          ia(current +1)  = 5*im*jm -2
          ja(current +1)  = 0
          jac(current +1) = 0.d0
        endif  
      else
        ja(current +1)  = m + valj * 5 + vali * jm*5
        jac(current +1) = - resd(i , j, e) / vol(i,j) 
        ! jac(current +1) = 1.d0 
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_withjn_dbyvol

subroutine testvector_partial(wd,m,l,k,gh,im,jm,istart,iend,jstart,jend)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in):: im,jm,gh,m,l,k,istart,iend,jstart,jend
  ! Returned objects ------------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(inout) :: wd
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j  
  stepj = 1+2*gh
  wd = 0.d0
  do j = jstart + k+1, jend+1, stepj
!DIR$ IVDEP  
  do i = istart + l+1, iend+1, stepj
    wd(i, j, m+1) = 1.d0
  enddo
  enddo
end subroutine testvector_partial

!! computejacobianfromjv_relaxed_withjn ONLY join in i-direction for two zones
subroutine computejacobianfromjv_relaxed_withjnandcheck(jac,ia,ja,resd,m,l,k,gh,im,jm,nbentry, coefdiag,mini,n)
  !
  implicit none
  ! variables for dimension -----------------------------------------
  integer, intent(in) :: im,jm,gh,nbentry,m,l,k,n
  real(8), intent(in) :: mini
  ! required arguments ----------------------------------------------
  real(8),dimension(1-gh:im+gh  ,1-gh:jm+gh  , 5   ),intent(in) :: resd
  real(8),dimension(im,jm),intent(in) :: coefdiag
  ! Returned objects ------------------------------------------------
  integer,dimension(nbentry),intent(inout) :: ia
  integer,dimension(nbentry),intent(inout) :: ja
  real(8),dimension(nbentry),intent(inout) :: jac
  ! Non-required arguments -------------------------------------------
  integer :: stepj,i,j,e,ii,jj, current, vali, valj
  stepj = 1+2*gh
  do e = 1,5
  ! do j = 1 + gh, jm
  do j = 1, jm
!DIR$ IVDEP  
  do i = 1, im
    !
    current = (i-1) + (j-1) * im + (e-1) * im*jm + k * im*jm*5 + l * im*jm*5*stepj + m * im*jm*5*stepj**2 + n*im*jm*5*5*stepj**2

    ia(current +1)  = e-1 + (j-1) * 5 + (i-1) * jm*5

    if (j <= k+1 + gh) then
      valj = k
    else
      valj = (j-gh - (k+1) + 2*gh) / stepj * stepj + k
    endif

    if (valj >= jm) then
      ia(current +1)  = 5*im*jm -2
      ja(current +1)  = 0
      jac(current +1) = 0.d0
    else  
      if (l <= gh) then
        if (i <= l+1 + gh) then
          vali = l 
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      else
        if (i <= l - gh) then
          vali = l
        else
          vali = (i - (l+1) + gh) / stepj * stepj + l
        endif
      endif
      if (n == 0) then
        if (vali >= im-2*gh) then
          vali = l
          if (ABS(jac(current +1)) < mini ) then
            ja(current +1)  = m + valj * 5 + vali * jm*5
            if (ia(current +1) == ja(current +1)) then
              jac(current +1) = coefdiag(i,j) - resd(i , j, e)
            else
              jac(current +1) = - resd(i , j, e)
            endif 
          else
            ia(current +1)  = 5*im*jm -2
            ja(current +1)  = 0
            jac(current +1) = 0.d0
          endif
        else
          if (ABS(jac(current +1)) < mini ) then
            ja(current +1)  = m + valj * 5 + vali * jm*5
            if (ia(current +1) == ja(current +1)) then
              jac(current +1) = coefdiag(i,j) - resd(i , j, e)
            else
              jac(current +1) = - resd(i , j, e)
            endif
          endif
        endif
      else
        if (vali <= 2*gh) then
          vali = im - MOD(im, stepj) + l
          if (vali > im-1) then
            vali = vali - stepj
          endif
          if (ABS(jac(current +1)) < mini ) then
            ja(current +1)  = m + valj * 5 + vali * jm*5
            if (ia(current +1) == ja(current +1)) then
              jac(current +1) = coefdiag(i,j) - resd(i , j, e)
            else
              jac(current +1) = - resd(i , j, e)
            endif
          endif
        else
          if (ABS(jac(current +1)) < mini ) then
            ja(current +1)  = m + valj * 5 + vali * jm*5
            if (ia(current +1) == ja(current +1)) then
              jac(current +1) = coefdiag(i,j) - resd(i , j, e)
            else
              jac(current +1) = - resd(i , j, e)  
            endif
          else
            ia(current +1)  = 5*im*jm -2
            ja(current +1)  = 0
            jac(current +1) = 0.d0    
          endif
        endif
      endif  
    endif    
    
  enddo
  enddo
  enddo
end subroutine computejacobianfromjv_relaxed_withjnandcheck
