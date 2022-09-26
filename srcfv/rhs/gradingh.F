  do i= 1, 2
  do h=1,gh
    grad_rho( :,1-h,:) =  TWO*grad_rho( :,2-h,:) - grad_rho( :,3-h,:)
    gradu(    :,1-h,:) =  TWO*gradu(    :,2-h,:) - gradu(    :,3-h,:)
    gradv(    :,1-h,:) =  TWO*gradv(    :,2-h,:) - gradv(    :,3-h,:)
    gradw(    :,1-h,:) =  TWO*gradw(    :,2-h,:) - gradw(    :,3-h,:)
    gradw(    :,1-h,:) =  TWO*gradw(    :,2-h,:) - gradw(    :,3-h,:)
    grad_temp(:,1-h,:) =  TWO*grad_temp(:,2-h,:) - grad_temp(:,3-h,:)
  
    grad_rho( :,jm+h,:) = TWO*grad_rho( :,jm-1+h,:) - grad_rho( :,jm-2+h,:)
    gradu(    :,jm+h,:) = TWO*gradu(    :,jm-1+h,:) - gradu(    :,jm-2+h,:)
    gradv(    :,jm+h,:) = TWO*gradv(    :,jm-1+h,:) - gradv(    :,jm-2+h,:)
    gradw(    :,jm+h,:) = TWO*gradw(    :,jm-1+h,:) - gradw(    :,jm-2+h,:)
    grad_temp(:,jm+h,:) = TWO*grad_temp(:,jm-1+h,:) - grad_temp(:,jm-2+h,:)
  
    grad_rho( 1-h,:,:) =  TWO*grad_rho( 2-h,:,:) - grad_rho( 3-h,:,:)
    gradu(    1-h,:,:) =  TWO*gradu(    2-h,:,:) - gradu(    3-h,:,:)
    gradv(    1-h,:,:) =  TWO*gradv(    2-h,:,:) - gradv(    3-h,:,:)
    gradw(    1-h,:,:) =  TWO*gradw(    2-h,:,:) - gradw(    3-h,:,:)
    grad_temp(1-h,:,:) =  TWO*grad_temp(2-h,:,:) - grad_temp(3-h,:,:)
  
    grad_rho( im+h,:,:) = TWO*grad_rho( im-1+h,:,:) - grad_rho( im-2+h,:,:)
    gradu(    im+h,:,:) = TWO*gradu(    im-1+h,:,:) - gradu(    im-2+h,:,:)
    gradv(    im+h,:,:) = TWO*gradv(    im-1+h,:,:) - gradv(    im-2+h,:,:)
    gradw(    im+h,:,:) = TWO*gradw(    im-1+h,:,:) - gradw(    im-2+h,:,:)
    grad_temp(im+h,:,:) = TWO*grad_temp(im-1+h,:,:) - grad_temp(im-2+h,:,:)
  enddo
  enddo
