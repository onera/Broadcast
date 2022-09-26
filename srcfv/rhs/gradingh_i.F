  do h=1,gh
    grad_rho( 1-h,j,:) =  TWO*grad_rho( 2-h,j,:) - grad_rho( 3-h,j,:)
    gradu(    1-h,j,:) =  TWO*gradu(    2-h,j,:) - gradu(    3-h,j,:)
    gradv(    1-h,j,:) =  TWO*gradv(    2-h,j,:) - gradv(    3-h,j,:)
    gradw(    1-h,j,:) =  TWO*gradw(    2-h,j,:) - gradw(    3-h,j,:)
    grad_temp(1-h,j,:) =  TWO*grad_temp(2-h,j,:) - grad_temp(3-h,j,:)
  
    grad_rho( im+h,j,:) = TWO*grad_rho( im-1+h,j,:) - grad_rho( im-2+h,j,:)
    gradu(    im+h,j,:) = TWO*gradu(    im-1+h,j,:) - gradu(    im-2+h,j,:)
    gradv(    im+h,j,:) = TWO*gradv(    im-1+h,j,:) - gradv(    im-2+h,j,:)
    gradw(    im+h,j,:) = TWO*gradw(    im-1+h,j,:) - gradw(    im-2+h,j,:)
    grad_temp(im+h,j,:) = TWO*grad_temp(im-1+h,j,:) - grad_temp(im-2+h,j,:)
  enddo
  
