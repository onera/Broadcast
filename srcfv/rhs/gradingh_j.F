  do h=1,gh
    grad_rho( i, 1-h,:) =  TWO*grad_rho( i,   2-h,:) - grad_rho( i,   3-h,:)
    gradu(    i, 1-h,:) =  TWO*gradu(    i,   2-h,:) - gradu(    i,   3-h,:)
    gradv(    i, 1-h,:) =  TWO*gradv(    i,   2-h,:) - gradv(    i,   3-h,:)
    gradw(    i, 1-h,:) =  TWO*gradw(    i,   2-h,:) - gradw(    i,   3-h,:)
    grad_temp(i, 1-h,:) =  TWO*grad_temp(i,   2-h,:) - grad_temp(i,   3-h,:)
  
    grad_rho( i,jm+h,:) =  TWO*grad_rho( i,jm-1+h,:) - grad_rho( i,jm-2+h,:)
    gradu(    i,jm+h,:) =  TWO*gradu(    i,jm-1+h,:) - gradu(    i,jm-2+h,:)
    gradv(    i,jm+h,:) =  TWO*gradv(    i,jm-1+h,:) - gradv(    i,jm-2+h,:)
    gradw(    i,jm+h,:) =  TWO*gradw(    i,jm-1+h,:) - gradw(    i,jm-2+h,:)
    grad_temp(i,jm+h,:) =  TWO*grad_temp(i,jm-1+h,:) - grad_temp(i,jm-2+h,:)
  enddo
  
