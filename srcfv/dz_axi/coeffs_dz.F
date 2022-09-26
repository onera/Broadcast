
    distm1    = ONE/abs( yc(i,j) - ym )

    coeffs(i,j,1,1) = ONE * distm1
    coeffs(i,j,2,1) = ONE * distm1
    coeffs(i,j,3,1) = ONE * distm1
    coeffs(i,j,4,1) = ONE * distm1
    coeffs(i,j,5,1) = ONE * distm1

    coeffs(i,j,1,2) = ZERO
    coeffs(i,j,2,2) = TWOTHIRD * gradmu(i,j,1) * distm1
    coeffs(i,j,3,2) = TWOTHIRD * gradmu(i,j,2) * distm1 + TWO * TWOTHIRD * mu(i,j) * distm1**2
    coeffs(i,j,4,2) = - gradmu(i,j,1) * distm1
    coeffs(i,j,5,2) = TWOTHIRD * ( gradmu(i,j,1) * velx(i,j) + mu(i,j) * gradu(i,j,1) ) * distm1

    coeffs(i,j,1,3) = ZERO
    coeffs(i,j,2,3) = TWOTHIRD * mu(i,j) * distm1
    coeffs(i,j,3,3) = TWOTHIRD * mu(i,j) * distm1
    coeffs(i,j,4,3) = - gradmu(i,j,2) * distm1 - mu(i,j) * distm1**2
    coeffs(i,j,5,3) = - 1.d0/3.d0 * mu(i,j) * velx(i,j) * distm1

    coeffs(i,j,1,4) = ZERO
    coeffs(i,j,2,4) = ZERO
    coeffs(i,j,3,4) = ZERO
    coeffs(i,j,4,4) = - mu(i,j) * distm1
    coeffs(i,j,5,4) = - ( gradmu(i,j,1) * velz(i,j) + mu(i,j) * gradw(i,j,1) ) * distm1

    coeffs(i,j,1,5) = ZERO
    coeffs(i,j,2,5) = ZERO
    coeffs(i,j,3,5) = ZERO
    coeffs(i,j,4,5) = ZERO
    coeffs(i,j,5,5) = - mu(i,j) * velz(i,j) * distm1

    coeffs(i,j,1,6) = ZERO
    coeffs(i,j,2,6) = ZERO
    coeffs(i,j,3,6) = ZERO
    coeffs(i,j,4,6) = ZERO
    coeffs(i,j,5,6) = TWOTHIRD * ( gradmu(i,j,2) * vely(i,j) + mu(i,j) * gradv(i,j,2) ) * distm1

    coeffs(i,j,1,7) = ZERO
    coeffs(i,j,2,7) = ZERO
    coeffs(i,j,3,7) = ZERO
    coeffs(i,j,4,7) = ZERO
    coeffs(i,j,5,7) = TWOTHIRD * mu(i,j) * vely(i,j) * distm1

    coeffs(i,j,1,8) = ZERO
    coeffs(i,j,2,8) = ZERO
    coeffs(i,j,3,8) = ZERO
    coeffs(i,j,4,8) = ZERO
    coeffs(i,j,5,8) = - ( gradmu(i,j,2) * velz(i,j) + mu(i,j) * gradw(i,j,2) ) * distm1

    coeffs(i,j,1,9) = ZERO
    coeffs(i,j,2,9) = ZERO
    coeffs(i,j,3,9) = ZERO
    coeffs(i,j,4,9) = ZERO
    coeffs(i,j,5,9) = - mu(i,j) * velz(i,j) * distm1

    coeffs(i,j,1,10) = ZERO
    coeffs(i,j,2,10) = ZERO
    coeffs(i,j,3,10) = ZERO
    coeffs(i,j,4,10) = ZERO
    coeffs(i,j,5,10) = - gradw(i,j,1) * distm1

    coeffs(i,j,1,11) = ZERO
    coeffs(i,j,2,11) = ZERO
    coeffs(i,j,3,11) = ZERO
    coeffs(i,j,4,11) = ZERO
    coeffs(i,j,5,11) = - ( gradw(i,j,2) - velz(i,j) * distm1 ) * distm1

    coeffs(i,j,1,12) = ZERO
    coeffs(i,j,2,12) = ZERO
    coeffs(i,j,3,12) = ZERO
    coeffs(i,j,4,12) = ZERO
    coeffs(i,j,5,12) = - mu(i,j) * vely(i,j) * distm1

    coeffs(i,j,1,13) = ZERO
    coeffs(i,j,2,13) = ZERO
    coeffs(i,j,3,13) = ZERO
    coeffs(i,j,4,13) = ZERO
    coeffs(i,j,5,13) = TWOTHIRD * ( gradu(i,j,1) + gradv(i,j,2) - TWO * vely(i,j) * distm1 ) * distm1

    coeffs(i,j,1,14) = ZERO
    coeffs(i,j,2,14) = ZERO
    coeffs(i,j,3,14) = ZERO
    coeffs(i,j,4,14) = ZERO
    coeffs(i,j,5,14) = TWOTHIRD * mu(i,j) * velz(i,j) * distm1

    coeffs(i,j,1,15) = ZERO
    coeffs(i,j,2,15) = ZERO
    coeffs(i,j,3,15) = ZERO
    coeffs(i,j,4,15) = ZERO
    coeffs(i,j,5,15) = ZERO

    coeffs(i,j,1,16) = ZERO
    coeffs(i,j,2,16) = ZERO
    coeffs(i,j,3,16) = ZERO
    coeffs(i,j,4,16) = ZERO
    coeffs(i,j,5,16) = ZERO


