

	coeffs(i,j,1,1) = ZERO
    coeffs(i,j,2,1) = - mu(i,j)
    coeffs(i,j,3,1) = - mu(i,j)
    coeffs(i,j,4,1) = - TWO * TWOTHIRD * mu(i,j)
    coeffs(i,j,5,1) = - mu(i,j) * cpprandtl

    coeffs(i,j,1,2) = ZERO
    coeffs(i,j,2,2) = ZERO
    coeffs(i,j,3,2) = ZERO
    coeffs(i,j,4,2) = ZERO
    coeffs(i,j,5,2) = - mu(i,j) * velx(i,j)

    coeffs(i,j,1,3) = ZERO
    coeffs(i,j,2,3) = ZERO
    coeffs(i,j,3,3) = ZERO
    coeffs(i,j,4,3) = ZERO
    coeffs(i,j,5,3) = - mu(i,j) * vely(i,j)

    coeffs(i,j,1,4) = ZERO
    coeffs(i,j,2,4) = ZERO
    coeffs(i,j,3,4) = ZERO
    coeffs(i,j,4,4) = ZERO
    coeffs(i,j,5,4) = - TWO * TWOTHIRD * mu(i,j) * velz(i,j)

