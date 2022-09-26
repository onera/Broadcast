
  dz2_out(i,j,1) = coeffs(i,j,1,1)  * func0d(i, j,  1)  + coeffs(i,j,1,2)  * func1d(i, j,  1) &
&                + coeffs(i,j,1,3)  * func2d(i, j,  1)  + coeffs(i,j,1,4)  * func3d(i, j,  1)

  dz2_out(i,j,2) = coeffs(i,j,2,1)  * func0d(i, j,  2)  + coeffs(i,j,2,2)  * func1d(i, j,  2) &
&                + coeffs(i,j,2,3)  * func2d(i, j,  2)  + coeffs(i,j,2,4)  * func3d(i, j,  2)

  dz2_out(i,j,3) = coeffs(i,j,3,1)  * func0d(i, j,  3)  + coeffs(i,j,3,2)  * func1d(i, j,  3) &
&                + coeffs(i,j,3,3)  * func2d(i, j,  3)  + coeffs(i,j,3,4)  * func3d(i, j,  3)

  dz2_out(i,j,4) = coeffs(i,j,4,1)  * func0d(i, j,  4)  + coeffs(i,j,4,2)  * func1d(i, j,  4) &
&                + coeffs(i,j,4,3)  * func2d(i, j,  4)  + coeffs(i,j,4,4)  * func3d(i, j,  4)

  dz2_out(i,j,5) = coeffs(i,j,5,1)  * func0d(i, j,  5)  + coeffs(i,j,5,2)  * func1d(i, j,  5) &
&                + coeffs(i,j,5,3)  * func2d(i, j,  5)  + coeffs(i,j,5,4)  * func3d(i, j,  5)


