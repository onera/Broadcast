
        r0 = HALF
        r1 = -3.d0/8.d0  !-FOURTH
        ! r2 = 0.2d0
        ! r3 = 1.3d0
        ! r4 = ONE

        !! Sponge layer upstream
        ! IF (Lmax - xc(i,j) .lt. ZERO) THEN
        !     ralpha = ZERO
        ! ELSE 
        !     ralpha = r2 * ((Lmax - xc(i,j))/Lmax)**r3
        ! ENDIF

        !! Sponge layer downstream
        IF (xc(i,j) - (xc(im,j) - Lmax) .lt. ZERO) THEN
            ralpha = ZERO
        ELSE 
            ralpha = r2 * ((xc(i,j) - (xc(im,j) - Lmax))/Lmax)**r3
        ENDIF

        wmoy(:) = w(i,j,:) - ralpha * (3.d0 * r0 * w(i,j,:) + r1 * (w(i+1,j,:) + w(i-1,j,:) + w(i,j+1,:) + w(i,j-1,:)) )

        source(:) = (r4 * ralpha * ( w(i,j,:) - wmoy(:) ) ) *vol(i,j)
        
        sourcear(i,j,:) = source(:)
