
        r2 = ONE
        
        !! Sponge layer upstream
        IF (Lmax - (xc(i,j)-xc(0,j)) .lt. ZERO) THEN
            ralpha = ZERO
        ELSE 
            ralpha = r2 * ((Lmax - (xc(i,j) -xc(0,j)))/Lmax)**r3
        ENDIF

        !! Sponge layer downstream
        ! IF (xc(i,j) - (xc(im,j) - Lmax) .lt. ZERO) THEN
        !     ralpha = ZERO
        ! ELSE 
        !     ralpha = r2 * ((xc(i,j) - (xc(im,j) - Lmax))/Lmax)**r3
        ! ENDIF

        source(:) = (r4 * ralpha * ( w(i,j,:) - wref(i,j,:) ) ) *vol(i,j)
        
        sourcear(i,j,:) = source(:)
