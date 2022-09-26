    
    ! Wiggle IVAN MARY pour enlver la dissipation dir i
    deltar(1:5) = (w(i-1,j,1:5) - w(i-2,j,1:5))*(w(i,j,1:5) - w(i-1,j,1:5))
    !
    deltal(1:5) = (w(i,j,1:5) - w(i-1,j,1:5))*(w(i+1,j,1:5) - w(i,j,1:5))
    !
    t(:) = sign(ONE,deltar(:)) + sign(ONE,deltal(:))
    !
    test = minval(t)
    !
    if(test.lt.-1e-12)then
        wiggle = ONE
    else
        wiggle = ZERO
    end if
    !
    eps4 = eps4*wiggle
    
    !
    
    
