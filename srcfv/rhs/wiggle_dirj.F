    
    ! Wiggle IVAN MARY pour enlver la dissipation dir i
    deltar(1:5) = (w(i,j-1,1:5) - w(i,j-2,1:5))*(w(i,j,1:5) - w(i,j-1,1:5))
    !
    deltal(1:5) = (w(i,j,1:5) - w(i,j-1,1:5))*(w(i,j+1,1:5) - w(i,j,1:5))
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
    



