    !2st direction
    
                 
    rhomr     = w(i,j,1)
    ur        = w(i,j,2)/rhomr
    vr        = w(i,j,3)/rhomr
    c2r       = gam*rgaz*tloc(i,j)
    !
    rhoml     = w(i,j-1,1)
    ul        = w(i,j-1,2)/rhoml
    vl        = w(i,j-1,3)/rhoml
    c2l       = gam*rgaz*tloc(i,j-1)
    !
    r         = sqrt( rhomr/rhoml)
    rr        = ONE/(ONE+r)
    omrr      = ONE-rr
    !
    u         =  ul*rr + ur*omrr
    v         =  vl*rr + vr*omrr
    !
    c2x       = c2l*rr + c2r*omrr
    nx2       = nx(i,j,2)*nx(i,j,2)+ny(i,j,2)*ny(i,j,2)

    vn2       = (u*u+v*v)
    mach2     = min(vn2/c2x,ONE)
    alpha     = HALF*( ONE - mach2 )
    
    cprim     = alpha*alpha*vn2+mach2*c2x
    !
    ab        = abs(nx(i,j,2)*u+ny(i,j,2)*v)
    sq        = sqrt(cprim*nx2)
    !
    rspec     = ab + sq            
                            
    
    !
