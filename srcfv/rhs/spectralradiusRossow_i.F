    !1st direction
    
                 
    rhomr     = w(i,j,1)
    ur        = w(i,j,2)/rhomr
    vr        = w(i,j,3)/rhomr
    c2r       = gam*rgaz*tloc(i,j)
    !
    rhoml     = w(i-1,j,1)
    ul        = w(i-1,j,2)/rhoml
    vl        = w(i-1,j,3)/rhoml
    c2l       = gam*rgaz*tloc(i-1,j)
    !
    r         = sqrt( rhomr/rhoml)
    rr        = ONE/(ONE+r)
    omrr      = ONE-rr
    !
    u         =  ul*rr + ur*omrr
    v         =  vl*rr + vr*omrr
    !
    c2x       = c2l*rr + c2r*omrr
    
    vn2       = (u*u+v*v)
    mach2     = min(vn2/c2x,ONE)
    alpha     = HALF*( ONE - mach2 )
    
    cprim     = alpha*alpha*vn2+mach2*c2x
    
    
    nx2       = nx(i,j,1)*nx(i,j,1)+ny(i,j,1)*ny(i,j,1)
    !
    ab        = abs(nx(i,j,1)*u+ny(i,j,1)*v)
    sq        = sqrt(cprim*nx2)
    !
    rspec     = ab + sq            
                            
    
    !
