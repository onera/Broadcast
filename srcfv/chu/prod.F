
a1 = gam1*gam*mach*mach*w(i,j,1)/tloc(i,j)
a2 = (HALF*(velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j))-cv*tloc(i,j))/w(i,j,1)

product(i,j,1) = HALF * vol(i,j) * (((velx(i,j)*velx(i,j)+vely(i,j)*vely(i,j))/w(i,j,1)+tloc(i,j)/(w(i,j,1)*gam*mach*mach) &
                                    +a1*a2*a2)*wf(i,j,1) &
                                    - velx(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,2)   &
                                    - vely(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,3)    &
                                    + ZERO*wf(i,j,4)   &
                                    + a1*a2/w(i,j,1)*wf(i,j,5) )

product(i,j,2) = HALF * vol(i,j) * (- velx(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,1) &
                                    + (ONE/w(i,j,1)+velx(i,j)*velx(i,j)*a1/(w(i,j,1)*w(i,j,1)))*wf(i,j,2)   &
                                    + velx(i,j)*vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,3)    &
                                    + ZERO*wf(i,j,4)   &
                                    - velx(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,5) )

product(i,j,3) = HALF * vol(i,j) * (- vely(i,j)*(ONE+a1*a2)/w(i,j,1)*wf(i,j,1) &
                                    + velx(i,j)*vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,2)    &
                                    + (ONE/w(i,j,1)+vely(i,j)*vely(i,j)*a1/(w(i,j,1)*w(i,j,1)))*wf(i,j,3)   & 
                                    + ZERO*wf(i,j,4)   &
                                    - vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,5) )

product(i,j,4) = HALF * vol(i,j) * (ZERO*wf(i,j,1) &
                                    + ZERO*wf(i,j,2)   &
                                    + ZERO*wf(i,j,3)    &
                                    + ONE/w(i,j,1)*wf(i,j,4)   &
                                    + ZERO*wf(i,j,5)  )

product(i,j,5) = HALF * vol(i,j) * (a1*a2/w(i,j,1)*wf(i,j,1) &
                                    - velx(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,2)   &
                                    - vely(i,j)*a1/(w(i,j,1)*w(i,j,1))*wf(i,j,3)   &
                                    + ZERO*wf(i,j,4)   &
                                    + a1/(w(i,j,1)*w(i,j,1))*wf(i,j,5) )