
    dist =  vol(i,j) / (sqrt((nx(i,j,1)*nx(i,j,1)+ny(i,j,1)*ny(i,j,1)))+ &
                        sqrt((nx(i,j,2)*nx(i,j,2)+ny(i,j,2)*ny(i,j,2)))  )
    

    cc    = sqrt(gam*pint*rom1)
    
    vitc  = sqrt(velxi*velxi+velyi*velyi) + cc
    
    dt_euler  = dist/vitc
    dt_ns     = 2.d0 * dist*dist*ror/(gampr*mu(i,j))
               ! 2.d0 is linked to the expression of gampr and gamprt
    
    
    dtm1      = cflm1/max(min(dt_euler,dt_ns),1.d-15) ! 1/dt
    dt        = max(min(dt_euler,dt_ns),1.d-15)*cfl ! 1/dt

    !write(30101,*) "i,j    =  ", i,j
    !write(30101,*) "dtm1   =  ", dtm1
    !write(30101,*) "distm1 =  ", distm1
    !write(30101,*) "cflm1  =  ", dtm1
    !write(30101,*) "mu     =  ", mu(i,j)

    
    
    
    
    
    
    
