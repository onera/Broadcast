
!$AD II-LOOP
  do j=1-gh, jm+gh
!$AD II-LOOP
!DIR$ IVDEP      
  do i=1-gh, im+gh
#include "phys/Primitives_var.F"
#include "phys/viscosity.F"
  enddo
  enddo  
    
  
