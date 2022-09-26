!$AD II-LOOP
      do j=1-gh,jm+gh
!$AD II-LOOP
!GCC$ IVDEP      
      do i=1-gh,im+gh
#include "phys/Primitives.F"
#include "phys/viscosity.F"
      enddo
      enddo
  
