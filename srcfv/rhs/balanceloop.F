
!$AD II-LOOP
  do j = 1 , jm
!$AD II-LOOP
!GCC$ IVDEP      
  do i = 1 , im
#include "rhs/balance.F"
  enddo
  enddo
      
