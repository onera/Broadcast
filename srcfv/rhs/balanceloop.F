
!$AD II-LOOP
  do j = 1 , jm
!$AD II-LOOP
!DIR$ IVDEP      
  do i = 1 , im
#include "rhs/balance.F"
  enddo
  enddo
      
