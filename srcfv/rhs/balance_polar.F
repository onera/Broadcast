!$AD II-LOOP
  do j = 1 , jm
!$AD II-LOOP
!DIR$ IVDEP      
  do i = 1 , im
#include "rhs/polarE.F"

residu(i  ,j,1) =   src1                              &
                  - ( hn(i+1,j  ,1,1) - hn(i,j,1,1) ) &
                  - ( hn(i  ,j+1,1,2) - hn(i,j,1,2) ) 
                  
residu(i  ,j,2) =   src2                              &
                  - ( hn(i+1,j  ,2,1) - hn(i,j,2,1) ) &
                  - ( hn(i  ,j+1,2,2) - hn(i,j,2,2) ) 

residu(i  ,j,3) =   src3                              &
                  - ( hn(i+1,j  ,3,1) - hn(i,j,3,1) ) &
                  - ( hn(i  ,j+1,3,2) - hn(i,j,3,2) )

residu(i  ,j,4) =   src4                              &
                  - ( hn(i+1,j  ,4,1) - hn(i,j,4,1) ) &
                  - ( hn(i  ,j+1,4,2) - hn(i,j,4,2) ) 

residu(i  ,j,5) =   src5                              &
                  - ( hn(i+1,j  ,5,1) - hn(i,j,5,1) ) &
                  - ( hn(i  ,j+1,5,2) - hn(i,j,5,2) ) 

  enddo
  enddo
      
