#include "geom/dxdy.F"

       gradu( i,j,1)  = dxm1*gui   + dxm2*guj
       gradv( i,j,1)  = dxm1*gvi   + dxm2*gvj
       gradw( i,j,1)  = dxm1*gwi   + dxm2*gwj 
       gradmu(i,j,1)  = dxm1*gmui  + dxm2*gmuj
                         
       gradu( i,j,2)  = dym1*gui   + dym2*guj
       gradv( i,j,2)  = dym1*gvi   + dym2*gvj
       gradw( i,j,2)  = dym1*gwi   + dym2*gwj 
       gradmu(i,j,2)  = dym1*gmui  + dym2*gmuj
                         
                         
