#include "rhs/spectralradius_j.F"
! correction Rossow JCP 2000
! #include "rhs/spectralradiusRossow_j.F"
#include "rhs/ducrosfordnc_j.F"
!   
!     if (eps4.gt.1.d-12) then  
! #include "rhs/wiggle_dirj.F"
!     endif
    
    !
    dissro2  = rspec * (eps2*diffro  + eps4*predro2 )
    dissrou2 = rspec * (eps2*diffrou + eps4*predrou2)
    dissrov2 = rspec * (eps2*diffrov + eps4*predrov2)
    dissrow2 = rspec * (eps2*diffrow + eps4*predrow2)
    dissroe2 = rspec * (eps2*diffroe + eps4*predroe2)
