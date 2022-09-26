#include "rhs/spectralradius_i.F"
! correction Rossow JCP 2000
! #include "rhs/spectralradiusRossow_i.F"
#include "rhs/ducrosfordnc_i.F"
! 
!     if (eps4.gt.1.d-12) then  
! #include "rhs/wiggle_diri.F"
!     endif
    
    dissro1  = rspec * (eps2*diffro  + eps4*predro1 )
    dissrou1 = rspec * (eps2*diffrou + eps4*predrou1)
    dissrov1 = rspec * (eps2*diffrov + eps4*predrov1)
    dissrow1 = rspec * (eps2*diffrow + eps4*predrow1)
    dissroe1 = rspec * (eps2*diffroe + eps4*predroe1)
    

