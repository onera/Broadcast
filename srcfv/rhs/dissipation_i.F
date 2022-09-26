#include "rhs/spectralradius_i.F"
! 
!     if (eps4.gt.1.d-12) then  
! #include "rhs/wiggle_diri.F"
!     endif
    
    dissro1  = rspec * eps4*predro1 
    dissrou1 = rspec * eps4*predrou1
    dissrov1 = rspec * eps4*predrov1
    dissrow1 = rspec * eps4*predrow1
    dissroe1 = rspec * eps4*predroe1
    

