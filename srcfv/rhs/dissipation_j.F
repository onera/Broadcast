#include "rhs/spectralradius_j.F"
!   
!     if (eps4.gt.1.d-12) then  
! #include "rhs/wiggle_dirj.F"
!     endif  
    !
    dissro2  = rspec * eps4*predro2 
    dissrou2 = rspec * eps4*predrou2
    dissrov2 = rspec * eps4*predrov2
    dissrow2 = rspec * eps4*predrow2
    dissroe2 = rspec * eps4*predroe2
