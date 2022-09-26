#include "rhs/spectralradius_i.F"
#include "rhs/spectralradius_j.F"
rspecbis = rspec
rspec = MAX(rspec, rspecbis)
sqbis = sq
sq = MAX(sq, sqbis)
#include "rhs/ducrosfordnc_maxSC.F"


    
    dissro1  = rspec * (eps2*diffro  + eps4*predro1 )
    dissrou1 = rspec * (eps2*diffrou + eps4*predrou1)
    dissrov1 = rspec * (eps2*diffrov + eps4*predrov1)
    dissrow1 = rspec * (eps2*diffrow + eps4*predrow1)
    dissroe1 = rspec * (eps2*diffroe + eps4*predroe1)
    
    diffro   = ( w(i,j,1) - w(i,j-1,1) )
    diffrou  = ( w(i,j,2) - w(i,j-1,2) )
    diffrov  = ( w(i,j,3) - w(i,j-1,3) )
    diffrow  = ( w(i,j,4) - w(i,j-1,4) )
    diffroe  = ( w(i,j,5) - w(i,j-1,5) )

    dissro2  = rspec * (eps2*diffro  + eps4*predro2 )
    dissrou2 = rspec * (eps2*diffrou + eps4*predrou2)
    dissrov2 = rspec * (eps2*diffrov + eps4*predrov2)
    dissrow2 = rspec * (eps2*diffrow + eps4*predrow2)
    dissroe2 = rspec * (eps2*diffroe + eps4*predroe2)