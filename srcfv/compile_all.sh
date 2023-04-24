#!/bin/bash

mpirun -np 1 python3 compile_lhs.py &
mpirun -np 1 python3 compile_borders.py &
mpirun -np 1 python3 compile_geom.py & 
mpirun -np 1 python3 compile_norm.py &
#mpirun -np 1 python3 compile_dz_master.py &
mpirun -np 1 python3 compile_dz.py &
#mpirun -np 1 python3 compile_tangentdz.py &
#mpirun -np 1 python3 compile_dz_axi_master.py &
mpirun -np 1 python3 compile_dz_axi.py &
mpirun -np 1 python3 compile_rhs.py &
mpirun -np 1 python3 compile_tangent.py &
mpirun -np 1 python3 compile_hessian.py &
mpirun -np 1 python3 compile_adj.py &
mpirun -np 1 python3 compile_chu.py &
mpirun -np 1 python3 compile_tangentchu.py 
#mpirun -np 1 python3 compile_costf.py &
#mpirun -np 1 python3 compile_tangentcostf.py &
#mpirun -np 1 python3 compile_poisson.py

