.. _srcfvfiles:

FORTRAN files
===========

* geom -> metrics computation.

* phys -> primitives and viscosity computations.

* borders -> boundary conditions.

* lhs -> basic matrix-free inversion for implicit time-stepping method.

* rhs -> spatial numerical scheme for residual computation.

* prepro -> preprocessed files for Algorithmic Diffenrentiation (AD).

* tangent -> linearised files computed by AD.

* adjoint -> linearised files computed by AD in backward mode.

* tangenttangentHess -> twice linearised files computed by AD (Hessian).

* dz -> transverse contributions of the Jacobian operator.

* tangentdz -> linearised transverse contributions.

* misc/ComputeJacobian.f90 -> functions to construct the Jacobian (test-vectors and index ordering).

