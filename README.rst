See BROADCAST full documentation at ...

#. Compilation in the root folder must be done to run the code.
#. 'srcfv' folder includes BROADCAST source -> Compilation inside this folder must be done to run the code.
#. 'misc' folder includes Jacobian construction and all PETSc functions -> Compilation inside this folder must be done to run the code.
#. 'ADFirstAidKit' folder to compute the adjoint through TAPENADE software -> unnecessary for BROADCAST.
#. 'SIM' folder includes the self-similar profile generator program -> useful to initialise boundary layer simulations.

Boundary Layer case files:
==========================

To compute the baseflow of a boundary layer
-------------------------------------------

* card_bl2d_fv_npz.py -> Python card to run BROADCAST code, call BROADCAST_npz.py
-> Must be fill in by the user (Numerical param., BROADCAST options...).

* BROADCAST_npz.py -> main program, call the functions compiled in 'srcfv' and 'misc' folders
-> Must be fill in by the user (Mesh, Adim. choice, BC location and type...).

* BROADCAST_split_func.py -> equal to BROADCAST.py, only different organisation.

* meshBL.py -> functions to build different cartesian stretched mesh.

* restart_init.py -> functions to read previous solution to restart simulation (only for cartesian rectangular grids).

* interpgrid.py -> 1st order interpolation functions (only for cartesian rectangular grids).

* compilef90.py -> compile initialisation.f90 and set_bnd.f90.

* initialisation.f90 -> approximated initialisation of blasius profile.

* set_bnd.f90 -> functions which fill the ref. state inside the BC.

For resolvent (global stability) analysis of a boundary layer
-------------------------------------------------------

* DzMatrix.py -> Compute the transverse (z-direction) contributions to the Jacobian for 3D stability.

* computeBLthickness -> Compute different boundary layer thicknesses, it can be used to restrict the resolvent (only for cartesian rectangular grids).

* resolvent_all.py -> main functions to perform 2D resolvent analysis, call add. functions inside misc/PETSc_func.py.

* resolvent_all3D_1block.py -> main functions to perform 3D resolvent analysis.

* resolvent_all3D_1block_control.py -> main functions to perform 3D resolvent analysis with a forcing constrained at the wall.

Cylinder case files:
==========================

To compute the baseflow of a cylinder
-------------------------------------------------------

* card_cyl2d.py -> Python card to run BROADCAST code for cylinder, call cylinder.py
-> Must be fill in by the user (Numerical param., BROADCAST options...).

* cylinder.py -> main program, call the functions compiled in 'srcfv' and 'misc' folders
-> Must be fill in by the user (Mesh, Adim. choice, BC location and type...).

* meshCyl.py -> functions to build O-mesh for cylinder.

For global stability analysis of a cylinder
-------------------------------------------------------

* DzMatrix_cyl.py -> Compute the transverse (z-direction) contributions to the Jacobian for 3D stability.

* biglobal_cyl.py -> main functions to perform 2D/3D (bi)global stability analysis, call add. functions inside MISC/PETSc_func.py.

To compute eigenvalue sensitivity of a cylinder
-------------------------------------------------------

* Product_Adjoint.py -> Normalise adjoint and direct modes with L2 norm or any other user's instructions.

* Hessian_cyl.py -> Compute the product of the Hessian operator with a user-provided mode.

* ProductHessian_cyl.py -> Apply the result returned by Hessian_cyl.py to the adjoint mode to compute the eigenvalue sensitivity.

To compute the coefficients of the Weakly Nonlinear Stability analysis (Sipp, Lebedev JFM 2007)
-------------------------------------------------------

* limitcycle_part_all.py -> compute the modes x22 & x20 and the coefficients \mu & \nu.

Guidelines for Weakly Nonlinear Stability coefficients computation:

#. Compute direct mode with biglobal_cyl.py.
#. Compute adjoint mode with biglobal_cyl.py.
#. Normalise direct & adjoint modes with Product_Adjoint.py.
#. Compute Hessian associated with the direct mode with Hessian_cyl.py.
#. Compute modes x22 & x20 and the coefficients \mu & \nu with limitcycle_part_all.py.



