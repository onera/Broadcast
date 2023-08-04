Methods
=======

Equations
------

:math:`\frac{dq}{dt} = R(q)`

with q the conservative state and R the residual of Navier-Stokes equations.

Numerical scheme
-----

Convective scheme
^^^^

The convective scheme implemented is the FE-MUSCL scheme, called DNC in the program files. This high-order scheme is available at order 3 (5 points stencil), 5 (7 points stencil), 7 (9 points stencil) and 9 (11 points stencil). Different versions of the scheme exist:

* :func:`flux_num_dnc3_nowall_2d` calls FE-MUSCL at order 3. It is the baseline scheme for cartesian equations.
* :func:`flux_num_dnc3_nowall_polar_2d` is the baseline scheme for axisymmetric equations.
* :func:`flux_num_dnc3_2d` is the scheme for cartesian equations with an adiabatic wall boundary at j=0 along all i-direction (flat plate case for instance). This scheme includes off-centered computation of gradients in order not to pick values far inside the wall. Results are very similar to those obtained with :func:`flux_num_dnc3_nowall_2d`.
* :func:`flux_num_dnc3_polar_2d` is the same as above for axisymmetric equations.


Viscous scheme
^^^^

The viscous scheme is a compact scheme at order 4 (5 point stencil) inside the domain and at order 2 (3 point stencil) near the wall boundary at j=0 if there is one.


Boundary conditions 
-----

* :func:`bc_general_2d`: Dirichlet BC (supersonic inflow BC). Possibility to prescribe different values towards the depth of the ghost cells (varying profile in the BC). Prescribed profile given with variable *field*. Input variables are (*w, loc, interf, field, im, jm*).

* :func:`bc_no_reflexion_2d`: Non-reflective BC. It prescribes a target value at the face center (same value within the depth of the ghost cells) computed from the characteristics through the state inside the domain and a prescribed state given with variable *wbd*. Input variables are (*w, wbd, loc, interf, nx, ny, gam, gh, im, jm*).

* :func:`bc_supandsubinlet_2d`: Mix of subsonic and supersonic inflow. If M<1, use a non-reflective BC (:func:`bc_no_reflexion_2d`) otherwise Dirichlet BC (:func:`bc_general_2d`). Input variables are (*w, loc, interf, field, nx, ny, gam, gh, im, jm*).

* :func:`bc_extrapolate_o2_2d`: Extrapolation BC at order 2 (supersonic outflow BC). Available at order 2, 3, 4, 5, 7 and 9. Input variables are (*w, loc, interf, im, jm, gh*).

* :func:`bc_symmetry_2d`: Symmetry BC. Input variables are (*w, loc, interf, nx, ny, gh, im, jm*).

* :func:`bc_antisymmetry_2d`: Anti-symmetry BC. Input variables are (*w, loc, interf, nx, ny, gh, im, jm*).

* :func:`jn_match_2d`: Join BC for periodic mesh (as the O-mesh of a cylinder for instance) or multi-block management. Copy the values given by the input *wd* into *wr*. Input variables are (*wr, prr, gh1r, gh2r, gh3r, gh4r, imr, jmr, wd, prd, gh1d, gh2d, gh3d, gh4d, imd, jmd, tr*).

* :func:`bc_wall_viscous_adia_2d`: Adiabatic viscous wall BC. Dirichlet BC for velocities :math:`u = v = 0`, Neumann BC for pressure with the assumption :math:`\frac{dp}{dn} = 0`. Input variables are (*w, loc, gam, interf, gh, im, jm*).

* :func:`bc_wall_viscous_iso_2d`: Constant isotherm viscous wall BC. Dirichlet BC for velocities :math:`u = v = 0`, Neumann BC for pressure with the assumption :math:`\frac{dp}{dn} = 0`. Prescribed constant wall temperature with variable *twall*. Input variables are (*w, twall, loc, gam, rgaz, interf, gh, im, jm*).

* :func:`bc_wall_viscous_iso_profile_2d`: Variable isotherm viscous wall BC. Dirichlet BC for velocities :math:`u = v = 0`, Neumann BC for pressure with the assumption :math:`\frac{dp}{dn} = 0`. Prescribed wall temperature profile with variable *twallprof*. Input variables are (*w, twallprof, loc, gam, rgaz, interf, gh, im, jm*).

* :func:`bc_wall_blow_profile_2d`: Adiabatic viscous wall BC with non-zero velocity in y-direction (equal to the wall-normal direction only if the wall is horizontal). Dirichlet BC for velocities :math:`u = 0`, :math:`v = velprof`, Neumann BC for pressure with the assumption :math:`\frac{dp}{dn} = 0`. Prescribed wall velocity profile in y-direction with variable *velprof*. Input variables are (*w, velprof, loc, gam, rgaz, interf, gh, im, jm*).


Inputs for the linearised boundary conditions are different: :ref:`linearisedbcinput`.


Linearised operators - Jacobian
-----

Exact linearisation of the residual is computed by the Algorithmic Differentiation tool. Then, the Jacobian is computed by series of test-vectors to fill in the different entries of the Jacobian without overlapping cross contributions. Test-vectors and indexing of matrix-vector products functions are inside *ComputeJacobian.f90*.

.. note::
   
   Opposite of the Jacobian is computed from the residual: :math:`A = - \frac{dR}{dq} \Rightarrow \frac{dq'}{dt} + Aq' = 0`

Time solvers
-----

Three (pseudo-)time solvers are available:

* *direct*: low-storage Runge-Kutta.
* *implicit*:  matrix-free implicit solver (similar to LU-SGS on approximated Jacobian).
* *fixed_point*: Newton solver.

