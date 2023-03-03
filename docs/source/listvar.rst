.. _listvar:


Variables
==========

Main variables of the solver:

Geometry:

* **im** (*integer*): number of cells in i-direction.
* **jm** (*integer*): number of cells in j-direction.
* **gh** (*integer*): number of ghost cells for the boundary conditions. **gh** depends on the order of the numerical scheme: **gh** = 2 (3th order), **gh** = 3 (5th order), **gh** = 4 (7th order) or **gh** = 5 (9th order).
* **xo** (*numpy array of size (im+1+2gh, jm+1+2gh)*): x coordinates of the cell nodes.
* **yo** (*numpy array of size (im+1+2gh, jm+1+2gh)*): y coordinates of the cell nodes.
* **xc** (*numpy array of size (im+2gh, jm+2gh)*): x coordinates of the cell centers.
* **yc** (*numpy array of size (im+2gh, jm+2gh)*): y coordinates of the cell centers.
* **vol** (*numpy array of size (im+2gh, jm+2gh)*): volumes of the cells.
* **volf** (*numpy array of size (im+2gh, jm+2gh, 2)*): inverse of the volume computed at the cell faces. 
* **nx** (*numpy array of size (im+1+2gh, jm+1+2gh, 2)*): normal in x-direction computed at the cell faces.
* **ny** (*numpy array of size (im+1+2gh, jm+1+2gh, 2)*): normal in y-direction computed at the cell faces.

Physics:

* **gam** (*real*): :math:`= \gamma` constant.
* **cp** (*real*): specific heat at constant pressure.
* **cv** (*real*): specific heat at constant volume.
* **rgaz** (*real*): specific gas constant.
* **muref** (*real*): Sutherland reference viscosity.
* **tref** (*real*): Sutherland reference tempearture.
* **cs** (*real*): Sutherland tempearture.
* **prandtl** (*real*): Prandtl number.
* **mach** (*real*): Mach number.

Solution:

* **w** (*numpy array of size (im+2gh, jm+2gh, 5)*): state. **w[:,:,0]** :math:`= \rho` is the density, **w[:,:,1]** :math:`= \rho u` is the momentum in x-direction, **w[:,:,2]** :math:`= \rho v` is the momentum in y-direction, **w[:,:,3]** :math:`= \rho w` is the momentum in z-direction and **w[:,:,4]** :math:`= \rho E` is the total energy.
* **res** (*numpy array of size (im+2gh, jm+2gh, 5)*): residual of Navier-Stokes equations.

Numerics:

* **k2** (*real*): low-order dissipation coefficient for shock sensor in the FE-MUSCL numerical scheme. Default value to 1.01.
* **k4** (*real*): high-order dissipation coefficient in the FE-MUSCL numerical scheme. Default value to 1.
* **sch** (*string*): numerical scheme to set. Only one implemented. **sch** = "dnc".
* **order** (*integer*): numerical scheme order. Availables options are 3, 5, 7 or 9.

Boundary conditions:

* **interf** (*numpy array of size (2, 2)*): range of the coordinates for one BC location. **interf[0,0]** = imin, **interf[0,1]** = jmin, **interf[1,0]** = imax and **interf[1,1]** = jmax.
* **loc** (*string*): location of the BC. Available values are "Ilo" (i=imin), "Ihi" (i=imax), "Jlo" (j=jmin) and "Jhi" (j=jmax).
* **lf** (*list of strings*): list of routines names for BC and numerical scheme.

Linearised operators - Jacobian:

* **IA** (*numpy array of size (number of non-zeros entries of the Jacobian)*): lists of row indices of the non-zeros entries of the Jacobian.
* **JA** (*numpy array of size (number of non-zeros entries of the Jacobian)*): lists of column indices of the non-zeros entries of the Jacobian.
* **Aij** (*numpy array of size (number of non-zeros entries of the Jacobian)*): lists of values of the non-zeros entries of the Jacobian.
