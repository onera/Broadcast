Tutorial: Boundary Layer
=====================


Base-flow
--------

To compute the base-flow of a boundary layer, two Python programs should be modified:

#. the card (choice of geometry, physics, solver methods...): **card_bl2d_fv_npz.py**.
#. the main program (choice of mesh, normalisation, BC...): **BROADCAST_npz.py**.

In this tutorial, we will come through each program step by step to setup a hypersonic boundary layer over an adiabatic flat plate.

When your case is setup, run the program to compute the base-flow:

.. code-block:: console

   $ mpirun -np 1 python card_bl2d_fv_npz.py

card_bl2d_fv_npz.py
^^^^^

First, inside **card_bl2d_fv_npz.py**, specify first the geometry of your domain: discretisation in i-direction (equal to x-direction only if the domain is horizontal) and j-direction, length, height and initial offset from origin:

.. code-block:: python

   im    = 300     # X discretization
   jm    = 150     # Y discretization 
   L     = 0.59    # FP length
   high  = 0.035   # FP high
   xini  = 0.0060  # beginning 

Select the folder *out_dir* to save your state solution and the *npz* file *treename* to save your geometry, state, residual, Jacobian.

.. code-block:: python

   out_dir = 'Wksp/'
   treename = 'tree300x150'

Select the number of iterations of your method and the CFL number (for any method). Select the frequency of printing/writing (for explicit/implicit methods only). For Newton method, by default it prints residual every iterations and writes the state at the last iteration. 

.. code-block:: python

   ite     = 10
   cfl     = 1.e10 #very large CFL for Newton method
   freqres = ite/2 # frequency to print residual defined by the integer, here 2
   freqsort= ite/1 # frequency to write solution defined by the integer, here 1

.. note::

   For the Newton method, start with a very large CFL as the scaling of this CFL is not normalised to 1. Try cfl = 1.e10.

Select the order of the convective scheme:

.. code-block:: python

   sch   = 'dnc'  # numerical scheme 
   order =   5    # formal FD order for dnc

Select the order of the extrapolation for the extrapolation BC (if used):

.. code-block:: python

   extraporder = 2   # extrapolation order for outflow

Select the equations/type of numerical scheme to use: *polar* for axisymmetric, *nowall* for general cartesian configuration and default one for an adiabatic flat plate:

.. code-block:: python

   if sch == 'dnc':
       routinesch = 'flux_num_dnc%i_2d' % order
       # routinesch = 'flux_num_dnc%i_nowall_2d' % order
       # routinesch = 'flux_num_dnc%i_nowall_polar_2d' % order 

Select the type of Boundary Conditions to use (location and application of the BC must be specified inside **BROADCAST_npz.py**):

.. code-block:: python

   routineout = 'bc_extrapolate_o%i_2d' % extraporder
   routinein  = 'bc_supandsubinlet_2d'
   # routinein  = 'bc_general_2d'
   routinenr  = 'bc_no_reflexion_2d'
   routinew   = 'bc_wall_viscous_adia_2d'


Select the solver method, *fixed_point* for Newton method:

.. code-block:: python

   compmode = 'fixed_point'  # computational mode in ['direct', 'impli', 'fixed_point']

Select your physical setup parameters with:

#. Mach number
#. Static free-stream temperature
#. Unit Reynold number

.. code-block:: python
   
   dphys['Mach']     = 4.5  
   dphys['T0']       = 288.  
   # dphys['P0']       = 728.312  
   dphys['Runit']    = 3.4e6

At the end of the file **card_bl2d_fv_npz.py**, the function :func:`bl2d_prepro` from **BROADCAST_npz.py** file initialises the geometry, the BC location and the normalisation. This function should also be modified by the user.

Then, the function :func:`bl2d_fromNPZtree` from **BROADCAST_npz.py** solves the configuration previously setup by :func:`bl2d_prepro`.

.. note::

   To restart a computation, comment the call to function :func:`bl2d_prepro` inside **card_bl2d_fv_npz.py**, otherwise you will repeat the pre-process.


BROADCAST_npz.py
^^^^^

Secondly, go inside the function :func:`bl2d_prepro` in **BROADCAST_npz.py**.

Specify the mesh in x-direction, the mesh is here uniform:

.. code-block:: python

   ## MESH in x-direction
   x  = _np.linspace(xini, xini+L , im+1)


Specify the mesh in y-direction, the mesh is splitted into two parts with different stretching if y < *deltaBL* or y > *deltaBL*:

.. code-block:: python

   ## MESH in y-direction
   Ny_in   = 80*jm//100 #number of points inside the BL    
   deltaBL = high/4     #height of the BL
   percent = 0.02       #growth factor increase inside the BL

   Ny_out  = jm - Ny_in 
   Nend    = high/deltaBL
   y_int   = mesh.bigeom_stretch_in(Ny_in, deltaBL, percent)
   y_out   = mesh.exp_stretch_out(Ny_out, deltaBL, percent, Nend)
   y       = _np.concatenate((y_int, y_out)) 

.. note::

   You can create your own mesh with an external meshing tool. For a cartesian rectangular mesh, import *x* and *y* grid point profiles as numpy arrays. Otherwise, import the full range of grid points as numpy arrays and store it inside the variables *x0* and *y0*.

Normalisation is performed with :math:`\rho_\infty`, :math:`U_\infty` and :math:`T_\infty`. In this example, normalisation of the length is performed with the unit Reynolds number.

.. code-block:: python

   ## Adim with ref length
   # Lref   = 8.e-2  
   # Muref  = Roref*Vref*Lref
   ## OR Adim with unit Reynolds
   Muref  = muinf
   Lref   = Muref/(Roref*Vref)
   ## OR no normalisation
   # Roref = 1.
   # Vref  = 1.
   # Tref  = 1.
   # Lref  = 1.
   # Muref = 1.

.. note::

   Dimensionalised data were provided inside **card_bl2d_fv_npz.py** because the normalisation is performed here. It is recommended to run the solver with normalisation as the operators are better conditionned for linear solvers. Resolvent and global stability analysis assumes normalised operators so normalisation is strongly recommended.

Specify the interfaces of the domain i.e. the location of the boundary conditions. Be careful, indexing is in FORTRAN (start at 1 for the first cell). Example for the inlet BC, it is along i=1, starts at the first bottom cell j=1 until the last top cell j=jm. 

.. code-block:: python

   # Ilo
   interf1      = _np.zeros((2,2), order='F')
   interf1[0,0] = 1  # imin
   interf1[0,1] = 1  # jmin
   interf1[1,0] = 1  # imax
   interf1[1,1] = jm # jmax 

A second example for bottom BC, it is along j=1, starts before the first left cell (located at i=1) inside the left ghost cells i=1-gh until the very last right ghost cells i=im+gh. 

.. code-block:: python

   # Jlo
   interf3      = _np.zeros((2,2), order='F')
   interf3[0,0] = 1-gh # imin 
   interf3[0,1] = 1  # jmin
   interf3[1,0] = im+gh # imax
   interf3[1,1] = 1  # jmax

.. note::

   Because the viscous fluxes are based on a compact stencil, boundary conditions should also be specified inside the ghost cells at the four corners of the domain. Notice the example of interf3 where the bottom boundary condition is applied from i=1-gh until i=im+gh. It results that boundary conditions should be applied in a good order. In this example, the inlet boundary condition should be applied before the bottom boundary condition.

Initialise the profiles for Dirichlet and non-reflection BC with the variables *field* and *wbd*. Be careful that they should be the same length as the corresponding interface. For instance, if *interf1* is the inlet BC where a Dirichlet is applied therefore the corresponding *field* has the length *jm* to match *interf1* range.

.. code-block:: python

   field = _np.zeros((jm, gh, 5), order = 'F') # profile for inlet, different values inside the ghost cells
   wbd   = _np.zeros((im+gh , 5), order = 'F') # profile for non-reflection top BC, value at the first ghost cell only

Compute and initialise all the state with a compressible self-similar solution:

.. code-block:: python

   road,uad,vad,Ead = blsim.BLprofile(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False)

Fill in the variables *field* and *wbd* for a boundary layer case with :func:`set_bndbl_2d`. Otherwise, these variables can be filled by the user with imported numpy arrays.

.. code-block:: python

   f_init.set_bndbl_2d(w, field, wbd, im)

Eventually, write all the setup inside a .npz file:

.. code-block:: python

   writeNPZ(filename, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, pinf=pinf)


Now, go inside the function :func:`bl2d_fromNPZtree` in **BROADCAST_npz.py**. Let's consider the Newton method solver:

.. code-block:: python

   elif compmode == 'fixed_point':

Apply the Boundary conditions before the computation of the residual (they should be applied in the good order):

.. code-block:: python

   # Boundary on state vector
   # finflow(w,'Ilo', interf1, field,im,jm)          
   finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
   fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
   foutflow(w,'Ihi', interf2, im, jm, gh)
   fwall(w,'Jlo', gam, interf3, gh, im, jm)

.. note::

   Be careful that the same boundary conditions should be applied three times in the code:

   #. BC on the state.
   #. Linearised BC to construct the Jacobian.
   #. Linearised BC to construct the 3D contributions of the Jacobian.

Compute the residual:

.. code-block:: python

   fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)

Then, the construction of the Jacobian follows an iterative procedure:

#. Definition of a test-vector with :func:`testvector`.
#. Apply the linearised BC.
#. Apply the linearised residual.
#. Indexing of the matrix-vector product to construct the Jacobian with :func:`computejacobianfromjv_relaxed`.

The Jacobian is constructed in a CSR PETSc format:

.. code-block:: python

   Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

Linear solver (LU-factorisation here) to invert the Jacobian :math:`A` is defined:

.. code-block:: python

   ksp  = pet.kspLUPetsc(Jacs)

Newton iteration is performed by solving :math:`\delta q = A^{-1} R(q)`.

.. code-block:: python

   ksp, dwtmp = pet.iterNewton(_np.ravel(res[gh:-gh,gh:-gh,:]), Jacs, ksp)

After convergence, the solution state is written at the cell center in a .dat file:

.. code-block:: python

   filename = out_dir + '/state_atcenter_ite%i.dat' % it
   __writestate_center(filename, im, jm, w, xc, yc, gh)

The Jacobian under the form of list of indices and values (and the solution state including the ghost cells) are written inside the setup .npz file:

.. code-block:: python

   fillNPZ(filename, w, res, IA, JA, Jacvol, gh)

In order to study three-dimensional periodic eigenmodes in resolvent or global stability analyses, the 3D contributions are also computed by the iterative procedure and stored in the same .npz file:

.. code-block:: python

   fillNPZ_3D(filename, IAdz, JAdz, Jacdz, IAdz2, JAdz2, Jacdz2)


Resolvent
--------

To compute the resolvent analysis of a boundary layer, go inside the file **resolvent_all3D_1block.py**.

Specify the equations where to apply forcing:

.. code-block:: python

   equations = [1, 2, 3]  #momentum
   # equations = [0, 1, 2, 3, 4]  #all equations

Compute Chu energy norm *Qq* and L2 norm *Qvol2*:

.. code-block:: python

   Qq = computeQ_Echu(w[gh:-gh,gh:-gh,:], vol[gh:-gh,gh:-gh], gam, mach)
   Qvol2, Qvol2inv = computeQ_L2(vol[gh:-gh,gh:-gh])

Restrict both forcing and response inside a sub-domain defined by xmin, xmax, ymin, ymax. The restriction here excludes the outlet (:math:`Re_x < 1.75 \times 10^6`) and top (half of the initial domain size) parts of the domain:

.. code-block:: python

   xmin = x[0,0]
   # xmin = 1.4e5

   xmax = x[-1,0]
   xmax = 1.75e6

   ymin = y[0,0]

   ymax = y[0,-1]
   ymax = y[0,-1]/2

Select the norm to apply for forcing and response and solve the eigenvalue problem. Be careful to be consistent with the choice of the equations where to apply forcing.

.. code-block:: python

   ##  Chu norm for response and L2 norm for forcing
   eigenvalue, eigenvector_forcing, eigenvector_response = resolvent(frequency, wavenumber, Jacsurvol, Qq, Qvol2, Qvol2inv, P, Dz, Dzz)
   ##  Chu norm for response and forcing
   # eigenvalue, eigenvector_forcing, eigenvector_response = resolvent(frequency, wavenumber, Jacsurvol, Qq, Qq,    Qvol2inv, P, Dz, Dzz)

Write the optimal gain, forcing and response:

.. code-block:: python

   filename = out_dir + "/eigenval.dat"
   if rank ==0: __writearray2(filename, eig, freq, wave)

   filename = out_dir + "/response_atcenter_eig_om{:.2}_be{:.2}_n{:d}_real.dat".format(freq, wave, k)
   __writestate_center_gh(filename, im, jm, w_response, x, y)
   f_opt = _np.reshape( _np.real(_np.array(eigenvector_forcing[k])), (im,jm,5))
   filename = out_dir + "/forcing_atcenter_eig_om{:.2}_be{:.2}_n{:d}_real.dat".format(freq, wave, k)
                  
Finally run the program with the four arguments:

#. Input .npz setup file.
#. Output folder.
#. Frequency :math:`\omega` (remember that this value is then multiplied by :math:`10^{-5}`).
#. Spanwise wavenumber :math:`\beta` (remember that this value is then multiplied by :math:`10^{-5}`).

For instance in the case of the first Mack mode:

.. code-block:: console

   $ mpirun -np 1 python resolvent_all3D_1block "tree300x150.npz" "./Wksp/firstmackmode" 3. 12.

.. note::

   The ansatz for the optimal response is :math:`q'=\check{q}e^{i(-\omega t + \beta z)}`. In this code, the resolvent operator writes :math:`\mathcal{R}=(i\omega I - A)` therefore from the definition of :math:`A`, one gets :math:`\check{q}=-\mathcal{R}P\check{f}`. Be aware that the output of **resolvent_all3D_1block.py** are :math:`\check{f}` and :math:`-\check{q}`.


Global stability analysis
--------

To compute the biglobal stability analysis of a boundary layer, go inside the file **biglobal.py**.

Select the number of eigenvalues to compute:

.. code-block:: python

   ## Number of eigenvalue to compute
   nev = 1  

Select options for the Krylov-Schur algorithm:

.. code-block:: python

   ## Maximum iterations of the Krylov-Schur method
   maxits = 30  
   ## Tolerance of the Krylov-Schur method
   tol = 1.e-5 

Select by commenting one of the line to solve the direct :math:`A\hat{q}=\lambda \hat{q}` or the adjoint :math:`A^*\hat{q}=\lambda \hat{q}` problem. 

.. code-block:: python

   ## OR manual shift-invert - It can compute direct & adjoint modes
   ## For direct mode
   A, ksp_A = pet.createShiftInvert(Jac3D, target)
   ## For adjoint mode
   # A, ksp_A = pet.createShiftInvert_Transpose(Jac3D, target)


Finally run the program with the four arguments:

#. Input .npz setup file.
#. Output folder.
#. Target eigenvalue, also called shift parameter :math:`s`. We look for the closest eigenvalue :math:`\lambda` to :math:`s`. It is a complex value. Real part is the growth rate and imaginary part gives the frequency.
#. Spanwise wavenumber :math:`\beta`.

For instance:

.. code-block:: console

   $ mpirun -np 1 python biglobal.py "tree300x150.npz" "./Wksp/globalmode" 0. 0.

.. note::

   The ansatz for the goblal mode is :math:`q'=\hat{q}e^{-\lambda t + i\beta z}` with :math:`\lambda=\sigma+i\omega`. From the definition of :math:`A`, one gets :math:`A\hat{q}=\lambda \hat{q}`. Therefore, unstable global modes have negative real part: :math:`\sigma < 0`.


