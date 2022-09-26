# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.


## Compute biGlobal stability analysis for SWBLI

import numpy as _np
from mpi4py import MPI
import misc.PETSc_func as pet
from petsc4py import PETSc

import resolvent_all as resol
import restart_init as ri
import timeit
import sys

if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    target =  0.7j  #0.7j   #Opposite frequency for adjoint
    nev = 1    #40  #1

    # wave = float(sys.argv[1])
    # # wavenumber = _np.array([wave])
    # wavenumber = wave

    maxits = 30  #30
    tol = 1.e-5  #1.e-5

    dir = './BASEFLOW_CYL/'


    out_dir = './Wksp/Cylinder/dnc_5/'


    file  = 'state_atcenter_630_y300_Re46p8'


    im = 630  # X discretization  #630
    jm = 300  # Y discretization  #300
    BLprof = _np.loadtxt(out_dir + file + '.dat',comments=('#','ZONE'),skiprows=3 )
    x   = _np.reshape(BLprof[:,0]  ,(im,jm), order='F')
    y   = _np.reshape(BLprof[:,1]  ,(im,jm), order='F')

    ## Read Jacsurvol
    viewer=PETSc.Viewer().createBinary(dir+'Jacsurvol', PETSc.Viewer.Mode.READ)
    Jacsurvol=PETSc.Mat()
    Jacsurvol.create(PETSc.COMM_WORLD)
    Jacsurvol.setType('mpiaij')
    Jacsurvol.load(viewer)
    Jacsurvol.assemble()

    ## Restriction in space of the Jacobian
    imin = 0

    imax = im-1

    jmax = jm-1
    print(jmax)
    # print x[-1,jmax]

    filename = out_dir + 'spectrum' + dir[10:-1]
    eigarray = _np.zeros(1)

    restrixJac = resol.restrix_cyl(im, jm, 5, x, y, imin, imax, jmax, square=False)
    Jacredux = restrixJac.transposeMatMult(Jacsurvol.matMult(restrixJac))

    Jac3D = Jacredux
    # Jac3D = -Jacredux


    ## For 3D
    # viewer=PETSc.Viewer().createBinary(dir+'Dz', PETSc.Viewer.Mode.READ)
    # Dz=PETSc.Mat()
    # Dz.create(PETSc.COMM_WORLD)
    # Dz.setType('mpiaij')
    # Dz.load(viewer)
    # Dz.assemble()
    # viewer=PETSc.Viewer().createBinary(dir+'Dz2', PETSc.Viewer.Mode.READ)
    # Dzz=PETSc.Mat()
    # Dzz.create(PETSc.COMM_WORLD)
    # Dzz.setType('mpiaij')
    # Dzz.load(viewer)
    # Dzz.assemble()
    # Dzredux = restrixJac.transposeMatMult(Dz.matMult(restrixJac))
    # Dzzredux = restrixJac.transposeMatMult(Dzz.matMult(restrixJac))
    # Jac3D = Jacredux - wavenumber**2 * Dzzredux + wavenumber * 1.j * Dzredux

    ksp = pet.kspLUPetsc(Jac3D)
    t1 = timeit.time.time()

    ## Shift-invert from Petsc - Only direct modes
    # vals,vecR,vecI,eps = pet.eigPetsc(Jac3D,ksp,target,nev, restrixJac, tol=tol, maxits=maxits)

    ## OR manual shift-invert - It can compute direct & adjoint modes
    ## For direct mode
    # A, ksp_A = pet.createShiftInvert(Jac3D, target)
    ## For adjoint mode
    A, ksp_A = pet.createShiftInvert_Transpose(Jac3D, target)

    Id = pet.createMatIDPetsc(Jac3D.getSize()[0], Jac3D.getSize()[1])
    eigshift, eigenvector, eps = pet.eigPetsc2(comm ,A , Id, Id ,nev=nev,tol=tol,maxits=maxits)
    vals = target + 1./_np.array(eigshift)
    vecRt = pet.computeEigenvector(comm, eigenvector, restrixJac)
    vecR = [vecRt]


    t2 = timeit.time.time()
    t_eig = t2 - t1

    if rank == 0:
        print('Time eig:', t_eig)
        print('Frequency n', 1)
        print(vals)
        print(_np.real(vals[0]))
        eigarray[0] = _np.real(vals[0])

    pet.saveSpectrumTecplot(eps,rank,filename)

    if rank == 0:
        for i in range(len(vecR)):
            w_response = _np.reshape( _np.real(vecR[i]), (im,jm,5))
            filename1 = out_dir + '/mode_atcenter_eig_n%i_real.dat' % i
            resol.__writestate_center_gh(filename1, im, jm, w_response, x, y)
            w_responsei = _np.reshape( _np.imag(vecR[i]), (im,jm,5))
            filename2 = out_dir + '/mode_atcenter_eig_n%i_imag.dat' % i
            resol.__writestate_center_gh(filename2, im, jm, w_responsei, x, y)
    comm.Barrier()

    if rank == 0:
        print(eigarray)    

