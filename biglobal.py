

## Compute biGlobal stability analysis

import numpy as _np
from mpi4py import MPI
import misc.PETSc_func as pet
from petsc4py import PETSc

import resolvent_all as resol
import restart_init as ri
import timeit
import sys
import os

if __name__ == '__main__':

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    treeBroadcast   = str(sys.argv[1])
    out_dir         = str(sys.argv[2])
    target          = float(sys.argv[3])
    wave            = float(sys.argv[4])

    os.system('mkdir -p %s' % out_dir)

    wavenumber = _np.array([wave])

    # target = -0.1  #0. #-0.1
    ## Number of eigenvalue to compute
    nev = 1  

    ## Maximum iterations of the Krylov-Schur method
    maxits = 30  
    ## Tolerance of the Krylov-Schur method
    tol = 1.e-5 

    dic  = _np.load(treeBroadcast)
    gh   = dic['gh']
    x    = dic['xc'][gh:-gh,gh:-gh]
    y    = dic['yc'][gh:-gh,gh:-gh]
    im   = dic['im']
    jm   = dic['jm']
    vol  = dic['vol']
    w    = dic['FlowSolutionEndOfRun']
    gam  = dic['Gamma']
    mach = dic['Mach']
    IA     = dic['IA']
    JA     = dic['JA']
    Jacvol = dic['Aij']
    IAdz   = dic['IAdz']
    JAdz   = dic['JAdz']
    Jacdz  = dic['Aijdz']
    IAdz2  = dic['IAdz2']
    JAdz2  = dic['JAdz2']
    Jacdz2 = dic['Aijdz2']

    Jacsurvol = pet.createMatPetscCSR(IA, JA, Jacvol, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
    Dz        = pet.createMatPetscCSR(IAdz, JAdz, Jacdz, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
    Dzz       = pet.createMatPetscCSR(IAdz2, JAdz2, Jacdz2, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

    ## Restriction in space of the Jacobian
    xmin = x[0,0]  

    xmax = x[-1,0]

    ymin = y[0,0]

    ymax = y[0,-1]

    restrixJac = resol.restrix(im, jm, 5, x, y, xmin, xmax, ymin, ymax, square=False)
    Jacredux = restrixJac.transposeMatMult(Jacsurvol.matMult(restrixJac))
    Dzredux = restrixJac.transposeMatMult(Dz.matMult(restrixJac))
    Dzzredux = restrixJac.transposeMatMult(Dzz.matMult(restrixJac))

    filename = out_dir + '/spectrum.dat'
    eigarray = _np.zeros(len(wavenumber))

    for k in range(len(wavenumber)):

        Jac3D = Jacredux - wavenumber[k]**2 * Dzzredux + wavenumber[k] * 1.j * Dzredux

        ksp = pet.kspLUPetsc(Jac3D)

        t1 = timeit.time.time()

        ## Shift-invert from Petsc - Only direct modes
        # vals,vecR,vecI,eps = pet.eigPetsc(Jac3D,ksp,target,nev, restrixJac, tol=tol, maxits=maxits)
        
        ## OR manual shift-invert - It can compute direct & adjoint modes
        ## For direct mode
        A, ksp_A = pet.createShiftInvert(Jac3D, target)
        ## For adjoint mode
        # A, ksp_A = pet.createShiftInvert_Transpose(Jac3D, target)

        Id = pet.createMatIDPetsc(Jac3D.getSize()[0], Jac3D.getSize()[1])
        eigshift, eigenvector, eps = pet.eigPetsc3(comm ,A , Id, Id ,nev=nev,tol=tol,maxits=maxits)
        vals = target + 1./_np.array(eigshift)
        vecRt = pet.computeEigenvector(comm, eigenvector, restrixJac)
        vecR = [vecRt]

        t2 = timeit.time.time()
        t_eig = t2 - t1

        if rank == 0:
            print('Time eig:', t_eig)
            print('Frequency n', k+1)
            print("Complex eigenvalues found:", vals)
            # print(_np.real(vals[0]))
            eigarray[k] = _np.real(vals[0])

        pet.saveSpectrumTecplot(eps,rank,filename)

        if rank == 0:
            # comm.Barrier()
            w_response = _np.reshape( _np.real(vecR[0]), (im,jm,5))
            filename1 = out_dir + '/mode_atcenter_eig_n%i_real.dat' % k
            resol.__writestate_center_gh(filename1, im, jm, w_response, x, y)
            # comm.Barrier()
            w_responsei = _np.reshape( _np.imag(vecR[0]), (im,jm,5))
            filename2 = out_dir + '/mode_atcenter_eig_n%i_imag.dat' % k
            resol.__writestate_center_gh(filename2, im, jm, w_responsei, x, y)
        comm.Barrier()

    if rank == 0:
        print(eigarray)    

