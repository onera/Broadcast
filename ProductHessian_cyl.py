# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

## Compute Product Hessian with optimal forcing

import numpy as _np

from petsc4py import PETSc
from mpi4py import MPI

import BROADCAST as toy
import restart_init as ri
import misc.PETSc_func as pet


dir   = 'Wksp/Cylinder'

dir2  = 'dnc_5'


dirout = './BASEFLOW_CYL/'


fileadj  = 'dnc_5/mode_630_y300_adjointL2'


viewer = PETSc.Viewer().createBinary(dirout+'H', PETSc.Viewer.Mode.READ)
# viewer = PETSc.Viewer().createBinary(dirout+'H2', PETSc.Viewer.Mode.READ)
H=PETSc.Mat()
H.create(PETSc.COMM_WORLD)
H.setType('mpiaij')
H.load(viewer)
H.assemble()

imold = 630 #630
jmold = 300 #300

fileadjt = './' + dir + '/' + fileadj
BLprof = _np.loadtxt(fileadjt+ '_real.dat',comments=('#','ZONE'),skiprows=3)
xc_tmp = _np.reshape(BLprof[:,0],(imold,jmold), order='F')
yc_tmp = _np.reshape(BLprof[:,1],(imold,jmold), order='F')
ro  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
rou = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rov = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
row = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roe = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
BLprof = _np.loadtxt(fileadjt+ '_imag.dat',comments=('#','ZONE'),skiprows=3)
roi  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
roui = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rovi = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
rowi = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roei = _np.reshape(BLprof[:,6],(imold,jmold), order='F')

im = imold
jm = jmold

adj = _np.zeros((im, jm, 5), order='F')
adj[:,:,0] = ro
adj[:,:,1] = rou 
adj[:,:,2] = rov 
adj[:,:,3] = row 
adj[:,:,4] = roe 
adji = _np.zeros((im, jm, 5), order='F')
adji[:,:,0] = roi
adji[:,:,1] = roui
adji[:,:,2] = rovi
adji[:,:,3] = rowi
adji[:,:,4] = roei
adjoint = _np.ravel(adj + 1.j * adji)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

a, b = H.getVecs()

rangeVec = b.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  b[k] = adjoint[k]
b.assemble()

b.conjugate()
H.multTranspose(b,a)
# H.mult(b,a)

#########
## With L2 norm
# viewer = PETSc.Viewer().createBinary(dirout+'Qvol2', PETSc.Viewer.Mode.READ)
# Qq=PETSc.Mat()
# Qq.create(PETSc.COMM_WORLD)
# Qq.setType('mpiaij')
# Qq.load(viewer)
# Qq.assemble()

# ksp = pet.kspLUPetsc(Qq)
# ksp.solve(a,b)
# a = b
##

#########
## Sensitivity to forcing
viewer = PETSc.Viewer().createBinary(dirout+'Qvol2', PETSc.Viewer.Mode.READ)
Qq=PETSc.Mat()
Qq.create(PETSc.COMM_WORLD)
Qq.setType('mpiaij')
Qq.load(viewer)
Qq.assemble()

viewer = PETSc.Viewer().createBinary(dirout+'Jacsurvol', PETSc.Viewer.Mode.READ)
Jac=PETSc.Mat()
Jac.create(PETSc.COMM_WORLD)
Jac.setType('mpiaij')
Jac.load(viewer)
Jac.assemble()

ksp = pet.kspLUPetsc(Qq)
Jac.transpose()
ksp2 = pet.kspLUPetsc(Jac)
ksp2.solve(a,b)
ksp.solve(-b,a)
# ksp.solve(b,a)
##

a.conjugate()
sol = pet.gatherVector2ArrayPetsc(a,MPI.COMM_WORLD,broadcast=True)

sensi = _np.reshape( _np.real(sol), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/sensitivity_atcenter_eig_n%i_real.dat' % 0
# filename1 = './' + dir + '/' + dir2 + '/sensitivityL2_atcenter_eig_n%i_real.dat' % 0
filename1 = './' + dir + '/' + dir2 + '/sensitivitytoforcing_atcenter_eig_n%i_real.dat' % 0
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)


sensi2 = _np.reshape( _np.imag(sol), (im,jm,5))
# sensi2 = _np.reshape( _np.imag(b), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/sensitivity_atcenter_eig_n%i_imag.dat' % 0
# filename1 = './' + dir + '/' + dir2 + '/sensitivityL2_atcenter_eig_n%i_imag.dat' % 0
filename1 = './' + dir + '/' + dir2 + '/sensitivitytoforcing_atcenter_eig_n%i_imag.dat' % 0
toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)

