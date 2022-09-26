# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

## Compute Limit-Cycle coefficient (Stuart-Landau equation) : 1st part

import numpy as _np

from petsc4py import PETSc
from mpi4py import MPI

import BROADCAST as toy
import restart_init as ri
import misc.PETSc_func as pet


dir   = 'Wksp/Cylinder'

dir2  = 'dnc_5'

dirout = './BASEFLOW_CYL/'

filemode = 'dnc_5/mode_630_y300Norm'

## omega at Re = 46.8 
omega = 0.7167814530347437 


imold = 630 #630
jmold = 300 #300
im = imold
jm = jmold

filemodet = './' + dir + '/' + filemode 
BLprof = _np.loadtxt(filemodet+ '_real.dat',comments=('#','ZONE'),skiprows=3)
xc_tmp = _np.reshape(BLprof[:,0],(imold,jmold), order='F')
yc_tmp = _np.reshape(BLprof[:,1],(imold,jmold), order='F')
ro  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
rou = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rov = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
row = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roe = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
BLprof = _np.loadtxt(filemodet+ '_imag.dat',comments=('#','ZONE'),skiprows=3)
roi  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
roui = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rovi = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
rowi = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roei = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
wmoder = _np.zeros((im, jm,5), order='F')
wmoder[:,:,0] = ro
wmoder[:,:,1] = rou 
wmoder[:,:,2] = rov 
wmoder[:,:,3] = row 
wmoder[:,:,4] = roe             
wmodei = _np.zeros((im, jm,5), order='F')
wmodei[:,:,0] = roi
wmodei[:,:,1] = roui
wmodei[:,:,2] = rovi
wmodei[:,:,3] = rowi
wmodei[:,:,4] = roei 
wmode = _np.ravel(wmoder + 1.j * wmodei)


viewer = PETSc.Viewer().createBinary(dirout+'Jacsurvol', PETSc.Viewer.Mode.READ)
# viewer = PETSc.Viewer().createBinary(dirout+'Jac', PETSc.Viewer.Mode.READ)
Jac=PETSc.Mat()
Jac.create(PETSc.COMM_WORLD)
Jac.setType('mpiaij')
Jac.load(viewer)
Jac.assemble()

# viewer = PETSc.Viewer().createBinary(dirout+'H', PETSc.Viewer.Mode.READ)
viewer = PETSc.Viewer().createBinary(dirout+'H2', PETSc.Viewer.Mode.READ)
H11=PETSc.Mat()
H11.create(PETSc.COMM_WORLD)
H11.setType('mpiaij')
H11.load(viewer)
H11.assemble()

viewer = PETSc.Viewer().createBinary(dirout+'Qvol2', PETSc.Viewer.Mode.READ)
Qq=PETSc.Mat()
Qq.create(PETSc.COMM_WORLD)
Qq.setType('mpiaij')
Qq.load(viewer)
Qq.assemble()

Id = pet.createMatIDPetsc(Jac.getSize()[0], Jac.getSize()[1])

LHS = 2.j*omega*Id + Jac

ksp = pet.kspLUPetsc(LHS)

a, b = H11.getVecs()
a22, c = H11.getVecs()
a20, d = H11.getVecs()

rangeVec = a.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  a[k] = wmode[k]
a.assemble()

H11.mult(a,b)
ksp.solve(-0.5*b,a22)

sol22 = pet.gatherVector2ArrayPetsc(a22,MPI.COMM_WORLD,broadcast=True)

sensi = _np.reshape( _np.real(sol22), (im,jm,5))
filename1 = './' + dir + '/' + dir2 + '/mode22_mesh9_real.dat' 
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)

sensi2 = _np.reshape( _np.imag(sol22), (im,jm,5))
filename1 = './' + dir + '/' + dir2 + '/mode22_mesh9_imag.dat' 
toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)

##################################################

ksp2 = pet.kspLUPetsc(Jac)

a.conjugate()
H11.mult(a,b)

breal = _np.real(pet.gatherVector2ArrayPetsc(b,MPI.COMM_WORLD,broadcast=True))

rangeVec = c.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  c[k] = breal[k]
c.assemble()

ksp2.solve(-c, a20)

sol20 = pet.gatherVector2ArrayPetsc(a20,MPI.COMM_WORLD,broadcast=True)

sensi = _np.reshape( _np.real(sol20), (im,jm,5))
filename1 = './' + dir + '/' + dir2 + '/mode20_mesh9_real.dat' 
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)

sensi2 = _np.reshape( _np.imag(sol20), (im,jm,5))
filename1 = './' + dir + '/' + dir2 + '/mode20_mesh9_imag.dat' 
toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)




