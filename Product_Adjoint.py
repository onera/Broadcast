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

filedir = 'dnc_5/mode_630_y300'

fileadj  = 'dnc_5/mode_630_y300_adjoint'


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

filedirt = './' + dir + '/' + filedir
BLprof = _np.loadtxt(filedirt+ '_real.dat',comments=('#','ZONE'),skiprows=3)
ro  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
rou = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rov = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
row = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roe = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
BLprof = _np.loadtxt(filedirt+ '_imag.dat',comments=('#','ZONE'),skiprows=3)
roi  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
roui = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rovi = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
rowi = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roei = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
dire = _np.zeros((im, jm, 5), order='F')
dire[:,:,0] = ro
dire[:,:,1] = rou 
dire[:,:,2] = rov 
dire[:,:,3] = row 
dire[:,:,4] = roe 
direi = _np.zeros((im, jm, 5), order='F')
direi[:,:,0] = roi
direi[:,:,1] = roui
direi[:,:,2] = rovi
direi[:,:,3] = rowi
direi[:,:,4] = roei
direct = _np.ravel(dire + 1.j * direi)
directbis = dire + 1.j * direi

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

viewer = PETSc.Viewer().createBinary(dirout+'Qvol2', PETSc.Viewer.Mode.READ)
Qq=PETSc.Mat()
Qq.create(PETSc.COMM_WORLD)
Qq.setType('mpiaij')
Qq.load(viewer)
Qq.assemble()

a, b = Qq.getVecs()
c, d = Qq.getVecs()

rangeVec = b.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  b[k] = adjoint[k]
b.assemble()

rangeVec = c.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  c[k] = direct[k]
c.assemble()

## Normalise the direct mode as Sipp, Lebedev JFM 2007

Xloc = 1.
# Yloc = 0.
indj = _np.searchsorted(xc_tmp[im/2,:], Xloc)
# print indj
valnow = 0.5*(directbis[im/2,indj,2] + directbis[im/2+1,indj,2])
vtarget = 0.4612
coefnorm = vtarget / valnow
print(coefnorm)
direct = direct*coefnorm
rangeVec = c.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  c[k] = direct[k]
c.assemble()
sensi = _np.reshape( _np.real(direct), (im,jm,5))
filename1 = './' + dir + '/' + filedir + 'Norm_real.dat'
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)
sensi2 = _np.reshape( _np.imag(direct), (im,jm,5))
filename1 = './' + dir + '/' + filedir + 'Norm_imag.dat'
toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)

####

# Qq.mult(b, a)

ksp = pet.kspLUPetsc(Qq)
ksp.solve(b,a)

## Normalisation of adjoint with Scalar Product

# Qq.mult(c, d)
# # print a.dot(d)
# e = a.dot(d)
# e.conjugate()
# a = a / e
# # print a.dot(d)

###

sol = pet.gatherVector2ArrayPetsc(a,MPI.COMM_WORLD,broadcast=True)

## Normalise adjoint mode given by user (as Sipp, Lebedev JFM 2007)

adjointbis = _np.reshape( _np.real(sol), (im,jm,5)) + 1.j*_np.reshape( _np.imag(sol), (im,jm,5))
Xloc = 1.
# Yloc = 0.
indj = _np.searchsorted(xc_tmp[im/2,:], Xloc)
# print indj
valnow = 0.5*(adjointbis[im/2,indj,2] + adjointbis[im/2+1,indj,2])
vtarget = 0.5
coefnorm = vtarget / valnow
print(coefnorm)
sol = sol*coefnorm

####

sensi = _np.reshape( _np.real(sol), (im,jm,5))
filename1 = './' + dir + '/' + fileadj + 'L2_real.dat'
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)

sensi2 = _np.reshape( _np.imag(sol), (im,jm,5))
filename1 = './' + dir + '/' + fileadj + 'L2_imag.dat'
toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)

