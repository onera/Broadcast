# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

## Compute Limit-Cycle coefficient (Stuart-Landau equation) : 1st part & 2nd part

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

fileadj  = 'dnc_5/mode_630_y300_adjointL2'


## omega at Re = 46.8 
omega = 0.7167814530347437 


imold = 630 #600 #630
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
Jac=PETSc.Mat()
Jac.create(PETSc.COMM_WORLD)
Jac.setType('mpiaij')
Jac.load(viewer)
Jac.assemble()

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
filename1 = './' + dir + '/' + dir2 + '/mode22_mesh7_real.dat' 
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)

sensi2 = _np.reshape( _np.imag(sol22), (im,jm,5))
filename1 = './' + dir + '/' + dir2 + '/mode22_mesh7_imag.dat' 
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

# ksp2.solve(c, a20)
ksp2.solve(-c, a20)

sol20 = pet.gatherVector2ArrayPetsc(a20,MPI.COMM_WORLD,broadcast=True)

sensi = _np.reshape( _np.real(sol20), (im,jm,5))
filename2 = './' + dir + '/' + dir2 + '/mode20_mesh7_real.dat' 
toy.__writestate_center_gh(filename2, im, jm, sensi, xc_tmp, yc_tmp)

sensi2 = _np.reshape( _np.imag(sol20), (im,jm,5))
filename2 = './' + dir + '/' + dir2 + '/mode20_mesh7_imag.dat' 
toy.__writestate_center_gh(filename2, im, jm, sensi2, xc_tmp, yc_tmp)


###################################################
###################################################
###################################################



fileadjt = './' + dir + '/' + fileadj
BLprof = _np.loadtxt(fileadjt+ '_real.dat',comments=('#','ZONE'),skiprows=3)
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

# filemodet = './' + dir + '/' + filemode22
filemodet = filename1[:-9]
BLprof = _np.loadtxt(filemodet+ '_real.dat',comments=('#','ZONE'),skiprows=3)
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
wmode22 = _np.ravel(wmoder + 1.j * wmodei)

# filemodet = './' + dir + '/' + filemode20
filemodet = filename2[:-9]
BLprof = _np.loadtxt(filemodet+ '_real.dat',comments=('#','ZONE'),skiprows=3)
ro  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
rou = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rov = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
row = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roe = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
wmoder = _np.zeros((im, jm,5), order='F')
wmoder[:,:,0] = ro
wmoder[:,:,1] = rou 
wmoder[:,:,2] = rov 
wmoder[:,:,3] = row 
wmoder[:,:,4] = roe             
wmode20 = _np.ravel(wmoder)


a, b = H11.getVecs()
c, d = H11.getVecs()
x1, x2 = H11.getVecs()
y1, y2 = H11.getVecs()
e, f = H11.getVecs()
g, h = H11.getVecs()
i, j = H11.getVecs()

rangeVec = a.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  a[k] = wmode[k]
a.assemble()

rangeVec = b.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  b[k] = adjoint[k]
b.assemble()

rangeVec = c.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  c[k] = wmode22[k]
c.assemble()

rangeVec = d.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  d[k] = wmode20[k]
d.assemble()


Qq.multTranspose(b, e)
print(a.dot(e))


IDE = pet.createMatIDPetsc(Qq.getSize()[0], Qq.getSize()[1])
b.conjugate()
IDE.setDiagonal(b,PETSc.InsertMode.INSERT_VALUES)
b.conjugate()


H11.mult(d,x2)
xt = x2
e.conjugate()
coef1 = xt.tDot(e) / a.tDot(e)
e.conjugate()

print('mu = ', coef1)


H11.conjugate()
H11.mult(c,y1)
yt = y1
e.conjugate()
coef2 = yt.tDot(e) / a.tDot(e)
e.conjugate()


print('nu = ', coef2)

print(coef1 + coef2)
print(_np.real(coef2)/_np.real(coef1))
print(_np.imag(coef2)/_np.imag(coef1))
print(1.*(_np.imag(coef1) + _np.imag(coef2)) / (_np.real(coef1) + _np.real(coef2)))


