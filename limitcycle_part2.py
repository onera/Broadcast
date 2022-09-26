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

fileadj  = 'dnc_5/mode_630_y300_adjointL2'

filemode22 = 'dnc_5/mode22'

filemode20 = 'dnc_5/mode20'


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

filemodet = './' + dir + '/' + filemode22
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

filemodet = './' + dir + '/' + filemode20
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

viewer = PETSc.Viewer().createBinary(dirout+'Qvol2', PETSc.Viewer.Mode.READ)
Qq=PETSc.Mat()
Qq.create(PETSc.COMM_WORLD)
Qq.setType('mpiaij')
Qq.load(viewer)
Qq.assemble()


# import resolvent_all as resol
# indstop = 15
# print xc_tmp[im/2,jm-1-indstop]
# restrixJac = resol.restrix_cyl(im, jm, 5, xc_tmp, yc_tmp, 0, im-1, jm-1-indstop, square=True)
# Qq = Qq.matMult(restrixJac)


# viewer = PETSc.Viewer().createBinary(dirout+'H', PETSc.Viewer.Mode.READ)
viewer = PETSc.Viewer().createBinary(dirout+'H2', PETSc.Viewer.Mode.READ)
H11=PETSc.Mat()
H11.create(PETSc.COMM_WORLD)
H11.setType('mpiaij')
H11.load(viewer)
H11.assemble()

viewer = PETSc.Viewer().createBinary(dirout+'Hmode22', PETSc.Viewer.Mode.READ)
H22=PETSc.Mat()
H22.create(PETSc.COMM_WORLD)
H22.setType('mpiaij')
H22.load(viewer)
H22.assemble()

viewer = PETSc.Viewer().createBinary(dirout+'Hmode20', PETSc.Viewer.Mode.READ)
H20=PETSc.Mat()
H20.create(PETSc.COMM_WORLD)
H20.setType('mpiaij')
H20.load(viewer)
H20.assemble()

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

H20.mult(a,x1)
H11.mult(d,x2)
print("----------------")
print(x1.dot(e))
print(x2.dot(e))
print("----------------")
## to plot mu
# H20.mult(a,g)
# H11.mult(d,h)
# IDE.mult(g,x1)
# IDE.mult(h,x2)

xt = x1 + x2
e.conjugate()
coef1 = 0.5 * xt.tDot(e) / a.tDot(e)
e.conjugate()

print('mu = ', coef1)


a.conjugate()
H11.conjugate()
H11.mult(c,y1)
H22.mult(a,y2)
a.conjugate()
print("----------------")
print(y1.dot(e))
print(y2.dot(e))
print("----------------")
## to plot nu
# a.conjugate()
# H11.mult(c,i)
# H22.mult(a,j)
# a.conjugate()
# IDE.mult(i,y1)
# IDE.mult(j,y2)

yt = y1 + y2
e.conjugate()
coef2 = 0.5 * yt.tDot(e) / a.tDot(e)
e.conjugate()


print('nu = ', coef2)

print(coef1 + coef2)
print(_np.real(coef2)/_np.real(coef1))
print(_np.imag(coef2)/_np.imag(coef1))
print(1.*(_np.imag(coef1) + _np.imag(coef2)) / (_np.real(coef1) + _np.real(coef2)))


# sol22 = pet.gatherVector2ArrayPetsc(x1,MPI.COMM_WORLD,broadcast=True)
# sensi = _np.reshape( _np.real(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/mu1_real.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)
# sensi2 = _np.reshape( _np.imag(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/mu1_imag.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)
# sol22 = pet.gatherVector2ArrayPetsc(x2,MPI.COMM_WORLD,broadcast=True)
# sensi = _np.reshape( _np.real(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/mu2_real.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)
# sensi2 = _np.reshape( _np.imag(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/mu2_imag.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)
# sol22 = pet.gatherVector2ArrayPetsc(y1,MPI.COMM_WORLD,broadcast=True)
# sensi = _np.reshape( _np.real(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/nu1_real.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)
# sensi2 = _np.reshape( _np.imag(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/nu1_imag.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)
# sol22 = pet.gatherVector2ArrayPetsc(y2,MPI.COMM_WORLD,broadcast=True)
# sensi = _np.reshape( _np.real(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/nu2_real.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)
# sensi2 = _np.reshape( _np.imag(sol22), (im,jm,5))
# filename1 = './' + dir + '/' + dir2 + '/nu2_imag.dat' 
# toy.__writestate_center_gh(filename1, im, jm, sensi2, xc_tmp, yc_tmp)




