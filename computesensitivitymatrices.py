#!/usr/bin/env python

'''
File: bl2d.py

Created on 21 july 2020

@author:       Cedric Content
@contact:      cedric.content@onera.fr
@organization: ONERA - DAAA

@summary:      This file is the main file of the program. It contains the
               routine "main" and other related routines.
'''

import BROADCAST_npz_sens as toy

import numpy as _np
# FORTRAN
import srcfv.f_geom    as f_geom
import srcfv.f_bnd     as f_bnd
import srcfv.f_sch     as f_sch
import srcfv.f_lhs     as f_lhs
import srcfv.f_lin     as f_lin
# import srcfv.f_adj     as f_adj
import srcfv.f_norm    as f_norm
import misc.f_misc     as f_misc

from petsc4py import PETSc
from mpi4py import MPI
import misc.PETSc_func as pet
import resolvent_all as resol

import os
import sys
import timeit

######################### CARTE ####################


out_dir = './Wksp/firstmackmode/'

treename = './Wksp/dnc_5/tree300x150'

## 2nd Mack mode
# fileforc = './Wksp/secondmackmode/forcing_atcenter_eig_om2.3e+01_be0.0_n0'
# fileresp = './Wksp/secondmackmode/response_atcenter_eig_om2.3e+01_be0.0_n0'

## 1st Mack mode
fileforc = './Wksp/firstmackmode/forcing_atcenter_eig_om3.0_be1.2e+01_n0'
fileresp = './Wksp/firstmackmode/response_atcenter_eig_om3.0_be1.2e+01_n0'

3Dmode = True  #True if oblique first Mack mode, False if 2D second Mack mode

## Frequency: omega mode for resolvent operator 
# omega = 23.e-5  #2nd Mack mode
omega = 3.e-5   #1st Mack mode

## Only for 3D modes, spanwise wavenumber (beta)
if 3Dmode:
  wave = 12.e-5  #1st Mack mode  

## Squarred optimal gain mu^2
# mu2 = 17945440.70685651**2  #for 2nd Mack mode (adiab)
mu2 = 1.15292235e+08**2  #for 1st Mack mode (adiab)

#################### MATRICES CONSTUCTION ########################

# Compute dR/dp
tinit = timeit.time.time()
toy.dRdp_fromNPZtree(treename, out_dir = out_dir)
tend = timeit.time.time()
tlaps = tend-tinit
print("Time Elapsed 0 =  ", tlaps)

## Compute dQchu/dq
tinit = timeit.time.time()
toy.dQchudq_fromNPZtree(treename, fileforc, fileresp, out_dir = out_dir)
tend = timeit.time.time()
tlaps = tend-tinit
print("Time Elapsed 1 =  ", tlaps)

## Compute H=dA/dq
tinit = timeit.time.time()
toy.dAdq_fromNPZtree(treename, fileresp, out_dir = out_dir)
tend = timeit.time.time()
tlaps = tend-tinit
print("Time Elapsed 2 =  ", tlaps)

## Only for 3D modes
if 3Dmode:
  tinit = timeit.time.time()
  toy.dAdq3D_fromNPZtree(treename, fileresp, out_dir = out_dir)
  tend = timeit.time.time()
  tlaps = tend-tinit
  print("Time Elapsed 2.5 =  ", tlaps)

# Compute \tilde{H}=dA/dp
tinit = timeit.time.time()
toy.dAdp_fromNPZtree(treename, fileresp, out_dir = out_dir)
tend = timeit.time.time()
tlaps = tend-tinit
print("Time Elapsed 3 =  ", tlaps)

## Only for 3D modes
if 3Dmode:
  tinit = timeit.time.time()
  toy.dAdp3D_fromNPZtree(treename, fileresp, out_dir = out_dir)
  tend = timeit.time.time()
  tlaps = tend-tinit
  print("Time Elapsed 3.5 =  ", tlaps)


#################### MATRIX-VECTOR PRODUCTs ########################

os.system('mkdir -p %s' % out_dir)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

dic = _np.load(treename+".npz")
im = dic['im']
jm = dic['jm']
gh = dic['gh']
xc = dic['xc']
yc = dic['yc']
vol = dic['vol']
gam = dic['Gamma']
mach = dic['Mach']
w = dic['FlowSolutionEndOfRun']

## Construct the Hessian dA/dq
IAh = dic['IAdAdq'] 
JAh = dic['JAdAdq'] 
Jach = dic['AijdAdq'] 
H = pet.createMatPetscCSR(IAh, JAh, Jach, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

## Only for 3D modes
if 3Dmode:
  IAdz = dic['IAdAdqdz'] 
  JAdz = dic['JAdAdqdz'] 
  Jacdz = dic['AijdAdqdz'] 
  IAdzz = dic['IAdAdqdz2'] 
  JAdzz = dic['JAdAdqdz2'] 
  Jacdzz = dic['AijdAdqdz2'] 
  Hdz = pet.createMatPetscCSR(IAdz, JAdz, Jacdz, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
  Hdzz = pet.createMatPetscCSR(IAdzz, JAdzz, Jacdzz, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
  H = H + wave*1.j*Hdz - wave**2*Hdzz

## Load optimal response
import restart_init as ri
xc_tmp, yc_tmp, ro, rou, rov, row, roe = ri.read_init(fileresp+ '_real.dat')
xc_tmp, yc_tmp, roi, roui, rovi, rowi, roei = ri.read_init(fileresp+ '_imag.dat')
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
response = _np.ravel(wmoder + 1.j * wmodei)
## Load optimal forcing
xc_tmp, yc_tmp, ro, rou, rov, row, roe = ri.read_init(fileforc+ '_real.dat')
xc_tmp, yc_tmp, roi, roui, rovi, rowi, roei = ri.read_init(fileforc+ '_imag.dat')
wforcr = _np.zeros((im, jm,5), order='F')
wforcr[:,:,0] = ro
wforcr[:,:,1] = rou 
wforcr[:,:,2] = rov 
wforcr[:,:,3] = row 
wforcr[:,:,4] = roe 
wforci = _np.zeros((im, jm,5), order='F')
wforci[:,:,0] = roi
wforci[:,:,1] = roui
wforci[:,:,2] = rovi
wforci[:,:,3] = rowi
wforci[:,:,4] = roei
forcing = _np.ravel(wforcr + 1.j * wforci)

## Create Identity matrix for resolvent operator
Id = pet.createMatIDPetsc(im*jm*5, im*jm*5)

## Construct the Jacobian A=dR/dq
IA = dic['IA'] 
JA = dic['JA'] 
Jac = dic['Aij'] 
Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
Jacs3D = Jacs

## Only for 3D modes
if 3Dmode:
  IAdz = dic['IAdz'] 
  JAdz = dic['JAdz'] 
  Jacdz = dic['Aijdz'] 
  IAdzz = dic['IAdz2'] 
  JAdzz = dic['JAdz2'] 
  Jacdzz = dic['Aijdz2'] 
  Jacsdz = pet.createMatPetscCSR(IAdz, JAdz, Jacdz, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
  Jacsdzz = pet.createMatPetscCSR(IAdzz, JAdzz, Jacdzz, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
  Jacs3D = Jacs - wave*1.j*Jacsdz - wave**2*Jacsdzz  #conjugated

kspA = pet.kspLUPetsc(Jacs)

## Construct dR/dp
IAdr = dic['IAdRdp'] 
JAdr = dic['JAdRdp'] 
Jacdr = dic['AijdRdp'] 
dRdp = pet.createMatPetscCSR(IAdr, JAdr, Jacdr, im*jm*5, im, im)

## Construct Chu energy matrix
Qchu = resol.computeQ_Echu(w[gh:-gh,gh:-gh,:], vol[gh:-gh,gh:-gh], gam, mach)
kspQ = pet.kspLUPetsc(Qchu)

## Construct derivative of Chu with response and forcing
IAdcqdq = dic['IAdChuqdq']
JAdcqdq = dic['JAdChuqdq']
Jacdcqdq = dic['AijdChuqdq']
IAdcfdq = dic['IAdChufdq']
JAdcfdq = dic['JAdChufdq']
Jacdcfdq = dic['AijdChufdq']
dChuq = pet.createMatPetscCSR(IAdcqdq, JAdcqdq, Jacdcqdq, im*jm*5, im*jm*5, 5)
dChuf = pet.createMatPetscCSR(IAdcfdq, JAdcfdq, Jacdcfdq, im*jm*5, im*jm*5, 5)

print("Matrices loaded, computation of sensitivity...", flush=True)

## Construct conjugate resolvent operator
Rconjinv = -1.j*omega*Id - Jacs3D
kspR = pet.kspLUPetsc(Rconjinv)

##Initialise vectors in PETSc format
a, b = H.getVecs()
c, d = H.getVecs()
e, f = H.getVecs()
g, h = dRdp.getVecs()

rangeVec = a.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  a[k] = response[k]
a.assemble()

rangeVec = b.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  b[k] = forcing[k]
b.assemble()

## Compute the first term with the Hessian for sensitivity to base-flow

## with P restriction
Qchu.mult(a,c)
kspR.solveTranspose(c,d)
d.conjugate()
H.multTranspose(d,e)
e.conjugate()

## without P restriction
# Qchu.mult(b,c)
# c.conjugate()
# H.multTranspose(c,e)
# e.conjugate()

## Compute both terms with Chu derivative for sensitivity to base-flow

a.conjugate()
dChuq.multTranspose(a,c)
a.conjugate()
c.conjugate()

b.conjugate()
dChuf.multTranspose(b,d)
b.conjugate()
d.conjugate()

## Sum the three terms to get sensitivity to base-flow

f = 2./mu2 * e + 1./mu2 * c - d   #with P restriction
# f = 2. * e + 1./mu2 * c - d   #without P restriction

## Normalise with Chu energy

kspQ.solve(f,c)

## Save sensitivity to base-flow (grad_q mu^2/(2mu^2))

sensb = pet.gatherVector2ArrayPetsc(c, MPI.COMM_WORLD, broadcast=True)
sensb = _np.reshape(_np.real(sensb), (im,jm,5))
filenamesensb = out_dir + 'sensitivitytobaseflow.dat'
# filenamesensb = out_dir + 'sensitivitytobaseflowNoRestrict.dat'
toy.__writestate_center_gh(filenamesensb, im, jm, sensb, xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])

print("...Sensitivity to base-flow written...", flush=True)

## Compute sensitivity to steady forcing

f.conjugate()
kspA.solveTranspose(f,d)
d.conjugate()

## Normalise with Chu energy

kspQ.solve(d,e)

## Save sensitivity to steady forcing (grad_f mu^2/(2mu^2))

sensf = pet.gatherVector2ArrayPetsc(e, MPI.COMM_WORLD, broadcast=True)
sensf = _np.reshape(_np.real(sensf), (im,jm,5))
filenamesensf = out_dir + 'sensitivitytosteadyforcing.dat'
# filenamesensf = out_dir + 'sensitivitytosteadyforcingNoRestrict.dat'
toy.__writestate_center_gh(filenamesensf, im, jm, sensf, xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])

print("...Sensitivity to steady forcing written...", flush=True)

###################################################################################################

## Compute the first term (var. from base-flow) for sensitivity to flow parameter p

dRdp.multTranspose(d,g)

## Save the first term for sensitivity to flow parameter p (grad_p mu^2/(2mu^2) from base-flow var.)

sensp = pet.gatherVector2ArrayPetsc(g, MPI.COMM_WORLD, broadcast=True)
filenamesensp = out_dir + 'sensitivitytowalltempbsf.dat'
# filenamesensp = out_dir + 'sensitivitytowalltempbsfNoRestrict.dat'
f_out = open(filenamesensp , 'w')
f_out.write('VARIABLES= "X" "Twall" \n')
for i in range(im):
  f_out.write(str(xc[gh+i,0]) + ' ' + str(_np.real(sensp[i])) + '\n')
f_out.close()


## Construct the second Hessian dA/dp
IAhp = dic['IAdAdp'] 
JAhp = dic['JAdAdp'] 
Jachp = dic['AijdAdp'] 
Hp = pet.createMatPetscCSR(IAhp, JAhp, Jachp, im*jm*5, im, im)

## Only for 3D modes
if 3Dmode:
  IAdz = dic['IAdAdpdz'] 
  JAdz = dic['JAdAdpdz'] 
  Jacdz = dic['AijdAdpdz'] 
  IAdzz = dic['IAdAdpdz2'] 
  JAdzz = dic['JAdAdpdz2'] 
  Jacdzz = dic['AijdAdpdz2'] 
  Hpdz = pet.createMatPetscCSR(IAdz, JAdz, Jacdz, im*jm*5, im, im)
  Hpdzz = pet.createMatPetscCSR(IAdzz, JAdzz, Jacdzz, im*jm*5, im, im)
  Hp = Hp - wave*1.j*Hpdz + wave**2*Hpdzz #good

## Compute the second term (var. from jacobian) for sensitivity to flow parameter p

k, l = Hp.getVecs()
m, n = Hp.getVecs()

## with P restriction
Qchu.mult(a,c)
kspR.solveTranspose(c,d)
d.conjugate()
# Hp.multTranspose(d,k)
Hp.multTranspose(-d,k)
k.conjugate()

## without P restriction
# Qchu.mult(b,c)
# c.conjugate()
# Hp.multTranspose(-c,k)
# k.conjugate()

m = 2./mu2 * k   #with P restriction 
# m = 2. * k   #without P restriction

## Save the second term for sensitivity to flow parameter p (grad_p mu^2/(2mu^2) from jacobian var.)

sensp2 = pet.gatherVector2ArrayPetsc(m, MPI.COMM_WORLD, broadcast=True)
filenamesensp2 = out_dir + 'sensitivitytowalltempjac.dat'
# filenamesensp2 = out_dir + 'sensitivitytowalltempjacNoRestrict.dat'
f_out = open(filenamesensp2 , 'w')
f_out.write('VARIABLES= "X" "Twall" \n')
for i in range(im):
  f_out.write(str(xc[gh+i,0]) + ' ' + str(_np.real(sensp2[i])) + '\n')
f_out.close()

print("...Sensitivity to flow parameter written.", flush=True)

