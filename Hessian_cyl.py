# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

## Compute the Hessian matrix

import BROADCAST as toy
import srcfv.f_geom    as f_geom
import srcfv.f_bnd     as f_bnd
import srcfv.f_lin     as f_lin
import misc.f_misc     as f_misc
import f_init
import SIM.BLprofiles_implicit as blsim
import misc.PETSc_func as pet

import srcfv.f_hess    as f_hess

import restart_init as ri

import numpy as _np
# import pickle
from petsc4py import PETSc
from mpi4py import MPI

################################################
## CARD ##

dir   = 'Wksp/Cylinder'

dir2  = 'dnc_5'

file = 'state_atcenter_630_y300_Re46p8'

dirout = './BASEFLOW_CYL/'

filemode = 'dnc_5/mode_630_y300Norm'
# filemode = 'dnc_5/mode20'
# filemode = 'dnc_5/mode22'

fileadj  = 'dnc_5/mode_630_y300_adjointL2'


dphys = dict()
dphys['Mach']    =  0.1   #0.1
dphys['T0']      =  288.  #288.
dphys['Runit']   =  46.8   #46.8

extraporder = 2

order = int(dir2[-1])
routinesch = 'flux_num_dnc%i_2d' % order
# routinesch = 'flux_num_dnc%i_nowall_2d' % order

routineout = 'bc_extrapolate_o%i_2d' % extraporder
routinein  = 'bc_supandsubinlet_2d'
# routinein  = 'bc_general_2d'
# routinenr  = 'bc_no_reflexion_2d'
routinenr  = 'bc_no_reflexion_forhessian_2d'
routinew   = 'bc_wall_viscous_adia_cyl_2d'
routinejn  = 'jn_match_2d'

##################################################
## PROGRAM ##

libbnd   = 'f_bnd'

finflow  = eval("%s.%s"    % (libbnd, routinein))
foutflow = eval("%s.%s"    % (libbnd, routineout))
fnoref   = eval("%s.%s"    % (libbnd, routinenr ))
fwall    = eval("%s.%s"    % (libbnd, routinew  ))
fjn      = eval("%s.%s"    % (libbnd, routinejn ))

libbnd   = 'f_lin'

routineout += '_d'
routinein  += '_d'
routinenr  += '_d'
routinew   += '_d'
routinesch += '_d'
routinejn  += '_d'

flininflow  = eval("%s.%s"    % (libbnd, routinein))
flinoutflow = eval("%s.%s"    % (libbnd, routineout))
flinnoref   = eval("%s.%s"    % (libbnd, routinenr ))
flinwall    = eval("%s.%s"    % (libbnd, routinew  ))
flinjn      = eval("%s.%s"    % (libbnd, routinejn))

libbnd   = 'f_hess'
libsch   = 'f_hess'

routineout += '_d'
routinein  += '_d'
routinenr  += '_d'
routinew   += '_d'
routinesch += '_d'
routinejn  += '_d'

flin2inflow  = eval("%s.%s"    % (libbnd, routinein))
# flin2outflow = eval("%s.%s"    % (libbnd, routineout))
flin2noref   = eval("%s.%s"    % (libbnd, routinenr ))
flin2wall    = eval("%s.%s"    % (libbnd, routinew  ))
flin2sch     = eval("%s.%s"    % (libsch, routinesch))
# flin2jn      = eval("%s.%s"    % (libbnd, routinejn))

dphys['gam']      =  1.4
dphys['cs']       =  110.4
dphys['Ts']       =  273.15      #273.15    #288
dphys['musuth']   =  1.716e-5  #1.716e-5  #1.711e-5
dphys['rgaz']     =  287.1
dphys['Prandtl']  =  0.72
k2 = 0.
k4 = 1.

gam      =  dphys['gam']
cs       =  dphys['cs']
tref     =  dphys['Ts']
muref    =  dphys['musuth']
rgaz     =  dphys['rgaz']
prandtl  =  dphys['Prandtl']
mach     =  dphys['Mach']
tinf     =  dphys['T0']
runit    =  dphys['Runit']

muinf   = toy.__comp_Sutherland(muref, tref, cs, tinf)
sound   = _np.sqrt(gam*rgaz*tinf)
uinf    = mach * sound
einf    = toy.__compute_tot_energy_inf(rgaz, gam, tinf, uinf)
rhoinf  = runit*muinf/uinf

dphys['mu0'] = muinf

cp = gam * rgaz /(gam-1.)
cv =       rgaz /(gam-1.)

Roref = rhoinf
Vref  = uinf
Tref  = tinf
Cvref = Vref**2/Tref
Rgpref = Cvref
Eref   = Vref**2

Lref   = 1.
Muref  = Roref*Vref*Lref
# Muref  = muinf
# Lref   = Muref/(Roref*Vref)

cp     = cp/Cvref
cv     = cv/Cvref
rgaz   = rgaz/Rgpref
tref  = tref/Tref
muref = muref/Muref
cs    = cs/Tref

uinf   = uinf/Vref
tinf   = tinf/Tref
rhoinf = rhoinf/Roref
einf   = einf/Eref
muinf  = muinf/Muref


filet = './' + dir + '/' + dir2 + '/' + file + '.dat'

im = 630   #600 #630
jm = 300   #300
BLprof = _np.loadtxt(filet,comments=('#','ZONE'),skiprows=3)
xc_tmp = _np.reshape(BLprof[:,0],(im,jm), order='F')
yc_tmp = _np.reshape(BLprof[:,1],(im,jm), order='F')

gh = (int(dir2[-1]) + 1) / 2

x0_tmp = _np.zeros((im+ 1, jm+ 1   ), order='F')
y0_tmp = _np.zeros((im+ 1, jm+ 1   ), order='F')

x0  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1   ), order='F')
y0  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1   ), order='F')
xc  = _np.zeros((im + 2*gh    , jm + 2*gh       ), order='F')
yc  = _np.zeros((im + 2*gh    , jm + 2*gh       ), order='F')
nx  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1, 2), order='F')
ny  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1, 2), order='F')
vol = _np.zeros((im + 2*gh    , jm + 2*gh       ), order='F')
volf= _np.zeros((im + 2*gh    , jm + 2*gh    , 2), order='F')
w   = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')

f_geom.bordersfromcenters_2d(x0_tmp,y0_tmp,xc_tmp,yc_tmp,im,jm)


## Important for accuracy: MESH to load
filet = './' + dir + '/' + dir2 + '/mesh_atnode.dat'
BLprof = _np.loadtxt(filet,comments=('#','ZONE'),skiprows=3)
x0_tmp = _np.reshape(BLprof[:,0],(im+1,jm+1), order='F')
y0_tmp = _np.reshape(BLprof[:,1],(im+1,jm+1), order='F')


jmin = 1-gh
jmax = jm+gh
prr1 = _np.array([[im+1  , jmin], [im+gh , jmax]], order='F')
prd1 = _np.array([[im-gh+1 , jmin], [im    , jmax]], order='F')
prr2 = _np.array([[ 1-gh , jmin], [ 0    , jmax]], order='F')
prd2 = _np.array([[ 1    , jmin], [ gh , jmax]], order='F')
tr1  = _np.array([1,2], order='F')
tr2  = _np.array([1,2], order='F')

# Compute Geometry
for j in range(jm+1):
    for i in range(im+1):
        x0[i+gh,j+gh] = x0_tmp[i,j]
        y0[i+gh,j+gh] = y0_tmp[i,j]    

f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)

wbd   = _np.zeros((im    , 5), order = 'F') # dummy ones in place of top domain state vector

f_bnd.jn_match_geom_2d(xc,prr1,gh,gh,gh,gh,im,jm,xc,prd2,gh,gh,gh,gh,im,jm,tr1)
f_bnd.jn_match_geom_2d(xc,prr2,gh,gh,gh,gh,im,jm,xc,prd1,gh,gh,gh,gh,im,jm,tr2)

f_bnd.jn_match_geom_2d(yc,prr1,gh,gh,gh,gh,im,jm,yc,prd2,gh,gh,gh,gh,im,jm,tr1)
f_bnd.jn_match_geom_2d(yc,prr2,gh,gh,gh,gh,im,jm,yc,prd1,gh,gh,gh,gh,im,jm,tr2)

f_bnd.jn_match_geom_2d(vol,prr1,gh,gh,gh,gh,im,jm,vol,prd2,gh,gh,gh,gh,im,jm,tr1)
f_bnd.jn_match_geom_2d(vol,prr2,gh,gh,gh,gh,im,jm,vol,prd1,gh,gh,gh,gh,im,jm,tr2)

nx[im+gh,:,1] = nx[gh,:,1]
ny[im+gh,:,1] = ny[gh,:,1]
nx[im+gh,:,0] = nx[gh,:,0]
ny[im+gh,:,0] = ny[gh,:,0]

prr1g = _np.array([[im+1+1  , jmin], [im+1+gh , jmax+1]], order='F')
prd1g = _np.array([[im+1-gh , jmin], [im    , jmax+1]], order='F')
prr2g = _np.array([[ 1-gh , jmin], [ 0    , jmax+1]], order='F')
prd2g = _np.array([[ 2    , jmin], [ 1+gh , jmax+1]], order='F')

fjn(nx,prr1g,gh,gh,gh,gh,im+1,jm+1,nx,prd2g,gh,gh,gh,gh,im+1,jm+1,tr1)
fjn(nx,prr2g,gh,gh,gh,gh,im+1,jm+1,nx,prd1g,gh,gh,gh,gh,im+1,jm+1,tr2)
fjn(ny,prr1g,gh,gh,gh,gh,im+1,jm+1,ny,prd2g,gh,gh,gh,gh,im+1,jm+1,tr1)
fjn(ny,prr2g,gh,gh,gh,gh,im+1,jm+1,ny,prd1g,gh,gh,gh,gh,im+1,jm+1,tr2)

for j in range(gh,jm+1+gh):
    for i in range(gh,im+1+gh):
        volf[i,j,0] = 2./(vol[i,j]+vol[i-1,j])
        volf[i,j,1] = 2./(vol[i,j]+vol[i,j-1])

### O mesh for the cylinder: #UNNECESSARY EXCEPT FOR PLOT
x0[:gh, :] = x0[im+1:im+1+gh,:]
y0[:gh, :] = y0[im+1:im+1+gh,:]
x0[im+1+gh:, :] = x0[gh:2*gh,:]
y0[im+1+gh:, :] = y0[gh:2*gh,:]

imold = im   #630
jmold = jm   #300
filet = './' + dir + '/' + dir2 + '/' + file + '.dat'
BLprof = _np.loadtxt(filet,comments=('#','ZONE'),skiprows=3)
Xin   = _np.reshape(BLprof[:,0],(imold,jmold), order='F')
Yin   = _np.reshape(BLprof[:,1],(imold,jmold), order='F')
roin  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
rouin = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
rovin = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
rowin = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
roein = _np.reshape(BLprof[:,6],(imold,jmold), order='F')
w[gh:-gh, gh:-gh, 0]     = roin
w[gh:-gh, gh:-gh, 1]     = rouin
w[gh:-gh, gh:-gh, 2]     = rovin
w[gh:-gh, gh:-gh, 3]     = rowin
w[gh:-gh, gh:-gh, 4]     = roein

wbd[:, 0]      = 1.
wbd[:, 1]      = 1.
wbd[:, 2]      = 0.
wbd[:, 3]      = 0.
wbd[:, 4]      = einf

#interfaces definitions (may be done at the begining)
# Ilo
interf1      = _np.zeros((2,2), order='F')
interf1[0,0] = 1  # imin
interf1[0,1] = 1-gh  # jmin
interf1[1,0] = 1  # imax
interf1[1,1] = jm+gh # jmax

# Ihi
interf2      = _np.zeros((2,2), order='F')
interf2[0,0] = im # imin
interf2[0,1] = 1-gh  # jmin
interf2[1,0] = im # imax
interf2[1,1] = jm+gh # jmax

# Jlo
interf3      =  _np.zeros((2,2), order='F')
interf3[0,0] = 1 # imin #1-gh #i_start-gh+1
interf3[0,1] = 1  # jmin
interf3[1,0] = im # imax
interf3[1,1] = 1  # jmax

# Jhi
interf4      =  _np.zeros((2,2), order='F')
interf4[0,0] = 1  # imin
interf4[0,1] = jm # jmin
interf4[1,0] = im # imax
interf4[1,1] = jm # jmax


# fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm) 
fnoref(w,w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm) 
# fwall(w,'Jlo', gam, interf3, gh, im, jm)
fwall(w,w,'Jlo', gam, interf3, gh, im, jm)
fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

filemodet = './' + dir + '/' + filemode 
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

wmoder = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
wmoder[gh:-gh,gh:-gh,0] = ro
wmoder[gh:-gh,gh:-gh,1] = rou 
wmoder[gh:-gh,gh:-gh,2] = rov 
wmoder[gh:-gh,gh:-gh,3] = row 
wmoder[gh:-gh,gh:-gh,4] = roe             
# flinnoref(w, wmoder,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
flinnoref(w, wmoder, w, wmoder,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
# flinwall(w, wmoder,'Jlo', gam, interf3, gh, im, jm)
flinwall(w, wmoder, w, wmoder,'Jlo', gam, interf3, gh, im, jm)
wmoder2 = _np.copy(wmoder, order='F')
flinjn(w, wmoder2, prr1, gh, gh, gh, gh, im, jm, w, wmoder, prd2, gh, gh, gh, gh, im, jm, tr1)
wmoder[-gh:,:,:] = wmoder2[-gh:,:,:]
wmoder2 = _np.copy(wmoder, order='F')
flinjn(w, wmoder2, prr2, gh, gh, gh, gh, im, jm, w, wmoder, prd1, gh, gh, gh, gh, im, jm, tr2)
wmoder[:gh,:,:] = wmoder2[:gh,:,:]
wmodei = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
wmodei[gh:-gh,gh:-gh,0] = roi
wmodei[gh:-gh,gh:-gh,1] = roui
wmodei[gh:-gh,gh:-gh,2] = rovi
wmodei[gh:-gh,gh:-gh,3] = rowi
wmodei[gh:-gh,gh:-gh,4] = roei           
# flinnoref(w, wmodei,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
flinnoref(w, wmodei, w, wmodei,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
# flinwall(w, wmodei,'Jlo', gam, interf3, gh, im, jm)
flinwall(w, wmodei, w, wmodei,'Jlo', gam, interf3, gh, im, jm)
wmodei2 = _np.copy(wmodei, order='F')
flinjn(w, wmodei2, prr1, gh, gh, gh, gh, im, jm, w, wmodei, prd2, gh, gh, gh, gh, im, jm, tr1)
wmodei[-gh:,:,:] = wmodei2[-gh:,:,:]
wmodei2 = _np.copy(wmodei, order='F')
flinjn(w, wmodei2, prr2, gh, gh, gh, gh, im, jm, w, wmodei, prd1, gh, gh, gh, gh, im, jm, tr2)
wmodei[:gh,:,:] = wmodei2[:gh,:,:]

## To check that the 2D baseflow was converged
import srcfv.f_norm    as f_norm
import srcfv.f_sch     as f_sch
fsch = eval("f_sch.flux_num_dnc5_2d")
res = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
norm, ninf = f_norm.compute_norml2inf(res ,im, jm, gh)
print(norm)
filename = './' + dir + '/' + dir2 + '/residualDz.dat'
toy.__writestate_center(filename, im, jm, res, xc, yc, gh)

wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
dhr  = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
dhi  = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
res  = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
resd = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5

Jach  = _np.zeros((nbentry), dtype=_np.complex128, order='F')
Jachr = _np.zeros((nbentry), order='F')
Jachi = _np.zeros((nbentry), order='F')
IAh   = _np.zeros((nbentry), dtype=_np.int32, order='F')
JAh   = _np.zeros((nbentry), dtype=_np.int32, order='F')

for m in range(5):
    for l in range(1 + 2*gh):
        for k in range(1 + 2*gh):
            wd *= 0.
            f_misc.testvector(wd,m,l,k,gh,im,jm)

            w[:gh,:,:]  = 0.
            w[:,:gh,:]  = 0.
            w[-gh:,:,:] = 0.
            w[:,-gh:,:] = 0.

            wdr = _np.copy(wd, order='F')               
            # flin2noref(w, wmoder, wmoder, wdr,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            flin2noref(w, w, wdr, w, wmoder, wdr,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            wd[:,-gh:,:] = wdr[:,-gh:,:]
            # fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm) 
            fnoref(w,w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)    
            wdr = _np.copy(wd, order='F')
            # flin2wall(w, wmoder, wmoder, wdr,'Jlo', gam, interf3, gh, im, jm)
            flin2wall(w, w, wdr, w, wmoder, wdr,'Jlo', gam, interf3, gh, im, jm)
            wd[:,:gh,:] = wdr[:,:gh,:]
            # fwall(w,'Jlo', gam, interf3, gh, im, jm)
            fwall(w,w,'Jlo', gam, interf3, gh, im, jm)
            fjn(wd,prr1,gh,gh,gh,gh,im,jm,wd,prd2,gh,gh,gh,gh,im,jm,tr1)
            fjn(wd,prr2,gh,gh,gh,gh,im,jm,wd,prd1,gh,gh,gh,gh,im,jm,tr2)
            fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
            fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)
            

            # flin2sch(res, resd, dhr, w, wd, wmoder, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
            flin2sch(res, resd, dhr, w, wmoder, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)

            wd *= 0.
            f_misc.testvector(wd,m,l,k,gh,im,jm)

            w[:gh,:,:]  = 0.
            w[:,:gh,:]  = 0.
            w[-gh:,:,:] = 0.
            w[:,-gh:,:] = 0.
 
            wdr = _np.copy(wd, order='F')              
            # flin2noref(w, wmodei, wmodei, wdr,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            flin2noref(w, w, wdr, w, wmodei, wdr,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            wd[:,-gh:,:] = wdr[:,-gh:,:]
            # fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)  
            fnoref(w,w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)     
            wdr = _np.copy(wd, order='F')
            # flin2wall(w, wmodei, wmodei, wdr,'Jlo', gam, interf3, gh, im, jm)
            flin2wall(w, w, wdr, w, wmodei, wdr,'Jlo', gam, interf3, gh, im, jm)
            wd[:,:gh,:] = wdr[:,:gh,:]
            # fwall(w,'Jlo', gam, interf3, gh, im, jm)
            fwall(w,w,'Jlo', gam, interf3, gh, im, jm)
            fjn(wd,prr1,gh,gh,gh,gh,im,jm,wd,prd2,gh,gh,gh,gh,im,jm,tr1)
            fjn(wd,prr2,gh,gh,gh,gh,im,jm,wd,prd1,gh,gh,gh,gh,im,jm,tr2)
            fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
            fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)
            

            # flin2sch(res, resd, dhi, w, wd, wmodei, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)  
            flin2sch(res, resd, dhi, w, wmodei, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)                               


            # f_misc.computejacobianfromjv_withjn(Jachr,IAh,JAh,dhr,m,l,k,gh,im,jm)
            # f_misc.computejacobianfromjv_withjn(Jachi,IAh,JAh,dhi,m,l,k,gh,im,jm)
            f_misc.computejacobianfromjv_withjn_dbyvol(Jachr,IAh,JAh,dhr,m,l,k,gh,im,jm,vol)
            f_misc.computejacobianfromjv_withjn_dbyvol(Jachi,IAh,JAh,dhi,m,l,k,gh,im,jm,vol)


Jach = Jachr + 1.j*Jachi

mini = 2.e-16
IAh, JAh, Jach = toy.remove_zero_jac(IAh, JAh, Jach, mini)

# print nbentry
nbentry = _np.shape(Jach)[0]
# print nbentry

# import scipy.sparse as sp
# import matplotlib.pyplot as plt
# Jacs = sp.csr_matrix((Jacdz2, (IAdz2, JAdz2)), shape=(im*jm*5, im*jm*5))
# plt.figure()
# plt.spy(Jacs)
# plt.show()

print("** Writing Hessian matrix **")

H = pet.createMatPetscCSR(IAh, JAh, Jach, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

# viewer = PETSc.Viewer().createBinary(dirout+'H', 'w')
viewer = PETSc.Viewer().createBinary(dirout+'H2', 'w')
# viewer = PETSc.Viewer().createBinary(dirout+'Hmode20', 'w')
# viewer = PETSc.Viewer().createBinary(dirout+'Hmode22', 'w')
viewer(H)


#########

# viewer = PETSc.Viewer().createBinary(dirout+'H', PETSc.Viewer.Mode.READ)
# H=PETSc.Mat()
# H.create(PETSc.COMM_WORLD)
# H.setType('mpiaij')
# H.load(viewer)
# H.assemble()

import misc.PETSc_func as pet

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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

a, b = H.getVecs()

rangeVec = b.getOwnershipRange()
for k in range(rangeVec[0],rangeVec[1]):
  b[k] = adjoint[k]
b.assemble()

b.conjugate()
H.multTranspose(b,a)

sol = pet.gatherVector2ArrayPetsc(a,MPI.COMM_WORLD,broadcast=True)

sensi = _np.reshape( _np.real(sol), (im,jm,5))
filename1 = './' + dir + '/' + dir2 + '/sensitivity_atcenter_eig_n%i_real.dat' % 0
toy.__writestate_center_gh(filename1, im, jm, sensi, xc_tmp, yc_tmp)

