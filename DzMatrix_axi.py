# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

## Compute the coupling matrices between 2 plans in z-direction (Dz and Dzz)
## For axisym.

import BROADCAST as toy
import srcfv.f_geom    as f_geom
import srcfv.f_bnd     as f_bnd
import srcfv.f_lin     as f_lin
import misc.f_misc     as f_misc
import f_init
import SIM.BLprofiles_implicit as blsim
import misc.PETSc_func as pet

import srcfv.f_dz_axi    as f_dz_axi

import restart_init as ri

import numpy as _np
# import pickle
from petsc4py import PETSc

################################################
## CARD ##

dir   = 'Wksp'

dir2  = 'dnc_5'


file  = 'state_atcenter_500_y150'
file  = 'state_atcenter_ite15'


dirout = './BASEFLOW_BL/'


dphys = dict()
dphys['Mach']    =  4.5   #4.5 
dphys['T0']      =  288.  #288.
dphys['Runit']   =  3.4e6 #3.4e6
ym               =  0.


extraporder = 2

routineout = 'bc_extrapolate_o%i_2d' % extraporder
routinein  = 'bc_supandsubinlet_2d'
# routinein  = 'bc_general_2d'
routinenr  = 'bc_no_reflexion_2d'
routinew   = 'bc_wall_viscous_adia_2d'

##################################################
## PROGRAM ##

libbnd   = 'f_bnd'

finflow  = eval("%s.%s"    % (libbnd, routinein))
foutflow = eval("%s.%s"    % (libbnd, routineout))
fnoref   = eval("%s.%s"    % (libbnd, routinenr ))
fwall    = eval("%s.%s"    % (libbnd, routinew  ))

libbnd   = 'f_lin'

routineout += '_d'
routinein  += '_d'
routinenr  += '_d'
routinew   += '_d'

flininflow  = eval("%s.%s"    % (libbnd, routinein))
flinoutflow = eval("%s.%s"    % (libbnd, routineout))
flinnoref   = eval("%s.%s"    % (libbnd, routinenr ))
flinwall    = eval("%s.%s"    % (libbnd, routinew  ))

dphys['gam']      =  1.4
dphys['cs']       =  110.4
dphys['Ts']       =  273.15      #273.15    #288
dphys['musuth']   =  1.716e-5  #1.716e-5  #1.711e-5
dphys['rgaz']     =  287.1
dphys['Prandtl']  =  0.72

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

Muref  = muinf
Lref   = Muref/(Roref*Vref)

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
xc_tmp, yc_tmp, ro, rou, rov, row, roe = ri.read_init(filet)

im = _np.shape(xc_tmp)[0]
jm = _np.shape(yc_tmp)[1]

ym *= 1./Lref

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

## More accurate if the user knows the minimal x and y points
xmin = 0.0048/Lref
ymin = ym
f_geom.bordersfromcenters_rectangular_2d(x0_tmp,y0_tmp,xc_tmp,yc_tmp,xmin,ymin)

for j in range(jm+1):
    for i in range(im+1):
        x0[i+gh,j+gh] = x0_tmp[i,j]
        y0[i+gh,j+gh] = y0_tmp[i,j]     

f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)

field = _np.zeros((jm, gh, 5), order = 'F') # dummy ones in place of blasius solution
wbd   = _np.zeros((im+gh    , 5), order = 'F') # dummy ones in place of top domain state vector

# import SIM.blasius_profiles as blasiussim
# road = _np.ones((im + 2*gh   , jm + gh     ), order='F')
# uad,vad = blasiussim.blasius_profiles(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False)
road,uad,vad,Ead = blsim.BLprofile(x0[:,:]*Lref, (y0[:,gh:]-ym)*Lref,mach, dphys, isplot=False, damped=False)

road = toy.centers_array(road)
uad  = toy.centers_array(uad)
vad  = toy.centers_array(vad)
Ead  = toy.centers_array(Ead)

w[:, gh:, 0]     = road[:,:]            * rhoinf
w[:, gh:, 1]     = road[:,:]*uad[:,:] * rhoinf * uinf
w[:, gh:, 2]     = road[:,:]*vad[:,:] * rhoinf * uinf
w[:, gh:, 4]     = road[:,:]*Ead[:,:] * rhoinf * einf

f_init.set_bndbl_2d(w, field, wbd, im)

w[gh:-gh, gh:-gh, 0]     = ro
w[gh:-gh, gh:-gh, 1]     = rou
w[gh:-gh, gh:-gh, 2]     = rov
w[gh:-gh, gh:-gh, 3]     = row
w[gh:-gh, gh:-gh, 4]     = roe

###
# import resolvent_all  as resol
# Qq = resol.computeQ_Ec(w[gh:-gh,gh:-gh,:], vol[gh:-gh,gh:-gh])
# viewer = PETSc.Viewer().createBinary(dirout+'QqEc', 'w')
# viewer(Qq)
# print 'Q Ec written'
###

#interfaces definitions (may be done at the begining)
# Ilo
interf1      = _np.zeros((2,2), order='F')
interf1[0,0] = 1  # imin
interf1[0,1] = 1  # jmin
interf1[1,0] = 1  # imax
interf1[1,1] = jm # jmax

# Ihi
interf2      = _np.zeros((2,2), order='F')
interf2[0,0] = im # imin
interf2[0,1] = 1  # jmin
interf2[1,0] = im # imax
interf2[1,1] = jm+gh # jmax

# Jlo
interf3      =  _np.zeros((2,2), order='F')
interf3[0,0] = 1-gh # imin #1-gh #i_start-gh+1
interf3[0,1] = 1  # jmin
interf3[1,0] = im+gh # imax
interf3[1,1] = 1  # jmax

# Jhi
interf4      =  _np.zeros((2,2), order='F')
interf4[0,0] = 1-gh  # imin
interf4[0,1] = jm # jmin
interf4[1,0] = im # imax
interf4[1,1] = jm # jmax

# foutflow(w,'Jhi', interf4, im, jm, gh)

# finflow(w,'Ilo', interf1, field,im,jm)    
finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm) 
fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)     
foutflow(w,'Ihi', interf2, im, jm, gh)
fwall(w,'Jlo', gam, interf3, gh, im, jm)

## To check that the 2D baseflow was converged
import srcfv.f_norm    as f_norm
import srcfv.f_sch     as f_sch
fsch = eval("f_sch.flux_num_dnc5_nowall_polar_2d")
res = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
k2 = 1.01
k4 = 1.
fsch(res, w, ym, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
norm, ninf = f_norm.compute_norml2inf(res ,im, jm, gh)
print(norm)
filename = './' + dir + '/' + dir2 + '/residualDz.dat'
toy.__writestate_center(filename, im, jm, res, xc, yc, gh)

wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
dz   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
dz2  = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5

Jacdz  = _np.zeros((nbentry), order='F')
IAdz   = _np.zeros((nbentry), dtype=_np.int32, order='F')
JAdz   = _np.zeros((nbentry), dtype=_np.int32, order='F')

Jacdz2 = _np.zeros((nbentry), order='F')
IAdz2  = _np.zeros((nbentry), dtype=_np.int32, order='F')
JAdz2  = _np.zeros((nbentry), dtype=_np.int32, order='F')

for m in range(5):
    for l in range(1 + 2*gh):
        for k in range(1 + 2*gh):
            wd *= 0.
            f_misc.testvector(wd,m,l,k,gh,im,jm)

            w[:gh,:,:]  = 0.
            w[:,:gh,:]  = 0.
            w[-gh:,:,:] = 0.
            w[:,-gh:,:] = 0.

            flininflow(w,wd,'Ilo',interf1,field,nx,ny,gam,im,jm)                
            flinnoref(w,wd,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            flinoutflow(w,wd,'Ihi', interf2, im, jm, gh)
            foutflow(w,'Ihi', interf2, im, jm, gh)
            flinwall(w,wd,'Jlo', gam, interf3, gh, im, jm)

            f_dz_axi.coeffs_5p_dz(dz, w, wd, ym, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, im, jm)
            # f_dz_axi.coeffs_9p_dz(dz, w, wd, ym, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, im, jm)
            f_dz_axi.coeffs_5p_dz2(dz2, w, wd, ym, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, im, jm)


            f_misc.computejacobianfromdz(Jacdz,IAdz,JAdz,dz,m,l,k,gh,im,jm)
            f_misc.computejacobianfromdz(Jacdz2,IAdz2,JAdz2,dz2,m,l,k,gh,im,jm)

mini = 2.e-16
IAdz, JAdz, Jacdz = toy.remove_zero_jac(IAdz, JAdz, Jacdz, mini)
IAdz2, JAdz2, Jacdz2 = toy.remove_zero_jac(IAdz2, JAdz2, Jacdz2, mini)

nbentry = _np.shape(Jacdz)[0]
# print nbentry

# import scipy.sparse as sp
# import matplotlib.pyplot as plt
# Jacs = sp.csr_matrix((Jacdz2, (IAdz2, JAdz2)), shape=(im*jm*5, im*jm*5))
# plt.figure()
# plt.spy(Jacs)
# plt.show()

# pickle.dump( [IAdz, JAdz, Jacdz], open( dirout + "AIJdz","wb") )
# pickle.dump( [IAdz2, JAdz2, Jacdz2], open( dirout + "AIJdz2","wb") )

print("** Writing 3D-component matrices **")

Dz = pet.createMatPetscCSR(IAdz, JAdz, Jacdz, im*jm*5, im*jm*5, 5*(2*gh+1)**2)
Dz2 = pet.createMatPetscCSR(IAdz2, JAdz2, Jacdz2, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

viewer = PETSc.Viewer().createBinary(dirout+'Dz', 'w')
# viewer = PETSc.Viewer().createBinary(dirout+'Dz_o8', 'w')
viewer(Dz)

viewer = PETSc.Viewer().createBinary(dirout+'Dz2', 'w')
viewer(Dz2)



