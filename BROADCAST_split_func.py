# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#!/usr/bin/env python

'''
File: Toy2D.py

Created on 21 january 2021

@author:       Cedric Content
@contact:      cedric.content@onera.fr
@organization: ONERA - DAAA

@summary:      This file is the main file of the program. It contains the
               routine "main" and other related routines.
'''
import srcfv.f_geom    as f_geom
import srcfv.f_bnd     as f_bnd
import srcfv.f_sch     as f_sch
import srcfv.f_lhs     as f_lhs
import srcfv.f_lin     as f_lin
import srcfv.f_adj     as f_adj
import srcfv.f_norm    as f_norm
# FROM A.POULAIN Thesis
import misc.f_misc     as f_misc
import misc.PETSc_func as pet
from mpi4py import MPI
import resolvent_all  as resol
import SIM
import SIM.BLprofiles_implicit as blsim
import f_init
import meshBL as mesh

import numpy as _np
import matplotlib.pyplot as plt

import os
import sys
import timeit



######################### Private functions ####################
def __writestate_node(filename, im, jm, w, x0, y0, gh) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "ro" "rou" "rov" "row" "roe" \n')
    f_out.write('ZONE I = ' + str(im) + ',  J = ' + str(jm) + '\n')
    for j in range(gh,jm+gh):
        for i in range(gh,im+gh):
            ro  = 0.25*(w[i-1,j-1,0] + w[i,j-1,0] + w[i-1,j,0] + w[i,j,0])
            rou = 0.25*(w[i-1,j-1,1] + w[i,j-1,1] + w[i-1,j,1] + w[i,j,1])
            rov = 0.25*(w[i-1,j-1,2] + w[i,j-1,2] + w[i-1,j,2] + w[i,j,2])
            row = 0.25*(w[i-1,j-1,3] + w[i,j-1,3] + w[i-1,j,3] + w[i,j,3])
            roe = 0.25*(w[i-1,j-1,4] + w[i,j-1,4] + w[i-1,j,4] + w[i,j,4])
            f_out.write(str(x0[i,j]) + ' ' + str(y0[i,j]) + ' ' +
                        str(ro)    + ' ' + str(rou)   + ' ' +
                        str(rov)   + ' ' + str(row)   + ' ' +
                        str(roe)   + '\n')
    f_out.close()

def __writestate_center(filename, im, jm, w, xc, yc, gh) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "ro" "rou" "rov" "row" "roe" \n')
    f_out.write('ZONE I = ' + str(im) + ',  J = ' + str(jm) + '\n')
    for j in range(gh,jm+gh):
        for i in range(gh,im+gh):
            ro  = w[i,j,0]
            rou = w[i,j,1]
            rov = w[i,j,2]
            row = w[i,j,3]
            roe = w[i,j,4]
            f_out.write(str(xc[i,j]) + ' ' + str(yc[i,j]) + ' ' +
                        str(ro)    + ' ' + str(rou)   + ' ' +
                        str(rov)   + ' ' + str(row)   + ' ' +
                        str(roe)   + '\n')
    f_out.close()

def __writestate_center_gh(filename, imloc, jmloc, w, xc, yc) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "ro" "rou" "rov" "row" "roe" \n')
    f_out.write('ZONE I = ' + str(imloc) + ',  J = ' + str(jmloc) + '\n')
    for j in range(jmloc):
        for i in range(imloc):
            ro  = w[i,j,0]
            rou = w[i,j,1]
            rov = w[i,j,2]
            row = w[i,j,3]
            roe = w[i,j,4]
            f_out.write(str(xc[i,j]) + ' ' + str(yc[i,j]) + ' ' +
                        str(ro)      + ' ' + str(rou)     + ' ' +
                        str(rov)     + ' ' + str(row)     + ' ' +
                        str(roe)     + '\n')
    f_out.close()

def __writeline(filename, imloc, w, xc,jloc) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "ro" "rou" "rov" "row" "roe" \n')
    f_out.write('ZONE I = ' + str(imloc)  + '\n')
    for i in range(imloc):
        ro  = w[i,0]
        rou = w[i,1]
        rov = w[i,2]
        row = w[i,3]
        roe = w[i,4]
        f_out.write(str(xc[i,jloc])  + ' ' +
                    str(ro)    + ' ' + str(rou)   + ' ' +
                    str(rov)   + ' ' + str(row)   + ' ' +
                    str(roe)   + '\n')
    f_out.close()

def __comp_Sutherland(propref, Ts, Cs, T):
    '''Dynamical viscosity / thermal conductivity from sutherland law'''
    return propref*_np.sqrt(T/Ts)*((1.+Cs/Ts)/(1.+Cs/T))

def __compute_tot_energy_inf(R_pg, gamma, t_inf, v_inf):
    '''Total energy E = R/(gamma-1)*Tinf+(uinf**2)/2'''
    return R_pg/(gamma-1.)*t_inf+0.5*v_inf*v_inf

def remove_zero_jac(IA, JA, Jac, mini=2.e-16):
    ''' Remove the zero components from the Jac list in order not to store any zero in the sparse matrix '''
    to_keep = _np.absolute(Jac) > mini
    Jac = Jac[to_keep,...]
    IA  = IA[to_keep,...]
    JA  = JA[to_keep,...]
    return IA, JA, Jac

def centers_array(A):
    ''' Compute the values of an array A at the centers'''
    return 0.25 * ( A[:-1,:-1] + A[1:,:-1] + A[:-1,1:] + A[1:,1:] )    

def read_input_dict_functions(lf):

    # get functions
    # lf = [finflow, foutflow, fnoref, fwall, fsch]
    finflow  = lf[0]
    foutflow = lf[1]
    fnoref   = lf[2]
    fwall    = lf[3]
    fsch     = lf[4]

    return finflow, foutflow, fnoref, fwall, fsch

def read_input_dict_geom(dgeom, dnum):
    ''' Compute/read only geometry '''
    im       = dgeom['im']
    jm       = dgeom['jm']
    L        = dgeom['length']
    high     = dgeom['high']
    xini     = dgeom['xini']
    sch      = dnum['sch']
    order    = dnum['order']
    if sch == 'dnc':
        gh = (order+1) / 2
    else:
        gh = (order-1) / 2 + 1

    return im, jm, gh, L, high, xini 

def read_input_dict_num(dnum):

    ite      = dnum['ite']
    cfl      = dnum['cfl']
    k2       = dnum['k2']
    k4       = dnum['k4']
    sch      = dnum['sch']
    order    = dnum['order']
    freqres  = dnum['freqres']
    freqsort = dnum['freqsort']

    return ite, cfl, k2, k4, sch, order, freqres, freqsort

def read_input_dict_phys(dphys):  

    # get physical constants
    gam      =  dphys['gam']
    cs       =  dphys['cs']
    tref     =  dphys['Ts']
    muref    =  dphys['musuth']
    rgaz     =  dphys['rgaz']
    prandtl  =  dphys['Prandtl']
    mach     =  dphys['Mach']
    tinf     =  dphys['T0']
    Lref     =  dphys['Lref']
    StateRef =  dphys['stateref'] 

    return gam, cs, tref, muref, rgaz, prandtl, mach, tinf, Lref, StateRef 

def compute_infq(muref, tref, cs, tinf, gam, rgaz, mach, StateRef, dphys):
    ''' compute the remaining far-field variables (those not included in the dictionnary dphys) '''
    muinf   = __comp_Sutherland(muref, tref, cs, tinf)
    sound   = _np.sqrt(gam*rgaz*tinf)
    uinf    = mach * sound
    einf    = __compute_tot_energy_inf(rgaz, gam, tinf, uinf)

    if StateRef == 'm0_p0_t0':
        pinf    =  dphys['P0']
        rhoinf  = pinf/(rgaz*tinf)
        runit   = rhoinf * uinf/muinf

    elif StateRef == 'm0_runit_t0':
        runit   =  dphys['Runit']
        rhoinf  = runit*muinf/uinf
        pinf    = rhoinf*rgaz*tinf

    # rey = runit * L
    cp = gam * rgaz /(gam-1.)
    cv =       rgaz /(gam-1.)  

    #for similitude sol
    dphys['Runit'] = runit
    dphys['mu0'] = muinf

    return muinf, sound, uinf, einf, pinf, rhoinf, runit, cp, cv, dphys

def construct_mesh(xini, L, high, im, jm):
    ''' Create the mesh '''
    x  = _np.linspace(xini, xini+L , im+1)
    ## MESH v2
    Ny_in   = 80*jm/100   #65%    #80%
    deltaBL = high/4      #high/7 #high/4 
    percent = 0.02 
    
    Ny_out  = jm - Ny_in 
    Nend    = high/deltaBL
    y_int   = mesh.bigeom_stretch_in(Ny_in, deltaBL, percent)
    y_out   = mesh.exp_stretch_out(Ny_out, deltaBL, percent, Nend)
    y       = _np.concatenate((y_int, y_out)) 

    return x, y

def initialize_fields(im, jm, gh):

    # Initialize all cfd fields
    x0  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1   ), order='F')
    y0  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1   ), order='F')
    xc  = _np.zeros((im + 2*gh    , jm + 2*gh       ), order='F')
    yc  = _np.zeros((im + 2*gh    , jm + 2*gh       ), order='F')
    nx  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1, 2), order='F')
    ny  = _np.zeros((im + 2*gh + 1, jm + 2*gh + 1, 2), order='F')
    vol = _np.zeros((im + 2*gh    , jm + 2*gh       ), order='F')
    volf= _np.zeros((im + 2*gh    , jm + 2*gh    , 2), order='F')
    w   = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
    res = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')

    return x0, y0, xc, yc, nx, ny, vol, volf, w, res

def compute_Timestep(L, im, Lref, y, mach, cfl):

    dx = L/((im-1)*Lref) # adim done after muinf A.Poulain
    dy = (y[1]-y[0])/Lref
    sound = 1./mach
    dt = cfl * min(dy,dx) / (sound+1.)

    return dt

def compute_Timestep_vadim(gh, yc, mach, cfl):

    dt  = cfl * (yc[gh,gh+1] - yc[gh,gh]) / (1./mach + 1.)
    # dt  = cfl * min(_np.amin(yc[gh+1:-gh]-yc[gh:-gh-1]),_np.amin(xc[gh+1:-gh]-xc[gh:-gh-1])) / (1./mach + 1.)

    return dt

def compute_GlobalTimestep_vadim(gh, x0, y0, mach, cfl):

    dt  = cfl * min(_np.amin(y0[gh+1:-gh]-y0[gh:-gh-1]),_np.amin(x0[gh+1:-gh]-x0[gh:-gh-1])) / (1./mach + 1.)

    return dt    

def compute_LocalTimestep_vadim(gh, x0, y0, mach, cfl):

    # dt  = cfl * _np.minimum(_np.abs(y0[:-1,1:]-y0[:-1,:-1]),_np.abs(x0[1:,:-1]-x0[:-1,:-1])) / (1./mach + 1.)
    dt  = cfl * _np.minimum( _np.sqrt((x0[:-1,1:]-x0[:-1,:-1])**2 + (y0[:-1,1:]-y0[:-1,:-1])**2), _np.sqrt((x0[1:,:-1]-x0[:-1,:-1])**2 + (y0[1:,:-1]-y0[:-1,:-1])**2) ) / (1./mach + 1.)

    return dt   

def compute_LocalTimestep_vadim_vol(gh, vol, mach, cfl):

    dt  = cfl * _np.sqrt(_np.abs(vol)) / (1./mach + 1.)

    return dt      

def compute_LocalTimestep_vadim_cyl(gh, x0, y0, mach, cfl):

    dt  = cfl * _np.maximum(_np.minimum( _np.sqrt((x0[:-1,1:]-x0[:-1,:-1])**2 + (y0[:-1,1:]-y0[:-1,:-1])**2), _np.sqrt((x0[1:,:-1]-x0[:-1,:-1])**2 + (y0[1:,:-1]-y0[:-1,:-1])**2) ), 1.e-15*_np.ones((_np.shape(x0)[0]-1,_np.shape(y0)[1]-1))) / (1./mach + 1.)

    return dt     

def computeAdimRVT(rhoinf, uinf, tinf, muinf, pinf, cp, cv, rgaz, einf, tref, muref, cs):
    ''' Adim (by RVT = rho, velo et temperature) '''

    Roref = rhoinf
    Vref  = uinf
    Tref  = tinf
    Pref  = Roref*Vref**2
    Cvref = Vref**2/Tref

    Eref   = Vref**2
    Rgpref = Cvref

    ## Adim with ref length
    # Lref   = 3.24e-3  #3.24e-3
    # Muref  = Roref*Vref*Lref
    ## OR Adim with unit Reynolds
    Muref  = muinf
    Lref   = Muref/(Roref*Vref)

    uinf   = uinf/Vref
    tinf   = tinf/Tref
    rhoinf = rhoinf/Roref
    # sound  = sound/Vref
    pinf   = pinf/Pref
    cp     = cp/Cvref
    cv     = cv/Cvref
    rgaz   = rgaz/Rgpref
    einf   = einf/Eref
    # sutherland
    tref  = tref/Tref
    muref = muref/Muref
    cs    = cs/Tref
    muinf = muinf/Muref

    return Lref, uinf, tinf, rhoinf, pinf, cp, cv, rgaz, einf, tref, muref, cs, muinf

def compute_geometryAdim(im, jm, gh, x, y, Lref, x0,y0,nx,ny,xc,yc,vol,volf):

    # Compute Geometry
    for i in range(im+1):
        x0[i+gh,:] = x[i]
    for j in range(jm+1):
        y0[:,j+gh] = y[j]

    # Adim Geom:
    x0 *= 1./Lref
    y0 *= 1./Lref

    f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)

    return x0, y0, nx, ny, xc, yc, vol, volf

def initialisationSelfSim(im, jm, gh, w, x0, y0, Lref, mach, dphys, rhoinf, uinf, einf, xc, yc, cv):
    ''' Initialize all the field + inlet and top BC with a self-similar solution '''

    state_adim    = _np.zeros(5)
    state_adim[0] = rhoinf
    state_adim[1] = rhoinf * uinf
    state_adim[2] = 0.
    state_adim[3] = 0.
    state_adim[4] = rhoinf * einf

    # Blasius for inlet

    field = _np.zeros((jm, gh, 5), order = 'F') # dummy ones in place of blasius solution
    wbd   = _np.zeros((im+gh    , 5), order = 'F') # dummy ones in place of top domain state vector

    # Initialization
    wbd[:, 0]      = state_adim[0]
    wbd[:, 1]      = state_adim[1]
    wbd[:, 2]      = state_adim[2]
    wbd[:, 3]      = state_adim[3]
    wbd[:, 4]      = state_adim[4]
    # print 'Wb shape at 2=', _np.shape(wbd)

    field[:, :, 0] = state_adim[0]
    field[:, :, 1] = state_adim[1]
    field[:, :, 2] = state_adim[2]
    field[:, :, 3] = state_adim[3]
    field[:, :, 4] = state_adim[4]

    # Initialize(field, w, )

    w[:, :, 0]     = state_adim[0]
    w[:, :, 1]     = state_adim[1]
    w[:, :, 2]     = state_adim[2]
    w[:, :, 3]     = state_adim[3]
    w[:, :, 4]     = state_adim[4]

    # Initialise from A.Poulain routine
    ## Compressible self-similar profile
    road,uad,vad,Ead = blsim.BLprofile(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False)
    ## OR Incompressible blasius profile 
    # # import SIM.blasius_profiles as blasiussim
    # # road = _np.ones((im + 2*gh   , jm + gh     ), order='F')
    # # uad,vad = blasiussim.blasius_profiles(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False)

    road = centers_array(road)
    uad  = centers_array(uad)
    vad  = centers_array(vad)
    Ead  = centers_array(Ead)

    w[:, gh:, 0]     = road[:,:]          * rhoinf
    w[:, gh:, 1]     = road[:,:]*uad[:,:] * rhoinf * uinf
    w[:, gh:, 2]     = road[:,:]*vad[:,:] * rhoinf * uinf
    w[:, gh:, 4]     = road[:,:]*Ead[:,:] * rhoinf * einf

    f_init.set_bndbl_2d(w, field, wbd, im)

    return w, wbd, field

def define_interf(im, jm, gh):
    ''' Define the interfaces of the BC '''

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

    return interf1, interf2, interf3, interf4

def read_lflin(lflin):

    flininflow  = lflin[0]
    flinoutflow = lflin[1]
    flinnoref   = lflin[2]
    flinwall    = lflin[3]
    flinsch     = lflin[4]
    # flinjac     = lflin[5]

    return   flininflow, flinoutflow, flinnoref, flinwall, flinsch

def updatefromBC(w, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, finflow=None, foutflow=None, fnoref=None, fwall=None, fjn=None, prr1=None, prr2=None, prd1=None, prd2=None, tr1=None, tr2=None, pinf=None):
    ''' Update the state w from the BC '''
    
    # finflow(w,'Ilo', interf1, field,im,jm)          
    finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
    fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
    # foutflow(w,'Ihi', interf2, im, jm, gh)
    foutflow(w,'Ihi',interf2,pinf,True,gam,nx,ny,im,jm,gh)
    fwall(w,'Jlo', gam, interf3, gh, im, jm)

    ## CYL
    # fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
    # fwall(w,'Jlo', gam, interf1, gh, im, jm)
    # fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
    # fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

    return w

def updatefromBC_cyl(w, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, finflow=None, foutflow=None, fnoref=None, fwall=None, fjn=None, prr1=None, prr2=None, prd1=None, prd2=None, tr1=None, tr2=None):
    ''' Update the state w from the BC '''

    ## CYL
    fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
    fwall(w,'Jlo', gam, interf1, gh, im, jm)
    fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
    fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

    return w    

def computeRes(sch, fsch, res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None):
    ''' Compute the opposite of the residual (not divided by the volume) '''
    if sch == 'dnc':
        fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
    else:
        fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

    return res

def computeRes_dbyvol(sch, fsch, res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None):
    ''' Compute the "classic" residual '''
    if sch == 'dnc':
        fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
    else:
        fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

    for m in range(_np.shape(res)[2]):
        res[gh:-gh,gh:-gh,m] = - res[gh:-gh,gh:-gh,m] / vol[gh:-gh,gh:-gh]

    return res 

def updatefromBC_lin(w, wd, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, flininflow=None, flinoutflow=None, flinnoref=None, flinwall=None, finflow=None, foutflow=None, fjn=None, prr1=None, prr2=None, prd1=None, prd2=None, tr1=None, tr2=None, pinf=None):
    ''' Update the state w from the BC '''

    # w[:gh,:,:]  = 0.
    # w[:,:gh,:]  = 0.
    # w[-gh:,:,:] = 0.
    # w[:,-gh:,:] = 0.

    flininflow(w,wd,'Ilo',interf1,field,nx,ny,gam,im,jm)                
    # finflow(w,'Ilo', interf1, field,im,jm) 
    flinnoref(w,wd,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
    # flinoutflow(w,wd,'Ihi', interf2, im, jm, gh)
    flinoutflow(w,wd,'Ihi',interf2,pinf,True,gam,nx,ny,im,jm,gh)
    # foutflow(w,'Ihi', interf2, im, jm, gh)
    flinwall(w,wd,'Jlo', gam, interf3, gh, im, jm)

    ## CYL
    # flinnoref(w,wd,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
    # flinwall(w,wd,'Jlo', gam, interf1, gh, im, jm)
    # fjn(wd,prr1,gh,gh,gh,gh,im,jm,wd,prd2,gh,gh,gh,gh,im,jm,tr1)
    # fjn(wd,prr2,gh,gh,gh,gh,im,jm,wd,prd1,gh,gh,gh,gh,im,jm,tr2)
    # fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
    # fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

    return w, wd  

def updatefromBC_adj(w, wd, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, flininflow=None, flinoutflow=None, flinnoref=None, flinwall=None, finflow=None, foutflow=None, fjn=None, prr1=None, prr2=None, prd1=None, prd2=None, tr1=None, tr2=None, pinf=None):
    ''' Update the state w from the BC '''

    flinwall(w,wd,'Jlo', gam, interf3, gh, im, jm) 
    flinoutflow(w,wd,'Ihi',interf2,pinf,True,gam,nx,ny,im,jm,gh)
    flinnoref(w,wd,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
    flininflow(w,wd,'Ilo',interf1,field,nx,ny,gam,im,jm)  
    
    return w, wd    

def updatefromBC_lin_cyl(w, wd, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, flininflow=None, flinoutflow=None, flinnoref=None, flinwall=None, finflow=None, foutflow=None, fjn=None, prr1=None, prr2=None, prd1=None, prd2=None, tr1=None, tr2=None):
    ''' Update the state w from the BC '''

    ## CYL
    flinnoref(w,wd,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
    flinwall(w,wd,'Jlo', gam, interf1, gh, im, jm)
    fjn(wd,prr1,gh,gh,gh,gh,im,jm,wd,prd2,gh,gh,gh,gh,im,jm,tr1)
    fjn(wd,prr2,gh,gh,gh,gh,im,jm,wd,prd1,gh,gh,gh,gh,im,jm,tr2)
    fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
    fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

    return w, wd      

def computeRes_lin(sch, flinsch, res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None):
    ''' Compute the opposite of the linearized residual (not divided by the volume) '''

    if 'dnc' in sch:      
        flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
    else:
        flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

    return resd   

def computeRes_dbyvol_lin(sch, flinsch, res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None):
    ''' Compute the linearized residual (divided by the volume) '''

    if 'dnc' in sch:      
        flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
    else:
        flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

    for m in range(_np.shape(resd)[2]):
        resd[gh:-gh,gh:-gh,m] = - resd[gh:-gh,gh:-gh,m] / vol[gh:-gh,gh:-gh]

    return resd       

def computeRes_dbyvol_adj(sch, flinsch, res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None):
    ''' Compute the adjoint residual (not divided by the volume) '''

    for m in range(_np.shape(resd)[2]):
        resd[gh:-gh,gh:-gh,m] = - resd[gh:-gh,gh:-gh,m] / vol[gh:-gh,gh:-gh]

    if 'dnc' in sch:      
        flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
    else:
        flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

    return wd

def computeCoeffDiagJac(dtm1, norm, norm0m1, ninf, ninf0m1, vol, gh):

    ## relaxed on diag:
    r = _np.max([norm[:3]*norm0m1[:3], ninf[:3]*ninf0m1[:3]])
    # r = _np.max([norm*norm0m1, ninf*ninf0m1])
    # r = _np.max([norm[2]*norm0m1[2], ninf[2]*ninf0m1[2]])
    cflm1 = r*dtm1
    print("1/cfl = ", cflm1)
    print(norm)
    coefdiag = cflm1 * vol[gh:-gh,gh:-gh]

    return coefdiag

def constructJacLists(w, res, coefdiag, im, jm, gh, wbd, field, nx, ny, gam, interf1, interf2, interf3, interf4, flininflow, flinoutflow, flinnoref, flinwall, sch, flinsch, finflow, foutflow, x0, y0, xc, yc, vol, volf, cp, cv, prandtl, rgaz, cs, muref, tref, k2, k4, it=None, ite=None, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None, pinf=None):

    # construct Jacobian
    wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
    resd = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

    nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
    Jac = _np.zeros((nbentry), order='F')
    IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
    JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')

    # timeinjac = timeit.time.time()

    for m in range(5):
        for l in range(1 + 2*gh):
            for k in range(1 + 2*gh):
                wd *= 0.
                f_misc.testvector(wd,m,l,k,gh,im,jm)

                w, wd = updatefromBC_lin(w, wd, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, flininflow, flinoutflow, flinnoref, flinwall, finflow, foutflow, pinf=pinf)

                resd = computeRes_lin(sch, flinsch, res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar, eps4ar, divu2ar, vort2ar)

                f_misc.computejacobianfromjv_relaxed(Jac,IA,JA,resd,m,l,k,gh,coefdiag)
                # if it == ite:
                #     f_misc.computejacobianfromjv(Jac,IA,JA,resd,m,l,k,gh,im,jm)
                # else:
                #     f_misc.computejacobianfromjv_relaxed(Jac,IA,JA,resd,m,l,k,gh,coefdiag)
                # f_misc.computejacobianfromjv_relaxed_dbyvol(Jac,IA,JA,resd,m,l,k,gh,coefdiag,vol)

    ## Remove the zero stored
    # timeinremove = timeit.time.time()
    # timeconstructjac += (timeinremove - timeinjac)/ite
    mini = 2.e-16
    IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
    nbentry = _np.shape(Jac)[0]
    # print(nbentry)
    # timeoutremove = timeit.time.time()
    # timeremove += (timeoutremove - timeinremove)/ite

    return IA, JA, Jac, nbentry

def constructJacobian(w, res, im, jm, gh, wbd, field, nx, ny, gam, interf1, interf2, interf3, interf4, sch, flinsch, x0, y0, xc, yc, vol, volf, cp, cv, prandtl, rgaz, cs, muref, tref, k2, k4, flininflow=None, flinoutflow=None, flinnoref=None, flinwall=None, finflow=None, foutflow=None, it=None, ite=None, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None, pinf=None, comm=MPI.COMM_WORLD):

    # construct Jacobian
    wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
    resd = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

    nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
    Jac = _np.zeros((nbentry), order='F')
    IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
    JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')

    # timeinjac = timeit.time.time()

    for m in range(5):
        for l in range(1 + 2*gh):
            for k in range(1 + 2*gh):
                wd *= 0.
                f_misc.testvector(wd,m,l,k,gh,im,jm)

                w, wd = updatefromBC_lin(w, wd, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, flininflow=flininflow, flinoutflow=flinoutflow, flinnoref=flinnoref, flinwall=flinwall, finflow=finflow, foutflow=foutflow, pinf=pinf)

                resd = computeRes_lin(sch, flinsch, res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar, eps4ar, divu2ar, vort2ar)

                f_misc.computejacobianfromjv_dbyvol(Jac,IA,JA,resd,m,l,k,gh,im,jm,vol)

    mini = 2.e-16
    IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
    nbentry = _np.shape(Jac)[0]
    # print(nbentry)

    Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2,comm=comm)
    # Jacs = pet.createMatPetscCSRNoMpi(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2,comm=comm)

    return Jacs  

def constructJacobian_cyl(w, res, im, jm, gh, wbd, field, nx, ny, gam, interf1, interf2, interf3, interf4, sch, flinsch, x0, y0, xc, yc, vol, volf, cp, cv, prandtl, rgaz, cs, muref, tref, k2, k4, flininflow=None, flinoutflow=None, flinnoref=None, flinwall=None, finflow=None, foutflow=None, fjn=None, prr1=None, prr2=None, prd1=None, prd2=None, tr1=None, tr2=None, it=None, ite=None, eps2ar=None, eps4ar=None, divu2ar=None, vort2ar=None):

    # construct Jacobian
    wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
    resd = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

    nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
    Jac = _np.zeros((nbentry), order='F')
    IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
    JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')

    # timeinjac = timeit.time.time()

    for m in range(5):
        for l in range(1 + 2*gh):
            for k in range(1 + 2*gh):
                wd *= 0.
                f_misc.testvector(wd,m,l,k,gh,im,jm)

                w, wd = updatefromBC_lin_cyl(w, wd, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, flinnoref=flinnoref, flinwall=flinwall, fjn=fjn, prr1=prr1, prr2=prr2, prd1=prd1, prd2=prd2, tr1=tr1, tr2=tr2)

                resd = computeRes_lin(sch, flinsch, res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar, eps4ar, divu2ar, vort2ar)

                f_misc.computejacobianfromjv_withjn_dbyvol(Jac,IA,JA,resd,m,l,k,gh,im,jm,vol)

    mini = 2.e-16
    IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
    nbentry = _np.shape(Jac)[0]
    print(nbentry)

    Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

    return Jacs      



######################### Private functions ####################

# solve monoblock Boundary Layer

def bl2d(dgeom = dict(), dphys = dict(), dnum = dict(), compmode = 'direct', lf = list(), lflin = list(), out_dir = 'totodir', isresol= False):
    '''
    exemple of monoblock use of 2DTOY
    to simulate 2D laminar Boundary Layer flow
    '''
    os.system('mkdir -p %s' % out_dir)

    finflow, foutflow, fnoref, fwall, fsch                          = read_input_dict_functions(lf)
    im, jm, gh, L, high, xini                                       = read_input_dict_geom(dgeom, dnum)
    ite, cfl, k2, k4, sch, order, freqres, freqsort                 = read_input_dict_num(dnum)
    gam, cs, tref, muref, rgaz, prandtl, mach, tinf, Lref, StateRef = read_input_dict_phys(dphys)

    interf1, interf2, interf3, interf4 = define_interf(im, jm, gh)

    muinf, sound, uinf, einf, pinf, rhoinf, runit, cp, cv, dphys = compute_infq(muref, tref, cs, tinf, gam, rgaz, mach, StateRef, dphys)

    if compmode == 'direct':
        rkcoefs = dnum['rkcoefs']
    elif compmode == 'fixed_point':
        lasolver = dnum['lasolver']
        if lasolver == 'gmres':
            tol = dnum['tol']

    x, y = construct_mesh(xini, L, high, im, jm)        

    x0, y0, xc, yc, nx, ny, vol, volf, w, res = initialize_fields(im, jm, gh)

    dt = compute_Timestep(L, im, Lref, y, mach, cfl)
    dtm1 = 1./dt

    print("============setup===============")
    print('scheme             = ', sch)
    print('order(or nb pts)   = ', order)
    print('dt                 = ', dt)
    print('StateRef           = ', StateRef)
    print("gam                = ", dphys['gam'])
    print("Ts                 = ", dphys['Ts'])
    print("cs                 = ", dphys['cs'])
    print("musuth             = ", dphys['musuth'])
    print("rgaz               = ", dphys['rgaz'])
    print("Prandtl            = ", dphys['Prandtl'])
    print("Mach               = ", dphys['Mach'])
    print("T0                 = ", dphys['T0'])
    print("Runit              = ", dphys['Runit'])
    print("Lref               = ", dphys['Lref'])
    print("============setup===============")
    print(" ")
    print(" ")

    print('mu_inf = ', muinf)
    print('u_inf = ', uinf)
    print('nu_inf = ', muinf/rhoinf)   

    Lref, uinf, tinf, rhoinf, pinf, cp, cv, rgaz, einf, tref, muref, cs, muinf = computeAdimRVT(rhoinf, uinf, tinf, muinf, pinf, cp, cv, rgaz, einf, tref, muref, cs)

    x0, y0, nx, ny, xc, yc, vol, volf = compute_geometryAdim(im, jm, gh, x, y, Lref, x0,y0,nx,ny,xc,yc,vol,volf)

    # StateAdim

    print('======StateAdim=========')
    print(' ')
    print('state_adim uinf    = ', uinf)
    print('state_adim tinf    = ', tinf)
    print('state_adim rhoinf  = ', rhoinf)
    print('state_adim sound   = ', sound)
    print('state_adim pinf    = ', pinf)
    print('state_adim cp      = ', cp)
    print('state_adim cv      = ', cv)
    print('state_adim rgaz    = ', rgaz)
    print('state_adim einf    = ', einf)
    print('state_adim tref    = ', tref)
    print('state_adim muref   = ', muref)
    print('state_adim cs      = ', cs)
    print('state_adim muinf   = ', muinf)
    print('state_adim runit   = ', runit)
    print('======StateAdim=========')
    print(' ')

    w, wbd, field = initialisationSelfSim(im, jm, gh, w, x0, y0, Lref, mach, dphys, rhoinf, uinf, einf, xc, yc, cv)

    # ######## Restart from a previous solution with exactly the same mesh

    # import restart_init as ri
    # filet = './Wksp/dnc_5/state_atcenter_ite10.dat'
    # Xin, Yin, roin, rouin, rovin, rowin, roein = ri.read_init(filet)

    # w[gh:-gh, gh:-gh, 0]     = roin
    # w[gh:-gh, gh:-gh, 1]     = rouin
    # w[gh:-gh, gh:-gh, 2]     = rovin
    # w[gh:-gh, gh:-gh, 3]     = rowin
    # w[gh:-gh, gh:-gh, 4]     = roein

    ## OR Restart from a previous solution with a different mesh (approximated interp. at 1st order, only for cartesian rectangular grids)
    # import interpgrid
    # w[gh:-gh, gh:-gh, 0] = interpgrid.interpgrid(Xin, Yin, roin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 1] = interpgrid.interpgrid(Xin, Yin, rouin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 2] = interpgrid.interpgrid(Xin, Yin, rovin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 3] = interpgrid.interpgrid(Xin, Yin, rowin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 4] = interpgrid.interpgrid(Xin, Yin, roein, xc[gh:-gh,:], yc[:,gh:-gh])


    filename = out_dir + '/initialisation_gh.dat'
    __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

    # Time Marching Loop
    if compmode == 'direct':
        freq = freqsort
        time = 0.
        wreal = w*1.
        denom = im*jm*freq*len(rkcoefs)
        timein0 = timeit.time.time()
        for it in range(1,ite+1):
            for rk in rkcoefs:
                # Boundary on state vector

                w = updatefromBC(w, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, finflow, foutflow, fnoref, fwall)

                # Compute spatial discretization

                res = computeRes(sch, fsch, res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm)
                
                # advance rk
                w[gh:-gh,gh:-gh,0] = wreal[gh:-gh,gh:-gh,0] + rk * dt * res[gh:-gh,gh:-gh,0] / vol[gh:-gh,gh:-gh]
                w[gh:-gh,gh:-gh,1] = wreal[gh:-gh,gh:-gh,1] + rk * dt * res[gh:-gh,gh:-gh,1] / vol[gh:-gh,gh:-gh]
                w[gh:-gh,gh:-gh,2] = wreal[gh:-gh,gh:-gh,2] + rk * dt * res[gh:-gh,gh:-gh,2] / vol[gh:-gh,gh:-gh]
                w[gh:-gh,gh:-gh,3] = wreal[gh:-gh,gh:-gh,3] + rk * dt * res[gh:-gh,gh:-gh,3] / vol[gh:-gh,gh:-gh]
                w[gh:-gh,gh:-gh,4] = wreal[gh:-gh,gh:-gh,4] + rk * dt * res[gh:-gh,gh:-gh,4] / vol[gh:-gh,gh:-gh]
            #Finalize time step
            if it == 1:
                norm0, nmoy0 = f_norm.compute_norml2(res ,im, jm, gh)
                for lala in range(5):
                    if (norm0[lala] <=3.e-16): norm0[lala] = 1.
            time += dt
            wreal = w * 1.
            if it%freqres == 0:
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm/norm0))
            # if it%freq == 0:
            #     timetosort = timeit.time.time()
            #     print 'Time in function = ', (timetosort- timein0) / denom
            #     norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
            #     print 'ite = %i , norm2(res) = %s' % (it, norm/norm0)
            #     # print 'write file'
            #     # usefull for plotting result
            #     fwall(w,'Jlo', gam, interf3, gh, im, jm)
            #     filename = out_dir + '/state_at_ite%i.dat' % it
            #     __writestate_node(filename, im, jm, w, x0, y0, gh)
            #     filename = out_dir + '/state_atcenter_ite%i.dat' % it
            #     __writestate_center(filename, im, jm, w, xc, yc, gh)
            #     timein0 = timeit.time.time()
            if it%freqsort == 0:
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm/norm0))
                # filename = out_dir + '/state_at_ite%i.dat' % it
                # __writestate_node(filename, im, jm, w, x0, y0, gh)
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
            if it == ite:
                filename = out_dir + '/state_atcentergh_ite%i.dat' % it
                __writestate_center_gh(filename, im, jm, w, xc, yc)

    elif compmode == 'impli':
        fimpli   = lf[-1]
        time = 0.
        dtcoef = 1.
        # wreal = w*1.
        dw = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
        
        w = updatefromBC(w, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, finflow, foutflow, fnoref, fwall)

        filename = out_dir + '/initialisation_gh.dat'
        __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

        lmax = 10  #4
        for it in range(1,ite+1):
            # if it > 4000: cfl = min(0.5 + (it-4000)*0.001 ,1.)
            # if it > 8000: cfl = min(3. + (it-8000)*0.001 ,10.)
            # Compute spatial discretization
            
            res = computeRes(sch, fsch, res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm)

            # implicit (MF)
            fimpli(dw,nx,ny,w,res,vol,volf,dtcoef,cfl,gam,rgaz,prandtl,lmax,gh,cv,cs,muref,tref,cs,im,jm)

            # advance BDF1
            w[gh:-gh,gh:-gh,0] += dw[gh:-gh,gh:-gh,0]
            w[gh:-gh,gh:-gh,1] += dw[gh:-gh,gh:-gh,1]
            w[gh:-gh,gh:-gh,2] += dw[gh:-gh,gh:-gh,2]
            w[gh:-gh,gh:-gh,3] += dw[gh:-gh,gh:-gh,3]
            w[gh:-gh,gh:-gh,4] += dw[gh:-gh,gh:-gh,4]

            # Boundary on state vector

            w = updatefromBC(w, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, finflow, foutflow, fnoref, fwall)

            #Finalize time step
            if it == 1:
                norm0, nmoy0 = f_norm.compute_norml2(res ,im, jm, gh)
                for lala in range(5):
                    if (norm0[lala] <=3.e-16): norm0[lala] = 1.
                impl, impl0 = f_norm.compute_norml2(dw ,im, jm, gh)
                print('ite = %i , norm2(imp) = %s' % (it, impl))

            if it%freqres == 0:
                print('cfl = ', cfl)
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm/norm0))
            if it%freqsort == 0:
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm/norm0))
            #     filename = out_dir + '/state_at_ite%i.dat' % it
            #     __writestate_node(filename, im, jm, w, x0, y0, gh)
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
            #     timein0 = timeit.time.time()
            if it == ite:
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
                filename = out_dir + '/state_atcentergh_ite%i.dat' % it
                __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)


    elif compmode == 'fixed_point':
        
        flininflow, flinoutflow, flinnoref, flinwall, flinsch = read_lflin(lflin)

        # timeconstructjac = 0.
        # timeremove  = 0.
        # timecoefdiag = 0.
        # timejaccsc   = 0.
        # timejacinv   = 0.
        cfl = 1.e4 #1.e1*4*10*10   # *1.e12    #1.e4 for hllc 
        dt  = compute_Timestep_vadim(gh, yc, mach, cfl)
        dtm1 = 1./dt
        for it in range(1,ite+1):

            # Boundary on state vector
            
            w = updatefromBC(w, wbd, field, nx, ny, gam, im, jm, gh, interf1, interf2, interf3, interf4, finflow, foutflow, fnoref, fwall)

            filename = out_dir + '/initialisation_gh.dat'
            __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

            # Compute spatial discretization
            eps2ar = _np.zeros((im+2*gh, jm+2*gh), order='F')
            eps4ar = _np.zeros((im+2*gh, jm+2*gh), order='F')
            divu2ar = _np.zeros((im+2*gh, jm+2*gh), order='F')
            vort2ar = _np.zeros((im+2*gh, jm+2*gh), order='F')
            
            res = computeRes(sch, fsch, res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, k2, k4, im, jm, eps2ar, eps4ar, divu2ar, vort2ar)

            # Newton residual:
            #norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
            norm, ninf = f_norm.compute_norml2inf(res ,im, jm, gh)
            if it == 1:
                for i in range(5):
                    norm[i] = max(norm[i],1.e-15)
                    ninf[i] = max(ninf[i],1.e-15)
                norm0m1 = 1./norm
                ninf0m1 = 1./ninf


            # construct Jacobian
            print('iter = ', it)
            coefdiag = computeCoeffDiagJac(dtm1, norm, norm0m1, ninf, ninf0m1, vol, gh)
            IA, JA, Jac, nbentry = constructJacLists(w, res, coefdiag, im, jm, gh, wbd, field, nx, ny, gam, interf1, interf2, interf3, interf4, flininflow, flinoutflow, flinnoref, flinwall, sch, flinsch, finflow, foutflow, x0, y0, xc, yc, vol, volf, cp, cv, prandtl, rgaz, cs, muref, tref, k2, k4, it, ite, eps2ar, eps4ar, divu2ar, vort2ar)
            
            # import scipy.sparse as sp
            # Jacs = sp.csr_matrix((Jac, (IA, JA)), shape=(im*jm*5, im*jm*5))
            Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

            # timeoutjaccsc = timeit.time.time()
            # timejaccsc += (timeoutjaccsc - timeoutremove)/ite

            # plt.figure()
            # plt.spy(Jacs)
            # plt.show()
            
            ksp  = pet.kspLUPetsc(Jacs)
            ksp, dwtmp = pet.iterNewton(_np.ravel(res[gh:-gh,gh:-gh,:]), Jacs, ksp)

            # Print RAM
            import psutil
            if it==1 or it==2:
                print(psutil.virtual_memory())
            # For parallel case, need to add MPI barrier before reshape
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            comm.Barrier()

            dwtmp = _np.real(dwtmp)
            dw = _np.reshape(dwtmp, (im,jm,5))

            # timeoutjacinv = timeit.time.time()
            # timejacinv += (timeoutjacinv - timeoutjaccsc)/ite
            
            # w[gh:-gh,gh:-gh,:] += dw
            ## Add a check before adding NaN
            if not _np.isnan(_np.sum(dw)):   
                w[gh:-gh,gh:-gh,:] += dw
                dw_old = dw
            else:
                w[gh:-gh,gh:-gh,:] -= dw_old
                cfl = cfl / 2  

            # if it%freqsort == 0:
            #     filename = out_dir + '/state_ite%i.dat' % it
            #     __writestate_node(filename, im, jm, w, x0, y0, gh)
            # if it%freqres == 0:
            #     nn = norm*norm0m1
            #     # print nn
            #     filename = out_dir + '/residual.dat'
            #     fout = open(filename , 'a')
            #     fout.write(str(it) + ' ' )
            #     for i in range(5):
            #         fout.write(str(nn[i]) + ' ')
            #     fout.write('\n')
            #     fout.close()
            if it == ite:
                filename = out_dir + '/fixedpoint.dat'
                __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
            # print 'Time: ', (timeconstructjac+timeremove+timecoefdiag+timejaccsc+timejacinv)*ite

        # print 'Time to construct Jacobian', timeconstructjac
        # print 'Time to remove zeros in Jacobian', timeremove
        # print 'Time to convert into csc = ', timejaccsc
        # print 'Time to invert = ', timejacinv
        # print 'Time Baseflow = ', (timeconstructjac+timeremove+timecoefdiag+timejaccsc+timejacinv)*ite        

        ### Resolvent : compute eigenvalues and eigenvectors

        if isresol:

            dir = './'
            dir = './BASEFLOW_BL/'

            os.system('mkdir -p %s' % dir)
            equations = [1, 2, 3] #Forcing on momentum equations
            print("** Writing matrices for resolvent **")
            resol.computeandwrite_PETSc(dir, gam, mach, vol, w, im, jm, gh, nbentry, Jac, IA, JA, equations)

