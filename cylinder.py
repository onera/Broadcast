# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#!/usr/bin/env python

'''
File: cylinder.py

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
import resolvent_all  as resol
# import SIM
# import SIM.BLprofiles_implicit as blsim
import f_init
import meshCyl as mesh
import meshBL as meshBL
import handleBC as handleBC

import numpy as _np
import matplotlib.pyplot as plt

import os
import sys
import timeit

from petsc4py import PETSc


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
    
def __writemesh(filename, im, jm, x, y) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" \n')
    f_out.write('ZONE I = ' + str(im) + ',  J = ' + str(jm) + '\n')
    for j in range(jm):
        for i in range(im):
            f_out.write(str(x[i,j]) + ' ' + str(y[i,j]) + '\n')
    f_out.close()

def __writenormals_node(filename, imloc, jmloc, nx, ny, x0, y0) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "nx_i" "nx_j" "ny_i" "ny_j" \n')
    f_out.write('ZONE I = ' + str(imloc) + ',  J = ' + str(jmloc) + '\n')
    for j in range(jmloc):
        for i in range(imloc):
            nx_i = nx[i,j,0]
            nx_j = nx[i,j,1]
            ny_i = ny[i,j,0]
            ny_j = ny[i,j,1]
            f_out.write(str(x0[i,j]) + ' ' + str(y0[i,j]) + ' ' +
                        str(nx_i)    + ' ' + str(nx_j)   + ' ' +
                        str(ny_i)   + ' ' + str(ny_j)   + '\n')
    f_out.close()    

def __writenormals_nodebis(filename, im, jm, nx, ny, x0, y0, gh) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "nx_i" "nx_j" "ny_i" "ny_j" \n')
    f_out.write('ZONE I = ' + str(im) + ',  J = ' + str(jm) + '\n')
    for j in range(gh, gh+jm):
        for i in range(gh, gh+im):
            nx_i = nx[i,j,0]
            nx_j = nx[i,j,1]
            ny_i = ny[i,j,0]
            ny_j = ny[i,j,1]
            f_out.write(str(x0[i,j]) + ' ' + str(y0[i,j]) + ' ' +
                        str(nx_i)    + ' ' + str(nx_j)   + ' ' +
                        str(ny_i)   + ' ' + str(ny_j)   + '\n')
    f_out.close()

def __writevolume(filename, im, jm, vol, volf, xc, yc, gh) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "vol" "volf_i" "volf_j" \n')
    f_out.write('ZONE I = ' + str(im) + ',  J = ' + str(jm+gh) + '\n')
    for j in range(gh, gh+jm+gh):
        for i in range(gh, gh+im):
            volfi = volf[i,j,0]
            volfj = volf[i,j,1]
            f_out.write(str(xc[i,j]) + ' ' + str(yc[i,j]) + ' ' +
                        str(vol[i,j])  + ' '   + str(volfi) + ' ' + 
                        str(volfj) + '\n')
    f_out.close() 

def __writevolumegh(filename, imloc, jmloc, vol, volf, xc, yc) :
    # print 'write file'
    f_out = open(filename , 'w')
    f_out.write('TITLE="state"\n')
    f_out.write('VARIABLES= "X" "Y" "vol" "volf_i" "volf_j" \n')
    f_out.write('ZONE I = ' + str(imloc) + ',  J = ' + str(jmloc) + '\n')
    for j in range(jmloc):
        for i in range(imloc):
            volfi = volf[i,j,0]
            volfj = volf[i,j,1]
            f_out.write(str(xc[i,j]) + ' ' + str(yc[i,j]) + ' ' +
                        str(vol[i,j])  + ' '   + str(volfi) + ' ' + 
                        str(volfj) + '\n')
    f_out.close()       

def __comp_Sutherland(propref, Ts, Cs, T):
    '''Dynamical viscosity / thermal conductivity from sutherland law'''
    return propref*_np.sqrt(T/Ts)*((1.+Cs/Ts)/(1.+Cs/T))

def __compute_tot_energy_inf(R_pg, gamma, t_inf, v_inf):
    '''Total energy E = R/(gamma-1)*Tinf+(uinf**2)/2'''
    return R_pg/(gamma-1.)*t_inf+0.5*v_inf*v_inf

def remove_zero_jac(IA, JA, Jac, mini=1.e-15):
    ''' Remove the zero components from the Jac list in order not to store any zero in the sparse matrix '''
    to_keep = _np.absolute(Jac) > mini
    Jac = Jac[to_keep,...]
    IA  = IA[to_keep,...]
    JA  = JA[to_keep,...]
    return IA, JA, Jac

def centers_array(A):
    ''' Compute the values of an array A at the centers'''
    return 0.25 * ( A[:-1,:-1] + A[1:,:-1] + A[:-1,1:] + A[1:,1:] )


######################### Private functions ####################

# solve monoblock  cylinder

def cyl2d(dgeom = dict(), dphys = dict(), dnum = dict(), compmode = 'direct', lf = list(), lflin = list(), out_dir = 'totodir', isresol= False):
    '''
    exemple of monoblock use of 2DTOY
    to simulate 2D circular cylinder
    '''
    os.system('mkdir -p %s' % out_dir)

    # get functions
    routinein  = lf[0]
    routineout = lf[1]
    routinenr  = lf[2]
    routinew   = lf[3]
    routinesch = lf[4]
    routinebw  = lf[5]
    routinejn  = lf[6]
    libbnd     = lf[7]
    libsch     = lf[8]

    finflow  = eval("%s.%s"    % (libbnd, routinein))
    foutflow = eval("%s.%s"    % (libbnd, routineout))
    fnoref   = eval("%s.%s"    % (libbnd, routinenr ))
    fwall    = eval("%s.%s"    % (libbnd, routinew  ))
    fsym     = eval("%s.%s"    % (libbnd, routinebw  ))
    fjn      = eval("%s.%s"    % (libbnd, routinejn  ))
    fsch     = eval("%s.%s"    % (libsch, routinesch))

    # Create mesh
    im       = dgeom['im']
    jm       = dgeom['jm']
    rint     = dgeom['ri']
    rext     = dgeom['rf']

    ite      = dnum['ite']
    cfl      = dnum['cfl']
    k2       = dnum['k2']
    k4       = dnum['k4']
    # alpha    = dnum['alpha']
    sch      = dnum['sch']
    order    = dnum['order']
    freqres  = dnum['freqres']
    freqsort = dnum['freqsort']
    gh       = dnum['gh']

    if compmode == 'direct':
        rkcoefs = dnum['rkcoefs']
    elif compmode == 'fixed_point':
        lasolver = dnum['lasolver']
        if lasolver == 'gmres':
            tol = dnum['tol']

    ## MESH v2
    Ny_in   = jm-15  #90%       #85%       #100%       #jm-4
    deltaBL = 50.-rint   #12.-rint  #10.-rint  #rext-rint  #50.-rint
    percent = 0.016        #0.01      #0.01      #0.016      #0.016
    
    Ny_out  = jm - Ny_in 
    high    = rext - rint
    Nend    = high/deltaBL
    r_int   = meshBL.bigeom_stretch_in(Ny_in, deltaBL, percent)
    r_out   = meshBL.exp_stretch_out(Ny_out, deltaBL, percent, Nend)
    r       = _np.concatenate((r_int, r_out)) + rint
    # r = meshBL.bigeom_stretch_in(Ny_in, deltaBL, percent) + rint

    x,y = mesh.cylinder2d_bisfull(r, 0., 0., im+1, jm+1)


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

    muinf   = __comp_Sutherland(muref, tref, cs, tinf)
    sound   = _np.sqrt(gam*rgaz*tinf)
    uinf    = mach * sound
    einf    = __compute_tot_energy_inf(rgaz, gam, tinf, uinf)

    dx = _np.abs(x[0,1]-x[0,0])/Lref # adim done after muinf A.Poulain
    sound = 1./dphys['Mach']
    dt = cfl * dx / (sound+1.)
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

    #for similitude sol
    dphys['mu0'] = muinf

    if StateRef == 'm0_p0_t0':
        pinf    =  dphys['P0']
        rhoinf  = pinf/(rgaz*tinf)
        runit   = rhoinf * uinf/muinf
        dphys['Runit'] = runit

    elif StateRef == 'm0_runit_t0':
        runit   =  dphys['Runit']
        rhoinf  = runit*muinf/uinf
        pinf    = rhoinf*rgaz*tinf

    print('mu_inf = ', muinf)
    print('u_inf = ', uinf)
    print('nu_inf = ', muinf/rhoinf)

    # Reynolds ------------------------------------------------------------
    rey     = runit * rint * 2
    print('Re_d =',  rey)

    cp = gam * rgaz /(gam-1.)
    cv =       rgaz /(gam-1.)

    # StateRef

    stateref    = _np.zeros(5)
    stateref[0] = rhoinf
    stateref[1] = rhoinf * uinf
    stateref[2] = 0.
    stateref[3] = 0.
    stateref[4] = rhoinf * einf

    # Adim (by RVT = rho, velo et temperature)

    Roref  = rhoinf
    Vref   = uinf
    Tref   = tinf
    Pref   = Roref*Vref**2
    Cvref  = Vref**2/Tref

    Eref   = Vref**2
    Rgpref = Cvref

    Lref   = 1.
    Muref  = Roref*Vref*Lref

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
    tref   = tref/Tref
    muref  = muref/Muref
    cs     = cs/Tref
    muinf  = muinf/Muref

    # StateAdim

    state_adim    = _np.zeros(5)
    state_adim[0] = rhoinf
    state_adim[1] = rhoinf * uinf
    state_adim[2] = 0.
    state_adim[3] = 0.
    state_adim[4] = rhoinf * einf


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

    #interfaces definitions (may be done at the begining)

    interfl = []
    interfdic = dict()
    interfloc = dict()
    BCdict = dict()

    # Jlo Wall
    interf1      =  _np.zeros((2,2), order='F')
    interf1[0,0] = 1     # imin
    interf1[0,1] = 1     # jmin
    interf1[1,0] = im    # imax
    interf1[1,1] = 1     # jmax
    interfl.append('interf1')
    interfdic['interf1'] = interf1
    interfloc['interf1'] = 'Jlo'

    # Jhi Noref
    interf2      =  _np.zeros((2,2), order='F')
    interf2[0,0] = 1     # imin
    interf2[0,1] = jm    # jmin
    interf2[1,0] = im    # imax
    interf2[1,1] = jm    # jmax
    interfl.append('interf2') 
    interfdic['interf2'] = interf2
    interfloc['interf2'] = 'Jhi'

    # Ilo 
    interf3      =  _np.zeros((2,2), order='F')
    interf3[0,0] = 1      # imin
    interf3[0,1] = 1-gh   # jmin
    interf3[1,0] = 1      # imax
    interf3[1,1] = jm+gh  # jmax
    interfl.append('interf3')
    interfdic['interf3'] = interf3
    interfloc['interf3'] = 'Ilo'

    # Ihi 
    interf4      =  _np.zeros((2,2), order='F')
    interf4[0,0] = im     # imin
    interf4[0,1] = 1-gh   # jmin
    interf4[1,0] = im     # imax
    interf4[1,1] = jm+gh  # jmax
    interfl.append('interf4')
    interfdic['interf4'] = interf4
    interfloc['interf4'] = 'Ihi'
    
    jmin = 1-gh
    jmax = jm+gh
    # pr = [[imin,jmin],[imax,jmax]]
    prr1 = _np.array([[im+1  , jmin], [im+gh , jmax]], order='F')
    prd1 = _np.array([[im-gh+1 , jmin], [im    , jmax]], order='F')
    prr2 = _np.array([[ 1-gh , jmin], [ 0    , jmax]], order='F')
    prd2 = _np.array([[ 1    , jmin], [ gh , jmax]], order='F')
    tr1  = _np.array([1,2], order='F')
    tr2  = _np.array([1,2], order='F')

    prl  = _np.array([[prr1, prd1],[prr2, prd2]]) 
    trl  = _np.array([[tr1],[tr2]]) 

    # Compute Geometry
    for j in range(jm+1):
        for i in range(im+1):
            x0[i+gh,j+gh] = x[i,j]
            y0[i+gh,j+gh] = y[i,j]
    
    f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)

    f_bnd.jn_match_geom_2d(xc,prr1,gh,gh,gh,gh,im,jm,xc,prd2,gh,gh,gh,gh,im,jm,tr1)
    f_bnd.jn_match_geom_2d(xc,prr2,gh,gh,gh,gh,im,jm,xc,prd1,gh,gh,gh,gh,im,jm,tr2)
    
    f_bnd.jn_match_geom_2d(yc,prr1,gh,gh,gh,gh,im,jm,yc,prd2,gh,gh,gh,gh,im,jm,tr1)
    f_bnd.jn_match_geom_2d(yc,prr2,gh,gh,gh,gh,im,jm,yc,prd1,gh,gh,gh,gh,im,jm,tr2)
    
    f_bnd.jn_match_geom_2d(vol,prr1,gh,gh,gh,gh,im,jm,vol,prd2,gh,gh,gh,gh,im,jm,tr1)
    f_bnd.jn_match_geom_2d(vol,prr2,gh,gh,gh,gh,im,jm,vol,prd1,gh,gh,gh,gh,im,jm,tr2)

    ## WRONG for an ODD number of points in x-direction
    # nx[im+gh,:,1] = nx[gh,:,1]
    # ny[im+gh,:,1] = ny[gh,:,1]
    # nx[im+gh,:,0] = nx[gh,:,0]
    # ny[im+gh,:,0] = ny[gh,:,0]
    
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

    wbd   = _np.zeros((im    , 5), order = 'F') # dummy zeros in place of top domain state vector
    field = _np.zeros((im, gh, 5), order = 'F')

    # Initialization
    wbd[:, 0]      = state_adim[0]
    wbd[:, 1]      = state_adim[1]
    wbd[:, 2]      = state_adim[2]
    wbd[:, 3]      = state_adim[3]
    wbd[:, 4]      = state_adim[4]

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


    ######## Restart from a previous solution with exactly the same mesh

    # import restart_init as ri
    # filet = './Wksp/Cylinder/dnc_5/state_atcenter_ite10.dat'
    # imold = 630   #630
    # jmold = 300   #300
    # # Xin, Yin, roin, rouin, rovin, rowin, roein = ri.read_init(filet)
    # BLprof = _np.loadtxt(filet,comments=('#','ZONE'),skiprows=3)
    # Xin   = _np.reshape(BLprof[:,0],(imold,jmold), order='F')
    # Yin   = _np.reshape(BLprof[:,1],(imold,jmold), order='F')
    # roin  = _np.reshape(BLprof[:,2],(imold,jmold), order='F')
    # rouin = _np.reshape(BLprof[:,3],(imold,jmold), order='F')
    # rovin = _np.reshape(BLprof[:,4],(imold,jmold), order='F')
    # rowin = _np.reshape(BLprof[:,5],(imold,jmold), order='F')
    # roein = _np.reshape(BLprof[:,6],(imold,jmold), order='F')

    # import scipy.interpolate as sci
    # finterp = sci.NearestNDInterpolator(_np.concatenate((_np.reshape(_np.ravel(Xin), (imold*jmold,1)), _np.reshape(_np.ravel(Yin), (imold*jmold,1))), axis=1), _np.ravel(roin))
    # roint = finterp(xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])
    # finterp = sci.NearestNDInterpolator(_np.concatenate((_np.reshape(_np.ravel(Xin), (imold*jmold,1)), _np.reshape(_np.ravel(Yin), (imold*jmold,1))), axis=1), _np.ravel(rouin))
    # rouint = finterp(xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])
    # finterp = sci.NearestNDInterpolator(_np.concatenate((_np.reshape(_np.ravel(Xin), (imold*jmold,1)), _np.reshape(_np.ravel(Yin), (imold*jmold,1))), axis=1), _np.ravel(rovin))
    # rovint = finterp(xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])
    # finterp = sci.NearestNDInterpolator(_np.concatenate((_np.reshape(_np.ravel(Xin), (imold*jmold,1)), _np.reshape(_np.ravel(Yin), (imold*jmold,1))), axis=1), _np.ravel(rowin))
    # rowint = finterp(xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])
    # finterp = sci.NearestNDInterpolator(_np.concatenate((_np.reshape(_np.ravel(Xin), (imold*jmold,1)), _np.reshape(_np.ravel(Yin), (imold*jmold,1))), axis=1), _np.ravel(roein))
    # roeint = finterp(xc[gh:-gh,gh:-gh], yc[gh:-gh,gh:-gh])

    # # w[gh:-gh, gh:-gh, 0] = _np.reshape(roint, (im,jm), order='C')
    # w[gh:-gh, gh:-gh, 1] = _np.reshape(rouint, (im,jm), order='C')
    # w[gh:-gh, gh:-gh, 2] = _np.reshape(rovint, (im,jm), order='C')
    # w[gh:-gh, gh:-gh, 3] = _np.reshape(rowint, (im,jm), order='C')
    # # w[gh:-gh, gh:-gh, 4] = _np.reshape(roeint, (im,jm), order='C')

    wlist = [w]
    BCtypel = [lf, lflin]
    wbddic = dict()
    fielddic = dict()
    twalldic = dict()
    prldic = dict()
    trldic = dict()

    BCdict['interf1'] = 'fwall'
    BCdict['interf2'] = 'fnoref'
    BCdict['interf3'] = 'fjn'
    BCdict['interf4'] = 'fjn'

    wbddic['interf2'] = wbd

    prldic['interf3'] = prl
    trldic['interf3'] = trl
    prldic['interf4'] = prl
    trldic['interf4'] = trl

    checkint = handleBC.checkinterf(interfl, interfdic, interfloc, im, jm, gh)
    if checkint != 'GoodJob':
        print(checkint)
        exit()
    interfl = handleBC.sortBC(interfl, interfdic, interfloc, im, jm)
    # print interfl

    filename = out_dir + '/initialisation.dat'
    # __writestate_center(filename, im, jm, w, xc, yc, gh)
    __writestate_node(filename, im+1, jm+1, w, x0, y0, gh)
    # filename = out_dir + '/mesh9_atnode.dat'
    # __writestate_node(filename, im+1, jm+1, w, x0, y0, gh)

    # Time Marching Loop
    if compmode == 'direct':
        freq = freqsort
        time = 0.
        wreal = w*1.
        denom = im*jm*freq*len(rkcoefs)
        timein0 = timeit.time.time()
        lm = max(im,jm)
        for it in range(1,ite+1):
            ##  print "iteration : ", it
            for rk in rkcoefs:
                # print "    rkstep = ", rk
                # Boundary on state vector

                fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                fwall(w,'Jlo', gam, interf1, gh, im, jm)
                fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
                fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

                # filename = out_dir + '/initialisation_gh.dat'
                # __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

                # Compute spatial discretization
                if sch in ['dnc', 'dnc_hll', 'dnc_hllc', 'rbc']:
                    # fwall needed for dissipation near bnd_wall
                    fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                elif 'ausm' in sch:
                    # fwall needed for dissipation near bnd_wall
                    fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, tRef,pRef,machRef, im, jm)
                elif sch == 'cdt':
                    fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, alpha, lm, im, jm)
                else:
                    fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)
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
            if it%freq == 0:
                timetosort = timeit.time.time()
                print('Time in function = ', (timetosort- timein0) / denom)
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm/norm0))
                # print 'write file'
                # usefull for plotting result
                fwall(w,'Jlo', gam, interf1, gh, im, jm)
                filename = out_dir + '/state_at_ite%i.dat' % it
                __writestate_node(filename, im, jm, w, x0, y0, gh)
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
                timein0 = timeit.time.time()
            if it == ite:
                filename = out_dir + '/state_atcentergh_ite%i.dat' % it
                __writestate_center_gh(filename, im, jm, w, xc, yc)

    elif compmode == 'impli':
        fimpli   = lf[-1]
        time = 0.
        dtcoef = 1.
        # wreal = w*1.
        dw = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
        # Boundary on state vector
        fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
        fwall(w,'Jlo', gam, interf1, gh, im, jm)
        fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
        fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

        filename = out_dir + '/initialisation_gh.dat'
        __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

        lm = max(im,jm)
        # ffiltilo(w,coef,coefbnd,im,jm,gh,cl)
        # TODO : Add time step computation!!
        lmax = 10  #8
        for it in range(1,ite+1):
            # if it > 4000: cfl = min(0.1 + (it-4000)*0.001 ,2.)
            # if it > 8000: cfl = min(3. + (it-8000)*0.001 ,10.)
            # Compute spatial discretization
            if sch in ['dnc', 'dnc_hll', 'dnc_hllc', 'rbc']:
                # fwall needed for dissipation near bnd_wall
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
            elif sch == 'cdt':
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, alpha, lm, im, jm)
            elif 'ausm' in sch:
                # fwall needed for dissipation near bnd_wall
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, tRef,pRef,machRef, im, jm)
            else:
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

            # implicit (MF)
            fimpli(dw,nx,ny,w,res,vol,volf,dtcoef,cfl,gam,rgaz,prandtl,lmax,gh,cv,cs,muref,tref,cs,im,jm)

            # advance BDF1
            w[gh:-gh,gh:-gh,0] += dw[gh:-gh,gh:-gh,0]
            w[gh:-gh,gh:-gh,1] += dw[gh:-gh,gh:-gh,1]
            w[gh:-gh,gh:-gh,2] += dw[gh:-gh,gh:-gh,2]
            w[gh:-gh,gh:-gh,3] += dw[gh:-gh,gh:-gh,3]
            w[gh:-gh,gh:-gh,4] += dw[gh:-gh,gh:-gh,4]

            # Boundary on state vector
            fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
            fwall(w,'Jlo', gam, interf1, gh, im, jm)
            fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
            fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)
            #Finalize time step
            if it == 1:
                norm0, nmoy0 = f_norm.compute_norml2(res ,im, jm, gh)
                for lala in range(5):
                    if (norm0[lala] <=3.e-16): norm0[lala] = 1.
                impl, impl0 = f_norm.compute_norml2(dw ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm0))
                print('ite = %i , norm2(imp) = %s' % (it, impl))

            if it%freqres == 0:
                print('cfl = ', cfl)
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                nn = norm/norm0
                print(" ")
                print("="*80)
                print('ite = %i , norm2(res) = %s' % (it, nn))
                print("="*80)
                # print nn
                filename = out_dir + '/residual.dat'
                fout = open(filename , 'a')
                fout.write(str(it) + ' ' )
                for i in range(5):
                    fout.write(str(nn[i]) + ' ')
                fout.write('\n')
                fout.close()
                maxnorm,imaxn,jmaxn = f_norm.compute_maxnorm(res ,im, jm, gh)
                print(" ")
                print("="*80)
                print("maxnorm = ", maxnorm)
                print("imax    = ", imaxn)
                print("jmax    = ", jmaxn)
                print("="*80)
                print(" ")

            if it%freqsort == 0:
                norm, nmoy = f_norm.compute_norml2(res ,im, jm, gh)
                print('ite = %i , norm2(res) = %s' % (it, norm/norm0))
                filename = out_dir + '/state_at_ite%i.dat' % it
                __writestate_node(filename, im+1, jm+1, w, x0, y0, gh)
                filename = out_dir + '/state_atcentergh_ite%i.dat' % it
                __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)
                # import sys
                # sys.exit(2)
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
                filename = out_dir + '/residual_at_ite%i.dat' % it
                # __writestate_node(filename, im+1, jm+1, res, x0, y0, gh)
                __writestate_center(filename, im, jm, res, xc, yc, gh)
                filename = out_dir + '/residualdbyvol_at_ite%i.dat' % it
                resdbyvol = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
                resdbyvol[:,:,0] = res[:,:,0]/vol
                resdbyvol[:,:,1] = res[:,:,1]/vol
                resdbyvol[:,:,2] = res[:,:,2]/vol
                resdbyvol[:,:,3] = res[:,:,3]/vol
                resdbyvol[:,:,4] = res[:,:,4]/vol
                # __writestate_node(filename, im+1, jm+1, resdbyvol, x0, y0, gh)
                __writestate_center(filename, im, jm, resdbyvol, xc, yc, gh)
                timein0 = timeit.time.time()
            if it == ite:
                filename = out_dir + '/state_atcentergh_ite%i.dat' % it
                __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)


    elif compmode == 'fixed_point':
        # get functions
        routinein  = lflin[0]
        routineout = lflin[1]
        routinenr  = lflin[2]
        routinew   = lflin[3]
        routinesch = lflin[4]
        routinebw  = lflin[5]
        routinejn  = lflin[6]
        libbnd     = lflin[7]
        libsch     = lflin[8]

        flinoutflow = eval("%s.%s"    % (libbnd, routineout))
        flininflow  = eval("%s.%s"    % (libbnd, routinein))
        flinnoref   = eval("%s.%s"    % (libbnd, routinenr ))
        flinwall    = eval("%s.%s"    % (libbnd, routinew  ))
        flinsym     = eval("%s.%s"    % (libbnd, routinebw ))
        flinjn      = eval("%s.%s"    % (libbnd, routinejn ))
        flinsch     = eval("%s.%s"    % (libsch, routinesch))

        timeconstructjac = 0.
        timeremove   = 0.
        timecoefdiag = 0.
        timejaccsc   = 0.
        timejacinv   = 0.

        wd   = _np.zeros((im+2*gh, jm+2*gh, 5), order='F')
        wlist.append(wd)

        cfl = 1.e21
        
        # for standard output:
        # dt *= 1./cfl

        lm = max(im,jm)

        for it in range(1,ite+1):
            dt  = cfl * _np.abs(xc[gh,gh] - xc[gh,gh+1]) / (1./mach + 1.)
            dtm1 = 1./dt

            # For AD
            wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
            resd = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

            # Boundary on state vector
            # fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
            # fwall(w,'Jlo', gam, interf1, gh, im, jm)
            # fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
            # fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

            modeBC = 0  #apply BC on w
            wlist[0] = w
            wlist = handleBC.applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic, fielddic, twalldic, prldic, trldic)
            w = wlist[0]

            filename = out_dir + '/initialisation_gh.dat'
            __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)
            # __writestate_node(filename, im+1, jm+1, w, x0, y0, gh)


            # Compute spatial discretization
            if sch in 'dnc':
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
            else:
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)
            # Newton residual:
            norm, ninf = f_norm.compute_norml2inf(res ,im, jm, gh)
            if it == 1:
                for i in range(5):
                    norm[i] = max(norm[i],1.e-15)
                    ninf[i] = max(ninf[i],1.e-15)
                norm0m1 = 1./norm
                ninf0m1 = 1./ninf

            # construct Jacobian

            nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
            Jac = _np.zeros((nbentry), order='F')
            IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
            JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')

            timeinjac = timeit.time.time()
            ## relaxed on diag:
            r = _np.max([norm[:3]*norm0m1[:3], ninf[:3]*ninf0m1[:3]])
            # r = _np.max([norm*norm0m1, ninf*ninf0m1])
            cflm1 = r*dtm1
            print('iter = ', it)
            print("1/cfl = ", cflm1)
            print(norm)
            coefdiag = cflm1 * vol[gh:-gh,gh:-gh]
                        
            Nzones = _np.array([])    
            if im%(2*gh+1) != 0:
                Nzones = _np.array([[0,0],[im/2/(2*gh+1)*(2*gh+1)-1,jm-1],[im/2/(2*gh+1)*(2*gh+1),0],[im-1,jm-1]])
            # Nzones = _np.array([[0,0],[im/2/(2*gh+1)*(2*gh+1)-1,jm-1],[im/2/(2*gh+1)*(2*gh+1),0],[im-1,jm-1]])    
            print(Nzones)
            if _np.shape(Nzones)[0] > 0:
                nbentry = (im*jm * (2*gh+1)*(2*gh+1) * 5*5)*_np.shape(Nzones)[0]/2
                Jac = _np.zeros((nbentry), order='F')
                IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
                JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
                for n in range(_np.shape(Nzones)[0]/2):
                    mini = 2.e-16
                    istart = Nzones[2*n,0]
                    iend = Nzones[2*n+1,0]
                    jstart = Nzones[2*n,1]
                    jend = Nzones[2*n+1,1]         
                    for m in range(5):
                        for l in range(1 + 2*gh):
                            for k in range(1 + 2*gh):
                                wd *= 0.
                                f_misc.testvector_partial(wd,m,l,k,gh,im,jm,istart,iend,jstart,jend)

                                # flinnoref(w,wd,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                                # flinwall(w,wd,'Jlo', gam, interf1, gh, im, jm)
                                # wdr = _np.copy(wd, order='F')
                                # flinjn(w, wdr, prr1, gh, gh, gh, gh, im, jm, w, wd, prd2, gh, gh, gh, gh, im, jm, tr1)
                                # wd[-gh:,:,:] = wdr[-gh:,:,:]
                                # wdr = _np.copy(wd, order='F')
                                # flinjn(w, wdr, prr2, gh, gh, gh, gh, im, jm, w, wd, prd1, gh, gh, gh, gh, im, jm, tr2)
                                # wd[:gh,:,:] = wdr[:gh,:,:]

                                modeBC = 1  #apply BC on wd
                                wlist[0] = w
                                wlist[1] = wd
                                wlist = handleBC.applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic, fielddic, twalldic, prldic, trldic)
                                w  = wlist[0]
                                wd = wlist[1]

                                if sch in 'dnc':
                                    flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                                else:
                                    flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

                                f_misc.computejacobianfromjv_relaxed_withjnandcheck(Jac,IA,JA,resd,m,l,k,gh,coefdiag,mini,n)
            else:
                for m in range(5):
                    for l in range(1 + 2*gh):
                        for k in range(1 + 2*gh):
                            wd *= 0.
                            f_misc.testvector(wd,m,l,k,gh,im,jm)
                            # f_misc.testvector_partial(wd,m,l,k,gh,im,jm,0,im-1,0,jm-1)

                            # w[   :gh,   :  ,:] = 0.
                            # w[-gh:  ,   :  ,:] = 0.
                            # w[   :  ,   :gh,:] = 0.
                            # w[   :  ,-gh:  ,:] = 0.

                            # flinnoref(w,wd,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                            # # fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                            # flinwall(w,wd,'Jlo', gam, interf1, gh, im, jm)
                            # # fwall(w,'Jlo', gam, interf1, gh, im, jm)
                            # wdr = _np.copy(wd, order='F')
                            # flinjn(w, wdr, prr1, gh, gh, gh, gh, im, jm, w, wd, prd2, gh, gh, gh, gh, im, jm, tr1)
                            # wd[-gh:,:,:] = wdr[-gh:,:,:]
                            # wdr = _np.copy(wd, order='F')
                            # flinjn(w, wdr, prr2, gh, gh, gh, gh, im, jm, w, wd, prd1, gh, gh, gh, gh, im, jm, tr2)
                            # wd[:gh,:,:] = wdr[:gh,:,:]
                            # fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
                            # fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

                            modeBC = 1  #apply BC on wd
                            wlist[0] = w
                            wlist[1] = wd
                            wlist = handleBC.applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic, fielddic, twalldic, prldic, trldic)
                            w  = wlist[0]
                            wd = wlist[1]

                            if sch in 'dnc':
                                flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                            else:
                                flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

                            f_misc.computejacobianfromjv_relaxed_withjn(Jac,IA,JA,resd,m,l,k,gh,coefdiag) 
                            # f_misc.computejacobianfromjv_relaxed_withjnandcheck(Jac,IA,JA,resd,m,l,k,gh,coefdiag,2.e-16,0)                

            ## Remove the zero stored
            timeinremove = timeit.time.time()
            timeconstructjac += (timeinremove - timeinjac)/ite
            mini = 2.e-16
            IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
            # print nbentry
            nbentry = _np.shape(Jac)[0]
            print(nbentry)
            timeoutremove = timeit.time.time()
            timeremove += (timeoutremove - timeinremove)/ite

            Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

            timeoutjaccsc = timeit.time.time()
            timejaccsc += (timeoutjaccsc - timeoutremove)/ite

            # plt.figure()
            # import scipy.sparse as sp
            # Jacs = sp.csc_matrix((Jac,(IA,JA)), shape=(im*jm*5,im*jm*5))
            # plt.spy(Jacs)
            # plt.show()

            # Iterative Linear Algebra Loop to solve Newton: Solve resd  = -res ==> J w_sol = - res
            if lasolver == 'gmres':
                nl2res = 1.
                print(" GMRES Not Yet implemented ")
                sys.exit(2)

            else:

                import psutil
                ksp  = pet.kspLUPetsc(Jacs)
                ksp, dwtmp = pet.iterNewton(_np.ravel(res[gh:-gh,gh:-gh,:]), Jacs, ksp)
                if it==1 or it==2:
                    print(psutil.virtual_memory())
                from mpi4py import MPI
                comm = MPI.COMM_WORLD
                comm.Barrier()
                dwtmp = _np.real(dwtmp)

                timeoutjacinv = timeit.time.time()
                timejacinv += (timeoutjacinv - timeoutjaccsc)/ite
                dw = _np.reshape(dwtmp, (im,jm,5))

            # w[gh:-gh,gh:-gh,:] += dw
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
            #     print 'iter = ', it
            #     print " "
            #     print "="*80
            #     print "norml2 = ", nn
            #     print "="*80
            #     print " "
            #     filename = out_dir + '/residual.dat'
            #     fout = open(filename , 'a')
            #     fout.write(str(it) + ' ' )
            #     for i in range(5):
            #         fout.write(str(nn[i]) + ' ')
            #     fout.write('\n')
            #     fout.close()
            #     maxnorm,imaxn,jmaxn = f_norm.compute_maxnorm(res ,im, jm, gh)
            #     print " "
            #     print "="*80
            #     print "maxnorm = ", maxnorm
            #     print "imax    = ", imaxn
            #     print "jmax    = ", jmaxn
            #     print "="*80
            #     print " "
                #
            # if it%freqsort == 0:
            #     filename = out_dir + '/state_at_ite%i.dat' % it
            #     __writestate_node(filename, im, jm, w, x0, y0, gh)
            if it == ite:
                filename = out_dir + '/fixedpoint.dat'
                __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)
                filename = out_dir + '/state_atcenter_ite%i.dat' % it
                __writestate_center(filename, im, jm, w, xc, yc, gh)
            print('Time: ', (timeconstructjac+timeremove+timecoefdiag+timejaccsc+timejacinv)*ite)

        # print 'Time to construct Jacobian', timeconstructjac
        # print 'Time to remove zeros in Jacobian', timeremove
        # print 'Time to convert into csc = ', timejaccsc
        # print 'Time to invert = ', timejacinv
        # print 'Time Baseflow = ', (timeconstructjac+timeremove+timecoefdiag+timejaccsc+timejacinv)*ite

        ### Resolvent : compute eigenvalues and eigenvectors

        # TODO : make a function for Resolvent
        if isresol:

            dir = './'
            dir = './BASEFLOW_CYL/'

            os.system('mkdir -p %s' % dir)
            equations = [1, 2, 3] #Forcing on momentum equations
            print("** Writing matrices for resolvent **")
            resol.computeandwrite_PETSc(dir, gam, mach, vol, w, im, jm, gh, nbentry, Jac, IA, JA, equations)


    elif compmode == 'adj':

        flininflow  = lflin[0]
        flinoutflow = lflin[1]
        flinnoref   = lflin[2]
        flinwall    = lflin[3]
        flinsch     = lflin[4]
        flinmirror  = lflin[5]
        flinjn      = lflin[6]
        wb = _np.zeros((im+2*gh, jm+2*gh, 5), order='F')
        resb = _np.zeros((im+2*gh, jm+2*gh, 5), order='F')
        # Boundary on state vector
        fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
        fwall(w,'Jlo', gam, interf1, gh, im, jm)
        fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
        fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

        # Compute spatial discretization
        if sch in ['dnc', 'dnc_hll', 'dnc_hllc', 'rbc']:
            fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
        else:
            fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)
        nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
        Jac = _np.zeros((nbentry), order='F')
        IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
        JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
        for m in range(5):
            for l in range(1 + 2*gh):
                for k in range(1 + 2*gh):
                    resb *= 0.
                    f_misc.testvector(resb,m,l,k,gh,im,jm)

                    # w[   :gh,   :  ,:] = 0.
                    # w[-gh:  ,   :  ,:] = 0.
                    # w[   :  ,   :gh,:] = 0.
                    # w[   :  ,-gh:  ,:] = 0.

                    # flinnoref(w,wb,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                    # fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                    # flinwall(w,wb,'Jlo', gam, interf1, gh, im, jm)
                    # fwall(w,'Jlo', gam, interf1, gh, im, jm)
                    # wbr = _np.copy(wb, order='F')
                    # flinjn(w, wbr, prr1, gh, gh, gh, gh, im, jm, w, wb, prd2, gh, gh, gh, gh, im, jm, tr1)
                    # wb[-gh:,:,:] = wbr[-gh:,:,:]
                    # wbr = _np.copy(wb, order='F')
                    # flinjn(w, wbr, prr2, gh, gh, gh, gh, im, jm, w, wb, prd1, gh, gh, gh, gh, im, jm, tr2)
                    # wb[:gh,:,:] = wbr[:gh,:,:]
                    # fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
                    # fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)

                    if sch in ['dnc', 'dnc_hll', 'dnc_hllc', 'rbc']:
                        flinsch(res, resb, w, wb, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                    else:
                        flinsch(res, resb, w, wb, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)

                    # flinnoref(w,wb,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                    # fnoref(w,wbd,'Jhi',interf2,nx,ny,gam,gh,im,jm)
                    # flinwall(w,wb,'Jlo', gam, interf1, gh, im, jm)
                    # fwall(w,'Jlo', gam, interf1, gh, im, jm)
                    # wbr = _np.copy(wb, order='F')
                    # flinjn(w, wbr, prr1, gh, gh, gh, gh, im, jm, w, wb, prd2, gh, gh, gh, gh, im, jm, tr1)
                    # wb[-gh:,:,:] = wbr[-gh:,:,:]
                    # wbr = _np.copy(wb, order='F')
                    # flinjn(w, wbr, prr2, gh, gh, gh, gh, im, jm, w, wb, prd1, gh, gh, gh, gh, im, jm, tr2)
                    # wb[:gh,:,:] = wbr[:gh,:,:]
                    # fjn(w,prr1,gh,gh,gh,gh,im,jm,w,prd2,gh,gh,gh,gh,im,jm,tr1)
                    # fjn(w,prr2,gh,gh,gh,gh,im,jm,w,prd1,gh,gh,gh,gh,im,jm,tr2)    
                        
                    # f_misc.computejacobianfromjv_withjn_dbyvol(Jac,IA,JA,wb,m,l,k,gh,im,jm, vol)
                    f_misc.computejacobianfromjv_withjn(Jac,IA,JA,wb,m,l,k,gh,im,jm)

        ## Remove the zero stored

        mini = 2.e-16
        IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
        # print nbentry
        nbentry = _np.shape(Jac)[0]
        print(nbentry)

        Jacs_adj = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

        print("** Writing adjoint **")

        dir = './'
        dir = './BASEFLOW_CYL/'

        viewer = PETSc.Viewer().createBinary(dir+'Jac_adj', 'w')
        viewer(Jacs_adj)




