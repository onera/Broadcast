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
import srcfv.f_geom   as f_geom
import srcfv.f_bnd     as f_bnd
import srcfv.f_sch     as f_sch
import srcfv.f_lhs     as f_lhs
import srcfv.f_lin     as f_lin
# import srcfv.f_adj     as f_adj
import srcfv.f_norm    as f_norm
# FROM A.POULAIN Thesis
import misc.f_misc     as f_misc
import misc.PETSc_func as pet
import resolvent_all  as resol
import SIM
import SIM.BLprofiles_implicit as blsim
import f_init
import meshBL as mesh
import handleBC as handleBC

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

######################### Private functions ####################

# solve monoblock Boundary Layer

def bl2d(dgeom = dict(), dphys = dict(), dnum = dict(), compmode = 'direct', lf = list(), lflin = list(), out_dir = 'totodir', isresol= False):
    '''
    exemple of monoblock use of 2DTOY
    to simulate 2D laminar Boundary Layer flow
    '''
    os.system('mkdir -p %s' % out_dir)

    # get functions
    # lf = [finflow, foutflow, fnoref, fwall, fsch]
    # finflow  = lf[0]
    # foutflow = lf[1]
    # fnoref   = lf[2]
    # fwall    = lf[3]
    # fsch     = lf[4]
    # fsym     = lf[5]

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
    L        = dgeom['length']
    high     = dgeom['high']
    xini     = dgeom['xini']
    ym       = dgeom['ym']

    ite      = dnum['ite']
    cfl      = dnum['cfl']
    k2       = dnum['k2']
    k4       = dnum['k4']
    sch      = dnum['sch']
    order    = dnum['order']
    freqres  = dnum['freqres']
    freqsort = dnum['freqsort']

    # Set ghost cells dimension
    if sch == 'dnc':
        gh = (order+1) // 2
    else:
        gh = (order-1) // 2 + 1 # +1 to avoid grad exchanges in multiblock configurations

    # gh = gh+2    

    if compmode == 'direct':
        rkcoefs = dnum['rkcoefs']
    elif compmode == 'fixed_point':
        lasolver = dnum['lasolver']
        if lasolver == 'gmres':
            tol = dnum['tol']


    x  = _np.linspace(xini, xini+L , im+1)
    ## MESH v2
    Ny_in   = 80*jm//100   #80%     
    deltaBL = high/4      #high/4 
    percent = 0.02 
    
    Ny_out  = jm - Ny_in 
    Nend    = high/deltaBL
    y_int   = mesh.bigeom_stretch_in(Ny_in, deltaBL, percent)
    y_out   = mesh.exp_stretch_out(Ny_out, deltaBL, percent, Nend)
    y       = _np.concatenate((y_int, y_out)) 
    # y = mesh.bigeom_stretch_in(Ny_in, deltaBL, percent)

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

    dx = L/((im-1)*Lref) # adim done after muinf A.Poulain
    dy = (y[1]-y[0])/Lref
    sound = 1./dphys['Mach']
    dt = cfl * min(dy,dx) / (sound+1.)
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
    rey     = runit * L

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

    Roref = rhoinf
    Vref  = uinf
    Tref  = tinf
    # Roref = 1.
    # Vref  = 1.
    # Tref  = 1.
    # Lref  = 1.
    Pref  = Roref*Vref**2
    Cvref = Vref**2/Tref

    Eref   = Vref**2
    Rgpref = Cvref

    ## Adim with ref length
    # Lref   = 8.e-2  #8.e-2
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

    # Compute Geometry
    for i in range(im+1):
        x0[i+gh,:] = x[i]
    for j in range(jm+1):
        y0[:,j+gh] = y[j]

    # Adim Geom:
    x0 *= 1./Lref
    y0 *= 1./Lref

    # ym *= 1./Lref
    # y0 += ym 

    # f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf)
    f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)

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
    # road,uad,vad,Ead = blsim.BLprofile(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False)
    road,uad,vad,Ead = blsim.BLprofile(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False, twall=4.395*Tref, isotherm=False)
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

    # w[:, :, 2]     = 0.

    f_init.set_bndbl_2d(w, field, wbd, im, jm, gh)


    ######## Restart from a previous solution with exactly the same mesh
    import restart_init as ri
    # filet = './Wksp/dnc_5/state_atcenter_ite14.dat'
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

    #interfaces definitions (may be done at the begining)
    interfl = []
    interfdic = dict()
    interfloc = dict()
    BCdict = dict()

    # Ilo
    interf1      = _np.zeros((2,2), order='F')
    interf1[0,0] = 1  # imin
    interf1[0,1] = 1  # jmin
    interf1[1,0] = 1  # imax
    interf1[1,1] = jm # jmax 
    interfl.append('interf1')
    interfdic['interf1'] = interf1
    interfloc['interf1'] = 'Ilo'

    # Ihi
    interf2      = _np.zeros((2,2), order='F')
    interf2[0,0] = im # imin
    interf2[0,1] = 1  # jmin
    interf2[1,0] = im # imax
    interf2[1,1] = jm+gh # jmax 
    interfl.append('interf2') 
    interfdic['interf2'] = interf2
    interfloc['interf2'] = 'Ihi'

    # Jlo
    interf3      =  _np.zeros((2,2), order='F')
    interf3[0,0] = 1-gh # imin 
    interf3[0,1] = 1  # jmin
    interf3[1,0] = im+gh # imax
    interf3[1,1] = 1  # jmax
    interfl.append('interf3')
    interfdic['interf3'] = interf3
    interfloc['interf3'] = 'Jlo'

    # Jhi
    interf4      =  _np.zeros((2,2), order='F')
    interf4[0,0] = 1-gh  # imin
    interf4[0,1] = jm # jmin
    interf4[1,0] = im # imax
    interf4[1,1] = jm # jmax
    interfl.append('interf4')
    interfdic['interf4'] = interf4
    interfloc['interf4'] = 'Jhi'

    # Jlow before plate Leading Edge
    # interf5      =  _np.zeros((2,2), order='F')
    # interf5[0,0] = 1-gh  # imin
    # interf5[0,1] = 1 # jmin
    # interf5[1,0] = 0 # imax
    # interf5[1,1] = 1 # jmax
    # interfl.append('interf5')
    # interfdic['interf5'] = interf5
    # interfloc['interf5'] = 'Jlo'

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
    

    wlist = [w]
    BCtypel = [lf, lflin]
    wbddic = dict()
    fielddic = dict()
    twalldic = dict()
    velwalldic = dict()
    prldic = dict()
    trldic = dict()

    BCdict['interf1'] = 'finflow'
    BCdict['interf2'] = 'foutflow'
    BCdict['interf3'] = 'fwall'
    BCdict['interf4'] = 'fnoref'
    # BCdict['interf5'] = 'fsym'

    wbddic['interf4'] = wbd
    fielddic['interf1'] = field

    # prldic['interf1'] = prl
    # trldic['interf1'] = trl
    # BCdict['interf1'] = 'fjn'

    checkint = handleBC.checkinterf(interfl, interfdic, interfloc, im, jm, gh)
    if checkint != 'GoodJob':
        print(checkint)
        exit()
    interfl = handleBC.sortBC(interfl, interfdic, interfloc, im, jm)

    # modeBC = 0  #apply BC on w
    # wlist = applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic, fielddic, twalldic, velwalldic, prldic, trldic)
    # w = wlist[0]

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
                
                finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
                # finflow(w,'Ilo', interf1, field,im,jm) 
                fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)         
                foutflow(w,'Ihi', interf2, im, jm, gh)
                fwall(w,'Jlo', gam, interf3, gh, im, jm)

                # Compute spatial discretization
                if sch == 'dnc':
                    # fwall needed for dissipation near bnd_wall
                    fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
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
        # Boundary on state vector
        finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
        # finflow(w,'Ilo', interf1, field,im,jm) 
        fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
        foutflow(w,'Ihi', interf2, im, jm, gh)
        fwall(w,'Jlo', gam, interf3, gh, im, jm)
        # fwall(w, tinf,'Jlo', gam, rgaz, interf3, gh, im, jm)

        filename = out_dir + '/initialisation_gh.dat'
        __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

        lmax = 10  #4
        for it in range(1,ite+1):
            # if it > 4000: cfl = min(0.5 + (it-4000)*0.001 ,1.)
            # if it > 8000: cfl = min(3. + (it-8000)*0.001 ,10.)
            # Compute spatial discretization
            if sch == 'dnc':
                # fwall needed for dissipation near bnd_wall
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
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
            finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
            # finflow(w,'Ilo', interf1, field,im,jm) 
            fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            foutflow(w,'Ihi', interf2, im, jm, gh)
            fwall(w,'Jlo', gam, interf3, gh, im, jm)
            # fwall(w, tinf,'Jlo', gam, rgaz, interf3, gh, im, jm)

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
        # get functions
        # lflin = [flininflow, flinoutflow, flinnoref, flinwall, flinsch]
        # flininflow  = lflin[0]
        # flinoutflow = lflin[1]
        # flinnoref   = lflin[2]
        # flinwall    = lflin[3]
        # flinsch     = lflin[4]
        # flinsym     = lflin[5]

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

        wd   = _np.zeros((im+2*gh, jm+2*gh, 5), order='F')

        wlist.append(wd)

        timeconstructjac = 0.
        timeremove  = 0.
        timecoefdiag = 0.
        timejaccsc   = 0.
        timejacinv   = 0.
        cfl = 1.e5
        dt  = cfl * (yc[gh,gh+1] - yc[gh,gh]) / (1./mach + 1.)
        dtm1 = 1./dt
        for it in range(1,ite+1):

            wd   = _np.zeros((im+2*gh, jm+2*gh, 5), order='F')
            resd = _np.zeros((im+2*gh, jm+2*gh, 5), order='F')
            nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
            Jac = _np.zeros((nbentry), order='F')
            IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
            JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')

            # Boundary on state vector
            # # finflow(w,'Ilo', interf1, field,im,jm)      
            # finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
            # fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            # foutflow(w,'Ihi', interf2, im, jm, gh)
            # fwall(w,'Jlo', gam, interf3, gh, im, jm)
            # # twall = 1.5*tinf
            # # fwall(w, twall,'Jlo', gam, rgaz, interf3, gh, im, jm)

            modeBC = 0  #apply BC on w
            wlist[0] = w
            wlist = handleBC.applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic, fielddic, twalldic, velwalldic, prldic, trldic)
            w = wlist[0]

            # fsym(w,'Jlo', interf5, gh, im, jm)
                        
            filename = out_dir + '/initialisation_gh.dat'
            __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

            # Compute spatial discretization
            if 'polar' in routinesch:
                fsch(res, w, ym, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
            elif 'dnc' in sch:    
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                # fsch(res, w, twall, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
            else:
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)
            norm, ninf = f_norm.compute_norml2inf(res ,im, jm, gh)
            if it == 1:
                for i in range(5):
                    norm[i] = max(norm[i],1.e-15)
                    ninf[i] = max(ninf[i],1.e-15)
                norm0m1 = 1./norm
                ninf0m1 = 1./ninf

            # Iterative Linear Algebra Loop to solve Newton: Solve resd  = -res ==> J w_sol = - res
            if lasolver == 'gmres':
                nl2res = 1.
                print(" GMRES Not Yet implemented ")
                sys.exit(2)

            else:
                # construct Jacobian
                timeinjac = timeit.time.time()
                ## relaxed on diag:
                r = _np.max([norm[:3]*norm0m1[:3], ninf[:3]*ninf0m1[:3]])
                cflm1 = r*dtm1
                print('iter = ', it)
                print("1/cfl = ", cflm1)
                print(norm)
                coefdiag = cflm1 * vol[gh:-gh,gh:-gh]
                for m in range(5):
                    for l in range(1 + 2*gh):
                        for k in range(1 + 2*gh):
                            wd *= 0.
                            f_misc.testvector(wd,m,l,k,gh,im,jm)

                            # w[:gh,:,:]  = 0.
                            # w[:,:gh,:]  = 0.
                            # w[-gh:,:,:] = 0.
                            # w[:,-gh:,:] = 0.

                            # flininflow(w,wd,'Ilo',interf1,field,nx,ny,gam,im,jm)                
                            # # finflow(w,'Ilo', interf1, field,im,jm) 
                            # flinnoref(w,wd,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
                            # flinoutflow(w,wd,'Ihi', interf2, im, jm, gh)
                            # foutflow(w,'Ihi', interf2, im, jm, gh)
                            # flinwall(w,wd,'Jlo', gam, interf3, gh, im, jm)
                            # # flinwall(w,wd,twall,'Jlo', gam,rgaz, interf3, gh, im, jm)

                            modeBC = 1  #apply BC on wd
                            wlist[0] = w
                            wlist[1] = wd
                            wlist = handleBC.applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic, fielddic, twalldic, velwalldic, prldic, trldic)
                            w  = wlist[0]
                            wd = wlist[1]

                            # flinsym(w, wd,'Jlo', interf5, gh, im, jm)
                            
                            if 'polar' in routinesch:
                                flinsch(res, resd, w, wd, ym, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)                            
                            elif 'dnc' in sch:      
                                flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm) 
                                # flinsch(res, resd, w, wd, twall, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)    
                            else:
                                flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)
                            
                            ## Finite Difference method (FD)
                            # epsilon = 1.e-6  #1.e-8 at order 1  #1.e-6 at order 2
                            # res1  = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
                            # w2 = w + epsilon*wd
                            # w2[:gh,:,:]  = 0.
                            # w2[:,:gh,:]  = 0.
                            # w2[-gh:,:,:] = 0.
                            # w2[:,-gh:,:] = 0.
                            # # finflow(w2,'Ilo', interf1, field,im,jm)
                            # finflow(w2,'Ilo',interf1,field,nx,ny,gam,im,jm)
                            # fnoref(w2,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)          
                            # foutflow(w2,'Ihi', interf2, im, jm, gh)
                            # fwall(w2,'Jlo', gam, interf3, gh, im, jm)
                            # if 'dnc' in sch: 
                            #     fsch(res1, w2, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)   
                            # else:
                            #     fsch(res1, w2, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)
                            # resd = (res1 - res) / epsilon
                            ## OR Finite Difference at order 2    
                            # resm1 = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
                            # wm = w - epsilon*wd
                            # wm[:gh,:,:]  = 0.
                            # wm[:,:gh,:]  = 0.
                            # wm[-gh:,:,:] = 0.
                            # wm[:,-gh:,:] = 0.
                            # # finflow(wm,'Ilo', interf1, field,im,jm)    
                            # finflow(wm,'Ilo',interf1,field,nx,ny,gam,im,jm)  
                            # fnoref(wm,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)    
                            # foutflow(wm,'Ihi', interf2, im, jm, gh)
                            # fwall(wm,'Jlo', gam, interf3, gh, im, jm)
                            # if 'dnc' in sch:
                            #     fsch(resm1, wm, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)       
                            # else:
                            #     fsch(resm1, wm, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, im, jm)               
                            # resd = (res1 - resm1) / (2*epsilon)

                            f_misc.computejacobianfromjv_relaxed(Jac,IA,JA,resd,m,l,k,gh,im,jm,coefdiag)
                            # if it == ite:
                            #     f_misc.computejacobianfromjv(Jac,IA,JA,resd,m,l,k,gh,im,jm)
                            # else:
                            #     f_misc.computejacobianfromjv_relaxed(Jac,IA,JA,resd,m,l,k,gh,im,jm,coefdiag)

                ## Remove the zero stored
                timeinremove = timeit.time.time()
                timeconstructjac += (timeinremove - timeinjac)/ite
                mini = 2.e-16
                IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
                nbentry = _np.shape(Jac)[0]
                print(nbentry)
                timeoutremove = timeit.time.time()
                timeremove += (timeoutremove - timeinremove)/ite

                # import scipy.sparse as sp
                # Jacs = sp.csr_matrix((Jac, (IA, JA)), shape=(im*jm*5, im*jm*5))
                Jacs = pet.createMatPetscCSR(IA, JA, Jac, im*jm*5, im*jm*5, 5*(2*gh+1)**2)

                timeoutjaccsc = timeit.time.time()
                timejaccsc += (timeoutjaccsc - timeoutremove)/ite

                # plt.figure()
                # plt.spy(Jacs)
                # plt.show()

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
            print('Time: ', (timeconstructjac+timeremove+timecoefdiag+timejaccsc+timejacinv)*ite)

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

