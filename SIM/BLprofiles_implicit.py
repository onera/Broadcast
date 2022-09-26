# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

## From the similarity functions and a grid, build the profiles for compressible BL
## only for flat plate and no pressure gradient

import numpy as _np
import matplotlib.pyplot as plt
import pickle
import os.path

import SIM.similarityBL as BL


def BLprofile(xc, yc, Mach, dphys, isotherm=False, prandtl=0.72, gam=1.4, isplot=False, damped=False, twall=1.):

    # Zero Pressure Gradient BL
    X = xc[:,0]
    Y = yc[0,:]
    
    infini  = 10.    # Considered as the infinite value (beyond the edge), maximal value computed
    eps     = 1.e-4 # Step for the functions f and g

    tinf    =  dphys['T0']
    g_wall  =  twall / ((1. + (gam -1.)/2. * Mach**2 )* tinf)

    if os.path.isfile('simfunctions_mach_%s_eps_%s_max_%s' % (Mach,eps,infini)):

        filein = open('simfunctions_mach_%s_eps_%s_max_%s' % (Mach,eps,infini),'rb')
        dict_pickle = pickle.load(filein)
        # print 'ok path'
        eta, f, df, g = dict_pickle['eta'], dict_pickle['f'], dict_pickle['df'], dict_pickle['g']

    else:

        tol     = 1.e-6
        itermax = 30
        eta, f, df, g = BL.similarity(Mach,g_wall,isotherm,prandtl,gam,eps,infini,tol,itermax)

    # print 'ok'
    # print _np.shape(eta)

    Nx = _np.shape(X)[0]
    Ny = _np.shape(Y)[0]

    u   = _np.zeros((Nx,Ny))
    v   = _np.zeros((Nx,Ny))
    rho = _np.zeros((Nx,Ny))
    E   = _np.zeros((Nx,Ny))

    etagrid = _np.zeros((Nx,Ny))

    cs      =  dphys['cs']
    tref    =  dphys['Ts']
    musuth  =  dphys['musuth']
    rgaz    =  dphys['rgaz']
    # tinf    =  dphys['T0']
    muinf   =  dphys['mu0']
    Runit   =  dphys['Runit']

    uinf    = Mach * _np.sqrt(gam * rgaz * tinf)
    rhoinf  = Runit * muinf / uinf
    cv      = rgaz / (gam - 1.)

    Twall   = (1. + (gam -1.)/2. * Mach**2 )* tinf * g[0]
    # Twall   = (1. + (gam -1.)/2. * Mach**2 )* tinf * 1.
    # print Twall/tinf

    Const   = musuth/muinf * tinf/tref * _np.sqrt(Twall/tref)* (tref + cs)/(Twall + cs)
    # Const   = _np.sqrt(tref/tinf) * (tinf + cs)/(tref + cs)
    # Const   = 1.
    # print Const

    # Hypo_rho = rhoinf / g ## Assumption for the density

    ## Other assumptions for the density
    # plt.plot(eta[:-2],Hypo_rho/rhoinf)
    RHSr     = g * (gam * cv * tinf ) + (g - df[:-1]**2.) * 0.5*uinf**2.
    Hypo_rho = rhoinf * (gam * cv * tinf) / RHSr
    # Hypo_rho = rhoinf * _np.ones(_np.shape(g))
    # Hypo_rho = rhoinf * ((gam * cv * tinf) / RHSr)**(1./(gam - 1.))
    # Hypo_rho = rhoinf / g
    # plt.plot(eta[:-2],Hypo_rho/rhoinf)
    # plt.show()
    invu2 = 1./uinf
    y=10.

    # print gam * cv * tinf
    # print 0.5*uinf**2.

    max_ind_eta = _np.shape(eta)[0]-1 - 2

    for i in range(Nx):

        # RHSr     = g * (gam * cv * tinf ) + (g - df[:-1]**2. - Const*muinf/(2.*rhoinf*uinf*X[i])*(eta[:-2]*df[:-1]-f[:-2]) ) * 0.5*uinf**2.
        # Hypo_rho = rhoinf * (gam * cv * tinf) / RHSr

        residu  = _np.sqrt(uinf / (2.*Const * rhoinf * muinf)) / X[i]**(0.5)

        u[i,0]        = 0.
        rho[i,0]      = Hypo_rho[0]
        v[i,0]        = 0.
        hi0           = (gam * cv * tinf + 0.5*uinf**2.) * g[0]
        E[i,0]        = hi0 / gam + (gam - 1.)/gam * 0.5*v[i,0]**2.

        eta_real      = 0.
        etagrid[i,0]  = eta_real
        ind_sol       = 0

        for j in range(1,Ny):

            eta_real = eta_real + residu * (Y[j] - Y[j-1]) * Hypo_rho[ind_sol]

            if eta_real >= eta[-3]:
                ind_sol = max_ind_eta
            else:
                ind_sol = _np.searchsorted(eta, eta_real)

            # print 'eta found: ', eta[ind_sol]

            u[i,j]   = uinf * df[ind_sol]
            rho[i,j] = Hypo_rho[ind_sol]
            v[i,j]   = _np.sqrt((Const * rhoinf * muinf * uinf)/(2.*X[i])) * 1./rho[i,j] * (eta[ind_sol] * df[ind_sol] - f[ind_sol])
            # v[i,j]   = _np.sqrt((Const * rhoinf * muinf * uinf)/(2.*X[i])) * 1./rhoinf * (eta[ind_sol] * df[ind_sol] - f[ind_sol])
            # Magical Tricks in hllc3
            # v[i,j]   = 2.*_np.pi*_np.sqrt((Const * rhoinf * muinf * uinf)/(2.*X[i])) * 1./rhoinf * (eta[ind_sol] * df[ind_sol] - f[ind_sol])
            # v[i,j]   = 6.*_np.sqrt((Const * rhoinf * muinf * uinf)/(2.*X[i])) * 1./rho[i,j] * (eta[ind_sol] * df[ind_sol] - f[ind_sol])
            hi       = (gam * cv * tinf + 0.5*uinf**2.) * g[ind_sol]
            if damped:
                etagrid[i,j] = eta_real
                utest = (u[i,j]-uinf)
                ytest = eta[ind_sol]
                if (utest*utest*invu2 < 0.95):
                    if (ytest < y): y = ytest
                    y2 = abs(y-ytest)
                    v[i,j] *= _np.exp(-y2**1.5)
            E[i,j]   = hi / gam + (gam - 1.)/gam * 0.5*(u[i,j]**2. + v[i,j]**2.)
            etagrid[i,j] = eta_real


    # for i in range(Nx-1):
    #     for j in range(1,Ny):
    #         v[i,j] = ( ( 2*rho[i,j]*u[i,j] - rho[i,j]*u[i+1,j] - rho[i+1,j]*u[i,j] ) / (X[i+1] - X[i]) * (Y[j] - Y[j-1]) + rho[i,j] * v[i,j-1] ) / ( 2*rho[i,j] - rho[i,j-1] )
    # v[-1,:] = v[-2,:]

    # plt.figure(10)
    # plt.title('eta')
    # # plt.contourf(X,Y,_np.transpose(etagrid),21)
    # CONT = plt.contour(X,Y,_np.transpose(etagrid),levels=_np.linspace(0.,10.,20))
    # plt.clabel(CONT, inline=1, fontsize=10)
    # # plt.colorbar()

    return rho/rhoinf, u/uinf, v/uinf, E/(cv*tinf + 0.5*uinf**2.)

if __name__ == '__main__':

    import meshBL as mesh

    Mach    = 4. #4.5
    prandtl = 0.72 #0.72
    gam     = 1.4

    Nx    = 100
    Ny    = 110

    high  = 0.03    #0.03
    L     = 0.4 #0.6
    x_ini = 0.01    #0.1

    Xtmp  = _np.linspace(0.,L,Nx)
    X     = Xtmp + x_ini
    Y1    = _np.linspace(0.,high,Ny)
    # Y        = _np.linspace(0.,high,Ny)
    Y     = mesh.stretch_tanh(Y1, 1., 3., 2.)

    dphys = dict()
    dphys['Ts']      = 273.15
    dphys['cs']      = 110.4
    dphys['musuth']  = 1.716e-5
    dphys['rgaz']    = 287.1
    dphys['mu0']     = 4.367e-6        #4.367e-6  #1.789e-5
    dphys['T0']      = 65.15    #65.15     #288.
    dphys['P0']      = 728.312    #728.312    #104001.
    # dphys['Runit']   = 5.8e6    #6.49e6   #5.8e6

    uinf    = Mach * _np.sqrt(gam * dphys['rgaz'] * dphys['T0'])

    # rhoinf  = dphys['Runit'] * dphys['mu0'] / uinf
    rhoinf  = dphys['P0'] / (dphys['rgaz']*dphys['T0'])
    dphys['Runit'] = uinf*rhoinf/dphys['mu0']
    print(dphys['Runit'])

    cv      = dphys['rgaz'] / (gam - 1.)

    rho, u, v, E = BLprofile(X, Y, Mach, dphys, prandtl)

    T = ((cv*dphys['T0'] + 0.5*uinf**2)*E - 0.5*uinf**2*(u**2 + v**2))/cv / dphys['T0']


    t_ss = []
    u_ss = []
    v_ss = []
    ro_ss= []
    var_ss=[]
    for i in range(0,Nx,1):
        for j in range(0,Ny,1):
            var_ss.append(Y[j]*_np.sqrt(rhoinf*uinf/(dphys['mu0']*X[i]))*_np.sqrt(1./2))
            # var_ss.append(Y[j]*_np.sqrt(rhoinf*uinf/(dphys['mu0']*X[i])))
            t_ss.append(T[i,j])
            u_ss.append(u[i,j])
            v_ss.append(v[i,j])
            ro_ss.append(rho[i,j])


    msize = 7.

    if isplot:
        plt.figure(1)
        plt.title('U')
        plt.xlabel('u/uinf')
        plt.ylabel(r"$y/x*\sqrt{{Re}_x/2}$")
        plt.plot(u_ss,var_ss,'.b',label='My software')
        plt.grid()
        plt.ylim(0.,15.)
        plt.xlim(0.,1.1)
        plt.legend()

        plt.figure(2)
        plt.title('T')
        plt.xlabel('T/Tinf')
        plt.ylabel(r"$y/x*\sqrt{{Re}_x/2}$")
        plt.plot(t_ss,var_ss,'.b',ms=msize,label='My software')
        plt.grid()
        plt.ylim(0.,15.)
        plt.xlim(0.9,4.5)
        plt.legend()


        # plt.figure(3)
        # plt.title('V')
        # plt.xlabel('v/uinf')
        # plt.ylabel(r"$y/x*\sqrt{{Re}_x/2}$")
        # plt.plot(v_ss,var_ss,'.b',ms=msize,label='My software')
        # plt.grid()
        # plt.ylim(0.,15.)
        plt.xlim(0.,1.1)
        # plt.legend()

        plt.figure(4)
        plt.title(r"$\rho$")
        plt.xlabel('rho/rhoinf')
        plt.ylabel(r"$y/x*\sqrt{{Re}_x/2}$")
        plt.plot(ro_ss,var_ss,'.b',ms=msize,label='My software')
        plt.grid()
        plt.ylim(0.,15.)
        # plt.xlim(0.,1.1)
        plt.legend()


        # plt.figure(1)
        # plt.title('U')
        # plt.contourf(X,Y,_np.transpose(u),11)
        # plt.colorbar()

        # plt.figure(2)
        # plt.title('v')
        # plt.contourf(X,Y,_np.transpose(v),11)
        # plt.colorbar()

        # plt.figure(3)
        # plt.title('rho')
        # plt.contourf(X,Y,_np.transpose(rho),11)
        # plt.colorbar()

        # plt.figure(4)
        # plt.title('E')
        # plt.contourf(X,Y,_np.transpose(E),11)
        # plt.colorbar()

        # plt.figure(5)
        # plt.title('U(y)')
        # plt.grid()
        # Nplot = 10
        # for k in range(Nplot):
            # indplot = k * Nx / Nplot
            # plt.plot(u[indplot,:],Y)

        # plt.figure(6)
        # plt.title('T')
        # plt.contourf(X,Y,_np.transpose(T),11)
        # plt.colorbar()



        plt.show()


