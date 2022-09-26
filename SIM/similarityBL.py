# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Compute the similarity solution for compressible BL
# Only for Flat Plate (no pressure gradient)

# Adiabatic Wall
# Minor changes required for Isotherm Wall

import numpy as _np
import matplotlib.pyplot as plt
import pickle


def similarity(Mach, g_wall, isotherm=False, prandtl=0.72, gam=1.4, eps=1.e-4, infini=10., tol=1.e-6, itermax=30):

    N_eta = int(infini/eps)
    eta = _np.linspace(0.,infini,N_eta)

    f   = _np.zeros((N_eta))
    g   = _np.zeros((N_eta-2))

    df  = _np.zeros((N_eta-1))

    ### Computation of the function f (used to compute u, v and the function g)

    f[0] = 0.	# v(0)=0
    f[1] = f[0] # f'(0)=0 <=> u(0)=0

    df_inf = (f[-1]-f[-2])/eps # f' on the infinite, f'(inf)=1 <=> u(inf)=Ue

    iter = 0

    while iter < itermax and abs(df_inf - 1.) > tol:

        if iter == 0:
            # First guess for f''(0)
            f2_cur = eps*eps * 0.3 # f[2] = eps**2 * f''(0)
        elif iter == 1:
            # Second guess
            f2_cur = eps**2. * 0.4
            f2_old = f[2]
        else:
            # Newton guess
            f2_cur = f[2] - (df_inf - 1.)/((df_inf-df_inf_old)/(f[2]-f2_old))
            f2_old = f[2]

        f[2] = f2_cur
        for k in range(3,_np.shape(f)[0]):
            # RHSf = - f[k-3] + 3*f[k-2] - 3*f[k-1] + eps * f[k-3] * (2*f[k-3] - 5*f[k-2] + 4*f[k-1])
            # f[k] = - RHSf / (1 - eps * f[k-3])
            ## Lower order but mathematically more accurate
            f[k] = - (- f[k-3] + 3.*f[k-2] - 3.*f[k-1] + eps * f[k-3] * (f[k-3] - 2.*f[k-2] + f[k-1]))

        df_inf_old = df_inf
        df_inf = (f[-1]-f[-2])/eps
        iter = iter + 1
        # print abs(df_inf - 1.)

    if iter == itermax:
        print(df_inf - 1.)
        print('No convergence in f computation after max iterations')
        return

    for k in range(_np.shape(df)[0]):
        df[k] = (f[k+1] - f[k])/eps

    # print iter
    # print "f''(0) = ", f[2]/eps**2

    ### Computation of the function g (used to compute rho and E)

    g_inf = 0. # g on the infinite, g(inf)=1 <=> hi(inf)=Hie
    Kmach = (1. - prandtl) * ((gam - 1.) * Mach**2.)/(1. + (gam - 1.)/2. * Mach**2.)

    iter = 0
    while iter < itermax and abs(g_inf - 1.) > tol:

        if iter == 0:
            # First guess for g(1)
            g1_cur = 1.
        elif iter == 1:
            # Second guess
            g1_cur = 0.9
            g1_old = g[1]
        else:
            # Newton guess
            g1_cur = g[1] - (g_inf - 1.)/((g_inf-g_inf_old)/(g[1]-g1_old))
            g1_old = g[1]

        g[1] = g1_cur

        if isotherm:
            g[0] = g_wall ## Isotherm wall <=> g(0)=g_wall
        else:
            g[0] = g[1] ## Adiabatic wall <=> g'(0)=0

        for k in range(2,_np.shape(g)[0]):
            # RHSg  = Kmach / eps**4. * ((- 1.833333*f[k-2] + 3.*f[k-1] - 1.5*f[k] + 0.333333*f[k+1])*(- f[k-2] + 3.*f[k-1] - 3.*f[k] + f[k+1]) + (2.*f[k-2] - 5.*f[k-1] + 4.*f[k] - f[k+1])**2.)
            ## Mathematically more accurate considering the orders
            RHSg  = Kmach / eps**4. * ((- 2.083333*f[k-2] + 4.*f[k-1] - 3.*f[k] + 1.333333*f[k+1] - 0.25*f[k+2])*(- 2.5*f[k-2] + 9.*f[k-1] - 12.*f[k] + 7.*f[k+1] - 1.5*f[k+2]) + (2.916667*f[k-2] - 8.666667*f[k-1] + 9.5*f[k] - 4.666667*f[k+1] + 0.916667*f[k+2])**2.)
            g[k]  = - g[k-2] + 2*g[k-1] - eps * prandtl * f[k-2] * (g[k-1] - g[k-2]) + eps**2. * RHSg

        g_inf_old = g_inf
        g_inf = g[-1]
        iter = iter + 1
        # print abs(g_inf - 1.)

    if iter == itermax:
        print(g_inf - 1.)
        print('No convergence in g computation after max iterations')
        return

    # print iter
    # print "g(0) = ", g[0]

    return eta, f, df, g

if __name__ == '__main__':
	Mach	= 4.5  #4.5
	prandtl = 0.72 #0.72
	gam		= 1.4

	tol 	= 1.e-6
	itermax = 30
	infini  = 10.	# Considered as the infinite value (beyond the edge), maximal value computed
	eps 	= 1.e-4 # Step for the functions f and g

	# eta, f, df, g = similarity(Mach,prandtl,gam,eps,infini,tol,itermax)
	eta, f, df, g = similarity(Mach,g_wall,prandtl,gam,eps,infini,tol,itermax)
	# fileout = open('simfunctions_mach_%s_eps_%s_max_%s' % (Mach,eps,infini),'wb')
	# pickle.dump({'eta' : eta, 'f' : f, 'df' : df, 'g' : g}, fileout)

	plt.figure(1)
	plt.title('f')
	plt.plot(eta,f)
	plt.figure(2)
	plt.title('df')
	plt.plot(eta[:-1],df)
	plt.figure(3)
	plt.title('g')
	plt.plot(eta[:-2],g)
	plt.figure(4)
	plt.title('eta*df - f')
	plt.plot(eta[:-1],eta[:-1]*df-f[:-1])
	plt.figure(5)
	plt.title('g - df**2')
	plt.plot(eta[:-2],g-df[:-1]**2)

	plt.show()



