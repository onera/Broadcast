# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

import numpy as _np


def distance(x1, x2):
    return _np.abs(x2 - x1)

def searchLeft(X, x0, start):
    ## search in a sorted 1D array

    Nx = _np.shape(X)[0]
    current = start
    while current <= Nx-1 and X[current] < x0:
        current = current + 1
    if current > Nx-1:
        return Nx-1    
    elif current == start:
        return None
    else: 
        return current - 1

def interpgrid(Xin, Yin, Zin, Xout, Yout):

    ## Structured cartesian rectangular grid & New grid borders included in the previous one (otherwise extrapolation at zero order of the border-points)
    ## Xin, Yin & Zin are (Nx,Ny) sorted numpy array
    ## Xout, Yout & Zout are (Nx2,Ny2) numpy array

    Nx = _np.shape(Xin)[0]
    Ny = _np.shape(Yin)[1]
    Nx2 = _np.shape(Xout)[0]
    Ny2 = _np.shape(Yout)[1]

    Ztmp = _np.zeros((Nx2,Ny))
    Zout = _np.zeros((Nx2,Ny2))
    
    ptstart = 0
    for i in range(Nx2):
        
        ptleft = searchLeft(Xin[:,0], Xout[i,0], ptstart)
        ptstart = ptleft 
        if ptleft == Nx-1:
            Ztmp[i,:] = Zin[Nx-1,:]
        elif ptleft == None:
            Ztmp[i,:] = Zin[0,:]
            ptstart = 0
        else:
            ptright = ptleft + 1
            dleft = distance(Xin[ptleft,0], Xout[i,0])
            dright = distance(Xin[ptright,0], Xout[i,0])
            ## Linear
            Ztmp[i,:] = (Zin[ptleft,:] * dright + Zin[ptright,:] * dleft) / (dleft + dright)

    ptstart = 0
    for j in range(Ny2):
        
        ptleft = searchLeft(Yin[0,:], Yout[0,j], ptstart)
        ptstart = ptleft 
        if ptleft == Ny-1:
            Zout[:,j] = Ztmp[:,Ny-1]
        elif ptleft == None:
            Zout[:,j] = Ztmp[:,0]   
            ptstart = 0 
        else:
            ptright = ptleft + 1
            dleft = distance(Yin[0,ptleft], Yout[0,j])
            dright = distance(Yin[0,ptright], Yout[0,j])
            ## Linear
            Zout[:,j] = (Ztmp[:,ptleft] * dright + Ztmp[:,ptright] * dleft) / (dleft + dright)        

    return Zout

def interpline(Xin, Yin, Xout):

    Nx = _np.shape(Xin)[0]
    Nx2 = _np.shape(Xout)[0]

    Yout = _np.zeros(Nx2)

    ptstart = 0
    for i in range(Nx2):
        
        ptleft = searchLeft(Xin[:], Xout[i], ptstart)
        ptstart = ptleft 
        if ptleft == Nx-1:
            Yout[i] = Yin[Nx-1]
        elif ptleft == None:
            Yout[i] = Yin[0]
            ptstart = 0
        else:
            ptright = ptleft + 1
            dleft = distance(Xin[ptleft], Xout[i,])
            dright = distance(Xin[ptright], Xout[i])
            ## Linear
            Yout[i] = (Yin[ptleft] * dright + Yin[ptright] * dleft) / (dleft + dright)

    return Yout

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import timeit

    Nx1 = 800
    Ny1 = 100

    x1  = _np.linspace(0., 1., Nx1)
    X1  = _np.zeros((Nx1,Ny1))
    for k in range(Ny1):
        X1[:,k] = x1

    y1  = _np.linspace(0.,0.1, Ny1)
    Y1  = _np.zeros((Nx1,Ny1))
    for k in range(Nx1):
        Y1[k,:] = y1
    
    print(_np.shape(X1))
    print(_np.shape(Y1))

    t0 = timeit.time.time()
    Z1 = _np.random.rand(Nx1,Ny1)
    t1 = timeit.time.time()
    print('Time random', t1 - t0)

    Nx2 = 900
    Ny2 = 100

    xini = X1[0,0]
    xend = X1[-1,0]
    yend = Y1[0,-1]

    x2  = _np.linspace(xini, xend, Nx2)
    X2  = _np.zeros((Nx2,Ny2))
    for k in range(Ny2):
        X2[:,k] = x2

    y2  = _np.linspace(0.,yend, Ny2)
    Y2  = _np.zeros((Nx2,Ny2))
    for k in range(Nx2):
        Y2[k,:] = y2


    t0 = timeit.time.time()
    Z2 = interpgrid(X1, Y1, Z1, X2, Y2)  
    t1 = timeit.time.time() 
    print('Time interp', t1 - t0) 


    plt.figure(1)
    plt.title('Original')
    plt.contourf(X1,Y1,Z1, 100)
    plt.colorbar()

    plt.figure(2)
    plt.title('Interpolated')
    plt.contourf(X2,Y2,Z2, 100)
    plt.colorbar()

    plt.show()



