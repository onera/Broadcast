# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import numpy as _np
import math as m
import matplotlib.pyplot as plt

def stretch_tanh(x0, a, b, c):
    '''
    From the input x0, this method use the function

               x0 |--> a/2. * tanh(b * (x0 - c))

    to strech the 1-D coordinate.

        c : location of the maximum gradian (in pourcentage)
        a :

    '''
    x1 = (x0 - x0[0]) / (x0[-1] - x0[0])
    x2 = a / 2. * _np.tanh(b * (x1 - c))
    x3 = (x2 - x2[0]) / (x2[-1] - x2[0])
    x4 = x3*(x0[-1] - x0[0]) + x0[0]
    return x4

def stretch_log(x0, d, e, f1, f2):
    '''
    ...
    '''
    # x0 = (x0 - x0[0]) / (x0[-1] - x0[0])
    num = _np.cosh(e * (x0 - f1))
    den = _np.cosh(e * (x0 - f2))
    x1 = (d - 1)/(2. * e)*_np.log(num / den)
    x1 = (x1-x1[0])/(x1[-1]-x1[0])
    return x1

def stretch_combined(x0, a, b, c, d, e, f1, f2):
    '''
    ...
    '''
    x1 = a / 2. * _np.tanh(b * (x0 - c))
    num = _np.cosh(e * (x0 - f1))
    den = _np.cosh(e * (x0 - f2))
    x2 = (d - 1)/(2. * e)*_np.log(num / den)
    x3 = x0 + x1 + x2
    x3 = x3 - x3[0]
    x3 = x3 / (x3[-1] - x3[0])
    return x3

def stretching_plot(r, x):
    plt.figure()
    p1, = plt.plot(r, r, 'k-+')
    p2, = plt.plot(r, x, 'r-+')
    plt.legend([p1, p2], ['r', 'x(r)'])
    plt.xlabel('r')
    # plt.show()


def exp_stretch_out(N_out, delta, percent, Nend):
    
    x_0 = _np.linspace(1.,N_out,N_out)
    percent_max = _np.log(Nend)/N_out
    # print 2*percent_max
    return delta*_np.exp((percent + (percent_max - percent) * x_0/N_out) * x_0) 


def bigeom_stretch_in(N_in, delta, percent):
    
    x0 = _np.linspace(0,N_in-1,N_in)
    b = 1. + percent
    a = 1. + b
    r1 = b
    y_min = delta * (b - 1.) / ( b**N_in - 1. )
    mu = - y_min / (b - 1.)
    xi = b / (b - 1.) * y_min
    # print r1, xi, mu
    print('Size 1st cell: ', y_min)
    return _np.concatenate((_np.zeros((1)), xi*r1**x0 + mu ))


def geom_stretch_out(N_out, delta, percent, Nend):

    x0 = _np.linspace(1.,N_out,N_out)
    percent_max = Nend ** (1./N_out) - 1.
    print(percent_max)
    return delta * ( 1 + percent_max ) **x0


def smooth_stretch_out(N_out, delta, percent, Nend):
    ## NOT WORKING WELL

    x0 = _np.linspace(1.,N_out,N_out)
    alphai = _np.ones(N_out) + percent
    step = 0.00001
    prod = _np.prod(alphai)
    while prod < Nend:
        for k in range(N_out):
            alphai[k] = 1. + percent + k*step
        prod = _np.prod(alphai)
        print(prod)    
        step += step
    y_out = _np.ones(N_out)*delta*alphai[0]   
    for k in range(1,N_out):
        y_out[k] = y_out[k-1]*alphai[k]    
    return y_out      


############################## MAIN ###########################
if __name__ == '__main__':

    L = 6   #0.50
    l = 79.5-40 #0.035 #0.150 #0.0012  #0.0752

    im = 600
    jm = 300   #240

    x  = _np.linspace(0., L, im)
    y  = _np.linspace(0., l, jm)
    a = 1.  #1.
    b = 1.5  #3.
    c = -1.  #2.
    y1 = stretch_tanh(y, a, b, c)
    x1 = stretch_tanh(x, a, b, c)

    # d = 2.
    # e = 1.
    # f1= 1.
    # f2= 0.
    # y1 = stretch_log(y, d, e, f1, f2)

    # stretching_plot(y, y1)
    stretching_plot(x, x1)

    x0 = _np.zeros((im,jm), order='F')
    y0 = _np.zeros((im,jm), order='F')

    for i in range(im):
        for j in range(jm):
            x0[i,j] = x[i]
            y0[i,j] = y1[j]

    # plt.figure()
    # plt.plot(x0,y0, 'r-')
    #plt.plot(z, sol, 'b--')

    height = y0[0,1:] - y0[0,:-1]

    # plt.figure()
    # plt.plot(y0[0,:-1], height,'b.')
    # plt.plot([y0[0,0],y0[0,-1]], [l/jm, l/jm],'k-')
    # plt.figure()
    # plt.plot(height,'b.')
    # plt.plot([0, jm-1], [l/jm, l/jm],'k-')

    height_percent = (height[1:] - height[:-1]) / height[:-1] * 100

    # plt.figure()
    # plt.plot(y0[0,:-2], height_percent,'g.')
    # # plt.figure()
    # # plt.plot(height_percent,'g.')

    print(_np.searchsorted(y0[0,:], y0[0,-1]/9))

    Ny_in = 90*jm/100
    
    delta = 11.5
    print(delta)
    percent = 0.016
    
    y_int = bigeom_stretch_in(Ny_in, delta, percent)

    Ny_out = jm - Ny_in
    N = l/delta

    # y_out = geom_stretch_out(Ny_out, delta, percent, N)
    y_out = exp_stretch_out(Ny_out, delta, percent, N)
    # y_out = smooth_stretch_out(Ny_out, delta, percent, N)

    heightin = (y_int[1:] - y_int[:-1])
    heightin_per = (heightin[1:] - heightin[:-1])/heightin[:-1] *100
    heightout= (y_out[1:] - y_out[:-1])
    heightout_per = (heightout[1:] - heightout[:-1])/heightout[:-1] *100
    
    print(_np.shape(y_int))
    print(_np.shape(y_out))
    y_t = _np.concatenate((y_int, y_out))
    # print _np.shape(y_t)
    # print y_t

    plt.figure()
    plt.title('y')
    plt.plot(y_int,'c.')
    plt.plot(_np.linspace(Ny_in+1,Ny_in+Ny_out,Ny_out),y_out,'m.')
    # plt.plot(y_t,'-k')

    plt.figure()
    plt.title('Height cell')
    plt.plot(heightin,'c.')
    plt.plot(_np.linspace(Ny_in,Ny_in+Ny_out-1,Ny_out-1),heightout,'m.')

    plt.figure()
    plt.title('Height cell increase in percent')
    plt.plot(heightin_per,'c.')
    plt.plot(_np.linspace(Ny_in-1,Ny_in+Ny_out-3,Ny_out-2),heightout_per,'m.')


    plt.show()

