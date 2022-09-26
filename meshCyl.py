# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import numpy as _np
import math as m
import matplotlib.pyplot as plt
import meshBL as mtools

def cylinder2d(ri=1.,rf=10.,c1=0.,c2=0.,imax=10,jmax=10):
    '''
    mesh generation of 2D cylinder

    '''
    im = imax
    jm = jmax

    x  = _np.zeros((im,jm), order = 'F')
    y  = _np.zeros((im,jm), order = 'F')
    t  = _np.linspace(0.,360.,im)
    r0 = _np.linspace(ri,rf,jm)
    a  = 1.
    b  = 3.     #3. #5.
    c  = 1.
    # xw  = mesh.stretch_tanh(x1, a, b, c)
    r = mtools.stretch_tanh(r0, a, b, c)



    for j in range(jm):
        for i in range(im):
            x[i,j] = c1+r[j]*_np.cos(t[i]*_np.pi/180.)
            y[i,j] = c2+r[j]*_np.sin(t[i]*_np.pi/180.)

    # __writemesh('totocyl.dat', im, jm, x, y)
    
    return x, y

def cylinder2d_bis(r=_np.zeros(10),c1=0.,c2=0.,imax=10,jmax=10):
    '''
    mesh generation of 2D cylinder with the array r given

    '''
    im = imax
    jm = jmax

    x  = _np.zeros((im,jm), order = 'F')
    y  = _np.zeros((im,jm), order = 'F')
    # t  = _np.linspace(0.,360.,im)  #Wrong dir
    # t  = _np.linspace(0.,180.,im)  #Half circle, wrong dir
    t  = _np.flipud(_np.linspace(0.,180.,im)) #Half circle top
    # t  = _np.flipud(_np.linspace(-180.,180.,im))
    # t  = _np.linspace(0.,-180.,im) #Half circle bottom

    for j in range(jm):
        for i in range(im):
            x[i,j] = c1+r[j]*_np.cos(t[i]*_np.pi/180.)
            y[i,j] = c2+r[j]*_np.sin(t[i]*_np.pi/180.)

    # __writemesh('totocyl.dat', im, jm, x, y)
    
    return x, y

def cylinder2d_bisfull(r=_np.zeros(10),c1=0.,c2=0.,imax=10,jmax=10):
    '''
    mesh generation of 2D cylinder with the array r given

    '''
    im = imax
    jm = jmax

    x  = _np.zeros((im,jm), order = 'F')
    y  = _np.zeros((im,jm), order = 'F')
    # t  = _np.flipud(_np.linspace(0.,180.,im)) #Half circle
    t  = _np.flipud(_np.linspace(-180.,180.,im))  #Connectivity upstream
    # t  = _np.flipud(_np.linspace(0.,360.,im))  #Connectivity downstream

    ## MESH 2
    # a  = 1.
    # b  = 4.     #3. #5.
    # c  = 0.5
    # t = mtools.stretch_tanh(t, a, b, c)


    ## MESH 3 - SYM
    # # t  = _np.flipud(_np.linspace(180.,360.,im/2+1))
    # t  = _np.flipud(_np.linspace(0.,180.,im/2+1))
    # a  = 1.
    # b  = 4.     #3. #5.
    # c  = 0.5
    # t1 = mtools.stretch_tanh(t, a, b, c)
    # # t  = _np.flipud(_np.linspace(0.,180.,im/2+1))
    # t  = _np.flipud(_np.linspace(180.,360.,im/2+1))
    # t2 = mtools.stretch_tanh(t, a, b, c)
    # t  = _np.concatenate((t1, t2[1:]))

    ## MESH 5 to 8 - ANTISYM - join upstream
    t  = _np.flipud(_np.linspace(0.,180.,im/2+1))
    a  = 1.
    b  = 1.5     #3. #5.
    c  = -1.
    t1 = mtools.stretch_tanh(t, a, b, c)
    t  = _np.linspace(180.,360.,im/2+1)
    t2 = mtools.stretch_tanh(t, a, b, c)
    t2 = _np.flipud(t2)
    t  = _np.concatenate((t1, t2[1:]))



    for j in range(jm):
        for i in range(im):
            x[i,j] = c1+r[j]*_np.cos(t[i]*_np.pi/180.)
            y[i,j] = c2+r[j]*_np.sin(t[i]*_np.pi/180.)

    # __writemesh('totocyl.dat', im, jm, x, y)
    
    return x, y    

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

if __name__ == '__main__':
    im = 100
    jm = 100
    x,y = cylinder2d(ri=0.1,rf=10.,c1=0.,c2=0.,imax=im,jmax=jm)
    __writemesh('totocyl.dat', im, jm, x, y)

