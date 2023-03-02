
## Initialise the state from a previous solution to restart

def read_init(file):
	''' Read from a previous solution file (only valid for cartesian rectangular mesh '''

	import numpy as _np

	BLprof = _np.loadtxt(file,comments=('#','ZONE'),skiprows=3 )

	X   = BLprof[:,0]
	Y   = BLprof[:,1]
	ro  = BLprof[:,2]
	rou = BLprof[:,3]
	rov = BLprof[:,4]
	row = BLprof[:,5]
	roe = BLprof[:,6]

	Ntot= _np.shape(X)[0]  
	# print Ntot
    
	count=1
	Y0=Y[0]
	while Y[count] == Y0:
		count+= 1
	im  = count
	jm  = Ntot//im

	# print('im = ', im)
	# print('jm = ', jm)

	X   = _np.reshape(X  ,(im,jm), order='F')
	Y   = _np.reshape(Y  ,(im,jm), order='F')
	ro  = _np.reshape(ro ,(im,jm), order='F')
	rou = _np.reshape(rou,(im,jm), order='F')
	rov = _np.reshape(rov,(im,jm), order='F')
	row = _np.reshape(row,(im,jm), order='F')
	roe = _np.reshape(roe,(im,jm), order='F')

	return X, Y, ro, rou, rov, row, roe


def interpolate_field_from_clicet(file ,x ,y ,R_pg ,gamma ,interptype='cubic'):
    '''Interpolate the conservative variables (except rhow) from 2D BL data from clicet file on the grid made by x and y'''

    from scipy.interpolate import griddata
    import numpy as _np

    BLprof = _np.loadtxt(file,comments=('#','ZONE'),skiprows=3 )

    X   = BLprof[:,0]
    Y   = BLprof[:,1]
    U   = BLprof[:,2]
    V   = BLprof[:,3]

    rho = BLprof[:,5]
    T   = BLprof[:,6]

    gridX, gridY = _np.meshgrid(x,y)

    Uinter   = griddata(_np.transpose([X,Y]), U, (gridX,gridY), interptype)
    Vinter   = griddata(_np.transpose([X,Y]), V, (gridX,gridY), interptype)
    rhointer = griddata(_np.transpose([X,Y]), rho, (gridX,gridY), interptype)
    Tinter   = griddata(_np.transpose([X,Y]), T, (gridX,gridY), interptype)

    return rhointer, rhointer*Uinter, rhointer*Vinter, rhointer*(R_pg/(gamma-1.)*Tinter+0.5*(Uinter**2+Vinter**2))


def interpolate_field_from_data(file ,x ,y ,interptype='cubic'):
    '''Interpolate the conservative variables (except rhow) from 2D BL data from previous solution file on the grid made by x and y'''

    from scipy.interpolate import griddata
    import numpy as _np
    
    BLprof = _np.loadtxt(file,comments=('#','ZONE'),skiprows=3 )

    X   = BLprof[:,0]
    Y   = BLprof[:,1]
    ro  = BLprof[:,2]
    rou = BLprof[:,3]
    rov = BLprof[:,4]
    row = BLprof[:,5]
    roe = BLprof[:,6]


    gridX, gridY = _np.meshgrid(x,y, sparse=True)

    rointer  = griddata(_np.transpose([X,Y]), ro , (gridX,gridY), interptype)
    roUinter = griddata(_np.transpose([X,Y]), rou, (gridX,gridY), interptype)
    roVinter = griddata(_np.transpose([X,Y]), rov, (gridX,gridY), interptype)
    roWinter = griddata(_np.transpose([X,Y]), row, (gridX,gridY), interptype)
    roEinter = griddata(_np.transpose([X,Y]), roe, (gridX,gridY), interptype)

    print(_np.shape(gridX))

    return gridX, gridY, rointer, roUinter, roVinter, roWinter, roEinter    


