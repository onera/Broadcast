
## Compute several BL thicknesses

import numpy as _np
import matplotlib.pyplot as plt

import restart_init as ri


def dichotomie(inflec, a, b, tol=1):

	while abs(b - a) > tol:
		c = (a + b) // 2
		if inflec[a] * inflec[c] <= 0:
			b = c
		else:
			a = c
	return b

def moving_avg(a, n):
	## offset moving average of length n on the array a (do not average the n last downstream values, start from the end and go upstream)
	Na = _np.shape(a)[0]
	sum = 0.
	sumarr = _np.zeros(Na)
	for i in range(n):
		sum += a[Na-1-i]
		sumarr[Na-1-i] = a[Na-1-i]
	for i in range(Na-n):
		sum += a[Na-n-1-i]
		sum -= a[Na-1-i]
		sumarr[Na-n-1-i] = sum/n
	return sumarr

def moving_avg_c(a, n):
	## centered moving average of length n on the array a (do not average the n last downstream values, start from the end and go upstream)
	Na = _np.shape(a)[0]
	sum = 0.
	sumarr = _np.zeros(Na)
	for i in range(n):
		sum += a[Na-1-i]
		sumarr[Na-1-i] = a[Na-1-i]
	for i in range(Na-n):
		sum += a[Na-n-1-i]
		sum -= a[Na-1-i]
		sumarr[Na-n//2-1-i] = sum/n
	return sumarr

def computeBLquant(filet):

	X, Y, ro, rou, rov, row, roe = ri.read_init(filet)

	u = rou/ro

	Nx = _np.shape(X)[0]
	Ny = _np.shape(X)[1]

	y_inf_ind = -1

	rho_inf = ro[:,y_inf_ind]
	# rho_inf = _np.ones(Nx)*rho_inf[-1]

	u_inf   = u[:,y_inf_ind]
	# u_inf = _np.ones(Nx)*u_inf[-1]

	BL_thick      = _np.zeros(Nx)
	BL_disp_thick = _np.zeros(Nx)
	BL_mom_thick  = _np.zeros(Nx)
	H12           = _np.zeros(Nx)
	test=_np.zeros((Nx,Ny))
	inflec        = _np.zeros((Nx,Ny-2))
	inflec_point  = _np.zeros(Nx)
	tmp           = _np.zeros((Nx,Ny-2))
	inflec2       = _np.zeros((Nx,Ny-2-1))
	inflec_point2 = _np.zeros(Nx)

	
	for k in range(Nx):

		## 99% BL thickness
		h = 0
		while u[k,h] < 0.99 * u_inf[k]:
			h = h+1
		BL_thick[k] = (Y[k,h] + Y[k,h-1])/2
		# print h
		
		## Displacement thickness
		for i in range(1,h+1):
			BL_disp_thick[k] = BL_disp_thick[k] + (Y[k,i] - Y[k,i-1])*( 1. - (rou[k,i] + rou[k,i-1])/2 / (rho_inf[k]*u_inf[k]) )
			test[k,i] = (Y[k,i] - Y[k,i-1])*( 1. - (rou[k,i] + rou[k,i-1])/2 / (rho_inf[k]*u_inf[k])) 

		## Momentum thickness
		for i in range(1,h+1):
			BL_mom_thick[k] = BL_mom_thick[k] + (Y[k,i] - Y[k,i-1])*( (rou[k,i] + rou[k,i-1])/2 / (rho_inf[k]*u_inf[k]) ) * ( 1. - (u[k,i] + u[k,i-1])/2 / u_inf[k] )

		## Shape factor
		H12[k] = BL_disp_thick[k] / BL_mom_thick[k]

		## Generalized inflection point
		loc = k
		for i in range(_np.shape(inflec)[1]):
			inflec[loc,i] = (ro[loc,i+1] - ro[loc,i]) / (Y[loc,i+1] - Y[loc,i]) * (u[loc,i+1] - u[loc,i]) / (Y[loc,i+1] - Y[loc,i]) + ro[loc,i] * (u[loc,i+2] - 2.*u[loc,i+1] + u[loc,i]) / (Y[loc,i+1] - Y[loc,i])**2
			tmp[loc,i] = ro[loc,i] * (u[loc,i+1] - u[loc,i]) / (Y[loc,i+1] - Y[loc,i])
		for i in range(_np.shape(inflec)[1]-1):	
			inflec2[k,i] = (tmp[loc,i+1] - tmp[loc,i]) / (Y[loc,i+1] - Y[loc,i])
		# if loc == Nx/2:
		# 	print dichotomie(inflec[loc,:], h/2, h)
		# inflec_point[k] = Y[loc, dichotomie(inflec[loc,:], h//2, h)]
		inflec_point[k] = Y[loc, dichotomie(inflec[loc,:], h//4, h)]
		inflec_point2[k] = Y[loc, dichotomie(inflec2[loc,:], h//4, h)]


	# plt.figure(5)
	# for i in range(0,Nx,20):
	# 	plt.plot(test[i,:],Y[i,:],'--o', label=' %i' %i)
	# plt.legend()
	# plt.figure(9)
	# plt.plot(Y[Nx/2,:-2], inflec[Nx/2,:])
	# plt.plot(Y[Nx/2,:-3], inflec2[Nx/2,:],'--r')
	# plt.grid()
	# plt.figure(4)
	# plt.plot(inflec[Nx/2,:])
	# plt.grid()

	return X[:,0], BL_thick, BL_disp_thick, BL_mom_thick, H12, inflec_point, inflec_point2

def computeBLquant2(xc,yc,w):

	X   = xc
	Y   = yc 
	ro  = w[:,:,0]
	rou = w[:,:,1]
	rov = w[:,:,2]
	row = w[:,:,3]
	roe = w[:,:,4]

	u = rou/ro

	Nx = _np.shape(X)[0]
	Ny = _np.shape(X)[1]

	y_inf_ind = -1

	rho_inf = ro[:,y_inf_ind]
	# rho_inf = _np.ones(Nx)*rho_inf[-1]

	u_inf   = u[:,y_inf_ind]
	# u_inf = _np.ones(Nx)*u_inf[-1]

	BL_thick      = _np.zeros(Nx)
	BL_disp_thick = _np.zeros(Nx)
	BL_mom_thick  = _np.zeros(Nx)
	H12           = _np.zeros(Nx)
	test=_np.zeros((Nx,Ny))
	inflec        = _np.zeros((Nx,Ny-2))
	inflec_point  = _np.zeros(Nx)
	tmp           = _np.zeros((Nx,Ny-2))
	inflec2       = _np.zeros((Nx,Ny-2-1))
	inflec_point2 = _np.zeros(Nx)

	derU = _np.gradient(u, Y[0,:], axis=1)
	gip = _np.gradient(ro*derU, Y[0,:], axis=1)


	for k in range(Nx):

		## 99% BL thickness
		h = 0
		while u[k,h] < 0.99 * u_inf[k]:
			h = h+1
		BL_thick[k] = (Y[k,h] + Y[k,h-1])/2
		# print h
		
		## Displacement thickness
		for i in range(1,h+1):
			BL_disp_thick[k] = BL_disp_thick[k] + (Y[k,i] - Y[k,i-1])*( 1. - (rou[k,i] + rou[k,i-1])/2 / (rho_inf[k]*u_inf[k]) )
			test[k,i] = (Y[k,i] - Y[k,i-1])*( 1. - (rou[k,i] + rou[k,i-1])/2 / (rho_inf[k]*u_inf[k])) 

		## Momentum thickness
		for i in range(1,h+1):
			BL_mom_thick[k] = BL_mom_thick[k] + (Y[k,i] - Y[k,i-1])*( (rou[k,i] + rou[k,i-1])/2 / (rho_inf[k]*u_inf[k]) ) * ( 1. - (u[k,i] + u[k,i-1])/2 / u_inf[k] )

		## Shape factor
		H12[k] = BL_disp_thick[k] / BL_mom_thick[k]

		## Generalized inflection point
		loc = k
		for i in range(_np.shape(inflec)[1]):
			inflec[loc,i] = (ro[loc,i+1] - ro[loc,i]) / (Y[loc,i+1] - Y[loc,i]) * (u[loc,i+1] - u[loc,i]) / (Y[loc,i+1] - Y[loc,i]) + ro[loc,i] * (u[loc,i+2] - 2.*u[loc,i+1] + u[loc,i]) / (Y[loc,i+1] - Y[loc,i])**2
			tmp[loc,i] = ro[loc,i] * (u[loc,i+1] - u[loc,i]) / (Y[loc,i+1] - Y[loc,i])
			# tmp[loc,i] = (u[loc,i+1] - u[loc,i]) / (Y[loc,i+1] - Y[loc,i])
		for i in range(_np.shape(inflec)[1]-1):	
			inflec2[k,i] = (tmp[loc,i+1] - tmp[loc,i]) / (Y[loc,i+1] - Y[loc,i])
		# if loc == Nx/2:
		# 	print dichotomie(inflec[loc,:], h/2, h)
		# inflec_point[k] = Y[loc, dichotomie(inflec[loc,:], h//2, h)]
		inflec_point[k] = Y[loc, dichotomie(inflec[loc,:], h//4, h)]
		# inflec_point2[k] = Y[loc, dichotomie(inflec2[loc,:], h//4, h)]
		inflec_point2[k] = Y[loc, dichotomie(gip[loc,:], h//4, h)]


	# plt.figure(5)
	# for i in range(0,Nx,20):
	# 	plt.plot(test[i,:],Y[i,:],'--o', label=' %i' %i)
	# plt.legend()
	# plt.figure(9)
	# plt.plot(Y[Nx/2,:-2], inflec[Nx/2,:])
	# plt.plot(Y[Nx/2,:-3], inflec2[Nx/2,:],'--r')
	# plt.grid()
	# plt.figure(4)
	# plt.plot(inflec[Nx/2,:])
	# plt.grid()

	return X[:,0], BL_thick, BL_disp_thick, BL_mom_thick, H12, inflec_point, inflec_point2

if __name__ == '__main__':


	dir   = 'Wksp'

	file  = 'hllc_3/state_atcenter_ite12'
	# file  = 'hllc_3/initialisation'
	file  = 'hllc_5/state_atcenter_mesh1500'
	file  = 'hllc_5/state_atcenter_mesh2500_y120_xini16e-4_longer'

	filet = './' + dir + '/' + file + '.dat'

	X, Y, ro, rou, rov, row, roe = ri.read_init(filet)

	Xplot, BL_thick, BL_disp_thick, BL_mom_thick, H12, inflec_point, inflec_point2 = computeBLquant(filet)


	print(' -------------------- ')
	print('Inlet delta 99 = ', BL_thick[0])
	print('Inlet delta * = ', BL_disp_thick[0])
	print('Inlet theta = ', BL_mom_thick[0])
	print('Inlet H12 = ', H12[0])
	print(' -------------------- ')
	print('Outlet delta 99 = ', BL_thick[-1])
	print('Outlet delta * = ', BL_disp_thick[-1])
	print('Outlet theta = ', BL_mom_thick[-1])
	print('Outlet H12 = ', H12[-1])


	Zoom = 4

	plt.figure(1)
	plt.title('u')
	plt.contourf(X,Y,rou/ro,11)
	# plt.contourf(X,Y,rou/ro,levels=list(_np.linspace(0.95,1.05,11)))
	plt.xlabel('x')
	plt.ylabel('y')
	plt.colorbar()
	plt.ylim(0., _np.amax(Y)/Zoom)
	plt.plot(Xplot, BL_thick,'k-', label='delta 99')
	plt.plot(Xplot, BL_disp_thick,'k--', label='delta*')
	plt.plot(Xplot, BL_mom_thick,'k.', label='theta')
	plt.plot(Xplot, inflec_point,'k-.', label='inflec. point')
	plt.legend()

	# plt.figure(2)
	# plt.title('H12')
	# plt.xlabel('x')
	# plt.plot(X[:,0], H12)
	# plt.grid()
	# # plt.ylim(0., 16.)

	plt.show()







