#!/usr/bin/env python

'''
File: bl2d.py

Created on 21 july 2020

@author:       Cedric Content
@contact:      cedric.content@onera.fr
@organization: ONERA - DAAA

@summary:      This file is the main file of the program. It contains the
               routine "main" and other related routines.
'''

import BROADCAST_npz_sens as toy

import numpy as _np
# FORTRAN
import srcfv.f_geom    as f_geom
import srcfv.f_bnd     as f_bnd
import srcfv.f_sch     as f_sch
import srcfv.f_lhs     as f_lhs
import srcfv.f_lin     as f_lin
# import srcfv.f_adj     as f_adj
import srcfv.f_norm    as f_norm
import misc.f_misc     as f_misc

import os
import sys
import timeit

######################### CARTE ####################
dgeom = dict()
dnum  = dict()
dphys = dict()

#GEOMETRY
im    = 300   # X discretization
jm    = 150  # Y discretization 
L     = 0.59  #0.59 # FP length #/2
high  = 0.035  #0.035 # FP high
xini  = 0.0060 #0.0024 0.0040 0.0060 # debut de la plaque


## stretching for outlet
# imst = 30  #10  #30  #30
# growthst = 1.1  #1.2  #1.05  #1.1
# dgeom['growth sp']  = growthst
# dgeom['im sp']      = imst


dgeom['im']      = im
dgeom['jm']      = jm
dgeom['length']  = L
dgeom['high']    = high
dgeom['xini']    = xini

isresolvant = False

out_dir = 'Wksp/'
treename = 'tree300x150_sens'

ite     = 10

cfl     = 1.e10 #very large CFL for Newton method
freqres = ite/2 # frequency to print residual defined by the integer, here 2
freqsort= ite/1 # frequency to write solution defined by the integer, here 1


sch   = 'dnc'   # numerical scheme 
order =   5     # formal FD order for dnc
k2 = 1.01       #1.01 
k4 = 1.         #1.

if xini < (order+1)/2*L/(im-1):
    print('Error : xini too small, go to negative x in ghost cells => xini should be higher than %s' %str((order+1)/2*L/(im-1)))
    exit()

extraporder = 2   # 2 #5   # extrapolation order for outflow

out_dir += '%s_%i' %(sch, order)
os.system('mkdir -p %s' % out_dir)

if sch == 'dnc':
    # routinesch = 'flux_num_dnc%i_2d' % order
    routinesch = 'flux_num_dnc%i_nowall_2d' % order
    # routinesch = 'flux_num_dnc%i_nowall_polar_2d' % order      

dnum['ite']      = ite
dnum['cfl']      = cfl
dnum['k2']       = k2
dnum['k4']       = k4
dnum['sch']      = sch
dnum['order']    = order
dnum['freqres']  = freqres
dnum['freqsort'] = freqsort

# Routines for BND:
routineout = 'bc_extrapolate_o%i_2d' % extraporder
routinein  = 'bc_supandsubinlet_2d'
# routinein  = 'bc_general_2d'
routinenr  = 'bc_no_reflexion_2d'
# routinew   = 'bc_wall_viscous_adia_2d'
# routinew   = 'bc_wall_viscous_iso_profile_2d'
routinew   = 'bc_wall_blow_profile_2d'

compmode = 'fixed_point'  # computational mode in ['direct', 'lin', 'adj', 'fixed_point']

dnum['lasolver'] = 'direct'    # linear algebra solver for fixed point in [direct, gmres]

dphys['gam']      = 1.4
dphys['Ts']       = 273.15  
dphys['cs']       = 110.4
dphys['musuth']   = 1.716e-5  
dphys['rgaz']     = 287.1
dphys['Prandtl']  = 0.72

dphys['Mach']     = 4.5  
dphys['T0']       = 288.  
# dphys['P0']       = 728.312  
dphys['Runit']    = 3.4e6

dphys['Lref']     = 1.
dphys['stateref'] = 'm0_runit_t0'
# StateRef = 'm0_p0_t0'


g1 = 1.0
g2 = 0.5  #  1./2.
g3 = 0.165919771368
g4 = 0.040919732041
g5 = 0.007555704391
g6 = 0.000891421261
rk6 = 1.
rk5 = g2
rk4 = g3/rk5
rk3 = g4/(rk4*rk5)
rk2 = g5/(rk3*rk4*rk5)
rk1 = g6/(rk2*rk3*rk4*rk5)
rkcoefs = [rk1, rk2, rk3, rk4, rk5, rk6]

dnum['rkcoefs'] = rkcoefs

#################### FACTORY ########################
lf    = list()
lflin = list()

if compmode == 'direct':

    libbnd = 'f_bnd'
    libsch = 'f_sch'

elif compmode == 'impli':

    libbnd = 'f_bnd'
    libsch = 'f_sch'
    liblhs = 'f_lhs'

elif compmode == 'lin':

    libbnd = 'f_lin'
    libsch = 'f_lin'

    routineout += '_d'
    routinenr  += '_d'
    routinew   += '_d'
    routinesch += '_d'

elif compmode == 'adj':

    libbnd = 'f_adj'
    libsch = 'f_adj'

    routineout += '_b'
    routinenr  += '_b'
    routinew   += '_b'
    routinesch += '_b'

elif compmode == 'fixed_point':

    libbnd = 'f_bnd'
    libsch = 'f_sch'


finflow  = eval("%s.%s"    % (libbnd, routinein))
foutflow = eval("%s.%s"    % (libbnd, routineout))
fnoref   = eval("%s.%s"    % (libbnd, routinenr ))
fwall    = eval("%s.%s"    % (libbnd, routinew  ))
fsch     = eval("%s.%s"    % (libsch, routinesch))

# lf = [finflow, foutflow, fnoref, fwall, fsch]
lf = [routinein, routineout, routinenr, routinew, routinesch, libbnd, libsch]

ffiltilo   = eval("f_bnd.bc_filteringilo_2d")
lf.append(ffiltilo)

if compmode == 'fixed_point':

    libbnd = 'f_lin'
    libsch = 'f_lin'

    routineout += '_d'
    routinein  += '_d'
    routinenr  += '_d'
    routinew   += '_d'
    routinesch += '_d'

    flinoutflow = eval("%s.%s"    % (libbnd, routineout))
    flininflow  = eval("%s.%s"    % (libbnd, routinein))
    flinnoref   = eval("%s.%s"    % (libbnd, routinenr ))
    flinwall    = eval("%s.%s"    % (libbnd, routinew  ))
    flinsch     = eval("%s.%s"    % (libsch, routinesch))

    # flininflow  = eval("f_bnd.bc_general_2d")

    # lflin = [flininflow, flinoutflow, flinnoref, flinwall, flinsch]
    lflin = [routinein, routineout, routinenr, routinew, routinesch, libbnd, libsch]


elif compmode == 'impli':
    lf.append(eval("f_lhs.impli_matrix_free_2d"))

tinit = timeit.time.time()
## Init the geometry and BC location
toy.bl2d_prepro(dgeom=dgeom, dphys=dphys, dnum=dnum , compmode=compmode, lf=lf, lflin=lflin, out_dir = out_dir, treename=treename, isresol= isresolvant)

## Solver
toy.bl2d_fromNPZtree(treename, ite, compmode = compmode, out_dir = out_dir, isresol= isresolvant)

tend = timeit.time.time()
tlaps = tend-tinit
denom = im*jm*ite*len(rkcoefs)
print("Time Elapsed           =  ", tlaps)
print("Time Elapsed/(cell*it) =  %s s  " % str(tlaps/denom))

#################### FACTORY ########################

