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
import srcfv.f_geom    as f_geom
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

import numpy as _np
import matplotlib.pyplot as plt

import srcfv.f_dz     as f_dz

import os
import sys
import timeit

import Converter.Internal as I
import Converter.PyTree as C
# import CGNS.MAP as cgm

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

def writeCGNS(filename, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, pinf=None):
    ''' Write the data in a CGNS tree '''    

    I.__FlowSolutionCenters__ = 'FlowSolution#Centers'
    w0 = _np.array(w, copy=True)
    res0 = _np.array(res, copy=True)

    t=I.newCGNSTree()
    b=I.newCGNSBase(name='Base',parent=t)
    # z=I.newZone(name="Zone1",parent=b,zsize=[[(im+1)*(jm+1),im*jm, 0]])
    z=I.newZone(name="Zone1",parent=b,zsize=[[(im+1)*(jm+1),im*jm,0]])
    Gc=I.newGridCoordinates(parent=z)
    I.newDataArray(name="CoordinateX",value=x0[:,:],parent=Gc)
    I.newDataArray(name="CoordinateY",value=y0[:,:],parent=Gc)
    # I.newDataArray(name="CoordinateZ",value=y0*0,parent=Gc)

    Geom = I.newDiscreteData(name="Geometry", parent=z)
    I.newDataArray(name="Volume",value=vol,parent=Geom)
    I.newDataArray(name="VolumeInv",value=volf,parent=Geom)
    I.newDataArray(name="NormalX",value=nx,parent=Geom)
    I.newDataArray(name="NormalY",value=ny,parent=Geom)
    I.newDataArray(name="DimX",value=im,parent=Geom)
    I.newDataArray(name="DimY",value=jm,parent=Geom)
    I.newDataArray(name="GhostCell",value=gh,parent=Geom)
    I.newDataArray(name="CoordinateCenterX",value=xc[:,:],parent=Geom)
    I.newDataArray(name="CoordinateCenterY",value=yc[:,:],parent=Geom)

    # Fs=I.newFlowSolution(name="FlowSolution#Init", gridLocation='CellCenter', parent=z)
    # # I.newGridLocation(value='CellCenter', parent=Fs)
    # I.newDataArray(name="Density",value=w0[gh:-gh,gh:-gh,0],parent=Fs)
    # I.newDataArray(name="MomentumX",value=w0[gh:-gh,gh:-gh,1],parent=Fs)
    # I.newDataArray(name="MomentumY",value=w0[gh:-gh,gh:-gh,2],parent=Fs)
    # I.newDataArray(name="MomentumZ",value=w0[gh:-gh,gh:-gh,3],parent=Fs)
    # I.newDataArray(name="EnergyStagnationDensity",value=w0[gh:-gh,gh:-gh,4],parent=Fs)
    # Res = I.newDiscreteData(name="Residual#Init", parent=z)
    # I.newGridLocation(value='CellCenter', parent=Res)
    # I.newDataArray(name="ResidualDensity",value=res0[gh:-gh,gh:-gh,0],parent=Res)
    # I.newDataArray(name="ResidualMomentumX",value=res0[gh:-gh,gh:-gh,1],parent=Res)
    # I.newDataArray(name="ResidualMomentumY",value=res0[gh:-gh,gh:-gh,2],parent=Res)
    # I.newDataArray(name="ResidualMomentumZ",value=res0[gh:-gh,gh:-gh,3],parent=Res)
    # I.newDataArray(name="ResidualEnergyStagnationDensity",value=res0[gh:-gh,gh:-gh,4],parent=Res)

    ## With Ghost cells
    Fs=I.newFlowSolution(name="FlowSolution#Init", gridLocation='CellCenter', parent=z)
    # I.newGridLocation(value='CellCenter', parent=Fs)
    I.newDataArray(name="Density",value=w0[:,:,0],parent=Fs)
    I.newDataArray(name="MomentumX",value=w0[:,:,1],parent=Fs)
    I.newDataArray(name="MomentumY",value=w0[:,:,2],parent=Fs)
    I.newDataArray(name="MomentumZ",value=w0[:,:,3],parent=Fs)
    I.newDataArray(name="EnergyStagnationDensity",value=w0[:,:,4],parent=Fs)
    Res = I.newDiscreteData(name="Residual#Init", parent=z)
    I.newGridLocation(value='CellCenter', parent=Res)
    I.newDataArray(name="ResidualDensity",value=res0[:,:,0],parent=Res)
    I.newDataArray(name="ResidualMomentumX",value=res0[:,:,1],parent=Res)
    I.newDataArray(name="ResidualMomentumY",value=res0[:,:,2],parent=Res)
    I.newDataArray(name="ResidualMomentumZ",value=res0[:,:,3],parent=Res)
    I.newDataArray(name="ResidualEnergyStagnationDensity",value=res0[:,:,4],parent=Res)

    routinein  = lf[0]
    routineout = lf[1]
    routinenr  = lf[2]
    routinew   = lf[3]
    routinesch = lf[4]
    libbnd     = lf[5]
    libsch     = lf[6]

    BC = I.newZoneBC(parent=z)
    BCin = I.newBC(name="Inflow",parent=BC, btype="FamilySpecified", family=libbnd + "." + routinein)
    I.newGridLocation(value='CellCenter', parent=BCin)
    I.newPointRange(name="PointRange",value=[interf1[0,0], interf1[1,0], interf1[0,1], interf1[1,1]],parent=BCin)
    I.newDataArray(name="Field",value=field,parent=BCin)
    BCout = I.newBC(name="Outflow",parent=BC, btype="FamilySpecified", family=libbnd + "." + routineout)
    I.newGridLocation(value='CellCenter', parent=BCout)
    I.newPointRange(name="PointRange",value=[interf2[0,0], interf2[1,0], interf2[0,1], interf2[1,1]],parent=BCout)

    I.newDataArray(name="Pref",value=pinf,parent=BCout)

    BCwall = I.newBC(name="Wall",parent=BC, btype="FamilySpecified", family=libbnd + "." + routinew)
    I.newGridLocation(value='CellCenter', parent=BCwall)
    I.newPointRange(name="PointRange",value=[interf3[0,0], interf3[1,0], interf3[0,1], interf3[1,1]],parent=BCwall)
    BCtop = I.newBC(name="NoRef",parent=BC, btype="FamilySpecified", family=libbnd + "." + routinenr)
    I.newGridLocation(value='CellCenter', parent=BCtop)
    I.newPointRange(name="PointRange",value=[interf4[0,0], interf4[1,0], interf4[0,1], interf4[1,1]],parent=BCtop)
    I.newDataArray(name="Wbd",value=wbd,parent=BCtop)

    Num = I.newIntegralData(name="NumericalScheme", parent=b)
    I.newDataArray(name="SchemeType",value=sch,parent=Num)
    I.newDataArray(name="Scheme",value=libsch + "." + routinesch,parent=Num)
    # I.newFamily(name=libsch + "." + routinesch,parent=Num)
    I.newDataArray(name="k2",value=k2,parent=Num)
    I.newDataArray(name="k4",value=k4,parent=Num)

    C._addState(b, 'Mach', mach)
    C._addState(b, 'Prandtl', prandtl)
    C._addState(b, 'Rgaz', rgaz)
    C._addState(b, 'Cp', cp)
    C._addState(b, 'Cv', cv)
    C._addState(b, 'Gamma', gam)
    C._addState(b, 'TemperatureSutherland', cs)
    C._addState(b, 'TemperatureRefSutherland', tref)
    C._addState(b, 'ViscosityRefSutherland', muref)

    C.convertPyTree2File(t, filename+'.hdf')

    return t

def writeNPZ(filename, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, cfl, freqres, freqsort, pinf=None):
    ''' Write the data in a NPZ file'''    

    w0 = _np.array(w, copy=True)
    res0 = _np.array(res, copy=True)

    routinein  = lf[0]
    routineout = lf[1]
    routinenr  = lf[2]
    routinew   = lf[3]
    routinesch = lf[4]
    libbnd     = lf[5]
    libsch     = lf[6]

    _np.savez(filename+'.npz',x0=x0,y0=y0,xc=xc,yc=yc,vol=vol,volf=volf,nx=nx,ny=ny,im=im,jm=jm,gh=gh,FlowSolutionInit=w0,ResidualInit=res0,interf1=interf1,interf2=interf2,interf3=interf3,interf4=interf4,BC1=libbnd + "." + routinein,BC2=libbnd + "." + routineout,BC3=libbnd + "." + routinew,BC4=libbnd + "." + routinenr,Wbd=wbd,Field=field,SchemeType=sch,Scheme=libsch + "." + routinesch,k2=k2,k4=k4,Mach=mach,Prandtl=prandtl,Rgaz=rgaz,Cp=cp,Cv=cv,Gamma=gam,TemperatureSutherland=cs,TemperatureRefSutherland=tref,ViscosityRefSutherland=muref,Pref=pinf,CFL=cfl,Freqres=freqres,Freqsort=freqsort)
    # print(_np.load(filename+'.npz')['BC4'])


def fillCGNS(filename, t, w, res, IA, JA, Jacvol, gh):
    ''' fill a given CGNS tree with current solution, residual and the lists (IA,JA,Aij) to build the Jacobian'''

    # I.__FlowSolutionCenters__ = 'FlowSolution#EndOfRun'

    z=I.getZones(t)[0]

    # Fs2=I.newFlowSolution(name="FlowSolution#EndOfRun", gridLocation='CellCenter', parent=z)
    # # I.newGridLocation(value='CellCenter', parent=Fs2)
    # I.newDataArray(name="Density",value=w[gh:-gh,gh:-gh,0],parent=Fs2)
    # I.newDataArray(name="MomentumX",value=w[gh:-gh,gh:-gh,1],parent=Fs2)
    # I.newDataArray(name="MomentumY",value=w[gh:-gh,gh:-gh,2],parent=Fs2)
    # I.newDataArray(name="MomentumZ",value=w[gh:-gh,gh:-gh,3],parent=Fs2)
    # I.newDataArray(name="EnergyStagnationDensity",value=w[gh:-gh,gh:-gh,4],parent=Fs2)
    # Res2 = I.newDiscreteData(name="Residual#EndOfRun", parent=z)
    # I.newGridLocation(value='CellCenter', parent=Res2)
    # I.newDataArray(name="ResidualDensity",value=res[gh:-gh,gh:-gh,0],parent=Res2)
    # I.newDataArray(name="ResidualMomentumX",value=res[gh:-gh,gh:-gh,1],parent=Res2)
    # I.newDataArray(name="ResidualMomentumY",value=res[gh:-gh,gh:-gh,2],parent=Res2)
    # I.newDataArray(name="ResidualMomentumZ",value=res[gh:-gh,gh:-gh,3],parent=Res2)
    # I.newDataArray(name="ResidualEnergyStagnationDensity",value=res[gh:-gh,gh:-gh,4],parent=Res2)

    ## With Ghost cells
    Fs2=I.newFlowSolution(name="FlowSolution#EndOfRun", gridLocation='CellCenter', parent=z)
    # I.newGridLocation(value='CellCenter', parent=Fs2)
    I.newDataArray(name="Density",value=w[:,:,0],parent=Fs2)
    I.newDataArray(name="MomentumX",value=w[:,:,1],parent=Fs2)
    I.newDataArray(name="MomentumY",value=w[:,:,2],parent=Fs2)
    I.newDataArray(name="MomentumZ",value=w[:,:,3],parent=Fs2)
    I.newDataArray(name="EnergyStagnationDensity",value=w[:,:,4],parent=Fs2)
    Res2 = I.newDiscreteData(name="Residual#EndOfRun", parent=z)
    I.newGridLocation(value='CellCenter', parent=Res2)
    I.newDataArray(name="ResidualDensity",value=res[:,:,0],parent=Res2)
    I.newDataArray(name="ResidualMomentumX",value=res[:,:,1],parent=Res2)
    I.newDataArray(name="ResidualMomentumY",value=res[:,:,2],parent=Res2)
    I.newDataArray(name="ResidualMomentumZ",value=res[:,:,3],parent=Res2)
    I.newDataArray(name="ResidualEnergyStagnationDensity",value=res[:,:,4],parent=Res2)

    Jac = I.newDiscreteData(name="JacobianListsCSR", parent=z)
    I.newGridLocation(value='CellCenter', parent=Jac)
    I.newDataArray(name="IA",value=IA,parent=Jac)
    I.newDataArray(name="JA",value=JA,parent=Jac)
    I.newDataArray(name="Aij",value=Jacvol,parent=Jac)

    C.convertPyTree2File(t, filename+'.hdf')

def fillNPZ(filename, w, res, IA, JA, Jacvol, gh):
    ''' fill a given NPZ file with current solution, residual and the lists (IA,JA,Aij) to build the Jacobian'''

    dic = dict(_np.load(filename+'.npz'))
    dic['ResidualEndOfRun'] = res
    dic['FlowSolutionEndOfRun'] = w
    dic['IA'] = IA
    dic['JA'] = JA
    dic['Aij'] = Jacvol
    _np.savez(filename+'.npz', **dic)


def fillNPZ_3D(filename, IAdz, JAdz, Jacdz, IAdz2, JAdz2, Jacdz2):
    ''' fill a given NPZ file with the lists (IA,JA,Aij) for Dz and Dz2 to build the 3D Jacobian'''

    dic = dict(_np.load(filename+'.npz'))
    dic['IAdz'] = IAdz
    dic['JAdz'] = JAdz
    dic['Aijdz'] = Jacdz
    dic['IAdz2'] = IAdz2
    dic['JAdz2'] = JAdz2
    dic['Aijdz2'] = Jacdz2
    _np.savez(filename+'.npz', **dic)

def readCGNStree(filename):
    ''' read a CGNS tree made to run BROADCAST code '''

    # T = cgm.load(filename+'.hdf')
    # treeBroadcast, cgnslinks, cgnspaths = T[0], T[1], T[2]
    treeBroadcast = C.convertFile2PyTree(filename+'.hdf')
    
    vol = I.getValue(I.getNodeFromName(treeBroadcast, "Volume"))
    volf = I.getValue(I.getNodeFromName(treeBroadcast, "VolumeInv"))
    nx = I.getValue(I.getNodeFromName(treeBroadcast, "NormalX"))
    ny = I.getValue(I.getNodeFromName(treeBroadcast, "NormalY"))
    im = I.getValue(I.getNodeFromName(treeBroadcast, "DimX"))
    jm = I.getValue(I.getNodeFromName(treeBroadcast, "DimY"))
    gh = I.getValue(I.getNodeFromName(treeBroadcast, "GhostCell"))
    x0 = I.getValue(I.getNodeFromName(treeBroadcast, "CoordinateX"))
    y0 = I.getValue(I.getNodeFromName(treeBroadcast, "CoordinateY"))
    xc = I.getValue(I.getNodeFromName(treeBroadcast, "CoordinateCenterX"))
    yc = I.getValue(I.getNodeFromName(treeBroadcast, "CoordinateCenterY"))
    wbd = I.getValue(I.getNodeFromName(treeBroadcast, "Wbd"))
    field = I.getValue(I.getNodeFromName(treeBroadcast, "Field"))

    pinf = I.getValue(I.getNodeFromName(treeBroadcast, "Pref"))

    sch = I.getValue(I.getNodeFromName(treeBroadcast, "SchemeType"))
    k2 = I.getValue(I.getNodeFromName(treeBroadcast, "k2"))
    k4 = I.getValue(I.getNodeFromName(treeBroadcast, "k4"))
    cp = I.getValue(I.getNodeFromName(treeBroadcast, "Cp"))
    cv = I.getValue(I.getNodeFromName(treeBroadcast, "Cv"))
    gam = I.getValue(I.getNodeFromName(treeBroadcast, "Gamma"))
    rgaz = I.getValue(I.getNodeFromName(treeBroadcast, "Rgaz"))
    cs = I.getValue(I.getNodeFromName(treeBroadcast, "TemperatureSutherland"))
    tref = I.getValue(I.getNodeFromName(treeBroadcast, "TemperatureRefSutherland"))
    muref = I.getValue(I.getNodeFromName(treeBroadcast, "ViscosityRefSutherland"))
    mach = I.getValue(I.getNodeFromName(treeBroadcast, "Mach"))
    prandtl = I.getValue(I.getNodeFromName(treeBroadcast, "Prandtl"))
    pathBCInflow = I.getPathsFromName(treeBroadcast,"Inflow", pyCGNSLike=True)[0]
    interf1 = _np.reshape(I.getValue(I.getNodeFromPath(treeBroadcast, pathBCInflow + "/PointRange")), (2,2), order='F')
    pathBCOutflow = I.getPathsFromName(treeBroadcast,"Outflow", pyCGNSLike=True)[0]
    interf2 = _np.reshape(I.getValue(I.getNodeFromPath(treeBroadcast, pathBCOutflow + "/PointRange")), (2,2), order='F')
    pathBCWall = I.getPathsFromName(treeBroadcast,"Wall", pyCGNSLike=True)[0]
    interf3 = _np.reshape(I.getValue(I.getNodeFromPath(treeBroadcast, pathBCWall + "/PointRange")), (2,2), order='F')
    pathBCNoRef = I.getPathsFromName(treeBroadcast,"NoRef", pyCGNSLike=True)[0]
    interf4 = _np.reshape(I.getValue(I.getNodeFromPath(treeBroadcast, pathBCNoRef + "/PointRange")), (2,2), order='F')
    libsch, routinesch = I.getValue(I.getNodeFromName(treeBroadcast, "Scheme")).split('.')
    libbnd, routinein  = I.getValue(I.getNodeFromPath(treeBroadcast, pathBCInflow + "/FamilyName")).split('.')
    libbnd, routineout = I.getValue(I.getNodeFromPath(treeBroadcast, pathBCOutflow + "/FamilyName")).split('.')
    libbnd, routinew   = I.getValue(I.getNodeFromPath(treeBroadcast, pathBCWall + "/FamilyName")).split('.')
    libbnd, routinenr  = I.getValue(I.getNodeFromPath(treeBroadcast, pathBCNoRef + "/FamilyName")).split('.')
    w   = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
    res = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
    if I.getPathsFromName(treeBroadcast,"FlowSolution#EndOfRun", pyCGNSLike=True) == []:
        pathw = I.getPathsFromName(treeBroadcast,"FlowSolution#Init", pyCGNSLike=True)[0]
        pathres = I.getPathsFromName(treeBroadcast,"Residual#Init", pyCGNSLike=True)[0]
    else:
        pathw = I.getPathsFromName(treeBroadcast,"FlowSolution#EndOfRun", pyCGNSLike=True)[0]
        pathres = I.getPathsFromName(treeBroadcast,"Residual#EndOfRun", pyCGNSLike=True)[0]
    w[:,:,0] = I.getValue(I.getNodeFromPath(treeBroadcast, pathw + "/Density"))
    w[:,:,1] = I.getValue(I.getNodeFromPath(treeBroadcast, pathw + "/MomentumX"))
    w[:,:,2] = I.getValue(I.getNodeFromPath(treeBroadcast, pathw + "/MomentumY"))
    w[:,:,3] = I.getValue(I.getNodeFromPath(treeBroadcast, pathw + "/MomentumZ"))
    w[:,:,4] = I.getValue(I.getNodeFromPath(treeBroadcast, pathw + "/EnergyStagnationDensity"))
    res[:,:,0] = I.getValue(I.getNodeFromPath(treeBroadcast, pathres + "/ResidualDensity"))
    res[:,:,1] = I.getValue(I.getNodeFromPath(treeBroadcast, pathres + "/ResidualMomentumX"))
    res[:,:,2] = I.getValue(I.getNodeFromPath(treeBroadcast, pathres + "/ResidualMomentumY"))
    res[:,:,3] = I.getValue(I.getNodeFromPath(treeBroadcast, pathres + "/ResidualMomentumZ"))
    res[:,:,4] = I.getValue(I.getNodeFromPath(treeBroadcast, pathres + "/ResidualEnergyStagnationDensity"))

    lf = [routinein, routineout, routinenr, routinew, routinesch, libbnd, libsch]
    
    # return treeBroadcast, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl
    return treeBroadcast, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, pinf

def readNPZtree(filename):
    ''' read a NPZ file made to run BROADCAST code '''

    dic = _np.load(filename+'.npz')
    im = dic['im']
    jm = dic['jm']
    gh = dic['gh']
    if 'FlowSolutionEndOfRun' in dic:
        w = dic['FlowSolutionEndOfRun']
        res = dic['ResidualEndOfRun']
        print('Solution read at EndOfRun')
    else:
        w = dic['FlowSolutionInit']
        res = dic['ResidualInit']
        print('Solution read at Init')
    x0 = dic['x0']
    y0 = dic['y0']
    xc = dic['xc']
    yc = dic['yc']
    vol = dic['vol']
    volf = dic['volf']
    nx = dic['nx']
    ny = dic['ny']
    field = dic['Field']
    wbd = dic['Wbd']
    sch = _np.array2string(dic['SchemeType'])[1:-1]
    k2 = dic['k2']
    k4 = dic['k4']
    cp = dic['Cp']
    cv = dic['Cv']
    gam = dic['Gamma']
    rgaz = dic['Rgaz']
    mach = dic['Mach']
    prandtl = dic['Prandtl']
    cs = dic['TemperatureSutherland']
    tref = dic['TemperatureRefSutherland']
    muref = dic['ViscosityRefSutherland']
    interf1 = dic['interf1']
    interf2 = dic['interf2']
    interf3 = dic['interf3']
    interf4 = dic['interf4']
    BC1 = dic['BC1']
    BC2 = dic['BC2']
    BC3 = dic['BC3']
    BC4 = dic['BC4']
    scheme = dic['Scheme']

    pinf = dic['Pref']
    cfl = dic['CFL']
    freqres = dic['Freqres']
    freqsort = dic['Freqsort']

    libsch, routinesch = _np.array2string(scheme)[1:-1].split('.')
    libbnd, routinein = _np.array2string(BC1)[1:-1].split('.')
    libbnd, routineout = _np.array2string(BC2)[1:-1].split('.')
    libbnd, routinew = _np.array2string(BC3)[1:-1].split('.')
    libbnd, routinenr = _np.array2string(BC4)[1:-1].split('.')

    lf = [routinein, routineout, routinenr, routinew, routinesch, libbnd, libsch]

    return im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, pinf, cfl, freqres, freqsort


######################### Private functions ####################

# solve monoblock Boundary Layer

def bl2d_prepro(dgeom = dict(), dphys = dict(), dnum = dict(), compmode = 'direct', lf = list(), lflin = list(), out_dir = 'totodir', treename='tree', isresol= False):
    '''
    exemple of monoblock use of 2DTOY
    to simulate 2D laminar Boundary Layer flow
    '''
    os.system('mkdir -p %s' % out_dir)

    # Create mesh
    im       = dgeom['im']
    jm       = dgeom['jm']
    L        = dgeom['length']
    high     = dgeom['high']
    xini     = dgeom['xini']

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


    ## MESH in x-direction
    x  = _np.linspace(xini, xini+L , im+1)
    ## With stretching
    # growthst = dgeom['growth sp']
    # imst = dgeom['im sp']
    # xst = mesh.stretch(imst, growthst, xini+L, L/im)
    # x  = _np.concatenate((x, xst)) 
    # im = im + imst

    ## MESH in y-direction
    Ny_in   = 80*jm//100 #number of points inside the BL  
    deltaBL = high/4     #height of the BL
    percent = 0.02       #growth factor increase inside the BL
    
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

    # Compute Geometry
    for i in range(im+1):
        x0[i+gh,:] = x[i]
    for j in range(jm+1):
        y0[:,j+gh] = y[j]
    # OR import your own mesh
    # x0 = ...

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

    ## Adim with ref length
    # Lref   = 8.e-2  
    # Muref  = Roref*Vref*Lref
    ## OR Adim with unit Reynolds
    Muref  = muinf
    Lref   = Muref/(Roref*Vref)
    ## OR no normalisation
    # Roref = 1.
    # Vref  = 1.
    # Tref  = 1.
    # Lref  = 1.
    # Muref = 1.

    Pref  = Roref*Vref**2
    Cvref = Vref**2/Tref
    Eref   = Vref**2
    Rgpref = Cvref

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

    # Adim Geom:
    x0 *= 1./Lref
    y0 *= 1./Lref

    f_geom.computegeom_2d(x0,y0,nx,ny,xc,yc,vol,volf,im,jm,gh)

    #interfaces definitions (may be done at the begining)
    # Ilo
    interf1      = _np.zeros((2,2), order='F')
    interf1[0,0] = 1  # imin
    interf1[0,1] = 1  # jmin
    interf1[1,0] = 1  # imax
    interf1[1,1] = jm # jmax 

    # Ihi
    interf2      = _np.zeros((2,2), order='F')
    interf2[0,0] = im # imin
    interf2[0,1] = 1  # jmin
    interf2[1,0] = im # imax
    interf2[1,1] = jm+gh # jmax  

    # Jlo
    interf3      = _np.zeros((2,2), order='F')
    interf3[0,0] = 1-gh # imin 
    interf3[0,1] = 1  # jmin
    interf3[1,0] = im+gh # imax
    interf3[1,1] = 1  # jmax

    # Jhi
    interf4      =  _np.zeros((2,2), order='F')
    interf4[0,0] = 1-gh  # imin
    interf4[0,1] = jm # jmin
    interf4[1,0] = im # imax
    interf4[1,1]= jm # jmax

    # Blasius for inlet
    
    field = _np.zeros((jm, gh, 5), order = 'F') # profile for inlet, different values inside the ghost cells
    wbd   = _np.zeros((im+gh , 5), order = 'F') # profile for non-reflection top BC, value at the first ghost cell only

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

    # Initialise from self-similar solution
    ## Compressible self-similar profile
    road,uad,vad,Ead = blsim.BLprofile(x0[:,:]*Lref, y0[:,gh:]*Lref,mach, dphys, isplot=False, damped=False)

    road = centers_array(road)
    uad  = centers_array(uad)
    vad  = centers_array(vad)
    Ead  = centers_array(Ead)

    w[:, gh:, 0]     = road[:,:]          * rhoinf
    w[:, gh:, 1]     = road[:,:]*uad[:,:] * rhoinf * uinf
    w[:, gh:, 2]     = road[:,:]*vad[:,:] * rhoinf * uinf
    w[:, gh:, 4]     = road[:,:]*Ead[:,:] * rhoinf * einf

    f_init.set_bndbl_2d(w, field, wbd, im)

    ######## Restart from a previous solution with exactly the same mesh
    import restart_init as ri
    filet = './Wksp/dnc_5/initialisation_gh.dat'
    # Xin, Yin, roin, rouin, rovin, rowin, roein = ri.read_init(filet)

    ## with exactly the same mesh
    # w[gh:-gh, gh:-gh, 0]     = roin
    # w[gh:-gh, gh:-gh, 1]     = rouin
    # w[gh:-gh, gh:-gh, 2]     = rovin
    # w[gh:-gh, gh:-gh, 3]     = rowin
    # w[gh:-gh, gh:-gh, 4]     = roein
    ## OR interpolate from a different mesh (only valid for cartesian rectangular grid)
    # import interpgrid
    # w[gh:-gh, gh:-gh, 0] = interpgrid.interpgrid(Xin, Yin, roin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 1] = interpgrid.interpgrid(Xin, Yin, rouin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 2] = interpgrid.interpgrid(Xin, Yin, rovin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 3] = interpgrid.interpgrid(Xin, Yin, rowin, xc[gh:-gh,:], yc[:,gh:-gh])
    # w[gh:-gh, gh:-gh, 4] = interpgrid.interpgrid(Xin, Yin, roein, xc[gh:-gh,:], yc[:,gh:-gh])
  
    filename = out_dir + '/initialisation_gh.dat'
    __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

    filename = out_dir + '/' + treename
    writeNPZ(filename, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, cfl, freqres, freqsort, pinf=pinf)


def bl2d_fromNPZtree(treename, ite = 10, compmode = 'fixed_point', out_dir = 'totodir', isresol= False):

    filename = out_dir + '/' + treename
    # tree, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl = readCGNStree(filename)
    # tree, im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, pinf = readCGNStree(filename)
    im, jm, gh, w, x0, y0, vol, volf, nx, ny, xc, yc, field, wbd, res, sch, k2, k4, interf1, interf2, interf3, interf4, lf, cp, cv, gam, cs, tref, muref, rgaz, mach, prandtl, pinf, cfl, freqres, freqsort = readNPZtree(filename)

    Lmax = 1.e6  #0.25e6  #0.10e6  #0.40e6  #1.e6
    r2 = 1.    #elsA 0.2  #0.1   #0.2   #1.
    r3 = 1.3   #elsA 1.3  #0.01   #1.3
    r4 = 30.     #els1 1.   #0.005  #0.01*100000    #0.02   #0.01?   #0.002  #14.  #30

    routinein  = lf[0]
    routineout = lf[1]
    routinenr  = lf[2]
    routinew   = lf[3]
    routinesch = lf[4]
    libbnd     = lf[5]
    libsch     = lf[6]

    finflow  = eval("%s.%s"    % (libbnd, routinein))
    foutflow = eval("%s.%s"    % (libbnd, routineout))
    fnoref   = eval("%s.%s"    % (libbnd, routinenr ))
    fwall    = eval("%s.%s"    % (libbnd, routinew  ))
    fsch     = eval("%s.%s"    % (libsch, routinesch))

    # Time Marching Loop
    if compmode == 'direct':
        # freqres = ite/2
        # freqsort= ite/1
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
        freq = freqsort
        dt  = cfl * (yc[gh,gh+1] - yc[gh,gh]) / (1./mach + 1.)
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

    ## Implicit
    elif compmode == 'impli':
        # freqres = ite/2
        # freqsort= ite/1
        # fimpli   = lf[-1]
        fimpli   = eval("f_lhs.impli_matrix_free_2d")
        time = 0.
        dtcoef = 1.
        dw = _np.zeros((im + 2*gh    , jm + 2*gh    , 5), order='F')
        # Boundary on state vector
        finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
        # finflow(w,'Ilo', interf1, field,im,jm)   
        fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)       
        foutflow(w,'Ihi', interf2, im, jm, gh)
        fwall(w,'Jlo', gam, interf3, gh, im, jm)

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

        routinein  += '_d'
        routineout += '_d'
        routinenr  += '_d'
        routinew   += '_d'
        routinesch += '_d'
        libbnd     = 'f_lin'
        libsch     = 'f_lin'

        flinoutflow = eval("%s.%s"    % (libbnd, routineout))
        flininflow  = eval("%s.%s"    % (libbnd, routinein))
        flinnoref   = eval("%s.%s"    % (libbnd, routinenr ))
        flinwall    = eval("%s.%s"    % (libbnd, routinew  ))
        flinsch     = eval("%s.%s"    % (libsch, routinesch))

        timeconstructjac = 0.
        timeremove  = 0.
        timecoefdiag = 0.
        timejaccsc   = 0.
        timejacinv   = 0.
        # cfl = 1.e10 
        dt  = cfl * (yc[gh,gh+1] - yc[gh,gh]) / (1./mach + 1.)
        dtm1 = 1./dt
        for it in range(1,ite+1):

            wd   = _np.zeros((im+2*gh, jm+2*gh), order='F')
            resd = _np.zeros((im+2*gh, jm+2*gh), order='F')

            # Boundary on state vector
            # finflow(w,'Ilo', interf1, field,im,jm)          
            finflow(w,'Ilo',interf1,field,nx,ny,gam,im,jm)
            fnoref(w,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
            foutflow(w,'Ihi', interf2, im, jm, gh)
            fwall(w,'Jlo', gam, interf3, gh, im, jm)
                        
            filename = out_dir + '/initialisation_gh.dat'
            __writestate_center_gh(filename, im+2*gh, jm+2*gh, w, xc, yc)

            # Compute spatial discretization
            # if sch == 'dnc':
            sourcear = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
            sourceard = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
            if 'dnc' in sch:    
                fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                # fsch(res, w, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm, Lmax, r2, r3, r4, sourcear)
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
            lasolver = 'direct'
            if lasolver == 'gmres':
                nl2res = 1.
                print(" GMRES Not Yet implemented ")
                sys.exit(2)

            else:
                # construct Jacobian
                wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
                resd = _np.zeros((im+2*gh, jm+2*gh,5), order='F')

                nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
                Jac = _np.zeros((nbentry), order='F')
                IA  = _np.zeros((nbentry), dtype=_np.int32, order='F')
                JA  = _np.zeros((nbentry), dtype=_np.int32, order='F')

                timeinjac = timeit.time.time()
                ## relaxed on diag:
                r = _np.max([norm[:3]*norm0m1[:3], ninf[:3]*ninf0m1[:3]])
                cflm1 = r*dtm1
                print('iter = ', it, flush=True)
                print("1/cfl = ", cflm1)
                print("Residual norm:", norm, flush=True)
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
            
                            flininflow(w,wd,'Ilo',interf1,field,nx,ny,gam,im,jm)                
                            flinnoref(w,wd,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
                            flinoutflow(w,wd,'Ihi', interf2, im, jm, gh)
                            flinwall(w,wd,'Jlo', gam, interf3, gh, im, jm)
                                                       
                            if 'dnc' in sch:      
                                flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm)
                                # flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm, Lmax, r2, r3, r4, sourcear, sourceard)
                                # flinsch(res, resd, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, k2, k4, im, jm, eps2ar, eps4ar, divu2ar, vort2ar)
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

                            f_misc.computejacobianfromjv_relaxed(Jac,IA,JA,resd,m,l,k,gh,coefdiag)
                            # if it == ite:
                            #     f_misc.computejacobianfromjv(Jac,IA,JA,resd,m,l,k,gh,im,jm)
                            # else:
                            #     f_misc.computejacobianfromjv_relaxed(Jac,IA,JA,resd,m,l,k,gh,coefdiag)

                ## Remove the zero stored
                timeinremove = timeit.time.time()
                timeconstructjac += (timeinremove - timeinjac)/ite
                mini = 2.e-16  #2.e-16
                IA, JA, Jac = remove_zero_jac(IA, JA, Jac, mini)
                nbentry = _np.shape(Jac)[0]
                print("Number of non-zero inside the Jacobian:", nbentry, flush=True)
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
                ksp.destroy()
                if it==1 or it==2 or it==3:
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

        print('Newton finished, compute Jacobian', flush=True)

        Jacvol = _np.zeros((nbentry), order='F')
        for k in range(nbentry):
            # Jacvol[k] = Jac[k]/vol[(IA[k]/(jm*5))+gh, ((JA[k]%(jm*5))/5)+gh] ##UNE GROSSE CONNERIE
            Jacvol[k] = Jac[k]/vol[(IA[k]//(jm*5))+gh, ((IA[k]%(jm*5))//5)+gh]
        Jacsurvol = pet.createMatPetscCSR(IA, JA, Jacvol, im*jm*5, im*jm*5, 5*(2*gh+1)**2)  

        print('filling NPZ', flush=True)
        filename = out_dir + '/' + treename
        fillNPZ(filename, w, res, IA, JA, Jacvol, gh)

        ### Compute Jacobian Dz and Dzz contributions

        print('Compute 3D contributions of the Jacobian', flush=True)

        wd   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
        dz   = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
        dz2  = _np.zeros((im+2*gh, jm+2*gh,5), order='F')
        nbentry = im*jm * (2*gh+1)*(2*gh+1) * 5*5
        Jacdz  = _np.zeros((nbentry), order='F')
        IAdz   = _np.zeros((nbentry), dtype=_np.int32, order='F')
        JAdz   = _np.zeros((nbentry), dtype=_np.int32, order='F')
        Jacdz2 = _np.zeros((nbentry), order='F')
        IAdz2  = _np.zeros((nbentry), dtype=_np.int32, order='F')
        JAdz2  = _np.zeros((nbentry), dtype=_np.int32, order='F')
        for m in range(5):
            for l in range(1 + 2*gh):
                for k in range(1 + 2*gh):
                    wd *= 0.
                    f_misc.testvector(wd,m,l,k,gh,im,jm)

                    flininflow(w,wd,'Ilo',interf1,field,nx,ny,gam,im,jm)                
                    flinnoref(w,wd,wbd,'Jhi',interf4,nx,ny,gam,gh,im,jm)
                    flinoutflow(w,wd,'Ihi', interf2, im, jm, gh)
                    flinwall(w,wd,'Jlo', gam, interf3, gh, im, jm)

                    f_dz.coeffs_5p_dz(dz, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, im, jm)
                    f_dz.coeffs_5p_dz2(dz2, w, wd, x0, y0, nx, ny, xc, yc, vol, volf, gh, cp, cv, prandtl, gam, rgaz, cs, muref, tref, cs, im, jm)

                    f_misc.computejacobianfromdz(Jacdz,IAdz,JAdz,dz,m,l,k,gh,im,jm)
                    f_misc.computejacobianfromdz(Jacdz2,IAdz2,JAdz2,dz2,m,l,k,gh,im,jm)

        IAdz, JAdz, Jacdz = remove_zero_jac(IAdz, JAdz, Jacdz, mini)
        IAdz2, JAdz2, Jacdz2 = remove_zero_jac(IAdz2, JAdz2, Jacdz2, mini)

        print('filling NPZ 3D', flush=True)
        fillNPZ_3D(filename, IAdz, JAdz, Jacdz, IAdz2, JAdz2, Jacdz2)
        
        ### Resolvent : compute eigenvalues and eigenvectors
