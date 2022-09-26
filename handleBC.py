# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#!/usr/bin/env python
## Functions to handle the BC in BROADCAST main

import numpy as _np
import srcfv.f_bnd     as f_bnd
import srcfv.f_lin     as f_lin
import srcfv.f_hess    as f_hess

def checkinterf(interfl, interfdic, interfloc, im, jm, gh):
    ''' Check if all the interfaces are defined especially inside the ghost cells '''

    summ = 0
    for k in range(len(interfl)):
        inter = interfl[k]
        if interfloc[inter] == 'Ilo':
            summ += interfdic[inter][1,1] - interfdic[inter][0,1] +1
        elif interfloc[inter] == 'Ihi':
            summ += interfdic[inter][1,1] - interfdic[inter][0,1] +1
        elif interfloc[inter] == 'Jlo':
            summ += interfdic[inter][1,0] - interfdic[inter][0,0] +1
        elif interfloc[inter] == 'Jhi':
            summ += interfdic[inter][1,0] - interfdic[inter][0,0] +1        

    if summ == 2*im + 2*jm + 4*gh:
        return 'GoodJob'
    elif summ > 2*im + 2*jm + 4*gh:
        return 'Error : Check definition of the interfaces, interfaces are overlapping, only one single BC can be applied to each cell'
    else:
        return 'Error : Check definition of the interfaces, BC must also be defined inside the ghost cells at the corners of the domain'

def sortBC(interfl, interfdic, interfloc, im, jm):
    ''' Sort the interfaces to apply the BC in the good order '''

    cornerdic = dict()
    k = 0
    for k in range(len(interfl)):
        inter = interfl[k]
        if interfloc[inter] == 'Ilo':
            if interfdic[inter][0,1] < 1:
                cornerdic['IloJlo'] = inter
            if interfdic[inter][1,1] > jm:
                cornerdic['IloJhi'] = inter
        elif interfloc[inter] == 'Ihi':
            if interfdic[inter][0,1] < 1:
                cornerdic['IhiJlo'] = inter
            if interfdic[inter][1,1] > jm:
                cornerdic['IhiJhi'] = inter     
        elif interfloc[inter] == 'Jlo':
            if interfdic[inter][0,0] < 1:
                cornerdic['IloJlo'] = inter
            if interfdic[inter][1,0] > im:
                cornerdic['IhiJlo'] = inter
        elif interfloc[inter] == 'Jhi':
            if interfdic[inter][0,0] < 1:
                cornerdic['IloJhi'] = inter
            if interfdic[inter][1,0] > im:
                cornerdic['IhiJhi'] = inter  
        if len(cornerdic) == 4:
            break  

    sortlist = []
    corner = cornerdic['IloJlo']
    sortlist.append(corner)
    for k in range(len(interfl)):
        inter = interfl[k]
        if interfdic[inter][0,0] == 1 and interfdic[inter][0,1] == 1:
            if interfloc[corner] == 'Ilo' and interfloc[inter] == 'Jlo':
                sortlist = [inter] + sortlist
                break
            elif interfloc[corner] == 'Jlo' and interfloc[inter] == 'Ilo':
                sortlist = [inter] + sortlist
                break    
    corner = cornerdic['IloJhi']
    if corner not in sortlist:
        sortlist.append(corner)
    for k in range(len(interfl)):
        inter = interfl[k]
        if interfdic[inter][0,0] == 1 and interfdic[inter][0,1] == jm:
            if inter not in sortlist:
                if interfloc[corner] == 'Ilo' and interfloc[inter] == 'Jhi':
                    sortlist = [inter] + sortlist
                    break
                elif interfloc[corner] == 'Jhi' and interfloc[inter] == 'Ilo':
                    sortlist = [inter] + sortlist
                    break
    corner = cornerdic['IhiJhi']
    if corner not in sortlist:
        sortlist.append(corner)
    for k in range(len(interfl)):
        inter = interfl[k]
        if interfdic[inter][0,0] == im and interfdic[inter][0,1] == jm:
            if inter not in sortlist:
                if interfloc[corner] == 'Ihi' and interfloc[inter] == 'Jhi':
                    sortlist = [inter] + sortlist
                    break
                elif interfloc[corner] == 'Jhi' and interfloc[inter] == 'Ihi':
                    sortlist = [inter] + sortlist
                    break
    corner = cornerdic['IhiJlo']
    if corner not in sortlist:
        sortlist.append(corner)
    for k in range(len(interfl)):
        inter = interfl[k]
        if interfdic[inter][0,0] == im and interfdic[inter][0,1] == 1:
            if interfloc[corner] == 'Ihi' and interfloc[inter] == 'Jlo':
                if inter not in sortlist:
                    sortlist = [inter] + sortlist
                    break
                else:
                    sortlist.remove(corner)
                    sortlist.append(corner)
                    break
            elif interfloc[corner] == 'Jlo' and interfloc[inter] == 'Ihi':
                if inter not in sortlist:
                    sortlist = [inter] + sortlist
                    break
                else:
                    sortlist.remove(corner)
                    sortlist.append(corner)
                    break
    for k in range(len(interfl)):
        if interfl[k] not in sortlist:
            sortlist.append(interfl[k])

    return sortlist


def applyBC(interfl, interfdic, interfloc, BCdict, modeBC, wlist, BCtypel, nx, ny, gam, rgaz, im, jm, gh, wbddic=None, fielddic=None, twalldic=None, velwalldic=None, prldic=None, trldic=None):
    ''' Apply the BC on w, wd, ... depending on modeBC '''
    ### modeBC=0 -> apply BC on w ; modeBC=1 -> apply linearised BC on wd ; ...

    lf = BCtypel[0]
    routinein  = lf[0]
    routineout = lf[1]
    routinenr  = lf[2]
    routinew   = lf[3]
    routinebw  = lf[5]
    routinejn  = lf[6]
    libbnd     = lf[7]

    finflow  = eval("%s.%s"    % (libbnd, routinein ))
    foutflow = eval("%s.%s"    % (libbnd, routineout))
    fnoref   = eval("%s.%s"    % (libbnd, routinenr ))
    fwall    = eval("%s.%s"    % (libbnd, routinew  ))
    fsym     = eval("%s.%s"    % (libbnd, routinebw ))
    fjn      = eval("%s.%s"    % (libbnd, routinejn ))

    if modeBC >= 1:
        wlist[0][:gh,:,:]  = 0.
        wlist[0][:,:gh,:]  = 0.
        wlist[0][-gh:,:,:] = 0.
        wlist[0][:,-gh:,:] = 0.

        if modeBC == 1:
            lflin = BCtypel[1]
            routinein  = lflin[0]
            routineout = lflin[1]
            routinenr  = lflin[2]
            routinew   = lflin[3]
            routinebw  = lflin[5]
            routinejn  = lflin[6]
            libbnd     = lflin[7]
            flinoutflow = eval("%s.%s"    % (libbnd, routineout))
            flininflow  = eval("%s.%s"    % (libbnd, routinein))
            flinnoref   = eval("%s.%s"    % (libbnd, routinenr ))
            flinwall    = eval("%s.%s"    % (libbnd, routinew  ))
            flinsym     = eval("%s.%s"    % (libbnd, routinebw ))
            flinjn      = eval("%s.%s"    % (libbnd, routinejn ))

        elif modeBC == 2:
            libbnd   = 'f_hess'
            routineout += '_d_d'
            routinein  += '_d_d'
            routinenr  += '_d_d'
            routinew   += '_d_d'
            routinebw  += '_d_d'
            routinejn  += '_d_d'
            flin2inflow  = eval("%s.%s"    % (libbnd, routinein))
            flin2noref   = eval("%s.%s"    % (libbnd, routinenr ))
            flin2wall    = eval("%s.%s"    % (libbnd, routinew  ))

    for k in range(len(interfl)):
        inter = interfl[k]
        if BCdict[inter] == 'finflow':
            if 'supandsub' in routinein:
                if modeBC == 1:
                    # flininflow(wlist[0],wlist[1],interfloc[inter],interfdic[inter],fielddic[inter],nx,ny,gam,im,jm) 
                    flininflow(wlist[0],wlist[1],interfloc[inter],interfdic[inter],fielddic[inter],nx,ny,gam,im,jm,gh) 
                elif modeBC == 2:
                    # flin2inflow(wlist[0],wlist[1],wlist[1],wlist[2],wlist[0],wlist[1],wlist[1],interfloc[inter],interfdic[inter],fielddic[inter],nx,ny,gam,im,jm)
                    flin2inflow(wlist[0],wlist[1],wlist[1],wlist[2],wlist[0],wlist[1],wlist[1],interfloc[inter],interfdic[inter],fielddic[inter],nx,ny,gam,im,jm,gh)
                # finflow(wlist[0],interfloc[inter],interfdic[inter],fielddic[inter],nx,ny,gam,im,jm) 
                finflow(wlist[0],interfloc[inter],interfdic[inter],fielddic[inter],nx,ny,gam,im,jm,gh) 
            else:      
                finflow(wlist[0],interfloc[inter],interfdic[inter],fielddic[inter],im,jm)
        elif BCdict[inter] == 'foutflow':
            if modeBC == 1:
                flinoutflow(wlist[0],wlist[1],interfloc[inter],interfdic[inter],im,jm,gh)
            foutflow(wlist[0],interfloc[inter],interfdic[inter],im,jm,gh)       
        elif BCdict[inter] == 'fwall':
            if inter in twalldic:
                if modeBC == 1:
                    flinwall(wlist[0],wlist[1],twalldic[inter],interfloc[inter],gam,rgaz,interfdic[inter],gh,im,jm)
                elif modeBC == 2:
                    flin2wall(wlist[0],wlist[0],wlist[2],wlist[0],wlist[1],wlist[2],twalldic[inter],interfloc[inter],gam,rgaz,interfdic[inter],gh,im,jm)
                    # flin2wall(wlist[0],wlist[1],wlist[2],_np.zeros((im+2*gh,jm+2*gh,5),order='F'),twalldic[inter],_np.zeros((_np.shape(twalldic[inter])[0]),order='F'),_np.zeros((_np.shape(twalldic[inter])[0]),order='F'),interfloc[inter],gam,rgaz,interfdic[inter],gh,im,jm)
                fwall(wlist[0],twalldic[inter],interfloc[inter],gam,rgaz,interfdic[inter],gh,im,jm)
            elif inter in velwalldic:
                if modeBC == 1:
                    flinwall(wlist[0],wlist[1],velwalldic[inter],interfloc[inter],gam,interfdic[inter],gh,im,jm)
                elif modeBC == 2:
                    flin2wall(wlist[0],wlist[0],wlist[2],wlist[0],wlist[1],wlist[2],velwalldic[inter],interfloc[inter],gam,interfdic[inter],gh,im,jm)
                    # flin2wall(wlist[0],wlist[1],wlist[2],_np.zeros((im+2*gh,jm+2*gh,5),order='F'),velwalldic[inter],_np.zeros((_np.shape(velwalldic[inter])[0]),order='F'),_np.zeros((_np.shape(velwalldic[inter])[0]),order='F'),interfloc[inter],gam,interfdic[inter],gh,im,jm)
                fwall(wlist[0],velwalldic[inter],interfloc[inter],gam,interfdic[inter],gh,im,jm)    
            else:
                if modeBC == 1:
                    flinwall(wlist[0],wlist[1],interfloc[inter],gam,interfdic[inter],gh,im,jm)
                elif modeBC == 2:
                    flin2wall(wlist[0],wlist[0],wlist[2],wlist[0],wlist[1],wlist[2],interfloc[inter],gam,rgaz,interfdic[inter],gh,im,jm)  
                fwall(wlist[0],interfloc[inter],gam,interfdic[inter],gh,im,jm)
        elif BCdict[inter] == 'fnoref':
            if modeBC == 1:
                flinnoref(wlist[0],wlist[1],wbddic[inter],interfloc[inter],interfdic[inter],nx,ny,gam,gh,im,jm) 
            elif modeBC == 2:
                flin2noref(wlist[0],wlist[0],wlist[2],wlist[0],wlist[1],wlist[2],wbddic[inter],interfloc[inter],interfdic[inter],nx,ny,gam,gh,im,jm) 
            fnoref(wlist[0],wbddic[inter],interfloc[inter],interfdic[inter],nx,ny,gam,gh,im,jm) 
        elif BCdict[inter] == 'fsym':
            if modeBC == 1:
                flinsym(wlist[0],wlist[1],interfloc[inter],interfdic[inter],nx,ny,gh,im,jm)
            fsym(wlist[0],interfloc[inter],interfdic[inter],nx,ny,gh,im,jm)  
        elif BCdict[inter] == 'fjn':
            if modeBC == 1:
                fjn(wlist[1],prldic[inter][0,0],gh,gh,gh,gh,im,jm,wlist[1],prldic[inter][1,1],gh,gh,gh,gh,im,jm,trldic[inter][0,0])
                fjn(wlist[1],prldic[inter][1,0],gh,gh,gh,gh,im,jm,wlist[1],prldic[inter][0,1],gh,gh,gh,gh,im,jm,trldic[inter][1,0])  
            elif modeBC == 2:
                fjn(wlist[2],prldic[inter][0,0],gh,gh,gh,gh,im,jm,wlist[2],prldic[inter][1,1],gh,gh,gh,gh,im,jm,trldic[inter][0,0])
                fjn(wlist[2],prldic[inter][1,0],gh,gh,gh,gh,im,jm,wlist[2],prldic[inter][0,1],gh,gh,gh,gh,im,jm,trldic[inter][1,0])  
            fjn(wlist[0],prldic[inter][0,0],gh,gh,gh,gh,im,jm,wlist[0],prldic[inter][1,1],gh,gh,gh,gh,im,jm,trldic[inter][0,0])
            fjn(wlist[0],prldic[inter][1,0],gh,gh,gh,gh,im,jm,wlist[0],prldic[inter][0,1],gh,gh,gh,gh,im,jm,trldic[inter][1,0])  

    return wlist 

