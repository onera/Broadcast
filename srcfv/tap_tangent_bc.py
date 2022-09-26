# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = 'prepro/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

tapmode = 'tangent/'
srcdir = './'
tapdir = srcdir + tapmode

os.system('mkdir -p %s' % tapdir)

for file in files:
    if 'rbc' in file:
        print('toto rbc niet %s' % file)
    elif ('bc' in file) :
        file1,ext = os.path.splitext(file)
        print(file1)
        print(ext)
        filein  = dir + file
        tapfile = tapdir + file1 + '_d' + ext
        if 'extrap' in file:
            for order in [2,3,4,5,7,9]:
                routine = 'bc_extrapolate_o%i_2d' % order
                exec_str = 'tapenade %s -head "%s(w)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
                print(exec_str)
                os.system(exec_str)
                os.system('mv %s %s' % (tapfile, tapdir + file1 + 'o%i_d%s' % (order, ext)))
        else:
            if 'sa' in file:
                if 'no' in file:
                    routine = 'bc_no_reflexion_sa_2d'
                elif 'viscous' in file:
                    routine = 'bc_wall_viscous_adia_sa_2d'
                elif 'sup' in file:
                    routine = 'bc_supandsubinlet_sa_2d'
            elif 'mirror' in file:
                routine = 'bc_mirror_sa_2d'
            else:
                if 'no' in file:
                    if 'hess' in file:
                        routine = 'bc_no_reflexion_forhessian_2d'
                    else:
                        routine = 'bc_no_reflexion_2d'
                elif 'viscous' in file:
                    if 'hess' in file:
                        routine = 'bc_wall_viscous_adia_forhessian_2d'
                    elif 'iso' in file:
                        routine = 'bc_wall_viscous_iso_2d'
                    else:
                       routine = 'bc_wall_viscous_adia_2d'
                elif 'sup' in file:
                    if 'hess' in file:
                        routine = 'bc_supandsubinlet_forhessian_2d'
                    else:
                        routine = 'bc_supandsubinlet_2d'
                elif 'mirror' in file:
                    if 'anti' in file:
                        routine = 'bc_mirror_anti_2d'
                    else:
                        routine = 'bc_mirror_2d'
                elif 'symmetry' in file:
                    if 'anti' in file:
                        routine = 'bc_antisymmetry_2d'
                    else:
                        routine = 'bc_symmetry_2d'
                elif 'pressure' in file:
                    routine = 'bc_pressure_2d'
                else:
                    routine = 'bc_general_2d'
            if 'hess' in file:            
                exec_str = 'tapenade %s -head "%s(wout)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            else:
                exec_str = 'tapenade %s -head "%s(w)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            print(exec_str)
            os.system(exec_str)
        print(tapfile)
        s1 += tapfile + ' '

# scheme
##  fcompiler = " --fcompiler='intelem' "
##  f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"#"--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
##  options = fcompiler + f90flags
##  tn = 'f_dirlin'
##  exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
##  os.system(exec_str)

