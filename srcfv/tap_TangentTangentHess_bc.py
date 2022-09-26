# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = 'tangent/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

tapmode = 'tangenttangentHess/'
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
        if (ext == '.f90') :
            filein  = dir + file
            tapfile = tapdir + file1 + '_d' + ext
            
            if 'reflexion_forhess' in file:
                routine = 'bc_no_reflexion_forhessian_2d_d'
                exec_str = 'tapenade %s -head "%s(woutd)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
           
                print(exec_str)
                os.system(exec_str)
                print(tapfile)
                s1 += tapfile + ' '
                
            elif 'no' in file:
                routine = 'bc_no_reflexion_2d_d'
                exec_str = 'tapenade %s -head "%s(wd)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
           
                print(exec_str)
                os.system(exec_str)
                print(tapfile)
                s1 += tapfile + ' '
                
            elif 'adia_forhess' in file:
                routine = 'BC_WALL_VISCOUS_ADIA_FORHESSIAN_2D_D'
                exec_str = 'tapenade %s -head "%s(woutd)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            
                print(exec_str)
                os.system(exec_str)
                print(tapfile)
                s1 += tapfile + ' '
                
            elif 'viscous' in file:
                routine = 'BC_WALL_VISCOUS_ADIA_2D_D'
                exec_str = 'tapenade %s -head "%s(wd)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            
                print(exec_str)
                os.system(exec_str)
                print(tapfile)
                s1 += tapfile + ' '
                
            elif 'supandsubinlet_forhess' in file:
                routine = 'BC_SUPANDSUBINLET_FORHESSIAN_2D_D'
                exec_str = 'tapenade %s -head "%s(woutd)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
                
                print(exec_str)
                os.system(exec_str)
                print(tapfile)
                s1 += tapfile + ' '

            elif 'sup' in file:
                routine = 'BC_SUPANDSUBINLET_2D_D'
                exec_str = 'tapenade %s -head "%s(wd0)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
                
                print(exec_str)
                os.system(exec_str)
                print(tapfile)
                s1 += tapfile + ' '

            # else:
            #     print 'NO 2ND derivative on %s' %file
            
        
# scheme
##  fcompiler = " --fcompiler='intelem' "
##  f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"#"--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
##  options = fcompiler + f90flags
##  tn = 'f_dirlin'
##  exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
##  os.system(exec_str)

