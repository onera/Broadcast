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
        tapfile = tapdir + file1 + '_dd' + ext
        flag = False
        if 'iso_profile' in file:
            routine = 'bc_wall_viscous_iso_profile_2d'
            exec_str = 'tapenade %s -head "%s(w)/(twallprof,w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            flag = True
        elif 'blow_profile' in file:
            routine = 'bc_wall_blow_profile_2d'
            exec_str = 'tapenade %s -head "%s(w)/(velprof,w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            flag = True
        if flag:
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

