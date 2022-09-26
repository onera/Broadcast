# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = 'prepro/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
tapdir = srcdir + 'tangent/'

os.system('mkdir -p %s' % tapdir)

for file in files:
    file1,ext = os.path.splitext(file)
    ##  print file1
    ##  print ext
    if 'rbc' in file:
        if (ext == '.f90') :
            filein  = dir + file
            tapfile = tapdir + file1 + '_2d_d' + ext
            exec_str = 'tapenade %s -head "%s_2d(residu)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
            os.system(exec_str)
            print(tapfile)
            s1 += tapfile + ' '

    elif 'bc' in file:
        print('%s has not been linearized ....' % file)
    elif 'geom' in file:
        print('%s has not been linearized ....' % file)
    elif 'jn' in file:
        print('%s has not been linearized ....' % file)
    elif 'sa' in file:
        print('%s has not been linearized ....' % file)
    else:
        if (ext == '.f90') :
            filein  = dir + file
            tapfile = tapdir + file1 + '_2d_d' + ext
            exec_str = 'tapenade %s -head "%s_2d(residu)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
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

