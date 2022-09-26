# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os
import string

# file to differentiate uT H v from tangent routines

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
    file1,ext = os.path.splitext(file)
    if 'bc' in file:
        print('%s has not been linearized ....' % file)
    elif 'geom' in file:
        print('%s has not been linearized ....' % file)
    elif 'jn' in file:
        print('%s has not been linearized ....' % file)
    else:
        if (ext == '.f90') :
            filein  = dir + file
            # file1str = string.split(file1,"_")
            file1str = file1.split("_")
            file2str = file1str[:-1]
            file2  = file2str[0]
            for s in file2str[1:]: file2 += '_' + s
            tapfile = tapmode + file2 + '_db'
            # print 'filein   = ', filein
            # print 'file1    = ', file1
            # print 'file2    = ', file2
            # print 'tapfile  = ', tapfile
            exec_str = 'tapenade %s -head "%s_2d_d(residud)/(w)" -d -optim diffliveness -optim statictape -O %s' % (filein, file2, tapdir)
            os.system(exec_str)
            s1 += tapfile + ' '



# scheme
##  fcompiler = " --fcompiler='intelem' "
##  f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"#"--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
##  options = fcompiler + f90flags
##  tn = 'f_sch'
##  exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
##  os.system(exec_str)

