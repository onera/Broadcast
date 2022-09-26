# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = 'prepro/'
lroutines = ['jn_match_2d', 'jn_match_grad_2d', 'jn_match_fx_2d']
lrgeo = ['jn_match_geom_2d', 'jn_match_face_geom_2d']
files = os.listdir(dir+'.')
print(files)
s1 = str()
s2 = str()

srcdir = './'
tapdir = srcdir + 'adjoint/'

os.system('mkdir -p %s' % tapdir)

for file in files:
    if (file == 'jn_match_geom.f90') :
        file1,ext = os.path.splitext(file)
        filein  = dir + file
        tapfile = tapdir + file1 + '_b%s' % (ext)
        for routine in lrgeo:
            exec_str = 'tapenade %s -head "%s(wr)/(wd)" -b -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            print(exec_str)
            os.system(exec_str)
            os.system('mv %s %s' % (tapfile, tapdir + '%s_b%s' % (routine, ext)))

    elif (file == 'jn_match.f90') :
        file1,ext = os.path.splitext(file)
        filein  = dir + file
        tapfile = tapdir + file1 + '_b%s' % (ext)
        for routine in lroutines:
            exec_str = 'tapenade %s -head "%s(wr)/(wd)" -b -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
            print(exec_str)
            os.system(exec_str)
            os.system('mv %s %s' % (tapfile, tapdir + '%s_b%s' % (routine, ext)))



# scheme
##  fcompiler = " --fcompiler='intelem' "
##  f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"#"--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
##  options = fcompiler + f90flags
##  tn = 'f_dirlin'
##  exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
##  os.system(exec_str)

