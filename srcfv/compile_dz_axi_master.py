# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

import os

#### Compile Matrix functions for Dz

dir = 'matrix_dz_axi/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
prepro_dir = srcdir + 'prepro_matrix_dz_axi'
os.system('mkdir -p %s' % prepro_dir)

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.F90') :
        filein  = dir + file
        fppfile = 'prepro_matrix_dz_axi/' + file1 + '.f90'
        exec_str = 'fpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        # exec_str = 'cpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        os.system(exec_str)

#### Linearize Matrix functions for Dz

dir = 'prepro_matrix_dz_axi/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
tapdir = srcdir + 'dz_axi/'

os.system('mkdir -p %s' % tapdir)

for file in files:
    file1,ext = os.path.splitext(file)
    ##  print file1
    ##  print ext
    if (ext == '.f90') :
        filein  = dir + file
        tapfile = tapdir + file1 + '_d' + ext
        exec_str = 'tapenade %s -head "%s(func0, func1, func2, func3, func4, func5, func6, func7, func8, func9, func10, func11, func12, func13, func14, func15)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
        #exec_str = 'tapenade %s -head "%s(func0)/(w)" -multi -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
        os.system(exec_str)
        print(tapfile)
        s1 += tapfile + ' '

#### Compile Matrix functions for Dzz

dir = 'matrix_dz2_axi/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
prepro_dir = srcdir + 'prepro_matrix_dz2_axi'
os.system('mkdir -p %s' % prepro_dir)

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.F90') :
        filein  = dir + file
        fppfile = 'prepro_matrix_dz2_axi/' + file1 + '.f90'
        exec_str = 'fpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        # exec_str = 'cpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        os.system(exec_str)

#### Linearize Matrix functions for Dzz

dir = 'prepro_matrix_dz2_axi/'
files = os.listdir(dir+'.')
# print files
s1 = str()
s2 = str()

srcdir = './'
tapdir = srcdir + 'dz_axi/'

os.system('mkdir -p %s' % tapdir)

for file in files:
    file1,ext = os.path.splitext(file)
    ##  print file1
    ##  print ext
    if (ext == '.f90') :
        filein  = dir + file
        tapfile = tapdir + file1 + '_d' + ext
        exec_str = 'tapenade %s -head "%s(func0, func1, func2, func3)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
        #exec_str = 'tapenade %s -head "%s(func0)/(w)" -multi -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
        os.system(exec_str)
        print(tapfile)
        s1 += tapfile + ' '

#### Compile the matrices built Dz and Dzz

dir = 'dz_axi/'
tn = 'f_dz_axi'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.F90') :
        filein  = dir + file
        fppfile = 'prepro_dz_axi/' + file1 + '.f90'
:x
        # exec_str = 'cpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        os.system(exec_str)
        s1 += filein + ' '

# scheme
fcompiler = " --fcompiler='intelem' "
# fcompiler = " --fcompiler='gnu95' "
if 'gnu95' in fcompiler:
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-flto -O3'"
    f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-Ofast'"
elif 'intel' in fcompiler :
    f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"
options = fcompiler + f90flags
exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
print(exec_str)
os.system(exec_str)        

