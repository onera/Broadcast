# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = 'prepro_matrix_dz2/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
tapdir = srcdir + 'dz/'

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


