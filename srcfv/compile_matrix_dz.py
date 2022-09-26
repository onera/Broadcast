# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = 'matrix_dz/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
prepro_dir = srcdir + 'prepro_matrix_dz'
os.system('mkdir -p %s' % prepro_dir)

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.F90') :
        filein  = dir + file
        fppfile = 'prepro_matrix_dz/' + file1 + '.f90'
        exec_str = 'fpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        os.system(exec_str)


