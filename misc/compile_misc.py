# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os

dir = './'
files = os.listdir(dir+'.')
#print files
s1 = str()

#srcdir = '/stck/ccontent/workspace/BL_2D/misc/'

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.f90') :
        filein  = dir + file

        s1 += filein + ' '
print(s1)

# scheme
fcompiler = " --fcompiler='intelem' "
# f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"
# elsA fortran compiler option---------------------
E_FORTRAN_FLAGS    = ' -fpp -w -w90 -w95 -cm -ip -axCORE-AVX512,CORE-AVX2 -fno-alias -no-vec-guard-write -qopt-streaming-stores=never -qopt-malloc-options=3 '+\
                     ' -fno-fnalias -convert big_endian -ip -funroll-loops -132 -fp-model strict -r8 -fma -qoverride-limits -unroll0 -O3'
# fcompiler = " --fcompiler='gnu95' "
if 'gnu95' in fcompiler:
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-flto -O3'"
    f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-Ofast'"
elif 'intel' in fcompiler :
    f90flags  = " --f90flags='%s' " % E_FORTRAN_FLAGS
# ---------------------elsA fortran compiler option
print(f90flags)
#f90flags  = "--f90flags='-mcmodel=medium -check bounds'"
options = fcompiler + f90flags
tn = 'f_misc'
exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
os.system(exec_str)

