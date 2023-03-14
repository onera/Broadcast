# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
import os
s1 = str()
s1 += 'initialisation.f90 '
s1 += 'set_bnd.f90 '
#
fcompiler = " --fcompiler='intelem' "
# fcompiler = " --fcompiler='gnu95' "
if 'gnu95' in fcompiler:
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3' "
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-flto -O3' "
    f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-Ofast' "
    # f90flags = "--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' "
elif 'intel' in fcompiler :
    f90flags = "--f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"
# f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"#"--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
options = fcompiler + f90flags
tn = 'f_init'
exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
os.system(exec_str)

