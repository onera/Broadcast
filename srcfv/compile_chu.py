import os

dir = 'chu/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
prepro_dir = srcdir + 'prepro_chu'
os.system('mkdir -p %s' % prepro_dir)

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.F90') :
        filein  = dir + file
        fppfile = 'prepro_chu/' + file1 + '.f90'
        exec_str = 'fpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        os.system(exec_str)

        s1 += fppfile + ' '

# scheme
fcompiler = " --fcompiler='intelem' "
# f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"
# elsA fortran compiler option---------------------
E_FORTRAN_FLAGS    = ' -fpp -w -w90 -w95 -cm -ip -axCORE-AVX512,CORE-AVX2 -fno-alias -no-vec-guard-write -qopt-streaming-stores=never -qopt-malloc-options=3 '+\
                     ' -fno-fnalias -convert big_endian -ip -funroll-loops -132 -fp-model strict -r8 -fma -qoverride-limits -unroll0 -O3'
f90flags  = " --f90flags='%s' " % E_FORTRAN_FLAGS
# ---------------------elsA fortran compiler option
print(f90flags)
#f90flags  = "--f90flags='-mcmodel=medium -check bounds'"
options = fcompiler + f90flags
tn = 'f_chu'
exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
os.system(exec_str)

