import os

dir = 'tangent_chu/'
tn = 'f_linchu'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'

for file in files:
    file1,ext = os.path.splitext(file)
    if (ext == '.f90') :
        filein  = dir + file
        #fppfile = 'prepro/' + file1 + '.f90'
        #exec_str = 'fpp -I%s -P -C %s %s' % (srcdir, filein, fppfile)
        #os.system(exec_str)
        s1 += filein + ' '

# scheme
fcompiler = " --fcompiler='intelem' "
f90flags  = " --f90flags='-mcmodel=medium -funroll-loops -r8' --opt='-ipo -O3'"#"--f90flags='-mcmodel=medium -funroll-loops -fdefault-real-8' --opt='-O3'"
options = fcompiler + f90flags
exec_str = 'f2py -m %s -c %s %s' % (tn, s1, options)
print(exec_str)
os.system(exec_str)

