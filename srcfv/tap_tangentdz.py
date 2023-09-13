import os

dir = 'prepro_dz/'
files = os.listdir(dir+'.')
#print files
s1 = str()
s2 = str()

srcdir = './'
tapdir = srcdir + 'tangentdz/'

os.system('mkdir -p %s' % tapdir)

for file in files:
    file1,ext = os.path.splitext(file)
    ##  print file1
    ##  print ext
    filein  = dir + file
    tapfile = tapdir + file1 + '_d' + ext
    flag = False
    if 'coeffs_5p_dz2' in file:
        routine = 'coeffs_5p_dz2'
        exec_str = 'tapenade %s -head "%s(dz2_out)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
        #exec_str = 'tapenade %s -head "%s(func0)/(w)" -multi -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
        flag = True
    elif 'coeffs_5p_dz' in file:
        routine = 'coeffs_5p_dz'
        exec_str = 'tapenade %s -head "%s(dz_out)/(w)" -d -optim diffliveness  -optim statictape -O %s' % (filein, routine, tapdir)
        #exec_str = 'tapenade %s -head "%s(func0)/(w)" -multi -optim diffliveness  -optim statictape -O %s' % (filein, file1, tapdir)
        flag = True
    if flag:
        print(exec_str)    
        os.system(exec_str)
        print(tapfile)
        s1 += tapfile + ' '


