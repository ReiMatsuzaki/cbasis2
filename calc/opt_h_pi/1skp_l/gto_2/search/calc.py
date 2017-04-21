import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## one GTO optimization value is 0.0105604,-0.0548991
with open("search.out", "w") as f:
    print_timestamp("start", f)

xos = [0.001, 0.0015, 0.002, 0.0025, 0.003,  0.005,  0.007,
       0.01,  0.015,  0.02,  0.0025, 0.03,   0.05,   0.07,
       0.1,   0.15,   0.2]
yos = [0.001, 0.0014, 0.002, 0.003, 0.005, 0.007,
       0.01,  0.014, 0.02,  0.03,  0.05,  0.07,
       0.1,   0.14,   0.2,   0.3,   0.5]
zos = [x-1.0j*y for (x,y) in product(xos, yos)]

f = open('search.csv', 'w')
f.write("re_0,im_0,re_1,im_1\n")
for (z0, z1) in combinations(zos, 2):
    basis_info = [('id',  True,  2, z0),
                  ('id',  True,  2, z1)]
    
    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
        w0 = 1.0,
        tol = pow(10.0, -5.0),
        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        outfile = 'res.out',
        print_level = 0)
    
    opt_res = res['w_res_list'][0][1]
    zopt1 = 0.0105604-0.0548991j
    z0 = opt_res.x[0]
    z1 = opt_res.x[1]
    eps = 0.00001
    if(z0.real > 0.0 and z0.imag < 0.0 and
       100.0 > z0.real and z0.imag > -100.0 and
       z1.real > 0.0 and z1.imag < 0.0 and
       100.0 > z1.real and z1.imag > -100.0 and
       opt_res.success):
       subprocess.call("cat res.out >> search.out", shell=True)
       f.write("{0},{1},{2},{3}\n".format(z0.real, z0.imag, z1.real, z1.imag))

f.close()
with open("search.out", "a") as f:
    print_timestamp("end", f)
