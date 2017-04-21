import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## one GTO optimization value is 0.0105604,-0.0548991
with open("search.out", "w") as f:
    print_timestamp("start", f)

yos = [0.0001, 0.00014, 0.0002, 0.0003, 0.0005, 0.0007,
       0.001,  0.0014,  0.002,  0.003,  0.005,  0.007,
       0.01,   0.014,   0.02,   0.03,   0.05,   0.07]

#zos = [x-1.0j*y for (x,y) in product(xos, yos)]

f = open('search.csv', 'w')
z0s = [0.01*1.8**n for n in range(15)]

for y0 in yos:
    basis_info = [('shift', True, 2, z0s, -1.0j*y0)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,

        target  = 'h_pi',
        w0      = 1.0,
        channel = '1s->kp',
        dipole  = 'length',        
        
        tol     = pow(10.0, -5.0),
        maxit   = 50,
        conv    = "grad",
        fdif    = 0.0001,
        grad    = False,
        hess    = False,

        outfile = 'res.out',
        print_level = 0)
    
    opt_res = res['w_res_list'][0][1]
    z0 = opt_res.x[0]
    eps = 0.00001
    if(z0.real > 0.0 and z0.imag < 0.0 and
       100.0 > z0.real and z0.imag > -100.0 and
       opt_res.success):
       subprocess.call("cat res.out >> search.out", shell=True)
       f.write("{0},{1}\n".format(z0.real, z0.imag))

f.close()
with open("search.out", "a") as f:
    print_timestamp("end", f)

