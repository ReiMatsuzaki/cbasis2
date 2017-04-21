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

ii = 1.0j
R0s = [1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]
t0s = [-4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4]
yos = [R * np.exp(-ii*t*np.pi/180.0) for R in R0s for t in t0s]

f = open('search.csv', 'w')
z0s = [0.005*1.8**n for n in range(15)]

for y0 in yos:
    print y0
    basis_info = [('scale', True,  2, z0s[:5],   y0),
                  ('scale', False, 2, z0s[5:], 1.0)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
        
        w0      = 1.0,
        target  = 'h_pi',
        channel = '1s->kp',
        dipole  = 'length',        

        tol  = pow(10.0, -5.0),
        grad = False,
        hess = False,
        fdif = 0.0001,

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

