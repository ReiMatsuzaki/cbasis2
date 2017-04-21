import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

#xos = [0.0001, 0.0002, 0.0003, 0.0005, 0.0007,
#       0.001,  0.002,  0.003,  0.005,  0.007]
#yos = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3]
zfs = [0.005*1.8**n for n in range(10, 15)]

x0s = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.01]
y0s = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.01]
zs = [x-1.0j*y for x in x0s for y in y0s]

ii = 1.0j
R0s = [1.6, 1.7, 1.8, 1.9, 2.0]
t0s = [-4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4]
rs = [R * np.exp(-ii*t*np.pi/180.0) for R in R0s for t in t0s]

zrs = [(z, r) for z in zs for r in rs]

subprocess.call("echo '' > search.out", shell=True)
#`echo > search.out`
for (z, r) in zrs:
    print z, r
    basis_info = [('geo', True,  2, 5,  z,   r),
                  ('geo', False, 2, 10, 0.1, 1.8)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,

        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',        
        w0 = 1.0,
        tol = pow(10.0, -5.0),
        grad=True,
        hess=True,
        
        outfile = 'res.out',
        print_level = 0)
    
    opt_res = res['w_res_list'][0][1]
    z0 = opt_res.x[0]
    z1 = opt_res.x[1]
    c0 = res['cs'][0][0]
    c1 = res['cs'][0][1]
    if(z0.real > 0.0 and z0.imag < 0.0 and
       opt_res.success):
        #        `cat res.out >> search.out`
        subprocess.call("cat res.out >> search.out", shell=True)
        

