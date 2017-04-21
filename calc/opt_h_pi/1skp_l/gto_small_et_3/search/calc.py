import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

with open("search.out", "w") as f:
    print_timestamp("start", f)

xos = [0.001, 0.003,  0.007,
       0.01,  0.03,   0.07,
       0.1,   0.3,    0.7,
       1.0,   3.0]

yos = [0.001, 0.003, 0.007,
       0.01,  0.03,  0.07,
       0.1,   0.3,   0.5]
zos = [x-1.0j*y for (x,y) in product(xos, yos)]

f = open('search.csv', 'w')
f.write("re_0,im_0,re_1,im_1\n")
for (z0, z1, z2) in combinations(zos, 3):
    basis_info = [('id',  True,  2, z0),
                  ('id',  True,  2, z1),
                  ('id',  True,  2, z2),
                  ('geo', False, 2, 5, 20.0, 1.0/1.8)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,

        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',        
        
        w0 = 1.0,
        
        tol = pow(10.0, -5.0),
        maxit = 50,
        conv = "grad",

        grad = True,
        hess = True,
        outfile = 'res.out',
        print_level = 0)
    
    opt_res = res['w_res_list'][0][1]
    z0 = opt_res.x[0]
    z1 = opt_res.x[1]
    z2 = opt_res.x[2]
    eps = 0.00001
    if(z0.real > 0.0 and z0.imag < 0.0 and
       100.0 > z0.real and z0.imag > -100.0 and
       z1.real > 0.0 and z1.imag < 0.0 and
       100.0 > z1.real and z1.imag > -100.0 and
       z2.real > 0.0 and z2.imag < 0.0 and
       100.0 > z2.real and z2.imag > -100.0 and
       abs(z0-z1) > eps and abs(z2-z1) > eps and abs(z2-z1) > eps and
       opt_res.success):
       subprocess.call("cat res.out >> search.out", shell=True)
       f.write("{0},{1},{2},{3}\n".format(z0.real, z0.imag, z1.real, z1.imag))

f.close()
with open("search.out", "a") as f:
    print_timestamp("end", f)
