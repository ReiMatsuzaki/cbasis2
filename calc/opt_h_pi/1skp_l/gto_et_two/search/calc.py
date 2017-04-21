import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *
#xos = [0.0001, 0.0002, 0.0003, 0.0005, 0.0007,
#       0.001,  0.002,  0.003,  0.005,  0.007]
#yos = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.3]
#xos = [0.001, 0.009]
#yos = [0.001, 0.03]
zos = [x-1.0j*y for (x,y) in product(xos, yos)]

subprocess.call("echo '' > search.out", shell=True)
#`echo > search.out`
for (z0, z1) in combinations(zos, 2):
    basis_info = [('id',  True,  2, z0),
                  ('id',  True,  2, z1),
                  ('geo', False, 2, 15, 0.01, 1.771600539669121)]

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
    z0 = opt_res.x[0]
    z1 = opt_res.x[1]
    c0 = res['cs'][0][0]
    c1 = res['cs'][0][1]
    if(z0.real > 0.0 and z0.imag < 0.0 and
       opt_res.success):
        #        `cat res.out >> search.out`
        subprocess.call("cat res.out >> search.out", shell=True)
        

