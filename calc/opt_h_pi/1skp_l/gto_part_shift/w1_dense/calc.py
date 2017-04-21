import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

numopt = 5

z0 = 0.005
z1 = 20.0
num_y0 = [(5, -0.004j),
          (10, -0.004j),
          (15, -0.004j),
          (20, 0.005-0.02j)]
for (num, y0) in num_y0:
    r = (z1/z0)**(1.0/(num-1))
    z0s = [z0*r**n for n in range(num)]
    basis_info = [('shift', True,  2, z0s[0:5], y0),
                  ('shift', False, 2, z0s[5:-1], 0.0)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,

        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        w0 = 1.0,

        tol = pow(10.0, -5.0),
        maxit = 50,
        conv = 'grad',
        fdif = 0.0001,
        grad = False,
        hess = False,

        print_level = 5,
        outfile = "res{0}.out".format(num),
        wf_outfile="psi{0}.csv".format(num),
        zeta_outfile="zeta{0}.csv".format(num),
        wf_rs = np.linspace(0, 40.0, 400))
    
    print res['w_res_list'][0][1].success
    y0 = res['w_res_list'][0][1].x[0]

    
