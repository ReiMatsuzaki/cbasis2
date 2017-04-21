import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## see ../search2

with open("search.out", "w") as f:
    print_timestamp("start", f)
    
#z0s = [0.005, 0.0048, 0.0046, 0.00458] + [0.00458 + n*0.00001 for n in range(10)]
z0s = [0.005 - n*0.0001 for n in range(40)]
y0s = [ 0.0024561470277-0.028652274653j,-0.0136819135026-0.0491311761337j]
for z0 in z0s:
    z0s = [z0*1.8**n for n in range(15)]
    basis_info = [('shift', True, 2, z0s[0:5], y0s[0]),
                  ('shift', True,2, z0s[5:-1], y0s[1])]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
    
        target  = 'h_pi',
        channel = '1s->kp',
        dipole  = 'length',
    
        w0      = 1.0,
        tol     = pow(10.0, -5.0),
        maxit   = 50,
        conv    = "grad",
        fdif    = 0.0001,
        grad    = False,
        hess    = False,
    
        outfile = 'res{0}.out'.format(n),
        wf_outfile="psi{0}.csv".format(n),
        wf_rs = np.linspace(0, 40, 400),
        print_level = 5)
    res_opt = res["w_res_list"][0][1]
    print z0, res_opt.success, res_opt.x
    y0s = res_opt.x
