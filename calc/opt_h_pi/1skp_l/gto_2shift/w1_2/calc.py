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
z0s = [0.005*1.8**n for n in range(15)]
with open("search.out", "w") as f:
    print_timestamp("start", f)

data = [(0, [0.00182904262031-0.0331253544721j, 0.0224517632247-0.0286633479976j]),
        (1, [0.00184604805831-0.0396484835939j, 0.0327547299986-0.0432288858114j])]

for (n, zs) in data:
    basis_info = [('shift', True, 2, z0s[0:2], zs[0]),
                  ('shift', True, 2, z0s[2:5], zs[1]),
                  ('shift', False,2, z0s[5:-1], 0.0)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
    
        target  = 'h_pi',
        channel = '1s->kp',
        dipole  = 'length',
    
        w0      = 1.0,
        tol     = pow(10.0, -6.0),
        maxit   = 50,
        conv    = "grad",
        fdif    = 0.0001,
        grad    = False,
        hess    = False,
    
        outfile = 'res{0}.out'.format(n),
        wf_outfile="psi{0}.csv".format(n),
        zeta_outfile="zeta{0}.csv".format(n),
        wf_rs = np.linspace(0, 40, 400),
        print_level = 5)

