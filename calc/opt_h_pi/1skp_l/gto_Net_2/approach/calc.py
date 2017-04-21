import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## initial guess is from ../N_increase/calc.py
z0 = 0.164466262866-0.126713448493j
z1 = 0.869882039309-0.52998454715j
ratio = 1.0/1.7716

f = open('opt.out', 'w')

for etz0 in [30000.0, 10000.0, 5000.0, 4000.0, 3500.0, 3300.0, 3200.0, 3180.0,
             3170.0, 3165.0]:
    
    basis_info = [('id',  True, 2, z0),
                  ('id',  True, 2, z1),
                  ('geo', False, 2, 15, etz0, ratio)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
        w0 = 1.0,
        tol = pow(10.0, -5.0),
        maxit = 100,
        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        print_level = 2,
        outfile = "res.out",
        wf_outfile = "wf.csv",
        wf_rs = np.linspace(0.0, 40.0, 200))

    opt_res = res['w_res_list'][0][1]
    z0 = opt_res.x[0]
    z1 = opt_res.x[1]
    print >>f, etz0, opt_res.success, z0, z1
    
f.close()    
