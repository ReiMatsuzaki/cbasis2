import sys
import numpy as np
import pandas as pd

cbasis_path = '../../../../../'
sys.path.append(cbasis_path + 'src_py/nnewton')
sys.path.append(cbasis_path + 'r1basis')
from r1basis import *
from opt_green import *

## see ../search_5gto
## see ../search5 etc

z0list = [0.001, 0.002, 0.005, 0.01]
z1 = 20.0
num = 15
y0 = 0.006985-0.032446838j
for z0 in z0list:
    r = (z1/z0)**(1.0/(num-1))
    basis_info = [('geo', False, 2, num, z0, r),
                  ('id' , True,  2, y0)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
        
        target = 'h_pi',
        w0 = 1.0,
        channel= '1s->kp',
        dipole = 'length',
             
        tol = pow(10.0, -5.0),
        maxit = 50,
        conv = 'grad',
        fdif = 0.0001,
        grad = False,
        hess = False,
             
        print_level = 5,
        outfile = "res{0}.out".format(z0),
        zeta_outfile="zeta{0}.csv".format(z0),
        wf_outfile="psi{0}.csv".format(z0),
        wf_rs = np.linspace(0, 40.0, 400))
    
    print res['w_res_list'][0][1].success
    y0 = res['w_res_list'][0][1].x[0]


