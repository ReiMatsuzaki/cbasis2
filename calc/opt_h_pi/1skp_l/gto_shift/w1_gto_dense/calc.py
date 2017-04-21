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
z0 = 0.005
z1 = 20.0
num_y0 = [(5,  0.0496525882029-0.0758265304904j),
          (10, 0.006985-0.032446838j),
          (15, 0.006985-0.032446838j),
          (20, 0.00255710899054-0.0213505048018j)]

for (num, y0) in num_y0:
    r = (z1/z0)**(1.0/(num-1))
    z0s = [z0 * r**n for n in range(num)]
    #z0s = [0.005*1.8**n for n in range(15)]

    basis_info = [('shift', True, 2, z0s, y0)]

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
        outfile = "res{0}.out".format(num),
        wf_outfile="psi{0}.csv".format(num),
        zeta_outfile="zeta{0}.csv".format(num),
        wf_rs = np.linspace(0, 40.0, 400))
    
    print res['w_res_list'][0][1].success
    y0 = res['w_res_list'][0][1].x[0]


