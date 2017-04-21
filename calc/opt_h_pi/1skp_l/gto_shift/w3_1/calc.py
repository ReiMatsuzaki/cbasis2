import sys
import numpy as np
import pandas as pd

cbasis_path = '../../../../../'
sys.path.append(cbasis_path + 'src_py/nnewton')
sys.path.append(cbasis_path + 'r1basis')
from r1basis import *
from opt_green import *

z0s = [0.005*1.8**n for n in range(15)]

## see ../search3/
y0_list = [0.006985-0.032446838j,
           0.00639157778689-0.0273146878907j,
           0.00460225081889-0.0432741394865j]

for (i, y0) in zip(range(len(y0_list)), y0_list):
    basis_info = [('shift', True, 2, z0s, y0)]
    opt_main(basis_type = 'GTO',
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
             outfile = "res{0}.out".format(i),
             wf_outfile="psi{0}.csv".format(i),
             wf_rs = np.linspace(0, 40.0, 400))

