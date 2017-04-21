import sys
import numpy as np
import pandas as pd

cbasis_path = '../../../../../'
sys.path.append(cbasis_path + 'src_py/nnewton')
sys.path.append(cbasis_path + 'r1basis')
from r1basis import *
from opt_green import *

## see ../search2/
z0s = [0.01*1.8**n for n in range(15)]
y0 = 0.00325743126863-0.0569381071243j

basis_info = [('shift', True, 2, z0s, y0)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,

         target = 'h_pi',
         w0 = 1.0,
         channel= '1s->kp',
         dipole = 'length',
         
         tol = pow(10.0, -7.0),
         maxit = 50,
         conv = 'grad',
         fdif = 0.0001,
         grad = False,
         hess = False,
         
         print_level = 5,
         outfile = "res.out",
         wf_outfile="psi.csv",
         wf_rs = np.linspace(0, 40.0, 400))

