import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## fixed even tempered basis was used with fitted GTOs in H2plus calculation.
## orbital exponent of cGTO are from ../search2/
## one gto optimization results is 0.0105604,-0.0548991
"""
z0z1_list = [[0.000133387-0.000115202j, 0.0105642-0.0548974j],
             [8.13222e-05-0.000123463j, 0.010563-0.0549j],
             [3.02395e-05,-0.000132463j, 0.0105612-0.0549013j],
             [0.0105622-0.0549025j, 4.41225e-05-0.000161345j],
             [0.00016083-0.000126075j, 0.0105655-0.0548966j],
"""

basis_info = [('log', True, 2, np.log(0.000133387-0.000115202j)),
              ('log', True, 2, np.log(0.0105642-0.0548974j)),
              ("geo", False, 2, 15, 0.01, 1.771600539669121)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,
         w0 = 1.0,
         tol = pow(10.0, -5.0),
         maxit = 50,
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 2,
         outfile = "res.out",
         wf_outfile = "wf.csv",
         zeta_outfile = "zeta.csv",
         wf_rs = np.linspace(0.0, 40.0, 200))
