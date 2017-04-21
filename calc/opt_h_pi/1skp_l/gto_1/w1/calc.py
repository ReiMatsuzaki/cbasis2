import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

opt_main(basis_type = 'GTO',
         basis_info = [('id', True, 2, 0.15-0.13j)],
         w0 = 1.0,
         tol = pow(10.0, -5.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 1,
         wf_outfile = "psi.csv",
         zeta_outfile="zeta.csv",
         wf_rs = np.linspace(0, 40, 400),
         outfile = "res.out")

