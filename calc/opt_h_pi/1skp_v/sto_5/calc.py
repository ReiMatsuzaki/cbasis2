import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../src_py/nnewton")
sys.path.append("../../../../r1basis")
from r1basis import *
from opt_green import *

opt_main(basis_type = 'STO',
         basis_info = [('id', True, 2, 0.65-0.43j)],
         w0 = 0.6,
         ws = list(np.linspace(0.6, 1.6, 11)),
         tol = pow(10.0, -8.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'velocity',
         print_level = 1,
         outfile = "res.out")

