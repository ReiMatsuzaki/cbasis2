import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

opt_main(basis_type = 'STO',
         basis_info = [("id", True, 2, 0.6-0.6j)],
         w0 = 1.0,
         ws = np.linspace(0.5, 1.5, 11),
         tol = pow(10.0, -5.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 1,
         outfile = "res.out")

