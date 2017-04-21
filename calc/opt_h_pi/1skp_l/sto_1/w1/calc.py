import sys
import numpy as np
import pandas as pd
sys.path.append("../../../../../src_py/driv_grid")
sys.path.append("../../../../../r1basis")
from driv_grid import solve_driv
from r1basis import *
from opt_green import *

channel = '1s->kp'
dipole = 'length'

opt_main(basis_type = 'STO',
         basis_info = [('id', True, 2, 0.6-0.6j)],
         w0 = 1.0,
         tol = pow(10.0, -5.0),
         target = 'h_pi',
         channel= channel,
         dipole = dipole,
         print_level = 2,
         outfile = "res.out",
         wf_outfile = "wf.csv",
         wf_rs = np.linspace(0, 40.0, 200))

