import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## initial guess is from opt-2GTO
z0 = 0.164834921195-0.127041608653j
z1 = 0.905433828185-0.58333509576j
basis_info = [('id',  True, 2, z0),
              ('id',  True, 2, z1),
              ('geo', False, 2, 5, 30.0, 1.0/1.7716)]

#              ("geo", False, 2, 2, 0.01, 1.771600539669121)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,
         w0 = 1.0,
         tol = pow(10.0, -8.0),
         maxit = 50,
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 2,
         outfile = "res.out",
         wf_outfile = "wf.csv",
         wf_rs = np.linspace(0.0, 40.0, 200))
