import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

z0s = [30.0/(1.7716**n) for n in range(15)]
y0 = 0.000939964489964-0.0340181892741j

basis_info = [('shift', True, 2, z0s, -1.0j*y0)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,
         w0 = 1.0,
         tol = pow(10.0, -8.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 1,
         outfile = "res.out",
         wf_outfile="psi.csv",
         wf_rs = np.linspace(0, 40.0, 400))

