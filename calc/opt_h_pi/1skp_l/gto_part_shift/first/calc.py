import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

z0s = [30.0/(2**n) for n in range(10)]
z1s = [30.0/(2**n) for n in range(10,15)]
print z0s
print z1s


#y0 = 0.000939964489964-0.0340181892741j

basis_info = [('shift', False, 2, z0s, 0.0),
              ('shift', True,  2, z1s, -1.0j*0.00001)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,
         w0 = 1.0,
         tol = pow(10.0, -5.0),
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 5,
         outfile = "res.out",
         wf_outfile="psi.csv",
         wf_rs = np.linspace(0, 40.0, 400))

