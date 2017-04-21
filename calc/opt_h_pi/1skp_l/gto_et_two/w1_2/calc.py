import sys
import numpy as np
import pandas as pd

## from main_opt_gto.org/Results/FOWF/N-Opt-GTO + M-ET-GTO/2opt-20et

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

"""
basis_info = [('id', True, 2, 0.0016664606357-0.0585549385367j),
              ('id', True, 2, 0.100486903923-0.0753439127634j),
              ("geo", False, 2, 10, 0.01, 2.5)]

basis_info = [('id', True, 2, 0.001-0.058j),
              ('id', True, 2, 0.100-0.075j),
              ("geo", False, 2, 10, 0.01, 2.5)]

basis_info = [('id', True, 2, 0.001-0.058j),
              ('id', True, 2, 0.100-0.075j),
              ("geo", False, 2, 12, 0.01, 2.3)]

basis_info = [('id', True, 2, 0.00112465837018-0.0563029012606j),
              ('id', True, 2, 0.090112316949-0.07608690597j),
              ("geo", False, 2, 12, 0.01, 2.3)]

basis_info = [('id', True, 2, 0.00112465837018-0.0563029012606j),
              ('id', True, 2, 0.090112316949-0.07608690597j),
              ("geo", False, 2, 15, 0.01, 2.2)]
"""
basis_info = [('id', True, 2, 0.000523842185913-0.0540785551876j),
              ('id', True, 2, 0.0823920345806-0.075252339163j),
              ("geo", False, 2, 15, 0.01, 2.2)]

opt_main(basis_type = 'GTO',
         basis_info = basis_info,
         w0 = 1.0,
         tol = pow(10.0, -8.0),
         maxit = 50,
         target = 'h_pi',
         channel= '1s->kp',
         dipole = 'length',
         print_level = 2,
         conv = "grad",
         outfile = "res.out",
         wf_outfile = "psi.csv",
         wf_rs = np.linspace(0.0, 40.0, 200))
