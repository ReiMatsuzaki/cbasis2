import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")

from r1basis import *
from opt_green import *

z0s = [0.005*1.8**n for n in range(15)]
y0list = [0.522564811269-1.0071321542j,
          0.389654320343-0.75372440293j,
          0.952648756905-1.81574855625j,
          0.7166392377-1.36257823221j]

for (i, y0) in zip(range(len(y0list)), y0list):
    opt_main(
        basis_info = [('scale', True, 2, z0s, y0)],
        basis_type = 'GTO',

        w0 = 1.0,
        tol = pow(10.0, -8.0),
        
        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        print_level = 1,
        
        outfile = "res{0}.out".format(i),
        wf_outfile="psi{0}.csv".format(i),
        zeta_outfile="zeta{0}.csv".format(i),
        wf_rs = np.linspace(0, 40.0, 400))

