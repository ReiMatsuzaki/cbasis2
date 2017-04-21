import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

#z0s = [30.0/(1.7716**n) for n in range(17)]
z0s = [30.0/(1.8**n) for n in range(15)]


for y0 in [0.00001, 0.00002, 0.00003, 0.00004, 0.00005,
           0.0001, 0.0002, 0.0003, 0.0004, 0.0005,
           0.001, 0.002, 0.003, 0.004, 0.005,
           0.01,  0.02,  0.03,  0.04,  0.05,
           0.1]:

    label = str(y0)[2:]
    basis_info = [('shift', False, 2, z0s, -1.0j*y0)]
    cbf_main(
        basis_type = 'GTO',
        basis_info = basis_info,
        w0 = 1.0,
        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        print_level = 1,
        outfile = "res{0}.out".format(label),
        wf_outfile="psi{0}.csv".format(label),
        wf_rs = np.linspace(0, 40.0, 400))

