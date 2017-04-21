import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *


z0s_list = [[0.0361960723459-0.0271313771956j, 0.148847470321-0.145461004611j],
            [0.164834921195-0.127041608653j,   0.905433828185-0.58333509576j]]
num = len(z0s_list)
for (n, z0s) in zip(range(num), z0s_list):
    res = opt_main(
        basis_type = 'GTO',
        basis_info = [('id', True, 2, z0s[0]), ('id', True, 2, z0s[1])],
        
        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        w0 = 1.0,

        tol = pow(10.0, -8.0),
        conv="grad",
        
        print_level = 5,
        outfile = "res{0}.out".format(n),
        zeta_outfile="zeta{0}.csv".format(n),
        wf_outfile="psi{0}.csv".format(n),
        wf_rs = np.linspace(0, 40.0, 400))


