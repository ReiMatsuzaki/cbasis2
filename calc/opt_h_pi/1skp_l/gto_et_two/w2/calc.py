import sys
import numpy as np
import pandas as pd

## from main_opt_gto.org/Results/FOWF/N-Opt-GTO + M-ET-GTO/2opt-20et

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

z0s_list = [
    [0.0016675881-0.058558696j,          0.1004093981-0.075346307j],
    [0.000523842185913-0.0540785551876j, 0.0823920345806-0.075252339163j],
    [0.00218864558775-0.0601584588628j,  0.103595012657-0.0762903003474j]]

num = len(z0s_list)
for (n, z0s) in zip(range(num), z0s_list):
    basis_info = [('id', True, 2, z0s[0]),
                  ('id', True, 2, z0s[1]),
                  ("geo", False, 2, 20, 2.5**(-12), 2.5)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
        w0 = 1.0,
        tol = pow(10.0, -8.0),
        maxit = 50,
        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        print_level = 5,
        conv = 'dx',
        outfile = "res{0}.out".format(n),
        wf_outfile = "psi{0}.csv".format(n),
        zeta_outfile = "zeta{0}.csv".format(n),
        wf_rs = np.linspace(0.0, 40.0, 200))
    print res["w_res_list"][0][1].success
