import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

subprocess.call("echo '' > search.out", shell=True)
#`echo > search.out`

#basis_info = [('geo', True, 2, 17, 0.01-0.0000001j, 1.771600539669121)]
basis_info = [('geo', True, 2, 10, 0.0022-0.005796j, 1.865+0.08j)]

res = opt_main(
    basis_type = 'GTO',
    basis_info = basis_info,
    w0 = 1.0,
    maxit = 100,
    tol = pow(10.0, -5.0),
    target = 'h_pi',
    channel= '1s->kp',
    dipole = 'length',
    outfile = 'res.out',
    wf_outfile = "psi.csv",
    wf_rs = np.linspace(0, 40, 400),
    print_level = 5)
