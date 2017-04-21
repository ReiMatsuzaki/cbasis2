import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

# one opt result is (0.0105604067821-0.0548990693108j)

#z0 = [0.01-0.05j, 0.199626942435-0.132318775251j]
#z0 = [0.01-0.055j, 0.002-0.002j]
#z0 = [0.005-0.02j, 0.002-0.002j]

basis_info = [('id',  False, 2, 0.005-0.02j),
              ('id',  True,  2, 0.005-0.005j),
              ("geo", False, 2, 15, 0.005, 1.8)]

res = opt_main(
    basis_type = 'GTO',
    basis_info = basis_info,
    
    target = 'h_pi',
    channel= '1s->kp',
    dipole = 'length',         
    w0 = 1.0,
    
    tol = pow(10.0, -5.0),
    maxit = 200,
    conv = "grad",
    
    print_level = 2,
    outfile = "res.out",
    wf_outfile = "psi.csv",
    zeta_outfile = "zeta.csv",
    wf_rs = np.linspace(0.0, 40.0, 200)
)


