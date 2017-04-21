import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

basis_info = [('id',  True, 2, 0.04-0.06j),
              ('id',  True, 2, 0.1-0.005j),
              ("geo", False, 2, 15, 0.005, 1.8)]
        
res = opt_main(
    basis_type = 'GTO',
    basis_info = basis_info,
    
    target = 'h_pi',
    channel= '1s->kp',
    dipole = 'length',         
    w0 = 1.0,
    
    tol = pow(10.0, -5.0),
    maxit = 50,
    conv = "grad",
    
    print_level = 2,
    outfile = "res.out"
    #         wf_outfile = "psi.csv",
    #         zeta_outfile = "zeta.csv",
    #    wf_rs = np.linspace(0.0, 40.0, 200)
)
