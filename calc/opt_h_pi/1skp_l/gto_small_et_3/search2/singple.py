import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

z0s_list = [[0.0230447241786-0.0497180498056j, 0.125250955338-0.0974602248987j]]
num = len(z0s_list)
for (n, z0s) in zip(range(num), z0s_list):
    basis_info = [('id',  True,  2, z0s[0]),
                  ('id',  True,  2, z0s[1]),
                  ('geo', False, 2, 5, 20.0, 1.0/2.5)]

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

        grad = True,
        hess = True,
        outfile = 'res{0}.out'.format(n),
        wf_outfile = "psi{0}.csv".format(n),
        wf_rs = np.linspace(0, 40.0, 401),
        print_level = 0)
    
    print res['w_res_list'][0][1].success
