import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## see ../search2
z0s = [0.005*1.8**n for n in range(15)]
with open("search.out", "w") as f:
    print_timestamp("start", f)

zs_list = [[ 0.000767660134519,-0.0168759624328,-0.140427371261,-0.131854426275],
           [ 0.0024561470277,-0.028652274653,-0.0136819135026,-0.0491311761337],
           [ 0.00384465965635,-0.0221516269296,-0.173656322395,-0.134760549518],
           [ 0.00571598742391,-0.0205920375373,-0.13918672031,-0.131833339129],
           [ 0.00782882844142,-0.0197790147774,-0.0725354502484,-0.0369446848077]]

for (n, zs) in zip(range(len(zs_list)), zs_list):
    basis_info = [('shift', True, 2, z0s[0:5], zs[0]+zs[1]*1.0j),
                  ('shift', True,2, z0s[5:-1], zs[2]+zs[3]*1.0j)]

    res = opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,
    
        target  = 'h_pi',
        channel = '1s->kp',
        dipole  = 'length',
    
        w0      = 1.0,
        tol     = pow(10.0, -6.0),
        maxit   = 50,
        conv    = "grad",
        fdif    = 0.0001,
        grad    = False,
        hess    = False,
    
        outfile = 'res{0}.out'.format(n),
        wf_outfile="psi{0}.csv".format(n),
        zeta_outfile="zeta{0}.csv".format(n),
        wf_rs = np.linspace(0, 40, 400),
        print_level = 5)

