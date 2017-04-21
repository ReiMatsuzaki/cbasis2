import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

numopt = 5
zs = [0.005*1.8**n for n in range(15)]
z0s = zs[0:numopt]
z1s = zs[numopt:]
y0list = [0.0121912652483,-0.0439556256182,
          0.00910004385327,-0.0308075005003,
          0.00731812496601,-0.0249479552054]

for (i, y0) in zip(range(3), y0list):

    basis_info = [('shift', True,  2, z0s, -1.0j*y0),
                  ('shift', False, 2, z1s, 0.0)]

    opt_main(
        basis_type = 'GTO',
        basis_info = basis_info,

        target = 'h_pi',
        channel= '1s->kp',
        dipole = 'length',
        w0 = 1.0,

        tol = pow(10.0, -5.0),
        maxit = 50,
        conv = 'grad',
        fdif = 0.0001,
        grad = False,
        hess = False,

         print_level = 5,
         outfile = "res{0}.out".format(i),
         wf_outfile="psi{0}.csv".format(i),
         wf_rs = np.linspace(0, 40.0, 400))

    
