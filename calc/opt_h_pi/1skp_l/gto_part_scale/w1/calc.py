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
y0list = [
    2.62213491278-5.75285174243j,
    0.838192881435-1.14147984663j,
    1.03788437358-1.64836990704j]
num = len(y0list)
for (i, y0) in zip(range(n), y0list):

    basis_info = [('scale', True,  2, z0s, y0),
                  ('scale', False, 2, z1s, 1.0)]

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
        grad=False,
        hess=False,

        print_level = 5,
        outfile = "res{0}.out".format(i),
        wf_outfile="psi{0}.csv".format(i),
        wf_rs = np.linspace(0, 40.0, 400))

    
