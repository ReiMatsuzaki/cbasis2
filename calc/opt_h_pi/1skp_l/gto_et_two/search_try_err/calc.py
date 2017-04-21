import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## one GTO optimization value is 0.0105604,-0.0548991

z0 = 0.002-0.1j
z1 = 0.05-0.54j

basis_info = [('id',  True,  2, z0),
              ('id',  True,  2, z1),
              ('geo', False, 2, 15, 0.01, 1.771600539669121)]

res = opt_main(
    basis_type = 'GTO',
    basis_info = basis_info,
    w0 = 1.0,
    tol = pow(10.0, -5.0),
    target = 'h_pi',
    channel= '1s->kp',
    dipole = 'length',
    print_level = 0)
    
