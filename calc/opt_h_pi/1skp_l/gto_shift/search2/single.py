import sys
import numpy as np
import pandas as pd
from itertools import product, combinations
import subprocess

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## one GTO optimization value is 0.0105604,-0.0548991


y0 = 0.00356251319821-0.032168804755j
f = open('search.csv', 'w')
z0s = [30.0/(1.8**n) for n in range(15)]
basis_info = [('shift', True, 2, z0s, -1.0j*y0)]
    
res = opt_main(
    basis_type = 'GTO',
    basis_info = basis_info,
    w0      = 1.0,
    tol     = pow(10.0, -5.0),
    target  = 'h_pi',
    channel = '1s->kp',
    dipole  = 'length',
    outfile = 'res.out',
    print_level = 0)
