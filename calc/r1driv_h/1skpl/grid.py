import sys
sys.path.append("../../../src_py/driv_grid")
from driv_grid import solve_driv
import numpy as np
import pandas as pd
#sys.path.append("../../../r1basis")
#from opt_green import *

#grid_csv = "grid.csv"
#channel = '1s->kp'
#dipole = 'length'
#h_pi = H_Photoionization(channel, dipole)
L = 1
w0 = 1.0
E0 = -0.5
n = 1000
h = 0.1

s = lambda r: 2.0*r*r*np.exp(-r)
v = lambda r: L*(L+1)/(2.0*r*r) -1.0/r
(rs, ys) = solve_driv(v, w0+E0, s, n, h)
df = pd.DataFrame([rs.real, -ys.real, -ys.imag]).T
df.columns = ["r", "re_y", "im_y"]
df.to_csv("grid.csv", index=False)

