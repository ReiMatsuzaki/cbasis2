import sys

import numpy as np
import pandas as pd
from scipy.linalg import solve

sys.path.append("../../../r1basis")
from r1basis import *

sys.path.append("../../../src_py/driv_grid")
from driv_grid import solve_driv

L = 1
ene = 0.5
n = 1000   ## grid num
h = 0.1    ## grid width
psi_csv = "psi.csv"
grid_csv= "grid.csv"
outname = "res.out"
f = open(outname, "w")

def zs_et_from_num(num, z0=0.01, z1=30.0):
    return np.array([z0 * (z1/z0)**(n*1.0/(num-1)) for n in range(num)])

def calc():
    ## optimized basis set
    b = STOs()
    b.add(2, 0.9965751177 -0.0013743026j)
    b.add(2, 1.0030366528 -0.2836728004j)
    b.add(2, 0.8462928140 -0.6952686244j)
    b.add(2, 0.4818046345 -1.0023929406j)
    b.add(2, 0.1412093744 -1.0662761427j)
    b.setup()
    
    ## velocity driven term
    driv = LC_STOs()
    driv.add(2.0, 2, 1.0)

    ## print claculation inputs
    f.write("calculation of hydrogen atom photoionization in length form\n")
    f.write("solve radial equation\n")
    f.write("(H-E)psi=driv\n")
    f.write("\n")
    f.write(">>> INPUTS >>> \n")
    f.write("L = {0}\n".format(L))
    f.write("energy = {0}\n".format(ene))
    f.write("driv = {0}\n".format(driv.str()))
    f.write("basis = \n")
    f.write("grid = ({0}, {1})".format(n, h))
    f.write(b.str())
        
    ## calculation
    lmat = ( b.calc_d2_mat()* (-0.5)
             + b.calc_rm_mat(-2)* (0.5*(L*(L+1)))
             + b.calc_rm_mat(-1) * (-1.0)
             + b.calc_rm_mat(0)  * (-ene))
    mvec = b.calc_vec(driv)
    cs = solve(lmat, mvec)
    alpha = np.dot(cs, mvec)

    ## calculate by CBF
    rs = np.linspace(0, 100, 300)
    ys = np.array(b.at_r(rs, cs))
    df = pd.DataFrame([rs, ys.real, ys.imag]).T
    df.columns = ["r", "re_y", "im_y"]
    df.to_csv(psi_csv, index=False)

    ## calculate by grid
    s = lambda r: driv.at_r([r])[0]
    v = lambda r: L*(L+1)/(2.0*r*r) -1.0/r
    (rs, ys) = solve_driv(v, ene, s, n, h)
    df = pd.DataFrame([rs.real, ys.real, ys.imag]).T
    df.columns = ["r", "re_y", "im_y"]
    df.to_csv(grid_csv, index=False)

    ## print outpus
    f.write("\n>>> RESULTS >>> \n")
    f.write("psi_csv = {0}\n".format(psi_csv))
    f.write("grid_csv = {0}\n".format(grid_csv))
    f.write("cs = \n")
    for c in cs:
        f.write("{0}\n".format(c))
    f.write("alpha = {0} {1}\n".format(alpha.real, alpha.imag))

calc()
f.close()
