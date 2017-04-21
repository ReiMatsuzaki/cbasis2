import os
import sys

import numpy as np
import pandas as pd

sys.path.append("../../../r1basis")
from r1basis import *

L = 1
ene = 0.5

def zs_et_from_num(num, z0=0.01, z1=30.0):
    return np.array([z0 * (z1/z0)**(n*1.0/(num-1)) for n in range(num)])

def calc():
    ## optimized basis set
    #basis_dir = os.path.expanduser('~/src/git/opt_cbf/lstsq_fit/calc/win_irr_coulomb_atan2/p_ene05/out/')
    #base = "5_30"    
    #cgto_csv = pd.read_csv(basis_dir+base + ".res.csv")
    #pns_c= cgto_csv["pn"].as_matrix()
    #zs_c = (cgto_csv["re_z"].as_matrix()
    #+ 1.0j * cgto_csv["im_z"].as_matrix())
    zs_c = [0.0812755955262-0.0613237786222j,
	    0.00497147387796-0.0113737972763j,
	    0.0323712622673-0.0451300037076j,
	    0.00317417887792-0.022582254987j,
	    0.0118391719646-0.0327847352576j]

    ## build basis functions
    zs = list(zs_et_from_num(15, 0.01, 30.0)) + list(zs_c)
    b = GTOs()
    b.add(2, zs)
    b.setup()
    
    ## velocity driven term
    driv = LC_STOs()
    driv.add(2.0, 1, 1.0)
    
    ## matrix and vector
    lmat = ( b.calc_rm_mat(0)  * ene
            + b.calc_d2_mat()* 0.5
            + b.calc_rm_mat(-2)* (-0.5*(-L*(L+1)))
            + b.calc_rm_mat(-1))
    mvec = ss.calc_vec(driv)
    cs = solve(lmat, mvec)
    alpha = np.dot(cs, mvec)
    
    print alpha

calc()
