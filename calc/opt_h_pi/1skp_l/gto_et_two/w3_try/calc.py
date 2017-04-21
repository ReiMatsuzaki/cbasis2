import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../src_py/nnewton")
sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

# one opt result is (0.0105604067821-0.0548990693108j)

#z0 = [0.01-0.05j, 0.199626942435-0.132318775251j]
#z0 = [0.01-0.055j, 0.002-0.002j]
#z0 = [0.005-0.02j, 0.002-0.002j]

fdata = open("data.dat", "w")
print >>fdata, "re", "im"
x0s = np.linspace(0.01, 0.1, 21)
y0s = np.linspace(-0.1, -0.01, 41)
i = 0
z0_y0 = 0.005-0.005j
for y0 in y0s:
    print "change"
    z0 = z0_y0
    for x0 in x0s:
        basis_info = [('id',  False, 2, x0+1.0j*y0),
                      ('id',  True,  2, z0),
                      ("geo", False, 2, 15, 0.005, 1.8)]
        
        res = opt_main(
            basis_type = 'GTO',
            basis_info = basis_info,

            target = 'h_pi',
            channel= '1s->kp',
            dipole = 'length',         
            w0 = 1.0,
            
            tol = pow(10.0, -5.0),
            maxit = 200,
            conv = "grad",
            
            print_level = 2,
            outfile = "res.out",
            #         wf_outfile = "psi.csv",
            #         zeta_outfile = "zeta.csv",
            #    wf_rs = np.linspace(0.0, 40.0, 200)
        )
        res_opt = res["w_res_list"][0][1]
        z0 = res_opt.x[0]
        if(x0==x0s[0]):
            z0_y0 = z0
            
        print i, 1.0*i/(len(x0s)*len(y0s)) * 100 ,"%", x0, y0
        if(not res_opt.success):
            print "convergence failed"
            exit(1)
            
        i = i+1
        print >>fdata, res_opt.val.real, res_opt.val.imag


