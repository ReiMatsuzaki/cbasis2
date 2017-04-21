import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as ip

def append_on_git(dirname):
    dn = os.environ['HOME'] + "/src/git/" + dirname
    sys.path.append(dn)

append_on_git("fowf_h_pi")
import fowf_h_pi as fowf

append_on_git("opt_cbf/lstsq_fit")
import gau_lstsq
import window_func

append_on_git("l2func/py_bind")
import l2func

os.chdir("out")

L = 1
pn = L+1
ks = [0.5, 1.0, 1.5]
z_arc = 1.0
x1 = 50.0
num_gto = 14

xmin = 5
xmax = 100

for k in ks:

    # numerical solution for FOWF
    grid = fowf.Grid(0.1, 2000)
    h_pi = fowf.HydrogenPI(k*k*0.5, 0, 1)
    driv_sys = fowf.DrivSys(grid)
    driv_sys.add_laplacian(h_pi.ene())
    driv_sys.add_mat_hydrogen_pi(h_pi)
    driv_sys.add_driv_1skp_length()
    wave_func = driv_sys.solve()
    res_data= fowf.calculate_fowf(grid, h_pi)

    # interpolation
    ip_re = ip.interp1d(res_data.T[0], res_data.T[1])
    ip_im = ip.interp1d(res_data.T[0], res_data.T[2])

    # define target function
    winfunc = lambda x: np.arctan(-z_arc*(x-x1))/np.pi+0.5
    target = lambda x: winfunc(x) * (ip_re(x) + 1.0j * ip_im(x))

    # for plot
    xs = np.linspace(0.01, 60.0, 400)
    ts = [target(x) for x in xs]
    trs = [t.real for t in ts]
    tis = [t.imag for t in ts]

    res = gau_lstsq.lstsq_from_region(num_gto, float(xmin), float(xmax), pn, 
                                      target)
    if(res):
        base = "{0}".format(k)
        
        # basis
        (cs, zs, other) = res
        fs = [l2func.GTO(c, pn, -z)  for (c, z) in zip(cs, zs)]
        l2func.write_func_list(fs, base + "_res.dat")

        # fitting plot
        fit_func = gau_lstsq.to_func(res)
        fs = [fit_func(x) for x in xs]
        frs = [f.real for f in fs]
        fis = [f.imag for f in fs]
        plt.plot(xs, frs)
        plt.plot(xs, fis)
        plt.plot(xs, trs, "k")
        plt.plot(xs, tis, "k--")
        plt.savefig(base + "_fit.png")
        plt.clf()
        
        # orbital exponents
        zr = [-z.real for z in zs]
        zi = [-z.imag for z in zs]
        plt.plot(zr, zi, "o")
        plt.savefig(base + "_exp.png")
        plt.clf()

    
