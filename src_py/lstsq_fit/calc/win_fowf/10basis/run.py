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
k = 1.0
z_arc = 1.0
num_gto = 10

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

xmin_list = [5, 7] # minimum of grid points
xw_list = [20, 30, 40] # location of window function
wg_list = [10, 20, 30] # window grid edge length
zarg_list = [4.0] # shape of window

for xmin in xmin_list:
    for zarg in zarg_list:
        for xw in xw_list:
            for wg in wg_list:
                
                xmax = xw + wg # location of grid edge

                # define target function
                winfunc = lambda x: np.arctan(-zarg*(x-xw))/np.pi+0.5
                target = lambda x: winfunc(x) * (ip_re(x) + 1.0j * ip_im(x))

                # for plot
                xs = np.linspace(0.01, 60.0, 400)
                ts = [target(x) for x in xs]
                trs = [t.real for t in ts]
                tis = [t.imag for t in ts]

                # lstsq results
                res = gau_lstsq.lstsq_from_region(num_gto, float(xmin), float(xmax), pn, target)
                
                if(res):
                    print xmin, zarg, xw, wg, ": OK"
                    base = "{0}_{1}_{2}_{3}".format(xmin, xw, wg, zarg)
            
                    # basis
                    (cs, zs, other) = res
                    fs = [l2func.GTO(c, pn, -z)  for (c, z) in zip(cs, zs)]
                    l2func.write_func_list(fs, base + "_res.dat")
                
                    # fitting plot
                    fit_func = gau_lstsq.to_func(res)
                    fs = [fit_func(x) for x in xs]
                    frs = [f.real for f in fs]
                    fis = [f.imag for f in fs]
                    plt.figure(figsize=(4, 3))
                    plt.plot(xs, frs)
                    plt.plot(xs, fis)
                    plt.plot(xs, trs, "k")
                    plt.plot(xs, tis, "k--")
                    plt.savefig(base + "_fit.png")
                    plt.cla()

                    # error plot
                    plt.plot(xs, [np.sqrt((fr-tr)**2 + (fi-ti)**2)
                                  for (fr,fi,tr,ti) in zip(frs,fis,trs,tis)])
                    plt.savefig(base + "_err.png")
                    plt.cla()

                    # orbital exponents
                    zr = [-z.real for z in zs]
                    zi = [-z.imag for z in zs]
                    plt.plot(zr, zi, "o")
                    plt.savefig(base + "_exp.png")
                    plt.cla()
                else:
                    print xmin, zarg, xw, wg, res, ": Failed"
                    
