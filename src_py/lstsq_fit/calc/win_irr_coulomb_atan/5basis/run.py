import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def append_on_git(dirname):
    dn = os.environ['HOME'] + "/src/git/" + dirname
    sys.path.append(dn)

append_on_git("opt_cbf/coulomb_func")
import coulomb

append_on_git("opt_cbf/lstsq_fit")
import gau_lstsq

append_on_git("l2func/py_bind")
import l2func

os.chdir("out1")

L = 1
pn = L+1
k = 1.0
z = 1.0
eta = -z/k

win_x1 = 50.0
z_arc = 1.0

numGto = 5

non_win = lambda x: coulomb.outgoing(eta, x, L, 0)
xs = np.linspace(0.01, 60.0, 400)

def lstsq_from_region(xmin, xmax, target, m):
    w = xmin*xmin
    num = int(xmax*xmax/w)
    try:
        return gau_lstsq.lstsq(m = m,
                               n = num,
                               w = w,
                               pn = pn,
                               f = target)
    except:
        return None

for xmin in [2, 3, 4, 5, 7]:
    for xmax in [50]:
        for win_x1 in [20, 30, 40]:

            # target
            winfunc = lambda x: np.arctan(-z_arc*(x-win_x1))/np.pi+0.5
            target = lambda x:winfunc(x)*non_win(x)
            ts = np.array([target(x) for x in xs])
            trs = [t.real for t in ts]
            tis = [t.imag for t in ts]
            
            res = lstsq_from_region(float(xmin), float(xmax), target, numGto)
            if(res==None):
                print(xmin, xmax, win_x1), "x"
            if(res):
                print(xmin, xmax, win_x1), "ok"
                base = "{0}_{1}_{2}".format(xmin, xmax, win_x1)
            
                # basis
                (cs, zs, other) = res
                fs = [l2func.GTO(c, pn, -z)  for (c, z) in zip(cs, zs)]
                l2func.write_func_list(fs, base + "_res.dat")

                fit_func = gau_lstsq.to_func(res)
                fs = [fit_func(x) for x in xs]
                frs = [f.real for f in fs]
                fis = [f.imag for f in fs]
                plt.figure(figsize=(4, 3))
                plt.plot(xs, frs)
                plt.plot(xs, fis)
                plt.plot(xs, trs, "k")
                plt.plot(xs, tis, "k--")
                plt.ylim(-1.1, 1.1)
                plt.savefig(base + "_fit.png")
                plt.clf()

                # orbital exponents
                zr = [-z.real for z in zs]
                zi = [-z.imag for z in zs]
                plt.plot(zr, zi, "o")
                plt.savefig(base + "_exp.png")
                plt.clf()
