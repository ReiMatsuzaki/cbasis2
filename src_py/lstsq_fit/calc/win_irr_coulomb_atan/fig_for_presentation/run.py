import os
import os.path
import sys
import numpy as np
import matplotlib.pyplot as plt

def append_on_git(dirname):
    dn = os.environ['HOME'] + "/src/git/" + dirname
    sys.path.append(dn)

append_on_git("opt_cbf/coulomb_func")
import coulomb

out_dir = "out"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

os.chdir(out_dir)

L = 1
pn = L+1
k = 1.0
z = 1.0
eta = -z/k
win_x1 = 30.0
z_arc = 1.0

non_win = lambda x: coulomb.outgoing(eta, x, L, 0)
winfunc = lambda x: np.arctan(-z_arc*(x-win_x1))/np.pi+0.5

xs = np.linspace(0.01, 50.0, 400)

ys = np.array([non_win(x) for x in xs])
plt.plot(xs, ys.real, "r", label="real part")
plt.plot(xs, ys.imag, "b", label="imaginary part")
plt.xlim(-2, 50.0)
plt.ylim(-2, 2)
plt.legend()
plt.savefig("outgoing_coulomb.png")
plt.clf()

ys = np.array([non_win(x)*winfunc(x) for x in xs])
plt.plot(xs, ys.real, "r", label="real part")
plt.plot(xs, ys.imag, "b", label="imaginary part")
ys_win = np.array([winfunc(x) for x in xs])
plt.plot(xs, ys_win, "k--", label="window function")
plt.xlim(-2, 50.0)
plt.ylim(-2, 2)
plt.legend()
plt.savefig("win_outgoing_coulomb.png")
plt.clf()
