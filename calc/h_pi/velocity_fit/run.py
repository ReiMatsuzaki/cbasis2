import sys
sys.path.append("../..")
from local import *

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import json


from symmolint import *
from r1gtoint import *
from minieigen import *
import grid_method
from cbasis import zs_from_rmax

import os.path

def plot_num():
    plt.plot(rs_g, ys_g.real, "k--", label="numerical")

def plot_psi(zs, zs2, lbl):
    gtos = R1GTOs()
    gtos.add(2, zs)
    gtos.add(2, zs2)
    gtos.normalize()
    cs = solve_alpha(driv, op, gtos)
    rs = rs_g
    ys = gtos.at_r(cs, rs)
    plt.plot(rs, ys.real, label=lbl)

def zs_et_from_num(num, z0=0.01, z1=30.0):
    return np.array([z0 * (z1/z0)**(n*1.0/(num-1)) for n in range(num)])

basis_dir = os.path.expanduser('~/src/git/opt_cbf/lstsq_fit/calc/win_irr_coulomb_atan2/p_ene05/out/')
base = "5_30"

## ==== Setting ====
setting_json = json.load(open(basis_dir + base + ".setting.json", "rb"))
ene = setting_json["E"]
L = setting_json["L"]

## ==== numerical solution ====
grid_fn = os.path.expanduser("~/src/git/opt_cbf/py_bind/grid_method/calc/h_velocity/h_velo.csv")
grid_csv = pd.read_csv(grid_fn)

## ==== Basis ======
cgto_csv = pd.read_csv(basis_dir + base + ".res.csv")
pns_c= cgto_csv["pn"].as_matrix()
zs_c = cgto_csv["re_z"].as_matrix() + 1.0j * cgto_csv["im_z"].as_matrix()
zs_et = zs_et_from_num(15, 0.01, 30.0)

gtos = R1GTOs()
for (pn,z) in zip(pns_c, zs_c):
    gtos.add(pn, z)
gtos.add(L+1, zs_et)
gtos.normalize()

## ==== Solve driven equation ====
sto_driv = R1STOs(); sto_driv.add(2.0, 1, 1.0);
driv = DrivSTO(sto_driv)
op = OpCoulomb(L, ene)

cs = solve_alpha(driv, op, gtos)
rs = np.linspace(0, 40, 300)
ys = gtos.at_r(cs, rs)


fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.set_title("(a) real part", loc="left")
ax1.set_xlim(0, 40)
ax1.set_xlabel(r"$r$/a.u.", fontsize=20)
ax1.plot(grid_csv["r"], grid_csv["re_y"], "k--", label="Grid")
ax1.plot(rs, ys.real, "r", label="GTO")
ax1.legend(bbox_to_anchor=(1.25,1),
           loc="upper right",
           borderaxespad=0)

ax2 = fig.add_subplot(212)
ax2.set_title("(b) imaginary part", loc="left")
ax2.set_xlim(0, 40)
ax2.set_xlabel(r"$r$/a.u.", fontsize=20)
ax2.plot(grid_csv["r"], grid_csv["im_y"], "k--", label="Grid")
ax2.plot(rs, ys.imag, "b", label="GTO")
ax2.legend(bbox_to_anchor=(1.25,1),
           loc="upper right",
           borderaxespad=0)

fig.tight_layout()
fig.subplots_adjust(right=0.78)
fig.savefig("fowf.png", dpi=50)
fig.savefig("fowf.eps")

