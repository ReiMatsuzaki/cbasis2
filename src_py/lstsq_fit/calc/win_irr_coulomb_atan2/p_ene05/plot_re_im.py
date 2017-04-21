import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import json
import sys
#sys.path.append("..")

## ==== Common ====
out_dir = "out/"
fig_dir = "fig/"
base = "5_30"

## ==== setting ====
#setting_json = json.load(open(out_dir + base + ".setting.json", "rb"))
#k = setting_json["k"]
#L = setting_json["L"]
#eta = setting_json["eta"]
#z_arc = setting_json["z_arc"]
#r0_win = setting_json["r0_win"]

## ==== target ====
t_csv = pd.read_csv(out_dir + "target.csv")
xs = t_csv["x"].as_matrix()
hs = t_csv["g"].as_matrix() + 1.0j * t_csv["f"].as_matrix()
whs = t_csv["wg"].as_matrix() + 1.0j * t_csv["wf"].as_matrix()

## ==== basis function expansion ====
res_csv = pd.read_csv(out_dir + base + ".res.csv")
cs = res_csv["re_c"].as_matrix() + 1.0j * res_csv["im_c"].as_matrix()
pns= res_csv["pn"].as_matrix()
zs = res_csv["re_z"].as_matrix() + 1.0j * res_csv["im_z"].as_matrix()
def func(r):
    return sum([c * r**pn * np.exp(-z*r*r)
                for (c, pn, z)
                in zip(cs, pns, zs)])
fs = np.array(map(func, xs))

## ==== Plot ====
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(xs,  fs.real, "r", label="GTO")
ax1.plot(xs_, hs.real, "k--", label=r"$g_{kl}(r)$")
ax1.plot(xs_, whs.real, "k-.", label=r"$w(r)g_{kl}(r)$")
#ax1.plot(xs[0:-1:5], ts_with_win.real[0:-1:5], "k-.", 
#         label=r"$w(r)g_{kl}(r)$")

ax1.set_title("(a)real part", loc="left")
ax1.set_xlim(0, +40.0)
ax1.set_ylim(-1.3, +1.6)
ax1.legend(bbox_to_anchor=(1.37,1),
           loc="upper right",
           borderaxespad=0)
plt.subplots_adjust(right=0.75)

ax2 = fig.add_subplot(212)
ax2.plot(xs, fs.imag, "b", label="GTO")
ax2.plot(xs, hs.imag, "k--", label=r"$f_{kl}(r)$")
ax2.plot(xs, whs.imag,"k-.", label=r"$w(r)f_{kl}(r)$")
ax2.set_title("(b)imaginary part", loc="left")
ax2.set_xlim(0, +40.0)
ax2.set_ylim(-1.3, +1.6)
ax2.legend(bbox_to_anchor=(1.37,1),
           loc="upper right",
           borderaxespad=0)

fig.savefig(fig_dir + "x0_5_r0_30_re_im.png", dpi=50)
fig.savefig(fig_dir + "x0_5_r0_30_re_im.eps")
