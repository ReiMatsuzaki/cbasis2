import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import sys
sys.path.append("..")
from common import *

## ==== common ====
out_dir = "out/"
fig_dir = "fig/"
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)
xs = np.linspace(0.01, 40.0, 200)

## ==== setting ====
setting_json = json.load(open(out_dir + "common.setting.json", "rb"))

keys = ["k", "L", "eta", "a_win0", "r0_win1", "a_win1"]
tmp = [setting_json[key] for key in keys]
(k, L, eta, a_win0, r0_win1, a_win1) = tmp

def plot_name(name, c):
    df = pd.read_csv(out_dir + name + ".res.csv")
    cs = df["re_c"].as_matrix() + 1.0j * df["im_c"].as_matrix()
    pns= df["pn"].as_matrix()
    zs = df["re_z"].as_matrix() + 1.0j * df["im_z"].as_matrix()
    def func(r):
        return sum([c * r**pn * np.exp(-z*r*r)
                    for (c, pn, z)
                    in zip(cs, pns, zs)])
    fs = np.array(map(func, xs))
    plt.plot(xs, fs.real, c, label=name)

plot_name("n_n", "r")
plot_name("n_y", "b")
plot_name("y_n", "rx")
plot_name("y_y", "bx")

non_win = lambda x: coulomb.outgoing(eta, k*x, L, 0)
ts = np.array(map(non_win, xs))
plt.plot(xs, ts.real, "k--", label="Coulomb")

plt.legend()
plt.xlim(0, 50.0)
plt.ylim(-1.3, 1.3)
plt.savefig(fig_dir + "fit.png", dpi=50)
