import matplotlib.pyplot as plt
import numpy as np
import os
#import csv
import pandas as pd
import json
import sys
sys.path.append("..")
from common import *

## ==== Common ====
out_dir = "out/"
fig_dir = "fig/"
xs = np.linspace(0.01, 60.0, 400)
setting_json = json.load(open(out_dir +  "common.setting.json", "rb"))

k = setting_json["k"]
L = setting_json["L"]
eta = setting_json["eta"]
target = lambda x: coulomb.outgoing(eta, k*x, L, 0)
ts = np.array(map(target, xs))

## ==== plot only real part ====
def plot_re_psi(base_list, filename):
    ## ==== Fitting target ====
    plt.plot(xs, ts.real, "k", label="Re[Target]")

    ## ==== Fitting results =====
    for base in base_list:
        res_csv = pd.read_csv(out_dir + base + ".res.csv")
        cs = res_csv["re_c"].as_matrix() + 1.0j * res_csv["im_c"].as_matrix()
        pns= res_csv["pn"].as_matrix()
        zs = res_csv["re_z"].as_matrix() + 1.0j * res_csv["im_z"].as_matrix()
        def func(r):
            return sum([c * r**pn * np.exp(-z*r*r)
                        for (c, pn, z)
                        in zip(cs, pns, zs)])
        fs = np.array(map(func, xs))
        plt.plot(xs, fs.real, label=base)
        
    ## ==== Plot ====
    plt.ylim(-1.3, 1.3)
    plt.legend()
    
    plt.savefig(fig_dir + filename + "_psi.png", dpi=50)
    plt.clf()

    ## ==== plot orbital exponents ====
    for base in base_list:
        res_csv = pd.read_csv(out_dir + base + ".res.csv")
        re_zs = res_csv["re_z"].as_matrix()
        im_zs = res_csv["im_z"].as_matrix()
        plt.plot(re_zs, im_zs, "o", label=base)

    plt.legend()
    plt.savefig(fig_dir + filename + "_oe.png", dpi=50)
    plt.clf()


if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)
plot_re_psi(["3_30", "5_30", "7_30", "9_30"], "r0_30")
plot_re_psi(["5_20", "5_30", "5_40"], "x0_5")
plot_re_psi(["7_20", "7_30", "7_40"], "x0_7")


