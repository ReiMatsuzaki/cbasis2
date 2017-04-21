import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
sys.path.append("../")
from ref import *

## ===== directory ====
out_dir = "out/"
fig_dir = "fig/"
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)

## ===== figure ====
fig_cs_sigu, ax_cs_sigu = plt.subplots()
fig_cs_piu , ax_cs_piu  = plt.subplots()
fig_cs_tot , ax_cs_tot  = plt.subplots()
fig_beta, ax_beta = plt.subplots()
    
name_list = ["p_p", "p_pf", "pf_p", "pf_pf"]
c_list    = ["r", "rx",    "b",  "bx"]

for (c, name) in zip(c_list, name_list):
    df = pd.read_csv(out_dir + name + ".res.csv")
    ene = df["energy"].as_matrix()
    cs_sigu = df["cs_sigu"]; ax_cs_sigu.plot(ene, cs_sigu, c, label=name)
    cs_piu = df["cs_piu"];   ax_cs_piu.plot( ene, cs_piu,  c, label=name)
    cs_tot = df["cs_total"]; ax_cs_tot.plot( ene, cs_tot,  c, label=name)
    beta = df["beta"];       ax_beta.plot(   ene, beta,    c, label=name)

ax_cs_sigu.legend()
ax_cs_sigu.plot(ref_ene, ref_cs_sigu, "ko", label="ref")
fig_cs_sigu.savefig(fig_dir + "cs_sigu.png", dpi=50)

ax_cs_piu.legend()
ax_cs_piu.plot(ref_ene, ref_cs_piu, "ko", label="ref")
fig_cs_piu.savefig( fig_dir + "cs_piu.png", dpi=50)

ax_cs_tot.legend()
ax_cs_tot.plot(ref_ene, ref_cs_tot, "ko", label="ref")
fig_cs_tot.savefig( fig_dir + "cs_tot.png", dpi=50)

ax_beta.legend()
ax_beta.plot(ref_ene, ref_beta, "ko", label="ref")
fig_beta.savefig(   fig_dir + "beta.png", dpi=50)


