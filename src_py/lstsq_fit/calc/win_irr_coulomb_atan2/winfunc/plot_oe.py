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

def plot_name(name, c):
    df = pd.read_csv(out_dir + name + ".res.csv")
    re_zs = df["re_z"].as_matrix()
    im_zs = df["im_z"].as_matrix()
    plt.plot(re_zs, im_zs, c, label=name)

plot_name("n_n", "ro")
plot_name("n_y", "bo")
plot_name("y_n", "rx")
plot_name("y_y", "bx")

plt.legend()
plt.savefig(fig_dir + "oe.png", dpi=50)
