import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

calc_csv = "wf.csv"
grid_csv = "../../grid.csv"

df_calc = pd.read_csv(calc_csv)
plt.plot(df_calc["r"], df_calc["re_y"], "r", label="Re[CBF]")
plt.plot(df_calc["r"], df_calc["im_y"], "b", label="Im[CBF]")

df_grid = pd.read_csv(grid_csv)
plt.plot(df_grid["r"], df_grid["re_y"], "k--", label="Re[Grid]")
plt.plot(df_grid["r"], df_grid["im_y"], "k-.", label="Im[Grid]")

plt.xlim(0, 40)
plt.legend(loc = "upper right")
plt.savefig("psi.png", dpi=50)

