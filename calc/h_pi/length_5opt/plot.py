import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

psi_csv = "psi.csv"
grid_csv = "grid.csv"

df_calc = pd.read_csv(psi_csv)
plt.plot(df_calc["r"], df_calc["re_y"], label="Re[CBF]")
plt.plot(df_calc["r"], df_calc["im_y"], label="Im[CBF]")

df_grid = pd.read_csv(grid_csv)
plt.plot(df_grid["r"], df_grid["re_y"], label="Re[Grid]")
plt.plot(df_grid["r"], df_grid["im_y"], label="Im[Grid]")

plt.xlim(0, 40)
plt.legend(loc = "upper right")
plt.title("psis")
plt.savefig("psi.png", dpi=50)



