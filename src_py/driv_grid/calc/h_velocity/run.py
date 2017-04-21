import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
sys.path.append("../../")
from driv_grid import solve_driv

## system
L = 1.0
ene = 0.5

## driven term
s = lambda r: 2.0 * r * np.exp(-r)

## potential 
v = lambda r: -1.0/r + L*(L+1)/(2.0*r*r)  

## grid
n = 1000
h = 0.1

## calculate and write
(rs, ys) = solve_driv(v, ene, s, n, h)
df = pd.DataFrame([rs.real, ys.real, ys.imag]).T
df.columns = ["r", "re_y", "im_y"]
df.to_csv("h_velo.csv", index=False)


df2 = pd.read_csv("h_velo.csv")
rs = df2["r"]
re_ys = df2["re_y"]
im_ys = df2["im_y"]
plt.plot(rs, re_ys, label="Re")
plt.plot(rs, im_ys, label="Im")
plt.xlim(0, 40)
plt.savefig("h_velo.png", dpi=50)
plt.clf()


