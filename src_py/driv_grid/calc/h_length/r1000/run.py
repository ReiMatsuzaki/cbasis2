import sys
import numpy as np
sys.path.append("../../../")
from driv_grid import solve_driv

## system
L = 1.0
ene = 0.5

## driven term
s = lambda r: 2.0 * r * r * np.exp(-r)

## potential 
v = lambda r: -1.0/r + L*(L+1)/(2.0*r*r)  

## grid points
n = 10000
h = 0.1

## calculate
(rs, ys) = solve_driv(v, ene, s, n, h)

## write
with open("h_length.csv", "w") as f:
    for (r, y) in zip(rs, ys):
        f.write("{0},{1},{2}\n".format(r, y.real, y.imag))


