import sys
import numpy as np
sys.path.append("../../")
from driv_grid import solve_driv

"""
solve driven equation appeared in 
.      F.Mrugala, J.Comput.Phys. 133*, 113 (1985)
"""

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

## compute dipole, written as M in the paper         
k = np.sqrt(2.0 * ene)
y = ys[-2]
r = rs[-2]
print "dip = ", abs(y*k / (4*np.exp(complex(0, k*r - np.pi*L/2.0 + 1.0/k * np.log(2*k*r) ))))

## compute dipole from matrix element of Green's operator
xs = np.array([(i+1)*h for i in range(n)])
ss = np.array([s(x) for x in xs])
alpha = np.dot(ys, ss) * h
k = np.sqrt(2.0 * ene)
print "alpha = ", alpha
print "dip = ", np.sqrt(alpha.imag * k / 2.0) / 2.0 


