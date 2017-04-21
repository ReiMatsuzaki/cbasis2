"""
compute initial state
"""
import sys
sys.path.append("../")
from common import *

g_i = get_gi_RM()
mat = calc_mat_complex(g_i, True)
t = mat.get_matrix("t", 0, 0)
v = mat.get_matrix("v", 0, 0)
h = t+v
s = mat.get_matrix("s", 0, 0)
tmp = ceig(h,s)
ei = tmp[0][0].real
ci = tmp[1].col(0)

f = open("res.out", "w")

print >>f, "sym:"
print >>f, sym

print >>f, "mole:"
print >>f, mole

print >>f, "gtos:"
print >>f, g_i.__repr__()

print >>f, "s(1,2) = ", s[1, 2]
print >>f, "s(29,29) = ", s[29, 29]
print >>f, "v(0,0) = ",  v[0, 0]
print >>f, "s(0,40) = ", s[0, 40]
print >>f, "t(0,40) = ", t[0, 40]
print >>f, "v(0,40) = ", v[0, 40]

print >>f, "s(40,40) = ", s[40, 40]
print >>f, "t(40,40) = ", t[40, 40]
print >>f, "v(40,40) = ", v[40, 40]

print >>f, "ei:"
print >>f, ei

print >>f, "ci:"
for c in ci:
    print >>f, c

f.close()
