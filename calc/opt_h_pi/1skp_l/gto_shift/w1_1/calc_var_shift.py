import sys
import numpy as np
import pandas as pd

sys.path.append("../../../../../r1basis")
from r1basis import *
from opt_green import *

## see old calculation in ~/src/git/opt_cbf/py_bind/calc/opt_gto/ret_shift/
# ==== ET-basis ====
def et_ox_from_max(r0, r1, n, num):
    """
    u(r) = r**n * exp(-zr**2)
    u'(r)= [nr**(n-1) -2zr**(n+1)] * exp(-zr**2)
    <=> n = 2zrr
    <=> z = n/(2* r_max**2)
    """
    z0 = n / (2.0 * r1*r1)
    z1 = n / (2.0 * r0*r0)
    ratio = (z1/z0)**(1.0/(num-1))
    return [z0*ratio**n for n in range(num)]


z0s = et_ox_from_max(1.0, 10.0, 2, 10)
print z0s[0], z0s[-1], z0s[1]/z0s[0], len(z0s)
yopt = 0.0037-0.039j

h_pi = H_Photoionization("1s->kp", "length")
basis_type = "GTO"
basis_info = [('shift', True, 2, z0s, yopt)]
(base_us, opt_list, var, y0s) = h_pi_read_info(basis_type, basis_info)
w = 1.0

v_w_ys = v_green_h_pi(h_pi, base_us, opt_list, var)
v_ys = v_w_ys(w)

re_data = []
im_data = []
#re_y0s = np.linspace(-0.001, 0.005, 60)
#im_y0s = np.linspace(-0.005, -0.001, 100)
"""
re_y0s = np.linspace(0.0037, 0.0038, 100)
im_y0s = np.linspace(-0.0395, -0.0385, 100)
"""
re_y0s = np.linspace(0.003, 0.004, 100)
im_y0s = np.linspace(-0.04, -0.03, 100)
for im_y0 in im_y0s:
    re_tmp = []
    im_tmp = []
    for re_y0 in re_y0s:
        y0 = re_y0 + 1.0j * im_y0
        v = v_ys([y0])
        re_tmp.append(v.real)
        im_tmp.append(v.imag)
    re_data.append(re_tmp)
    im_data.append(im_tmp)
    
df = pd.DataFrame(re_data)
df.to_csv('re_alpha.csv', index=False)
df = pd.DataFrame(im_data)
df.to_csv('im_alpha.csv', index=False)

re_z_f = open('re_z.csv', 'w')
for re_y0 in re_y0s:
    print >>re_z_f, re_y0

im_z_f = open('im_z.csv', 'w')
for im_y0 in im_y0s:
    print >>im_z_f, im_y0
