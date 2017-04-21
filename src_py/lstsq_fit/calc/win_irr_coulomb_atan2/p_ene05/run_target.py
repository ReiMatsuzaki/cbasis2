import sys
sys.path.append("..")
import json
from common import *

## ==== Common ====
out_dir = "out/"
xs = np.linspace(0.01, 60.0, 400)
base = "5_30"

## ==== setting ====
setting_json = json.load(open(out_dir + base + ".setting.json", "rb"))
k = setting_json["k"]
L = setting_json["L"]
eta = setting_json["eta"]
z_arc = setting_json["z_arc"]
r0_win = setting_json["r0_win"]

## ==== Fitting Target ====
h_func = lambda x: coulomb.outgoing(eta, k*x, L, 0)
win_func = lambda x: np.arctan(-z_arc * (x-r0_win))/np.pi + 0.5
wh_func = lambda x: h_func(x) * win_func(x)

hs = np.array(map(h_func, xs))
whs = np.array(map(wh_func, xs))

with open(out_dir + "target.csv", 'w') as f:
    f.write("x,g,f,wg,wf\n")
    for (x, h, wh) in zip(xs, hs, whs):
        f.write("{0}, {1}, {2}, {3}, {4}\n".format(x, h.real, h.imag, wh.real, wh.imag))

