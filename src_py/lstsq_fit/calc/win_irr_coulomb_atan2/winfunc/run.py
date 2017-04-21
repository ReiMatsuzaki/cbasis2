import sys
sys.path.append("..")
from common import *
import json

dir_name = "out"
if not os.path.exists(dir_name):
    os.mkdir(dir_name)

os.chdir(dir_name)

xmin_list = [3, 4, 5, 6, 7, 8, 9, 10]
r0_win_list = [20, 30, 40]

L = 1
pn = L+1
ene = 0.5
k = np.sqrt(ene*2.0)
z = 1.0
eta = -z/k
num_gto = 5
xmin = 5.0

#r0_win0 = 1.0
a_win0  = 1.0
r0_win1 = 30.0
a_win1  = 1.0

## number of fitting grid 
num = 100

## json
out_dict = {
    "L": L,
    "E": ene,
    "k": k,
    "Z": z,
    "eta": eta,
    "num_grid": num_gto,
    "a_win0": a_win0,
    "r0_win1": r0_win1,
    "a_win1": a_win1
}
json.dump(out_dict, open("common.setting.json", "w"), indent=2)

## ==== prepation of calculation ====
non_win = lambda x: coulomb.outgoing(eta, k*x, L, 0)
#s = lambda x: np.exp( - a*(x-x0))
#win0 = lambda x: (x**(2*L+1) + s(x)) / (1.0 + s(x))
win0 = lambda x: 1.0 - np.exp(-a_win0 * x)
win1= lambda x: np.arctan(-a_win1*(x-r0_win1))/np.pi+0.5

b2c = lambda b: "y" if b else "n"
w = xmin * xmin

def fit_write(name, target):
    ## fitting
    res = lstsq_from_region(num_gto, num, w, pn, target)

    ## print and write results
    if res:
        print name, ": ok"
        ## write to file
        (cs, zs, other) = res
        f = open(name + ".res.csv", "w")
        f.write("re_c,im_c,pn,re_z,im_z\n")
        for (c, z) in zip(cs, zs):
            f.write("{0},{1},{2},{3},{4}\n".format(c.real,c.imag,pn,-z.real,-z.imag))
        f.close()
    else:
        print name, ": x"    

fit_write("y_y", lambda x: non_win(x) * win0(x) * win1(x))
fit_write("n_y", lambda x: non_win(x) * win1(x))
fit_write("y_n", lambda x: non_win(x) * win0(x))
fit_write("n_n", non_win)

