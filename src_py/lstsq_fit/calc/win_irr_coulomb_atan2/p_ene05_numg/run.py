import sys
sys.path.append("..")
from common import *
import json

dir_name = "out"
if not os.path.exists(dir_name):
    os.mkdir(dir_name)

os.chdir(dir_name)

xmin_list = [5, 7, 9]
r0_win_list = [30]
num_gto_list = [2, 3, 4, 5]

L = 1
pn = L+1
ene = 0.5
k = np.sqrt(ene*2.0)
z = 1.0
eta = -z/k
z_arc = 1.0

## number of fitting grid 
num = 100

## json
out_dict = {
    "L": L,
    "E": ene,
    "k": k,
    "Z": z,
    "eta": eta,
    "z_arc": z_arc,
}
json.dump(out_dict, open("common.setting.json", "w"), indent=2)

## prepation of calculation
non_win = lambda x: coulomb.outgoing(eta, k*x, L, 0)

## start calculation
for xmin in xmin_list:
    for r0_win in r0_win_list:
        for num_gto in num_gto_list:
            # target
            winfunc = lambda x: np.arctan(-z_arc*(x-r0_win))/np.pi+0.5
            target = lambda x:winfunc(x)*non_win(x)

            w = xmin * xmin
            res = lstsq_from_region(num_gto, num, w, pn, target)
            base = "{0}_{1}_{2}".format(xmin, r0_win, num_gto)
            if(res==None):
                print base, "x"
            if(res):
                print base, "ok"
                (cs, zs, other) = res
                out_dict["xmin_grid"] =  xmin
                out_dict["r0_win"] = r0_win
                
                
                filename = base + ".setting.json"
                f = open(filename, "w")
                json.dump(out_dict, f, indent=2)
                f.close()
            
                filename = base + ".res.csv"
                f = open(filename, "w")
                f.write("re_c,im_c,pn,re_z,im_z\n")
                for (c, z) in zip(cs, zs):
                    f.write("{0},{1},{2},{3},{4}\n".format(c.real, c.imag, pn, -z.real, -z.imag))
                f.close()
                
