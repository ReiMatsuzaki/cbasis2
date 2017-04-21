from scipy.linalg import eigh
import sys
sys.path.append("../")
from common import *


g_i = get_gi_RM()
mat  = calc_mat_complex(g_i, False)
s = np.array(mat.get_matrix("s", 0, 0)).real
(val, vec) = eigh(s)

out_dir = "out/"
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

f = open("out/eigval.csv", "w")
f.write("eigval\n")
for v in val:
    f.write("{0}\n".format(v))
f.close()

## see out/eigval.csv file
min_index_list = [0, 1, 2]
for min_index in min_index_list:
    min_vec = vec[:,min_index]
    f = open("out/eigvec{0}.csv".format(min_index), "w")
    f.write("re,im,abs\n")
    for (idx,v) in zip(range(len(min_vec)), min_vec):
        f.write("{0},{1}\n".format(idx,v))
    f.close()

## seeing out/min_eigvec.csv, most big component of basis is 
## eigvec0 : 2, 23  (10+6+4+ 3)
## eigvec1 : 0, 25  (10+6+4+ 5)
## eigvec2 : 1, 24  (10+6+4+ 4)


