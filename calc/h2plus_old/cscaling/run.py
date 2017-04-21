## from 2016_8_31.ipynb
import sys
from itertools import product
sys.path.append("../")
from common import *

enes = [0.001] + list(np.linspace(0.1,3.0, 30))

dir_name = "out"
if not os.path.exists(dir_name):
    os.mkdir(dir_name)
os.chdir(dir_name)
              
## ==== Initial state ====
g_i = get_gi_RM()
(ei, ci) = solve_init(g_i)

## ==== Final state ====
ls0 = [1,3]
ls1 = [1,3]
l2c = {1:"p", 3:"f", 5:"h"}

sub_s = (SubSymGTOs()
	 .xyz((0,0,+R0/2.0))
	 .xyz((0,0,-R0/2.0))
	 .ns( (0,0,0))
	 .rds(Reduction(sym.irrep_z, [[1], [-1]]))
	 .zeta([3,9]))

num_et = 15
theta_list = [20, 30, 40, 50]
for theta in theta_list:
    name = str(theta)
    
    ## -- basis --
    zs_real = et_zeta(num_et, 0.002, 30.0)
    zs = [np.exp(-1.0j * theta * np.pi / 180.0) * z for z in zs_real]
    g_psi1 = get_gpsi1_simple(   ls1,  zs, [sub_s])
    g_psi0_L = get_gpsi0_L_simple(ls0, zs)
    g_chi_L = get_gchi_L(ls0, 1.0)
    solver = PhotoIonizationCoulombPlus(2.0, 1, ls0)
    solver.precalc_1ele(g_i, g_chi_L, g_psi1, g_psi0_L)

    ## calculate cross sections and beta
    cs_sigu = [solver.calc_one(ene-ei, True, "sigu")[0] for ene in enes]
    cs_beta=  np.array([solver.calc_one(ene-ei, True, "total") for ene in enes])
    cs_total = list(cs_beta.T[0])
    beta = list(cs_beta.T[1])
    cs_piu = [a-b for (a,b) in zip(cs_total, cs_sigu)]
    
    df = pd.DataFrame([enes, cs_sigu, cs_piu, cs_total, beta]).T
    df.columns = ["energy", "cs_sigu", "cs_piu", "cs_total", "beta"]
    df.to_csv(name + ".res.csv", index=False)    
