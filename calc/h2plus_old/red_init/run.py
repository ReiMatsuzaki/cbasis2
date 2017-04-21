import sys
from itertools import product
sys.path.append("../")
from common import *

dir_name = "out"
if not os.path.exists(dir_name):
    os.mkdir(dir_name)
os.chdir(dir_name)

## ==== Final state ====
def get_sub_s(zs):
    sub_s = (SubSymGTOs()
	     .xyz((0,0,+R0/2.0))
	     .xyz((0,0,-R0/2.0))
	     .ns( (0,0,0))
	     .rds(Reduction(sym.irrep_z, [[1], [-1]]))
	     .zeta(zs))
    return sub_s

## ==== setting basis set ====
ls0 = [1,3]
ls1 = [1,3]
enes = [0.001] + list(np.linspace(0.1,3.0, 30))
rzeta = et_zeta(15, 0.01, 30.0)
g_psi0_L = get_gpsi0_L(ls0, rzeta)
g_chi_L = get_gchi_L(ls0, 1.0)

for red_num in [0,1,2,3]:
    name = str(red_num)
    
    ## -- Initial state --
    g_i = get_gi_RM(red_num)
    (ei, ci) = solve_init(g_i)
    
    ## -- setting basis set --
    g_psi1 = get_gpsi1(   ls1, rzeta, [get_sub_s([3,9])])    

    solver = PhotoIonizationCoulombPlus(2.0, 1, ls0)
    solver.precalc_1ele(g_i, g_chi_L, g_psi1, g_psi0_L)
    
    ## -- calculate cross sections and beta --
    cs_sigu = [solver.calc_one(ene-ei, True, "sigu")[0] for ene in enes]
    beta=     [solver.calc_one(ene-ei, True, "total")[1] for ene in enes]
    df = pd.DataFrame([enes, cs_sigu, beta]).T
    df.columns = ["energy", "cs_sigu", "beta"]
    df.to_csv(name + ".res.csv", index=False)


