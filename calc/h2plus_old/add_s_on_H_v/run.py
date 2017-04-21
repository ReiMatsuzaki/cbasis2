import sys
from itertools import product
sys.path.append("../")
from common import *

dir_name = "out"
if not os.path.exists(dir_name):
    os.mkdir(dir_name)


os.chdir(dir_name)

## ==== Setting ====
rs = np.linspace(0.0, 100.0, 200)

## ==== Initial state ====
g_i = get_gi_RM()
(ei, ci) = solve_init(g_i)

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


for (name, zs) in [("3_9", [3,9]),
                   ("3", [3]),
                   ("null", None)]:
    print name
    
    ## setting basis set
    print "setting"
    if zs == None:
        g_psi1 = get_gpsi1(   ls1, rzeta)
    else:
        g_psi1 = get_gpsi1(   ls1, rzeta, [get_sub_s(zs)])    

    solver = PhotoIonizationCoulombPlus(2.0, 1, ls0, "velocity")
    solver.precalc_1ele(g_i, g_chi_L, g_psi1, g_psi0_L)
    
    ## calculate cross sections and beta
    print "calc"
    ws = [ene-ei for ene in enes]
    cs_sigu = [solver.calc_one(w, True, "sigu")[0]   for w in ws]
    cs_sigu_ima = [solver.calc_cs_imalpha(w, "sigu") for w in ws]
    beta=     [solver.calc_one(w, True, "total")[1]  for w in ws]
    df = pd.DataFrame([enes, ws, cs_sigu, beta, cs_sigu_ima]).T
    df.columns = ["energy", "w", "cs_sigu", "beta", "cs_sigu_ima"]
    df.to_csv(name + ".res.csv", index=False)

    ## calculate psi at E=1.0
    print "psi"
    ene = 1.0
    w = ene - ei    
    c1 = solver.calc_psi1(w, 0)
    
    (ys, dy1s) = g_psi1.at_r_ylm_cpp(1, 0, sym.irrep_z, c1, rs)
    y1s = np.array(ys)
    (ys, dy1s) = g_psi1.at_r_ylm_cpp(3, 0, sym.irrep_z, c1, rs)
    y3s = np.array(ys)
    (ys, dy1s) = g_psi1.at_r_ylm_cpp(5, 0, sym.irrep_z, c1, rs)
    y5s = np.array(ys)
    
    df_psi1 = pd.DataFrame([rs,
                            y1s.real, y1s.imag,
                            y3s.real, y3s.imag,
                            y5s.real, y5s.imag]).T
    
    df_psi1.columns = ["r", "re_1", "im_1", "re_3", "im_3", "re_5", "im_5"]
    df_psi1.to_csv(name + ".psi.csv", index=False)
    
