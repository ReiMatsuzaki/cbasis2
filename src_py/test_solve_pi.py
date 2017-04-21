import sys
sys.path.append("coulomb/")

import unittest
import random
from part_wave import *
from solve_pi import *
from scipy import integrate
from sympy.mpmath import coulombf
import time

from coulomb import coulomb_phase

def calc_ip_numint(f, g, rs, method=0):
    """
    f : continuum wave function
    g : L2 function which is not normalized
    """
    fgs = np.array([f(r) * g(r) for r in rs])
    ggs = np.array([g(r) * g(r) for r in rs])
    if method == 0:
	int_fgs = integrate.simps(fgs, rs)
	int_ggs = integrate.simps(ggs, rs)
    else:
	int_fgs = integrate.romb(fgs, rs[1]-rs[0])
	int_ggs = integrate.romb(ggs, rs[1]-rs[0])
    return int_fgs / np.sqrt(int_ggs) 

def calc_ip_numint_2(f, g1, g2, c1, c2, rs):
    """
    f : continuum wave function
    g1, g2: L2 function not normalized
    c1, c2: coefficient of g1 and g2
    """
    
    g1g1 = np.array([ g1(r) * g1(r) for r in rs])
    n1 = 1.0/np.sqrt( integrate.simps(g1g1, rs))
    g2g2 = np.array([ g2(r) * g2(r) for r in rs])
    n2 = 1.0/np.sqrt( integrate.simps(g2g2, rs))
    
    fg1g2 = np.array([f(r) * (c1 * n1 * g1(r) + c2 * n2 * g2(r)) for r in rs])
    return integrate.simps(fg1g2, rs)


def get_gi_RM(sym, R0):
    """ see T.Rescigno and C.McCurdy(1985) PRA 31, 624
    """
    mult = lambda x: 3.0*x
    zeta_s_h = map(mult,
		   [0.0285649, 0.0812406, 0.190537, 0.463925, 1.202518,
		    3.379649, 10.60720, 38.65163, 173.5822, 1170.498])
    zeta_p_h = map(mult,
		   [0.015442, 0.035652, 0.085676,
		    0.227763, 0.710128, 3.009711])
    zeta_d_h = [15.0, 5.0, 5.0/3.0, 5.0/9.0]
    zeta_s_cen = [14.4/(3.0**n) for n in range(12)]
    zeta_d_cen = [10.0/(3.0**n) for n in range(9)]
    g_i = (SymGTOs.create().sym(sym)
	   .sub(SubSymGTOs()
		.xyz((0,0,+R0/2.0)).xyz((0,0,-R0/2.0))
		.ns( (0,0,0))
		.rds(Reduction(sym.irrep_s, [[1], [1]]))
		.zeta(zeta_s_h))
	   .sub(SubSymGTOs()
		.xyz((0,0,+R0/2.0)).xyz((0,0,-R0/2.0))
		.ns( (0,0,1))
		.rds((Reduction(sym.irrep_s, [[1], [-1]])))
		.zeta(zeta_p_h))
	   .sub(SubSymGTOs()
		.xyz((0,0,+R0/2.0)).xyz((0,0,-R0/2.0))
		.ns( (0,0,2))
		.rds(Reduction(sym.irrep_s, [[1], [1]]))
		.zeta(zeta_d_h))
	   .sub(SubSymGTOs()
		.xyz((0,0,0))
		.ns( (0,0,0))
		.rds(Reduction(sym.irrep_s, [[1]]))
		.zeta(zeta_s_cen))
           .sub(SubSymGTOs()
		.xyz((0,0,0))
		.ns( (0,0,2))
		.rds(Reduction(sym.irrep_s, [[1]]))
		.zeta(zeta_d_cen))
	   .atom((0,0,0+R0/2.0), 1)
	   .atom((0,0,0-R0/2.0), 1)
	   .setup())
    return g_i

def get_gi(sym, R0):
    zeta_s_h = [0.0285649, 0.0812406, 0.190537, 0.463925, 1.202518, 3.379649, 10.60720, 38.65163, 173.5822, 1170.498]
    zeta_p_h = [0.015442, 0.035652, 0.085676, 0.227763, 0.710128, 3.009711]


    sub_s_H = (SubSymGTOs()
	       .xyz((0, 0, +R0/2.0)).xyz((0, 0, -R0/2.0))
	       .ns( (0, 0, 0))
	       .rds(Reduction(sym.irrep_s, [[1], [1]]))
	       .zeta(zeta_s_h))
    sub_p_H = (SubSymGTOs()
	       .xyz((0, 0, +R0/2.0)).xyz((0, 0, -R0/2.0))
	       .ns( (0, 0, 1))
	       .rds(Reduction(sym.irrep_s, [[1],[-1]]))
	       .zeta(zeta_p_h))
    
    g_i = (SymGTOs.create().sym(sym)
           .sub(sub_s_H).sub(sub_p_H)
           .atom((0,0,+R0/2.0),1)
           .atom((0,0,-R0/2.0),1)
           .setup())
    return g_i

def get_gi_small(sym, R0):
    zeta_s_h = [0.0285649, 0.0812406]
    zeta_p_h = [0.015442]

    sub_s_H = (SubSymGTOs()
	       .xyz((0, 0, +R0/2.0)).xyz((0, 0, -R0/2.0))
	       .ns( (0, 0, 0))
	       .rds(Reduction(sym.irrep_s, [[1], [1]]))
	       .zeta(zeta_s_h))
    sub_p_H = (SubSymGTOs()
	       .xyz((0, 0, +R0/2.0)).xyz((0, 0, -R0/2.0))
	       .ns( (0, 0, 1))
	       .rds(Reduction(sym.irrep_s, [[1],[-1]]))
	       .zeta(zeta_p_h))
    
    g_i = (SymGTOs.create().sym(sym)
           .sub(sub_s_H).sub(sub_p_H)
           .atom((0,0,+R0/2.0),1)
           .atom((0,0,-R0/2.0),1)
           .setup())
    return g_i

def get_basic_basis(sym, R0, ls):
    zeta_p_et = [0.01,0.0177160054,0.0313856847, 0.0556028960, 0.0985061205, 0.174513496,
	       0.309168204, 0.547722558, 0.970345578, 1.71906475, 3.04549604, 5.39540243,
	       9.55849785, 16.9338400, 30.0000000]

    ## below basis is from /Users/rei/src/git/opt_cbf/lstsq_fit/calc/win_irr_coulomb_atan/5basis
    zeta_list = {}
    zeta_p = [0.0812755955262-0.0613237786222j,
              0.00497147387796-0.0113737972763j,
              0.0323712622673-0.0451300037076j,
              0.00317417887792-0.022582254987j,
              0.0118391719646-0.0327847352576j] + zeta_p_et
    
    ## below basis is 4/14 fitted-GTO + ET-GTO. It is not efficient basis
    ## but is usefull for compare in this test
    zeta_p = [ 0.00256226-0.01559939j, 0.00389597-0.02240632j,
	       0.01229986-0.03080238j,0.03010506-0.04378147j,
	       0.01, 0.0177160054, 0.0313856847, 0.0556028960, 0.0985061205, 0.174513496,
	       0.309168204, 0.547722558, 0.970345578, 1.71906475, 3.04549604, 5.39540243,
	       9.55849785, 16.9338400, 30.0000000]
    zeta_f = [
        0.114846838205-0.0528053215279j,
        0.0502495991466-0.0443217040212j,
        0.00291922532665-0.0190617777251j,
        0.0223369214731-0.0358454508219j,
        0.0083110955684-0.0268941896615j] + zeta_p_et

    zeta_list[1] = zeta_p
    zeta_list[3] = zeta_f
    zeta_list[5] = zeta_p

    """
    zeta_p = [ 0.00256226-0.01559939j, 0.00389597-0.02240632j,
	       0.01229986-0.03080238j,0.03010506-0.04378147j,
	       0.01, 0.0177160054, 0.0313856847, 0.0556028960, 0.0985061205, 0.174513496,
	       0.309168204, 0.547722558, 0.970345578, 1.71906475, 3.04549604, 5.39540243,
	       9.55849785, 16.9338400, 30.0000000]
    """
    
    g_psi1 = (SymGTOs.create().sym(sym)
              .atom((0,0,+R0/2.0), 1).atom((0,0,-R0/2.0), 1))
    
    for L in ls:
        for (Mc, irrep) in zip([-1,0,1],
                               [sym.irrep_y, sym.irrep_z, sym.irrep_x]):
            sub_p_Cen_xyz = sub_solid_sh(L, Mc, zeta_list[L], irrep)
            g_psi1.sub(sub_p_Cen_xyz)

    g_psi1.setup()
    
    g_psi0_L = {}
    g_chi_L = {}
    for L in ls:
        g_psi0_L[L] = solid_part_waves(sym, L, [-1,0,1], zeta_list[L])
        g_psi0_L[L].atom((0,0,-R0/2.0), 1).atom((0,0,+R0/2.0), 1).setup()
        g_chi_L[L] = solid_part_waves(sym,  L, [-1,0,1], [1])
        
    return (g_psi1, g_psi0_L, g_chi_L)
    
def get_basic_basis_small(sym, R0, ls):
    ## below basis is from /Users/rei/src/git/opt_cbf/lstsq_fit/calc/win_irr_coulomb_atan/5basis
    zeta_list = {}
    zeta_p = [0.0812755955262-0.0613237786222j,
              0.00497147387796-0.0113737972763j,
              0.0323712622673-0.0451300037076j,
              0.00317417887792-0.022582254987j,
              0.0118391719646-0.0327847352576j]

    zeta_f = [
        0.114846838205-0.0528053215279j,
        0.0502495991466-0.0443217040212j,
        0.00291922532665-0.0190617777251j,
        0.0223369214731-0.0358454508219j,
        0.0083110955684-0.0268941896615j]

    zeta_list[1] = zeta_p
    zeta_list[3] = zeta_f
    zeta_list[5] = zeta_p

    g_psi1 = (SymGTOs.create().sym(sym)
              .atom((0,0,+R0/2.0), 1).atom((0,0,-R0/2.0), 1))
    
    for L in ls:
        for (Mc, irrep) in zip([-1,0,1],
                               [sym.irrep_y, sym.irrep_z, sym.irrep_x]):
            sub_p_Cen_xyz = sub_solid_sh(L, Mc, zeta_list[L], irrep)
            g_psi1.sub(sub_p_Cen_xyz)

    g_psi1.setup()
    
    g_psi0_L = {}
    g_chi_L = {}
    for L in ls:
        g_psi0_L[L] = solid_part_waves(sym, L, [-1,0,1], zeta_list[L])
        g_psi0_L[L].atom((0,0,-R0/2.0), 1).atom((0,0,+R0/2.0), 1).setup()
        g_chi_L[L] = solid_part_waves(sym,  L, [-1,0,1], [1])
        
    return (g_psi1, g_psi0_L, g_chi_L)
                          
class Test_pi_class(unittest.TestCase):
    def setUp(self):
        pass

    def _test_p(self):
        sym = D2h()
        R0 = 2.0
        g_i = get_gi(sym, R0)
        (g_psi1, g_psi0_L, g_chi_L) = get_basic_basis(sym, R0, [1])
        (ei, ci) = solve_init(g_i)

        solver = PhotoIonizationCoulombPlus(2.0, 1)
        solver.precalc_1ele(g_i, g_chi_L, g_psi1, g_psi0_L, [1])

        tmp = solver.calc_one(0.5-ei, [1], False, "total")
        (cs1, beta_1) = tmp
        tmp = solver.calc_one(0.5-ei, [1], True, "total")
        (cs2, beta_2) = tmp
        
        tmp = solver.calc_one(0.5-ei, [1], False, "sigu")
        (cs1_sig, dum) = tmp
        tmp = solver.calc_one(0.5-ei, [1], True, "sigu")
        (cs2_sig, dum) = tmp
        
        cs1_pi = cs1 - cs1_sig
        cs2_pi = cs2 - cs2_sig
        
        #
        # the reference values are taken from 2016/7/2016_7_21.ipynb
        #
        self.assertAlmostEqual(0.0419107791263, cs1)
        self.assertAlmostEqual(0.247993940969,  cs2)
        
        self.assertAlmostEqual(0.020047948729647776, cs1_pi)
        self.assertAlmostEqual(0.02186283039664566, cs1_sig)
        
        self.assertAlmostEqual(0.2479935960314486, cs2_pi)
        self.assertAlmostEqual(3.4493719728780974 * 10.0 ** (-7.0), cs2_sig)        

    def _test_pf_phase(self):
        sym = D2h()
        R0 = 2.0
        g_i = get_gi_small(sym, R0)
        (g_psi1, g_psi0_L, g_chi_L) = get_basic_basis_small(sym, R0, [1,3])
        (ei, ci) = solve_init(g_i)

        solver = PhotoIonizationCoulombPlus(2.0, 1)
        solver.precalc_1ele(g_i, g_chi_L, g_psi1, g_psi0_L, [1,3])
        
        ps_sigu = solver.phase_shift_easy(1.0-ei, 1, "sigu")
        ps_piu = solver.phase_shift_easy(1.0-ei, 1, "piu")
        
        #
        # these reference values is from M.Brosolo, P.Decleva(1992)
        #
        print "calc_pahse:", ps_sigu/np.pi, ps_piu/np.pi
        print "reference_phase:", 0.12205, -0.14923
        #self.assertAlmostEqual(0.12205*np.pi, ps_sigu)
        #self.assertAlmostEqual(-0.14923*np.pi, ps_piu)
        
        solver.precalc_phase_shift_1ele(g_chi_L, g_psi0_L, [1,3])
        t_mat = solver.t_matrix(0.4, [1,3], "sigu")
        s_mat = t_mat + np.identity(2)
        print "---s---"
        print s_mat
        print s_mat.conj().transpose() * s_mat
        print "-------"
        k_mat = t2k(t_mat)
        print t_mat
        print k_mat

    def _test_pf(self):
        sym = D2h()
        R0 = 2.0
        g_i = get_gi(sym, R0)
        (g_psi1, g_psi0_L, g_chi_L) = get_basic_basis(sym, R0, [1,3])
        (ei, ci) = solve_init(g_i)

        solver = PhotoIonizationCoulombPlus(2.0, 1)
        solver.precalc_1ele(g_i, g_chi_L, g_psi1, g_psi0_L, [1,3])

        tmp = solver.calc_one(0.5-ei, [1], False, "total")
        (cs1, beta_1) = tmp
        tmp = solver.calc_one(0.5-ei, [1], True, "total")
        (cs2, beta_2) = tmp
        
        tmp = solver.calc_one(0.5-ei, [1], False, "sigu")
        (cs1_sig, dum) = tmp
        
        tmp = solver.calc_one(0.5-ei, [1], True, "sigu")
        (cs2_sig, dum) = tmp
        
        cs1_pi = cs1 - cs1_sig
        cs2_pi = cs2 - cs2_sig
        
        #
        # the reference values are taken from 2016/7/2016_7_21.ipynb
        #
        self.assertAlmostEqual(0.0419107791263, cs1)
        self.assertAlmostEqual(0.247993940969,  cs2)
        
        self.assertAlmostEqual(0.020047948729647776, cs1_pi)
        self.assertAlmostEqual(0.02186283039664566, cs1_sig)
        
        self.assertAlmostEqual(0.2479935960314486, cs2_pi)
        self.assertAlmostEqual(3.4493719728780974 * 10.0 ** (-7.0), cs2_sig)                

    def test_p_he(self):
        
        sym = D2h()

        zeta0 = [0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694]
        g_i_he = (SymGTOs.create().sym(sym)
	          .sub(SubSymGTOs()
	               .xyz((0,0,0))
	               .ns( (0,0,0))
	               .rds(Reduction(sym.irrep_s, [[1]]))
	               .zeta(zeta0))
	          .atom((0,0,0), 2)
	          .setup())
        
        zeta1 = [0.09154356-0.24865707j]

        g_psi0_he_L = {1: solid_part_waves(sym, 1, [-1,0,1], zeta1).atom((0,0,0), 2)}
        g_chi_he_L =  {1: solid_part_waves(sym, 1, [-1,0,1], [1])}
        
        solver_he_p = PhotoIonizationCoulombPlus(1.0, 2, [1])

        eri_method = ERI_method().use_symmetry(1).coef_R_memo(1)

        solver_he_p.precalc_stex(g_i_he, g_chi_he_L,
			         g_psi0_he_L[1], g_psi0_he_L,
                                 eri_method)

        w0 = 0.91794
        (cs1, beta) = solver_he_p.calc_one(w0, False, "sigu")
        (cs2, beta) = solver_he_p.calc_one(w0, True, "sigu")
        cs_im_a = solver_he_p.calc_cs_imalpha(w0, "total")

        solver_he_p = PhotoIonizationCoulombPlus(1.0, 2, [1], "velocity")
        solver_he_p.precalc_stex(g_i_he, g_chi_he_L,
			         g_psi0_he_L[1], g_psi0_he_L,
			         eri_method)
        w0 = 0.91794
        (cs1_v, beta) = solver_he_p.calc_one(w0, False, "sigu")
        (cs2_v, beta) = solver_he_p.calc_one(w0, True, "sigu")
        cs_ima_v = solver_he_p.calc_cs_imalpha(w0, "total")

        print 7.1066699, 5.5353368
        print cs1, cs2, -2.0*cs_im_a
        print cs1_v, cs2_v, -2.0*cs_ima_v
        
    def _test_eri(self):
        
        ## this test is written for check segmentation error 
        
        eri_method = ERI_method().use_symmetry(1)
        R0 = 1.4
        ls = [1]
        (sym, g_i, g_psi1, g_psi0_LM, g_chi_LM) = get_basic_basis(R0, ls)
        eri = calc_ERI_complex(g_psi1,eri_method)
        return 0

    def _test_calc_jk(self):
        
        ## this test is written for check segmentation error
        R0 = 1.4
        ls = [1]
        (sym, g_i, g_psi1, g_psi0_LM, g_chi_LM) = get_basic_basis(R0, ls)

        eri_method = ERI_method().use_symmetry(1)
        eri_method = ERI_method()
        mat_i = calc_mat_complex(g_i, True)
        eri_i = calc_ERI_complex(g_i, eri_method)
        mo = calc_RHF(sym, mat_i, eri_i, 2, 20, 0.00001, 0)

        g_i = (SymGTOs.create()
               .sym(sym)
               .sub(sub_solid_sh(0, 0, [1.2, 2.3], sym.irrep_s))
               .atom((0,0,-R0/2.0), 1)
               .atom((0,0,+R0/2.0), 1)
               .setup())
        
        g_psi0 = (SymGTOs.create()
	          .sym(sym)
	          .sub(sub_solid_sh(1, -1, [1.1, 1.5, 2.1], sym.irrep_y))
	          .atom((0,0,+R0/2.0), 1)
	          .atom((0,0,-R0/2.0), 1)
	          .setup())

        g_psi1 = g_psi0.clone()
        eri_J = calc_ERI(g_psi0, g_psi1, g_i, g_i,    eri_method)
        ci = me.VectorXc.Zero(2)
        ci[0] = 1.0; ci[1] = 1.2
        J_C = calc_J(eri_J, ci, 0, g_psi0, g_psi1)
        eri_K = calc_ERI(g_psi0, g_i,    g_i, g_psi1, eri_method)
        K_C = calc_K(eri_K, mo.C[0,0].col(0), 0, g_psi0, g_psi1)
        ir = sym.irrep_y
        jk_C_y = J_C[ir, ir] + K_C[ir, ir]
        print jk_C_y[1, 2]
        
    def _test_H2(self):
        print "H2"

        R0 = 1.4
        ls = [1]
        sym = D2h()
        g_i = get_gi(sym, R0)
        
        eri_method = ERI_method().use_symmetry(1).coef_R_memo(1)
        
        mat = calc_mat_complex(g_i, True)
        eri = calc_ERI_complex(g_i, eri_method)
        mo = calc_RHF(g_i.sym_group, mat, eri, 2, 20, 0.00001, 0)
        print "mo.energy:", mo.energy
        
        (g_psi1, g_psi0_L, g_chi_L) = get_basic_basis(sym, R0, ls)
        solver = PhotoIonizationCoulombPlus(1.0, 2)

        t0 = time.clock(); print 0
        solver.precalc_stex(g_i, g_chi_L, g_psi1, g_psi0_L, ls, eri_method)
        t1 = time.clock(); print 1
        au2ev = 27.2114 

        wev =  43.6
        w = wev / au2ev
        (cs, beta) = solver.calc_one(w, ls, False, "total")
        print wev, cs, beta

        wev =  70.8
        w = wev / au2ev
        (cs, beta) = solver.calc_one(w, ls, False, "total")
        print wev, cs, beta

        wev =  125.2
        w = wev / au2ev
        (cs, beta) = solver.calc_one(w, ls, False, "total")
        print wev, cs, beta

        t2 = time.clock(); print 2

        print t1-t0, t2-t1
        
        """
        results 2016/8/25
        43.6 0.739582832143 1.97846871416
        70.8 0.103428651007 1.96888856731
        125.2 0.000682559486628 1.66990961502
        
        time ~~ 4min

        reference:
        43.6 1.133 1.919
        70.8 0.232 1.917
        125.2 0.0323 1.91

        results 2016/8/26
        43.6 0.964913783581 1.98146303085
        70.8 0.20289113545 1.97025519216
        125.2 0.0248550006626 1.96179809097

        time = 46.388012 0.028387
        """
        

class Test_coulomb_cck(unittest.TestCase):
    def setUp(self):
        print "c"
        self.zeta_p = [ 0.00256226-0.01559939j, 0.00389597-0.02240632j,
	           0.01229986-0.03080238j,0.03010506-0.04378147j,
	           0.01, 0.0177160054, 0.0313856847, 0.0556028960,
                   0.0985061205, 0.174513496, 0.309168204, 0.547722558,
	           0.970345578, 1.71906475, 3.04549604, 5.39540243,
	           9.55849785, 16.9338400, 30.0000000]

    def _test_complex(self):
        zeta_p = self.zeta_p
        z = 1.0
        k = 1.0
        ene = 0.5 * k * k
        zeta1 = 0.6-0.3j; c1 = 0.4
        zeta2 = 1.2-0.2j; c2 = 0.6
        
        # how to calculate
        # general case:
        # <Im y | s> = <(y - y*)/2i | s> = <y-y* | s> / (-2i)
        # s is real function:
        # = (<y*|s> - <y|s> )/2i = Im[<y*|s>]
        
        for L in [1, 2, 3, 5]:
            for M in [-1, 0, +1]:
                g_psi0 = (SymGTOs.create()
                          .sym(C1())
                          .sub(sub_solid_sh(L, M, zeta_p, 0))
                          .setup())
                g_a = (SymGTOs.create()
                       .sym(C1())
                       .sub(sub_solid_sh(L, M, [1.0], 0))
                       .setup())
                g_b = (SymGTOs.create()
                       .sym(C1())
                       .sub(sub_solid_sh(L, M, [zeta1, zeta2], 0))
                       .setup())
                
                c_npsi0 = calc_coulomb_cck(g_psi0, g_a, 0, z, ene)
                s_psi0_b_C = calc_mat(g_psi0,        g_b, False)["s"][0, 0]
                s_psi0_b_H = calc_mat(g_psi0.conj(), g_b, False)["s"][0, 0]
                c_b = np.array([c1, c2])
                calc_ip = (+np.dot(c_npsi0.conj(), s_psi0_b_H * c_b)
                           -np.dot(c_npsi0,        s_psi0_b_C * c_b)) / (-2.0j)
                rs = np.linspace(0, 10, num=100)
                ref_ip = calc_ip_numint_2(lambda r: coulombf(L, -z/k, r*k),
                                          lambda r: r**(L+1) * np.exp(-zeta1 * r * r),
                                          lambda r: r**(L+1) * np.exp(-zeta2 * r * r),
                                          c1, c2, rs)
                print L, M, calc_ip, ref_ip
                self.assertAlmostEqual(calc_ip, ref_ip, places=3)

    def _test_complex_v(self):
        zeta_p = self.zeta_p
        z = 1.0
        k = 1.0
        ene = 0.5 * k * k
        zeta1 = 0.6-0.3j; c1 = 0.4
        zeta2 = 1.2-0.2j; c2 = 0.6
        
        # how to calculate
        # general case:
        # <Im y | s> = <(y - y*)/2i | s> = <y-y* | s> / (-2i)
        # s is real function:
        # = (<y*|s> - <y|s> )/2i = Im[<y*|s>]
        
        for L in [5]:
            for M in [-1, 0, +1]:
                g_psi0 = (SymGTOs.create()
                          .sym(C1())
                          .sub(sub_solid_sh(L, M, zeta_p, 0))
                          .atom((0,0,0), 1)
                          .setup())
                g_a = (SymGTOs.create()
                       .sym(C1())
                       .sub(sub_solid_sh(L, M, [1.0], 0))
                       .atom((0,0,0), 1)
                       .setup())
                g_b = (SymGTOs.create()
                       .sym(C1())
                       .sub(sub_solid_sh(L, M, [zeta1, zeta2], 0))
                       .atom((0,0,0), z)
                       .setup())
                
                c_npsi0 = calc_coulomb_cck(g_psi0, g_a, 0, z, ene)
                v_psi0_b_C = calc_mat(g_psi0,        g_b, True)["v"][0, 0]
                v_psi0_b_H = calc_mat(g_psi0.conj(), g_b, True)["v"][0, 0]
                c_b = np.array([c1, c2])
                calc_ip = (+np.dot(c_npsi0.conj(), v_psi0_b_H * c_b)
                           -np.dot(c_npsi0,        v_psi0_b_C * c_b)) / (-2.0j)
                rs = np.linspace(0, 10, num=100)
                eps = 0.000000001
                ref_ip = calc_ip_numint_2(lambda r: -z/(r+eps) * coulombf(L, -z/k, r*k),
                                          lambda r: r**(L+1) * np.exp(-zeta1 * r * r),
                                          lambda r: r**(L+1) * np.exp(-zeta2 * r * r),
                                          c1, c2, rs)
                print L, M, calc_ip, ref_ip
                self.assertAlmostEqual(calc_ip, ref_ip, places=4)

    def _test_real(self):
        zeta_p = self.zeta_p
        z = 1.0
        k = 1.0
        ene = 0.5 * k * k
        zeta_b = 1.1        

        print ""
        for L in [1, 2, 3]:
            for M in [-1,0,1]:
                g_psi0 = (SymGTOs.create()
                          .sym(C1())
                          .sub(sub_solid_sh(L, M, zeta_p, 0))
                          .setup())
                g_a = (SymGTOs.create()
                       .sym(C1())
                       .sub(sub_solid_sh(L, M, [1.0], 0))
                       .setup())
                g_b = (SymGTOs.create()
                       .sym(C1())
                       .sub(sub_solid_sh(L, M, [zeta_b], 0))
                       .setup())

                c_npsi0 = calc_coulomb_cck(g_psi0, g_a, 0, z, ene)
                s_psi0_b = calc_mat(g_psi0, g_b, False)["s"][0, 0].col(0)
                calc_ip = np.dot(c_npsi0, s_psi0_b).imag
                rs = np.linspace(0, 10, num=100)
                ref_ip  = calc_ip_numint(lambda r: coulombf(L, -z/k, r*k),
                                         lambda r: r**(L+1) * np.exp(-zeta_b * r * r),
                                         rs)

                msg = """
                (L, M) = ({0}, {1}) 
                calc = {2}
                ref  = {3}
                err  = {4}
                """.format(L, M, calc_ip, ref_ip, abs(calc_ip-ref_ip))
                print L, M, calc_ip, ref_ip, abs(calc_ip-ref_ip)
                self.assertAlmostEqual(calc_ip, ref_ip, places=3, msg=msg)


if __name__ == '__main__':
    unittest.main()    
    

