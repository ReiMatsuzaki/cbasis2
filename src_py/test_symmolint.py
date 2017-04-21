import numpy as np
from numpy import pi, exp, sqrt
from symmolint import *

import unittest
import minieigen as me
import scipy.linalg as la

class Test_JK(unittest.TestCase):
    def setUp(self):
        pass

    def _test_jk0(self):
        ## this test is for checking segmentation error
        xatom = (0, 0, 0)
        sym = Cs()
        Ap = sym.get_irrep("A'")
        App= sym.get_irrep("A''")
        zeta0 = [0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694]
        zeta1 = zeta0

        gtos = (SymGTOs.create()
                .sym(sym)
                .sub(SubSymGTOs()
                     .xyz(xatom).ns((0, 0, 0))
                     .rds(Reduction(Ap,  [[1]]))
                     .zeta(zeta0))
                .sub(SubSymGTOs()
                     .xyz(xatom).ns((0, 0, 1))
                     .rds(Reduction(App, [[1]]))
                     .zeta(zeta1))
                .atom(xatom, 2.0)
                .setup())

        mat_set = calc_mat_complex(gtos, True)
        eri = calc_ERI_complex(gtos, ERI_method().use_symmetry(1))

        w = 1.0
        # mo = calc_RHF(sym, mat_set, eri, 2, 20, 0.00001,  0)
        # jk = calc_JK(eri, mo.C, 0, 0, 1.0, 1.0)
        C = BMat()
        c00 = me.MatrixXc.Zero(len(zeta0), len(zeta0))
        C.set_matrix(0, 0, c00)
        c11 = me.MatrixXc.Zero(len(zeta1), len(zeta1))
        C.set_matrix(1, 1, c11)
        jk = calc_JK(eri, C, 0, 0, 1.0, 1.0)
    
    def _test_jk(self):
        pass
    """
        sym = Cs()
        g_i = SymGTOs.create().sym(sym)
        g_0 = SymGTOs.create().sym(sym)
        g_1 = SymGTOs.create().sym(sym)
        g_f = SymGTOs.create().sym(sym)
        
        zeta_i = [0.4, 1.0]
        sub_i = (SubSymGTOs()
                 .xyz((0,0,0))
                 .ns((0,0,0))
                 .zeta(zeta_i)
                 .rds(Reduction(0, [[1]])))
        g_i.sub(sub_i)
        g_f.sub(sub_i)
        
        zeta_0 = [0.1, 0.2]
        sub_0 = (SubSymGTOs()
                 .xyz((0,0,0))
                 .ns((0,0,1))
                 .zeta(zeta_0)
                 .rds(Reduction(1, [[1]])))
        g_0.sub(sub_0)
        g_f.sub(sub_0)

        zeta_1 = [0.03, 1.1]
        sub_1 = (SubSymGTOs()
                 .xyz((0,0,0))
                 .ns((0,0,1))
                 .zeta(zeta_1)
                 .rds(Reduction(1, [[1]])))
        g_1.sub(sub_1)
        g_f.sub(sub_1)        

        g_i.setup()
        g_0.setup()
        g_1.setup()
        g_f.setup()
        
        C = BMat()
        c00 = me.MatrixXc.Zero(2, 2)
        c00[0,0] = 1.0; c00[0,1] = 1.2
        c00[1,0] = 0.2; c00[1,1] = 0.3
        C.set_matrix(0, 0, c00)
        c11 = me.MatrixXc.Zero(4, 4)
        c11[0,0] = 1.0; c11[0,1] = 1.2
        c11[1,0] = 0.2; c11[1,1] = 0.3
        C.set_matrix(1, 1, c11)
        
        method = ERI_method()

        coef_J = 1.1
        coef_K = 1.2

        eri_f = calc_ERI_complex(g_f, method)
        J_f = calc_JK(eri_f, C, 0, 0, coef_J, coef_K)

        eri_J = calc_ERI(g_0, g_1, g_i, g_i, method)
        J_01 = calc_J(eri_J, C[0,0].col(0), 0, g_0, g_1, coef_J)

        eri_K = calc_ERI(g_0, g_i, g_i, g_1, method)
        K_01 = calc_K(eri_K, C[0,0].col(0), 0, g_0, g_1, coef_K)

        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(J_f[1,1][i,j+2],
                                       J_01[1,1][i,j] +
                                       K_01[1,1][i,j],
                                       msg = "{0}, {1}".format(i,j))
    """

class Test_minieigen(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_print(self):
        v = me.VectorXc([5, 6])
        x = me.MatrixXc([[1, 2], [3, 4]])
        xv = x*v

class Test_ceig(unittest.TestCase):
    def setUp(self):
        pass

    def test_ceig(self):
        h = me.MatrixXc(
            [[1.0+0.1j,  0.2,   0.1-0.1j ],
             [0.2,       0.3,   0.4      ],
             [0.1-0.1j,  0.4,   1.0-0.2j ]])
        s = me.MatrixXc(
            [[1.0,       0.2,   0.2+0.00j ],
             [0.2,       1.0,   0.1       ],
             [0.2+0.00j, 0.1,   1.0       ]])
        (eigs, eigvec) = ceig(h, s)
        for i in range(h.cols()):
            lexp = h*eigvec.col(i)
            rexp = s*eigvec.col(i)*eigs[i]
            self.assertAlmostEqual(0.0, abs(lexp-rexp))

        for i in range(eigs.rows()-1):
            self.assertTrue(eigs[i].real < eigs[i+1].real)

    def test_ceig_canonical(self):
        eps = 0.000001
        h2 = me.MatrixXc([[1.0, 0.5],
                          [0.5, 1.5]])
        s2 = me.MatrixXc([[1.4, 0.2],
                          [0.2, 0.4]])
        h3 = me.MatrixXc([[1.0,     0.5,     0.5+eps],
                          [0.5,     1.5,     1.5+eps],
                          [0.5+eps, 1.5+eps, 1.5+eps]])
        s3 = me.MatrixXc([[1.4,     0.2,     0.2+eps],
                          [0.2,     0.4,     0.4+eps],
                          [0.2+eps, 0.4+eps, 0.4+eps]])
        (eigs, eigvec) = ceig_canonical(h3, s3, 2)

        for i in range(2):
            lexp = h3 * eigvec.col(i)
            rexp = s3 * eigvec.col(i)*eigs[i]
            for j in range(2):
                msg = "{0}, {1}\n lexp={2}\nrexp={3}".format(i,j,lexp[j],rexp[j])
                self.assertAlmostEqual(lexp[j], rexp[j], msg=msg, places=5)

            
class Test_BMatSet(unittest.TestCase):
    def setUp(self):
        pass

    def test_getset(self):
        cs = Cs()
        print type(cs)
        print cs.get_irrep("A'")
        bmat = BMatSet(cs.order)
        s00 = me.MatrixXc.Zero(4, 4)
        s00[0, 1] = 1.1; s00[1, 0] = 1.2
        Ap = cs.get_irrep("A'")
        bmat.set_matrix("s", Ap, Ap, s00)
        s00_get1 = bmat.get_matrix("s", Ap, Ap)
        s00_get2 = bmat.get_matrix("s", Ap, Ap)
        s00_get1[0, 1] = 2.1
        self.assertAlmostEqual(1.2, s00_get2[1, 0])
        self.assertAlmostEqual(2.1, s00_get2[0, 1])
        
        
class Test_SymMolInt(unittest.TestCase):

    def setUp(self):
        pass
    
    def test_symmetry_group(self):
        sym = Cs()
        self.assertEqual(2, sym.order)
        self.assertEqual("Cs", sym.name)

    def test_reduction_sets(self):
        """
        coef_iat_ipn = me.MatrixXc.Zero(2, 3)
        coef_iat_ipn[0, 1] = 1.1
        coef_iat_ipn[1, 0] = 1.2
        """
        Ap = Cs().get_irrep("A'")
        # App = Cs().get_irrep("A''")
        coef_iat_ipn = [[0.0, 1.1], [1.2, 0.0]]
        rds1 = Reduction(Ap, coef_iat_ipn)
        self.assertEqual(Ap, rds1.irrep)
        self.assertAlmostEqual(1.1, rds1.coef_iat_ipn()[0, 1])
        self.assertAlmostEqual(1.2, rds1.coef_iat_ipn()[1, 0])

    def test_sub_ms(self):
        sym = Cs()
        atom = Atom("A", 1.0, [(0,0,0)])
        for L in range(1,4):
            sub = sub_solid_sh(sym, atom, L, [-1,0,1], [1.1-0.2j])

    def test_at_r_center(self):
        sym = Cs()
        Ap = sym.get_irrep("A'")
        z = 1.3
        c = 1.2
        r = 2.5
        
        ## -- see support/normalized_gto.py --
        expo=exp(-r*r*z)
        y_ref = {}
        y_ref[0]=c*2*2**(0.75)/(pi**(0.25)*sqrt(z**(-1.5)))*r*expo
        y_ref[1]=c*4*2**(0.75)*sqrt(3)/(3*pi**(0.25)*sqrt(z**(-2.5)))*r**2*expo
        y_ref[2]=c*8*sqrt(15.0)*2**(0.75)/(15*pi**(0.25)*sqrt(z**(-3.5)))*r**3*expo
        y_ref[3]=c*16*sqrt(105.0)*2**(0.75)/(105*pi**(0.25)*sqrt(z**(-4.5)))*r**4*expo

        mole = Molecule(sym, [Atom("cen", 0, [(0, 0, 0)])])
        cen  = mole.atom("cen")
        for L in [0,1,2,3]:
            M = 0
            gtos = (SymGTOs(mole)
                    .sub(sub_solid_sh(sym, cen, L, M, [z]))
                    .setup())
            irrep = sym.irrep_s if L==0 else sym.irrep_z
            (y0, dy0) = gtos.at_r_ylm(L, M, irrep, [c], [r])
            msg = "L={0}\ncalc={1}\nref={2}".format(L, y0[0], y_ref[L])
            self.assertAlmostEqual(y_ref[L], y0[0], msg=msg)
        
    def test_clone(self):
        sym = Cs()
        Ap = sym.get_irrep("A'")
        zs = [2.0**n-0.1j for n in range(-2,2)]
        mole = Molecule(sym, [Atom("A", 1.1, [(0, 0, 0)])])
        gtos = (SymGTOs(mole)
                .sub(SubSymGTOs(sym, mole.atom("A"))
                     .ns( (0, 0, 0))
                     .rds(Reduction(Ap, [[1]]))
                     .zeta(zs))
                .setup())
        gtos_copy = gtos.clone()
        gtos_cc   = gtos.conj()
        mat = calc_mat_complex(gtos, False)
        mat_copy = calc_mat_complex(gtos, False)
        mat_cc = calc_mat_complex(gtos, False)
        self.assertAlmostEqual(mat.get_matrix("v", 0, 0)[0, 1], 
                               mat_copy.get_matrix("v", 0, 0)[0, 1])
        self.assertAlmostEqual(mat_cc.get_matrix("v", 0, 0)[0, 1].conjugate(), 
                               mat_copy.get_matrix("v", 0, 0)[0, 1])

        
    def test_calc_v(self):
        xyz = (0.3, 0.4, 0.1)
        zs = [2.0**n-0.1j for n in range(-2, 2)]
        sym = Cs()
        mole = Molecule(sym,
                        [Atom("A", 0.0, [(0,0.2,0)]),
                         Atom("B", 1.1, [(0,0,0)])])
        gtos = (SymGTOs(mole)
                .sub(SubSymGTOs(sym, mole.atom("A"))
                     .ns( (1, 0, 0))
                     .rds(Reduction(0, [[1]]))
                     .zeta(zs))
                .setup())
        tmp_v1 = calc_mat_hermite(gtos, True)
        tmptmp_v1 = tmp_v1["v"]
        v1 = tmptmp_v1[0, 0]
        v2 = calc_mat(gtos.conj(), gtos, True)["v"][0, 0]
        self.assertAlmostEqual(v1[1, 2], v2[1, 2])
        """
        gtos = (SymGTOs()
                .sym(Cs())
                .sub(SubSymGTOs()
                     .xyz((0, 0.2, 0))
                     .ns( (1, 0, 0))
                     .rds(Reduction(0, [[1]]))
                     .zeta(zs))
                .atom(xyz, q)
                .setup())
        gtos2 = (SymGTOs.create()
                 .sym(Cs())
                 .sub(SubSymGTOs()
                      .xyz((0, 0.2, 0))
                      .ns( (1, 0, 0))
                      .rds(Reduction(0, [[1]]))
                      .zeta(zs))
                 .atom((0,0,0), 0)
                 .setup())
        v1 = calc_mat_hermite(gtos, True)["v"][0, 0]
        v2 = calc_v_mat(gtos2.conj(), gtos2, xyz, q)[0, 0]
        """
                                                    


class Test_H_atom(unittest.TestCase):

    def setUp(self):
        xatom = (0, 0, 0)
        sym = Cs()
        mole = Molecule(sym, [Atom("H", 1.0).add(0, 0, 0)])
        hatom = mole.atom("H")
        Ap = sym.get_irrep("A'")
        App= sym.get_irrep("A''")
        zeta = [2.0**n for n in range(-10, 10)]
        
        self.gtos = (SymGTOs(mole)
                     .sub(sub_solid_sh(sym, hatom, 0, 0, zeta))
                     .sub(sub_solid_sh(sym, hatom, 1, 0, zeta)))
        mat = calc_mat_complex(self.gtos, True)
        
        s = mat.get_matrix("s", Ap, Ap)
        t = mat.get_matrix("t", Ap, Ap)
        v = mat.get_matrix("v", Ap, Ap)
        h = t + v
        (self.eigs, self.eigvecs) = ceig(h, s)
     
    def test_energy(self):
        
        for i in range(0, 4):
            n = i + 1
            ene = -0.5/(n*n);
            self.assertTrue(abs(ene-self.eigs[i])   < 0.0001);

    def test_wavefunc(self):
        c = self.eigvecs.col(0)
        c = self.gtos.correct_sign(0, 0, 0, c)
        rs = [1.1]
        Ap = 0
        (ys_calc, ds_calc) = self.gtos.at_r_ylm(0, 0, Ap, c, rs)
        ys_refs = [2.0*r*np.exp(-r) for r in rs]
        self.assertTrue(abs(ys_calc[0] - ys_refs[0]) < 0.0001)

    def test_2p(self):
        App= Cs().get_irrep("A''")
        mat = calc_mat_complex(self.gtos, True)
        s = mat.get_matrix("s", App, App)
        t = mat.get_matrix("t", App, App)
        v = mat.get_matrix("v", App, App)
        (eigs, eigenvecs) = ceig(t+v, s)

        self.assertTrue(abs(eigs[0]+0.125) < 0.000001)

        c = self.gtos.correct_sign(1, 0, App, eigenvecs.col(0))
        rs = [1.1]
        (ys_calc, ds_calc) = self.gtos.at_r_ylm(1, 0, App, c, rs)

        ys_refs = [(1.0/(2.0*np.sqrt(6.0)))*r*r*np.exp(-0.5*r) for r in rs]
        self.assertTrue(abs(ys_calc[0] - ys_refs[0]) < 0.0001)        
        

class Test_H2_plus(unittest.TestCase):

    def setUp(self):
        self.sym = Cs()
        
        mole = Molecule(self.sym,
                        [Atom("H", 1.0, [(0,0,1), (0,0,-1)])])
        Ap = self.sym.get_irrep("A'")
        self.gtos = SymGTOs(mole,
                            [("H",
                              [(0,0,0),
                               (0,0,1)],
                              [Reduction(Ap, [[1,0],
                                              [1,0]]),
                               Reduction(Ap, [[0,+1],
                                              [0,-1]])],
                              [2.0**n for n in range(-10,10)])])

    def test_energy(self):
        """
        reference energy is from 
        calc/ccolumbus/h2/theta_traj/lie_2s3p_4s6dcstong/out/theta10_10g/out
        in rcclsc
        """
        self.gtos.setup()

        mat = calc_mat_complex(self.gtos, True)
        Ap = self.sym.get_irrep("A'")
        s = mat.get_matrix("s", Ap, Ap)
        t = mat.get_matrix("t", Ap, Ap)
        v = mat.get_matrix("v", Ap, Ap)
        h = t + v
        (eigs, eigvecs) = ceig(h, s)
        self.assertAlmostEqual(-1.1026342144949, eigs[0], places=3,
                               msg = "eigs[0] = {0}".format(eigs[0]))

class Test_H_photoionization(unittest.TestCase):
    def setUp(self):
        self.calc_full()
    
    def calc_full(self):
        xatom = (0, 0, 0)
        sym = Cs()
        Ap = sym.get_irrep("A'")
        App= sym.get_irrep("A''")
        zeta0 = [2.0**n for n in range(-10, 10)]
        zeta1 = [2.0**n-0.025j for n in range(-15, 5)]
        mole = Molecule(sym, [Atom("H", 1.0, [(0,0,0)])])
        hatom = mole.atom("H")
        self.gtos = (SymGTOs(mole)
                     .sub(sub_solid_sh(sym, hatom, 0, 0, zeta0))
                     .sub(sub_solid_sh(sym, hatom, 1, 0, zeta1)))
        self.gtos.setup()
        mat = calc_mat_complex(self.gtos, True)
        h0 = mat.get_matrix("t", Ap, Ap) + mat.get_matrix("v", Ap, Ap)
        s0 = mat.get_matrix("s", Ap, Ap)
        (eigs0, eigvecs0) = ceig(h0, s0)
        ene0 = eigs0[0]
        c0 = self.gtos.correct_sign(0, 0, Ap, eigvecs0.col(0))
        self.assertAlmostEqual(ene0, -0.5, places=4)
        
        h1 = mat.get_matrix("t", App, App) + mat.get_matrix("v", App, App)
        s1 = mat.get_matrix("s", App, App)

        ## length form
        z10 = mat.get_matrix("z", App, Ap)
        w = 1.0
        m1 = z10 * c0
        c1 = la.solve(h1 - (w+ene0)*s1, m1)
        self.alpha_full = (m1 * c1).sum()
        
        ## velocity form
        dz10 = mat.get_matrix("dz", App, Ap)
        w = 1.0
        m1 = dz10 * c0
        c1 = la.solve(h1-(w+ene0)*s1, m1)
        self.alpha_full_v = (m1 * c1).sum()

    def calc_part(self):
        xatom = (0, 0, 0)
        sym = Cs()
        Ap = sym.get_irrep("A'")
        App= sym.get_irrep("A''")
        zeta0 = [2.0**n for n in range(-10, 10)]
        zeta1 = [2.0**n-0.02j for n in range(-15, 5)]
        self.gtos0 = (SymGTOs.create()
                      .sym(sym)
                      .sub(SubSymGTOs()
                           .xyz(xatom)
                           .ns((0, 0, 0))
                           .rds(Reduction(Ap, [[1]]))
                           .zeta(zeta0))
                      .atom(xatom, 1.0))
        self.gtos1 = (SymGTOs.create()
                      .sym(sym)
                      .sub(SubSymGTOs()
                           .xyz(xatom)
                           .ns((0, 0, 0))
                           .rds(Reduction(App, [[1]]))
                           .zeta(zeta1))
                      .atom(xatom, 1.0)
                      .setup())

        mat0 = calc_mat_complex(self.gtos0, True)
        mat1 = calc_mat_complex(self.gtos1, True)
        mat10= calc_mat(self.gtos1, self.gtos0, False)
        h0 = mat0.get_matrix("t", Ap, Ap) + mat0.get_matrix("v", Ap, Ap)
        s0 = mat0.get_matrix("s", Ap, Ap)
        (eigs0, eigvecs) = ceig(h0, s0)
        ene0 = eigs0[0]
        c0 = eigvecs.col(0)
        h1 = mat1.get_matrix("t", App, App) + mat1.get_matrix("v", App, App)
        s1 = mat1.get_matrix("s", App, App)
        z10 = mat10.get_matrix("z", App, Ap)
        m1 = z10 * c0
        w = 1.0
        c1 = la.solve(h1 - (w+ene0)*s1, m1)
        self.alpha_part = (m1 * c1).sum()
        

    def test_1skp(self):
        # see calc/stoh/1skp/l_5 in rcclsc
        alpha_ref = -1.88562800720386+0.362705406693342j
        self.assertAlmostEqual(self.alpha_full, alpha_ref, places=3)
        #self.assertAlmostEqual(self.alpha_part, self.alpha_part)
        
        # see calc/stoh/1skp/l_5 in rcclsc (E=0.5)
        alpha_ref_v = 0.385628007210641-0.362705406667124j
        self.assertAlmostEqual(self.alpha_full_v,
                               -alpha_ref_v,
                               places=3)


class Test_He(unittest.TestCase):

    def setUp(self):
        xatom = (0, 0, 0)
        self.sym = Cs()
        sym = self.sym
        mole = Molecule(sym)
        hatom = Atom("H", 1.0).xyz(0, 0, 0)
        mole.atom(hatom)
        Ap = sym.get_irrep("A'")
        App= sym.get_irrep("A''")
        zeta0 = [0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694]
        zeta1 = zeta0
        
        self.gtos = (SymGTOs(mole)
                     .sub(SubSymGTOs(sym, hatom)
                          .ns((0, 0, 0))
                          .rds(Reduction(Ap,  [[1]]))
                          .zeta(zeta0))
                     .sub(SubSymGTOs(sym, hatom)
                          .ns((0, 0, 1))
                          .rds(Reduction(App, [[1]]))
                          .zeta(zeta1))
                     .setup())
        
    def _test_ERI(self):
        pass
        """
        gtos = self.gtos
        sym = self.sym
        mat_set = calc_mat_complex(gtos, True)
        eri = calc_ERI_complex(gtos, ERI_method().use_symmetry(1))

        w = 1.0
        mo = calc_RHF(sym, mat_set, eri, 2, 20, 0.00001,  0)
        self.assertAlmostEqual(-2.8617, mo.energy, 4)
        # h_se = calc_SEHamiltonian(mo, eri, 0, 0)
        # alpha = calc_alpha(mo, mat_set, 0, 0, h_se, w)
        # cs = pi_total_crosssection(alpha, w, 2)

        self.assertEqual([1, 0], mo.num_occ())
        self.assertAlmostEqual(-0.91795, mo.eigs[0][0], 4)

        eri.write("eri.bin")
        eri2 = ERI_read("eri.bin")

        # ==== check for JK ====
        jk = calc_JK(eri, mo.C, 0, 0, 1.0, 1.0)
        h0 = calc_SEHamiltonian(mo, eri, 0, 0)[1, 1]
        h1 = (mat_set.get_matrix("t", 1, 1) +
              mat_set.get_matrix("v", 1, 1) +
              jk[1, 1])
        """
        
if __name__ == '__main__':
    unittest.main()
        
