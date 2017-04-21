import unittest
import random
from part_wave import *
from symmolint import *

class Test_s2y(unittest.TestCase):
    def setUp(self):
        pass

    def test_d(self):
        s = {}
        s[-2] = 1.0
        s[-1] = 1.1
        s[0]  = 1.2
        s[1]  = 1.3
        s[2]  = 1.4

        y = s2y(s, 2)
        self.assertAlmostEqual(y[-2], 1.0/np.sqrt(2.0) * (s[2]-1.0j*s[-2]))
        self.assertAlmostEqual(y[-1], 1.0/np.sqrt(2.0) * (s[1]-1.0j*s[-1]))
        self.assertAlmostEqual(y[0], s[0])
        self.assertAlmostEqual(y[+1], -1.0/np.sqrt(2.0) * (s[1]+1.0j*s[-1]))
        self.assertAlmostEqual(y[+2], 1.0/np.sqrt(2.0) * (s[2]+1.0j*s[-2]))

        self.assertAlmostEqual(s[2],  (-1)**2 * y[2].real * np.sqrt(2))
        self.assertAlmostEqual(s[-2], (-1)**2 * y[2].imag * np.sqrt(2))
        self.assertAlmostEqual(s[1],  (-1)**1 * y[1].real * np.sqrt(2))
        self.assertAlmostEqual(s[-1], (-1)**1 * y[1].imag * np.sqrt(2))
    
    def test_sp(self):

        ss0 = {}
        ss0[0,0] = 1.0
        ss0[0,1] = 1.1
        ss0[0,-1] = 1.2

        yy0_ref = {}
        yy0_ref[0,0] = ss0[0,0]
        yy0_ref[0,1]  = -np.sqrt(0.5) * (ss0[0,1]  + 1.0j * ss0[0,-1])
        yy0_ref[0,-1] = +np.sqrt(0.5) * (ss0[0,1]  - 1.0j * ss0[0,-1]) 
        
        yy0_calc = ss0_2_yy0(ss0, 0, 1)

        for m in [-1,0,1]:
            self.assertAlmostEqual(yy0_ref[0,m], yy0_calc[0,m])

    def test_ps(self):

        ss0 = {}
        ss0[-1,0] = 1.0
        ss0[0,0] = 1.1
        ss0[+1,0] = 1.2

        yy0_ref = {}
        yy0_ref[0,0] = ss0[0,0]
        yy0_ref[1,0]  = -np.sqrt(0.5) * (ss0[1,0]  - 1.0j * ss0[-1,0])
        yy0_ref[-1,0] = +np.sqrt(0.5) * (ss0[1,0]  + 1.0j * ss0[-1,0]) 
        
        yy0_calc = ss0_2_yy0(ss0, 1, 0)

        for m in [-1,0,1]:
            self.assertAlmostEqual(yy0_ref[m,0], yy0_calc[m,0])

    def test_exception(self):
        ss0 = {}
        ss0[1, 0] = 1.1
        ss0[0, 0] = 0.1
        ss0[-1,0] = 0.2
        self.assertRaises(Exception, lambda: ss0_2_yy0(ss0, 2, 0))
        self.assertRaises(Exception, lambda: ss0_2_yy0(ss0, 1, 1))

    def test_compare_with_old(self):
        L = 3
        Lp = 1
        Ms = range(-L, L+1)
        Mps = range(-Lp, Lp+1)
        ss0 = { (M, Mp): random.random() for M in Ms for Mp in Mps}
        yy0_new = ss0_2_yy0(ss0, L, Lp)
        yy0_old = ss0_2_yy0_old(ss0, L, Lp)
        for M in Ms:
            for Mp in Mps:
                msg =  "key : ({0}, {1}). \n".format(M, Mp)
                msg += "a : {0}\n".format(yy0_new[M, Mp])
                msg += "b : {0}\n".format(yy0_old[M, Mp])
                self.assertAlmostEqual(yy0_new[M, Mp],
                                       yy0_old[M, Mp],
                                       msg = msg)


class Test_PartialWave(unittest.TestCase):
    ## copied to test_solve_driv.py
    def setUp(self):
        pass

    def test_orthonormality(self):
        zeta = [1.1-0.3j]
        g = (SymGTOs.create().sym(Cs())
             .sub(sub_solid_sh(0, 0, zeta, 0))
             .sub(sub_solid_sh(1,-1, zeta, 0))
             .sub(sub_solid_sh(1, 0, zeta, 0))
             .sub(sub_solid_sh(1, 1, zeta, 0))
             .sub(sub_solid_sh(2,-2, zeta, 0))
             .sub(sub_solid_sh(2,-1, zeta, 0))
             .sub(sub_solid_sh(2, 0, zeta, 0))
             .sub(sub_solid_sh(2,+1, zeta, 0))
             .sub(sub_solid_sh(2,+2, zeta, 0))
             .sub(sub_solid_sh(3,-1, zeta, 0))
             .sub(sub_solid_sh(3, 0, zeta, 0))
             .sub(sub_solid_sh(3,+1, zeta, 0))
             .sub(sub_solid_sh(5,-1, zeta, 0)) #12
             .sub(sub_solid_sh(5, 0, zeta, 0)) #13
             .sub(sub_solid_sh(5,+1, zeta, 0)) #14
             .setup())
        s = calc_mat_complex(g, False)["s"][0, 0]
        n = s.rows()
        for i in range(0, n):
            self.assertAlmostEqual(1.0, s[i, i])
            for j in range(0, i):
                msg = "(i,j)=({0},{1})".format(i, j)
                self.assertAlmostEqual(0.0, s[i, j], msg=msg)

    def test_solid_part_waves(self):
        zeta = [1.1-0.3j]
        sym = D2h()

        for L in [1,2,3, 5]:
            gs = solid_part_waves(D2h(), L, [-1,0,1], zeta)
            for (irrep, M) in [(sym.irrep_y, -1),
                               (sym.irrep_z, 0),
                               (sym.irrep_x, 1)]:
                g = (SymGTOs.create().sym(sym)
                     .sub(sub_solid_sh(L, M, zeta, irrep))
                     .setup())
                s = calc_mat(g, gs, False)["s"][irrep, irrep]
                msg = "(L,M) = ({0},{1})".format(L, M)
                self.assertAlmostEqual(1.0, s[0,0], msg=msg)

    def test_exception(self):
        zeta = [1.1-0.3j]
        self.assertRaises(lambda: solid_part_waves(Cs(), 1, [0], zeta))
        self.assertRaises(lambda: solid_part_waves(D2h(), 1, [2], zeta))

                        

if __name__ == '__main__':
    unittest.main()
