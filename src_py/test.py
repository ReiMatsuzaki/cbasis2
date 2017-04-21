import os
import sys
import unittest
from l2func import *

# sys.path.append("/Users/rei/src/git/opt_cbf/py_bind/nnewton")
# from nnewton import *

class TestL2Func(unittest.TestCase):
    def setUp(self):
        pass

    def test_func(self):
        s1 = STO(0.1, 2, 0.3)
        s2 = STO(0.1, 2, 0.3)
        g1 = GTO(0.1, 2, 0.3)
        g2 = STO(0.1, 2, 0.3)

        psi = s1 * 0.2 + (0.1 * g1 - 0.3 * g2) * 0.1
        x = 0.7
        self.assertAlmostEqual(psi.at(x),
                               s1.at(x) * 0.2 + (0.1*g1.at(x)-0.3*g2.at(x)) * 0.1)

        print s1
        print g1

    def test_CIP_d2(self):
        g1 = GTO(0.2, 2, 0.3)
        g2 = GTO(0.3, 4, 0.1)
        self.assertAlmostEqual(cip(g1, dr(2), g2),
                               cip(g2, dr(2), g1))

    def test_CIP_op(self):
        s1 = STO(1.0, 2, 1.1)
        g1 = GTO(1.1, 2, 1.2)
        self.assertAlmostEqual(cip(s1, dr(2) + 0.2*rm(2), g1),
                               cip(s1, D2(), g1)+ 0.2 * cip(s1, Rm(2), g1))

    def test_CIP(self):
        s1 = STO(0.1, 2, 0.3) 
        g1 = GTO(0.2, 1, 0.2)

        s = s1 + 0.2 * s1
        g = g1
        #op = 0.5 * D2() + (Rm(1) + 0.1 * Rm(2)) * 0.2
        op = 0.5*dr(2) + (rm(1) + 0.1*rm(2))*0.2

        print 0.4 * (s1 + 0.2*g1)

        self.assertAlmostEqual(cip(s, op, g),
                               0.5 * 1.2 * cip(s1, dr(2), g) +
                               0.2 * 1.2 * cip(s1, rm(1), g) +
                               0.2 * 1.2 * 0.1 * cip(s1, rm(2), g))

    def test_CIP_sss(self):
        s1 = STO(0.1, 2, 0.3)
        s2 = STO(0.2, 3, 0.4)
        s3 = STO(0.3, 1, 0.2)
        s12= STO(0.02, 5, 0.7)
        self.assertAlmostEqual(cip(s1, s2, s3), cip(s12, s3))

    def test_ricatti_functions(self):
        k = 1.1
        r = 1.2
        rh1 = rhankel(k, 0)
        self.assertAlmostEqual(np.exp(1.0j*k*r), rh1.at(r))

        rb1 = rbessel(k, 0)
        self.assertAlmostEqual(np.sin(k*r), rb1.at(r))
        
    def test_op_apply(self):
        g2 = FuncZero() + GTO(0.2, 1, 0.2) + 0.2 * GTO(0.1, 3, 0.1)
        g1 = GTO(0.2, 2, 0.3)
        self.assertAlmostEqual(cip(g1, dr(1), g2),
                               cip(g1, op_apply(dr(1), g2)))

    def test_add_sto(self):
        s1 = STO(1.0, 2, 1.1)
        s2 = STO(0.2, 3, 0.2)
        g1 = GTO(1.1, 2, 1.2)
        
        self.assertAlmostEqual(cip(FuncAdd(s1, s2), g1),
                               cip(s1, g1) + cip(s2, g1))
        self.assertAlmostEqual(cip(s1 + s2, g1),
                               cip(s1, g1) + cip(s2, g1))

    def test_cnorm(self):
        s = 0.2 * STO(0.1, 2, 0.3) + 0.2*(0.3 * GTO(0.2, 3, 0.1) + STO(0.4, 3, 0.2))
        n = cnormalize(s)
        self.assertAlmostEqual(1.0, cnorm2(n))

    def test_linear_combination(self):
        us = [STO(1.0, 2, n*0.3+0.2) for n in range(3)]
        cs = [(n+1)*0.1 for n in range(3)]
        
        expe = cs[0]*us[0] + cs[1]*us[1] + cs[2]*us[2]
        calc = linear_combination(cs, us)

        r0 = 0.7
        self.assertAlmostEqual(expe.at(r0), calc.at(r0))

    def test_sum_func(self):
        fs = [GTO(1.0, n, 1.0) for n in [2,3,4]]
        sum_fs = sum_func(fs)
        
        r0 = 0.7
        self.assertAlmostEqual(fs[0].at(r0) + fs[1].at(r0) + fs[2].at(r0),
                               sum_fs.at(r0))

class TestCutExpFunc(unittest.TestCase):
    def setUp(self):
        pass

    def test_cip(self):
        s1 = CutSTO(1.2, 2, 2.5, 10.0)
        s2 = CutSTO(1.1, 3, 1.5, 10.0)
        sol = 0.0386718749998404;
        self.assertAlmostEqual(cip(s2, s1), cip(s1, s2))
        self.assertAlmostEqual(sol, cip(s1, s2))
        
        r0 = 4.0
        cut_s = CutSTO(1.1, 2, 0.2, r0)
        s =        STO(1.1, 2, 0.2)
        x = 2.0
        self.assertAlmostEqual(s.at(x), cut_s.at(x))
        self.assertAlmostEqual(s.at(r0), cut_s.at(r0))
        
        x = 6.0
        self.assertAlmostEqual(0.0, cut_s.at(x))

    def test_cip_op(self):
        s1 = CutSTO(1.2, 2, 2.5, 10.0)
        s2 = CutSTO(1.3, 3, 0.1, 10.0)

        r2_s1 = CutSTO(1.2, 4, 2.5, 10.0)

        self.assertAlmostEqual(cip(s2, r2_s1),
                               cip(s2, rm(2), s1))

    def test_cut_op(self):
        r0 = 10.0
        cs1 = CutSTO(1.2, 2, 2.5, r0)
        cs2 = CutSTO(1.3, 3, 0.1, r0)
        s1 = STO(1.2, 2, 2.5)
        s2 = STO(1.3, 3, 0.1)
        self.assertAlmostEqual(cip(cs1, cs2),
                               cip(op_apply(Cut(r0), s1), 
                                   op_apply(Cut(r0), s2)))

class TestHatom(unittest.TestCase):
    def setUp(self):
        pass

    def test_eigen(self):
        hatom = HAtom(1.0)
        f10 = hatom.eigenfunc(1, 0)
        f20 = hatom.eigenfunc(2, 0)
        f21 = hatom.eigenfunc(2, 1)
        h0 = hatom.hop(0)
        h1 = hatom.hop(1)
        print h0
        self.assertAlmostEqual(1.0, cnorm2(f10))
        self.assertAlmostEqual(1.0, cnorm2(f20))
        self.assertAlmostEqual(1.0, cnorm2(f21))
        self.assertAlmostEqual(0.0, cip(f10, f20))
        self.assertAlmostEqual(-0.5, cip(f10, h0, f10))
        self.assertAlmostEqual(-0.0, cip(f10, h0, f20))
        self.assertAlmostEqual(-0.0, cip(f20, h0, f10))
        self.assertAlmostEqual(-0.125, cip(f21, h1, f21))
        
    def test_length(self):
        f01 = HAtom(1.0).length(1, 0, 1) # 1s->kp
        ref = STO(2.0, 2, 1.0)
        self.assertAlmostEqual(f01.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).length(2, 1, 0) # 2p->ks
        ref = STO(1.0/(2.0*np.sqrt(6.0)), 3, 0.5)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).length(2, 1, 2) # 2p->kd
        ref = STO(1.0/(2.0*np.sqrt(6.0)), 3, 0.5)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).length(3, 2, 1) # 3d->kp
        ref = STO(4.0/(81.0*np.sqrt(30.0)), 4, 1.0/3.0)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).length(3, 2, 3) # 3d->kf
        ref = STO(2.0/81.0*np.sqrt(2.0/15.0), 4, 1.0/3.0)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))

    def test_velocity(self):
        f01 = HAtom(1.0).velocity(1, 0, 1) # 1s->kp
        ref = STO(-2.0, 1, 1.0)
        self.assertAlmostEqual(f01.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).velocity(2, 1, 0) # 2p->ks
        ref = STO(3.0/(2.0*np.sqrt(6.0)), 1, 0.5) + STO(-1.0/(4.0*np.sqrt(6.0)), 2, 0.5)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))
        
        calc = HAtom(1.0).velocity(2, 1, 2) # 2p->kd
        ref = STO(-1.0/(4.0*np.sqrt(6.0)), 2, 0.5)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).velocity(3, 2, 1) # 3d->kp
        ref = STO(2.0*np.sqrt(30.0)/243.0, 2, 1.0/3.0)+STO(-2.0*np.sqrt(30.0)/3645.0, 3, 1.0/3.0)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))

        calc = HAtom(1.0).velocity(3, 2, 3) # 3d->kf
        ref = STO(-2.0/243.0*np.sqrt(2.0/15.0), 3, 1.0/3.0)
        self.assertAlmostEqual(calc.at(2.2), ref.at(2.2))
        
class TestDBasis(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_d0(self):
        z0 = 1.0-0.3j
        s0 = d0_basis(STO, 2, z0)
        self.assertAlmostEqual(1.0, cnorm(s0))

        g0 = d0_basis(GTO, 2, z0)
        self.assertAlmostEqual(1.0, cnorm(g0))

    def test_d1(self):
        r0 = 4.0
        z0 = 0.6-0.3j
        h = 0.0001

        def d0_at(z):
            return d0_basis(STO, 2, z).at(r0)

        
        expe = (d0_at(z0+h) - d0_at(z0-h) + 1.0j*d0_at(z0-1.0j*h) - 1.0j*d0_at(z0+1.0j*h))/(4*h)
        calc = d1_basis(STO, 2, z0).at(r0)
        print "expe, calc: ", expe, calc
        self.assertAlmostEqual(expe, calc)


if __name__ == '__main__':
    unittest.main()



