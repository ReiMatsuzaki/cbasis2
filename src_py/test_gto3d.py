import os
import sys
import unittest
from l2func import *
import numpy as np
import scipy.linalg as la

class TestGtoSet(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_create(self):
        basis_set = SphericalGTOSet()
        xyz = (0.0, 0.0, 0.0)
        r0 = 10.0
        for n in range(-10, 10):
            z = 2.0**n
            basis_set.add_one_basis(0, 0, xyz, z)
            for L in [0]:
                basis_set.add_basis(L, (0.0, 0.0, +r0), z)
                basis_set.add_basis(L, (0.0, 0.0, -r0), z)
                basis_set.add_basis(L, (0.0, +r0, 0.0), z)
                basis_set.add_basis(L, (0.0, -r0, 0.0), z)
                basis_set.add_basis(L, (+r0, 0.0, 0.0), z)
                basis_set.add_basis(L, (-r0, 0.0, 0.0), z)

        smat = basis_set.s_mat()
        for e in sorted(abs(la.eig(smat)[0]))[0:10]:
            print e

        print "-----"
        hmat = basis_set.t_mat() + basis_set.v_mat(1.0, xyz)
        zmat = basis_set.xyz_mat((0, 0, 1))
        for e in  sorted(la.eig(hmat, smat)[0].real)[0:10]:
            print e
        
        """
        basis_set.add_basis(0, (0.1, 0.2, 0.3), 1.1-0.3j)
        basis_set.add_basis(1, (0.1, 0.2, 0.3), 1.1-0.3j)
        smat = basis_set.s_mat()
        tmat = basis_set.t_mat()
        vmat = basis_set.v_mat(1.1, (0.1, 0.2, 0.3))
        dip = basis_set.xyz_mat((0, 0, 1))
        print smat
        print vmat 
        print tmat
        print dip
        """
"""

class TestCartGTO(unitttest.TestCase):
    def setUp(self):
        pass

    def test_at(self):
        x0 = 0.1; y0 = 0.2; z0 = 0.3;
        n = 1; m = 2; l = 3;
        zeta = 0.5-0.1j
        g1 = CartGTO(1.2,
                     (n, m, l),
                     (x0, y0, z0),
                     zeta)
        x = 0.2; y = 0.33; z = 0.44;
        dx = x-x0; dy = y-y0; dz = z-z0; 
        self.assertAlmostEqual(1.2 * dx**n * dy**m * dz**l *
                               np.exp(-zeta*dx*dx*dy*dy*dz*dz),
                               g1.at(x, y, z))

    def test_S(self):
        c1 = 1.1
        c2 = 1.2
        g1 = cartGTO(c1, (2, 1, 1), (0.0, 0.0, 0.0),      0.2+0.3j)
        g2 = cartGTO(c2, (2, 1, 2), (0.0, 0.2-0.1j, 0.3), 0.2-0.7j)
        s = cip(g1, g2)
        self.assertAlmostEqual(c1*c2*(-14.3092-3.67256j))

    def test_T(self):
        c1 = 1.1
        c2 = 1.2
        g1 = cartGTO(c1, (0, 0, 0), (0.1, 0.2,    0.3), 1.1+0.2j)
        g2 = cartGTO(c2, (0, 0, 0), (0.1, -0.22, -0.3), 1.3-0.5j)
        s = cip(g1, KE(), g2)
        self.assertAlmostEqual(c1*c2*(-14.3092-3.67256j))

    def test_V(self):
        c1 = 0.33
        c2 = 0.11
        q  = 1.1
        g1 = cartGTO(c1, (2, 1, 0), (0.1, 0.2-0.03j, 0.3), 1.1+0.2j)
        g2 = cartGTO(c2, (0, 1, 3), (0.1, -0.22, -0.3-0.7j), 1.3+0.02j)
        se = -0.0111712963403-0.0039848461450j
        s  = cip(g1, na(q, (-0.1, -0.1+0.1j, 0.3)), g2)

class TestSphericalGTO(unitttest.TestCase):
    def setUp(self):
        pass

    def test_overlap(self):
        gtos = SphericalGTOSet()
        gtos.add_one_basis(1, 0.1, 0.2, 0.3, 1.1)




class TestSphericalGTO(unittest.TestCase):
    def setUp(self):
        pass


    def test_innert_product(self):


        self.assertAlmostEqual(0.88910293518825-0.5004029971544j,
                               overlap(1.1+0.2j,
                                       0, 0, 0,
		                       0.1,  0.2,   0.3,
		                       1.3+0.5j,
                                       0, 0, 0,
		                       0.1, -0.22, -0.3))

        self.assertAlmostEqual(1.422908003563741695-0.476445156336799588j,
                               kinetic(1.1+0.2j,
                                       0, 0, 0,
		                       0.1,  0.2,   0.3,
		                       1.3+0.5j,
                                       0, 0, 0,
		                       0.1, -0.22, -0.3))

        self.assertAlmostEqual(-0.0111712963403-0.0039848461450j,
                               nuclear_attraction(0.1,  C(0.2, 0.03), 0.3,
                                                  1.0, 
			                          2, 1, 0, C(1.1,0.2),
			                          0.1, -0.22, C(-0.3,-0.7),
                                                  1.0,
			                          0, 1, 3, C(1.3,0.02),
			                          -0.1, C(-0.1,0.1), 0.3))


        
    def test_array3(self):
        xs = tuple_c3(1.1, 1.2, 1.3)
        self.assertAlmostEqual(1.1, xs.x)
        self.assertAlmostEqual(1.2, xs.y)
        self.assertAlmostEqual(1.3, xs.z)

    def test_create(self):
        g1 = SphericalGTO(2, 1, tuple_c3(0.1, 0.2, 0.3), 0.2-0.1j)
        self.assertEqual(2, g1.l)
        self.assertEqual(1, g1.m)
        self.assertAlmostEqual(0.2-0.1j, g1.zeta)
        self.assertAlmostEqual(0.1, g1.xyz.x)

"""
if __name__ == '__main__':
    unittest.main()
