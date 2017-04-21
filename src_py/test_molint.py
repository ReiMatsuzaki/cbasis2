import os 
import sys
import unittest
from l2func_molint import *
import numpy as np
import scipy.linalg as la

class Test_sym_eig(unittest.TestCase):
    def setUp(self):
        pass

    def test_cnorm(self):
        x = 1.1-0.2j
        y = 0.2+0.1;
        xs = np.array([x, y])
        self.assertAlmostEqual(np.sqrt(x*x+y*y), cnorm(xs))

    def test_inv_sqrt(self):
        s = np.array([[1.0,       0.2, 0.2+0.1j],
                      [0.2,       1.0, 0.1],
                      [0.2+0.1j, 0.1, 1.0]])
        s2 = matrix_inv_sqrt(s)
        sss = np.dot(s2, np.dot(s2, s))

        for i in range(3):
            for j in range(3):
                if i==j:
                    self.assertAlmostEqual(1.0, sss[i, i])
                else:
                    self.assertAlmostEqual(0.0, sss[i, j])

    def test_sym_eig(self):
        """
        h = np.array([[1.0-0.3j, 0.2,      0.1+0.1j],
                      [0.2,      0.3+0.1j, 0.4],
                      [0.1+0.1j,      0.4,      1.0-0.2j]])
        """
        h = np.array([[1.0+0.1j, 0.2, 0.1-0.1j],
                      [0.2,      0.3, 0.4],
                      [0.1-0.1j, 0.4, 1.0-0.2j]])
        s = np.array([[1.0,       0.2, 0.2+0.00j],
                      [0.2,       1.0, 0.1],
                      [0.2+0.00j, 0.1, 1.0]])


        (eigs, eigvecs) = sym_eig(h, s)
        for i in range(3):
            xs = np.dot(h, eigvecs.T[i])
            ys = eigs[i]*np.dot(s, eigvecs.T[i])
            for (x, y) in zip(xs, ys):
                self.assertAlmostEqual(x, y)

        for i in range(3):
            for j in range(3):
                sij = np.dot(eigvecs.T[i], np.dot(s, eigvecs.T[j]))
                if(i == j):
                    self.assertAlmostEqual(1.0, sij)
                else:
                    self.assertAlmostEqual(0.0, sij)

        (eigs, eigvecs) = c2_eig(h, s)
        for i in range(3):
            xs = np.dot(h, eigvecs.T[i])
            ys = eigs[i]*np.dot(s, eigvecs.T[i])
            for (x, y) in zip(xs, ys):
                self.assertAlmostEqual(x, y)
        for i in range(3):
            for j in range(3):
                sij = np.dot(eigvecs.T[i], np.dot(s, eigvecs.T[j]))
                if(i == j):
                    self.assertAlmostEqual(1.0, sij)
                else:
                    self.assertAlmostEqual(0.0, sij)

        
class TestMolint(unittest.TestCase):
    def setUp(self):
        pass

    def test_first(self):
        gtos = GTOs()
        for n in range(-5, 5):
            zeta = 2.0**n
            gtos.add_gtos(0, (0.0, 0.0, 0.0), zeta)
        gtos.add_atom(1.0, (0.0, 0.0, 0.0))
        res = gtos.calc().to_dict()
        """
        s = res.get("s")
        t = res.get("t")
        v = res.get("v")
        z = res.get("z")
        h = v + t
        """
        print "sym_eig"
        h = res["v"] + res["t"]
        s = res["s"]
        (eigs, eigvecs) = sym_eig(h, s)
        print eigs

    def test_at_r(self):
        gtos = GTOs()
        gtos.add_gtos(0, (0.0, 0.0, 0.0), 0.1+0.01j)
        gtos.add_gtos(0, (0.0, 0.0, 0.0), 0.1+0.02j)

        vs = gtos.at_r_ylm(0, 0, rs=[1.1, 1.2], cs=[1.1, 1.3])

        # -- below values are copied from test_gto3d.cpp --
        #(1.1,0), (1.05578093454658,0.101631683045226)
        #(1.2,0), (1.12598919239445,0.104354649962767)

        self.assertAlmostEqual(1.05578093454658+0.101631683045226j, vs[0])
        self.assertAlmostEqual(1.12598919239445+0.104354649962767j, vs[1])

    def test_zmat_other(self):
        gtos1 = GTOs()
        gtos2 = GTOs()
        gtosf = GTOs()

        for i in range(2):
            xyz = (0.1-0.1j, 1.0, 0.2)
            gtos1.add_gtos(0, xyz, 1.3+i)
            gtosf.add_gtos(0, xyz, 1.3+i)

        for L in range(0, 3):
            xyz = (0.0, 0.0, 0.0)
            gtos2.add_gtos(L, xyz, 1.2)
            gtosf.add_gtos(L, xyz, 1.2)

        res12 = gtos1.calc_zmat_other(gtos2).to_dict()
        resf  = gtosf.calc().to_dict()
        self.assertEqual(gtos1.size_basis(), res12['z'].shape[0])
        self.assertEqual(gtos2.size_basis(), res12['z'].shape[1])
        for i in range(gtos1.size_basis()):
            for j in range(gtos2.size_basis()):
                self.assertAlmostEqual(resf['z'][i][2+j],
                                       res12['z'][i][j])

    def test_show(self):
        gtos = GTOs()
        gtos.add_gtos(0, (1.0-0.1j, 0.0, 0.2), 1.1-0.2j)
        gtos.add_gtos(1, (1.0-0.1j, 0.0, 0.2), 1.1)
        gtos.add_gtos(2, (1.0-0.1j, 0.0, 0.2), 123.0)
        gtos.add_gtos(0, (0.0, 2.3-12.0j, 0.2), 123.0)
        gtos.show()

if __name__ == '__main__':
    unittest.main()
