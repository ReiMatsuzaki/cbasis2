import unittest
from nnewton import *
import numpy as np

def func(k, xs):
    x = xs[0]
    a = 1.5
    return np.cos(a*(x-k))

def func_grad_hess(k, xs):
    x = xs[0]
    a = 1.5
    return (np.cos(a*(x-k)),
            np.array([-a*np.sin(a*(x-k))]),
            np.array([[-a*a*np.cos(a*(x-k))]]))

def func_keyarg(**keyargs):
    print keyargs
    print 'h' in keyargs

class TestNNewton(unittest.TestCase):
    def setUp(self):
        pass

    def test_pd(self):
        x = [1.1]
        (val, grad, hess) = func_grad_hess(0.3, x)
        h = 0.001
        target = lambda xs:func(0.3, xs)
        self.assertAlmostEqual(grad[0], num_pd(target, x, h, 0, "c1"))

    def test_pd12(self):
        x = [1.1]
        (val, grad0, hess0) = func_grad_hess(0.3, x)
        h = 0.001
        target = lambda xs:func(0.3, xs)
        (grad1, hess1) = num_grad_hess(target, x, h, "r1")
        self.assertAlmostEqual(grad0[0], grad1[0], places=5)
        self.assertAlmostEqual(hess0[0,0], hess1[0,0], places=5)
        (grad1, hess1) = num_grad_hess(target, x, h, "c1")
        self.assertAlmostEqual(grad0[0],   grad1[0]  , places=10)
        self.assertAlmostEqual(hess0[0,0], hess1[0,0], places=9)

    def test_pd_scale(self):
        x = [1.1]
        (val, grad, hess) = func_grad_hess(0.3, x)
        h = 0.001
        target = lambda xs:func(0.3, xs)
        calc = num_pd(target, x, lambda az:0.001*az, 0, "c1")
        self.assertAlmostEqual(grad[0], calc)

    def test_nnewton(self):
        k0 = 0.2
        target = lambda xs: func(k0, xs)
        (convq, xs_o, val, grad) = nnewton(target, [0.11])
        self.assertTrue(convq)
        self.assertAlmostEqual(k0, xs_o[0])

    def test_nnewton_analytic(self):
        k0 = 0.2
        target = lambda xs:func_grad_hess(k0, xs)

        self.assertTrue(is_reasonable_val_grad_hess(target,
                                                    [0.1],
                                                    method="c1",
                                                    h=0.001))

        (convq, xs_o, val, grad) = nnewton(target, [0.1], func="vgh", show_lvl=1)
        self.assertTrue(convq)
        self.assertAlmostEqual(k0, xs_o[0])

    def test_nnewton_params(self):
        ks = [-0.1, 0.0, 0.1]
        res = nnewton_params(func, ks, (0.0, [0.03]))
        (k, xs_o, val) = res[0]
        self.assertAlmostEqual(k, -0.1)
        self.assertAlmostEqual(k, xs_o[0])
        (k, xs_o, val) = res[1]
        self.assertAlmostEqual(k, 0.0)
        self.assertAlmostEqual(k, xs_o[0])
        (k, xs_o, val) = res[2]
        self.assertAlmostEqual(k, 0.1)
        self.assertAlmostEqual(k, xs_o[0])

    def test_keyargs(self):
        func_keyarg(h=1.0, x=1.1)

if __name__ == '__main__':
    unittest.main()
