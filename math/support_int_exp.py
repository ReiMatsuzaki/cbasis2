"""
from numpy import exp
from scipy.integrate import quad
print quad(lambda x: x**3 * exp(-1.1*x),   0, 15.0)
print quad(lambda x: x**3 * exp(-1.1*x*x), 0, 15.0)
print quad(lambda x: x**2 * exp(-1.1*x*x), -15.0, 15.0)
print quad(lambda x: x**3 * exp(-1.1j*x+1.3*x*x), 0, 15.0)
print quad(lambda x: x**3 * exp(-1.1j*x+1.3*x*x), -15.0, 15.0)
"""

from sympy import *
x = Symbol('x')
n = Symbol('n')
z = Symbol('z')
k = Symbol('k')
print integrate(x**3 * exp(-1.1*x), (x,0,oo))
print integrate(x**3 * exp(-1.1*x*x), (x,0,oo))
print N(integrate(x**2 * exp(-1.1*x*x), (x,-oo,oo)))
inte = simplify(integrate(x**n * exp(-k*x - z*x*x), (x,0,oo), conds='none'))
print N(inte.subs({n:3, k:1.1, z:1.3-0.2j}))
inte =  integrate(x**3 * exp( -1.1j*x - 1.2*x*x), (x,-oo,oo))
print N(inte)

"""
4.09808073219042
0.413223140495868
0.768167471819406
0.0873211906359305 + 0.0197914200245872*I
0.e-21 - 0.599365436693823*I
"""

