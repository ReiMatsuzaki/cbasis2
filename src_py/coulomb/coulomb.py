"""
compute coulomb function.

eta = zA zB / k

2015/12/3 R.Matsuzaki
"""

import coulomb_bind
import sympy
import sympy.mpmath as symmath
import cmath
import numpy as np

def coulomb(eta, x, l_f, k=0):
    """ compute regular and irregular coulomb function and these derivative with error.
    F_L(eta, x), F_L'(eta, x), G_L(eta, x), G_L'(eta, x)
    
    Input
    _____
    eta : scalar
          first argument of Coulomb functions(zAzB/k)
    x : scalar
          second argument of Coulomb functions(rho=kr)
    l_f : scalar
          angular quantum number.
    k : scalar
          normaly 0.
    Returns
    -------
    f  : F_L(eta, x)
    g  : G_L(eta, x)
    fp : F_L'(eta, x)
    gp : G_L'(eta, x)
    err:
    """
    if( x < 0.00000001):
        raise Exception("x must be positive real number")
    func = coulomb_bind.coulomb(eta, x, l_f, k)
    return (func.f, func.g, func.fp, func.gp, func.err)

def asym_arg(eta, x, l_f):
    sigma_l = cmath.phase(complex(symmath.gamma(l_f + 1 + 1.0j*eta)))
    return x - eta * np.log(2.0*x) -0.5*l_f*np.pi + sigma_l

def outgoing(eta, x, l_f, k=0):
    """ compute irregular coulomb function with outgoing boundary condition.
    
    Input
    _____
    eta : scalar
          first argument of Coulomb functions (zAzB/k)
    x : scalar
          second argument of Coulom functions (rho=kr)
    l_f : scalar
          angular quantum number
    k : scalar
          normaly 0
    Returns
    -------
    out: complex(f, g)
    """
    (f, g, fp, gp, err) = coulomb(eta, x, l_f, k)
    return g + 1.0j * f

def coulomb_phase(L, eta):
    """
    Gives the Coulomb phase shift.
    -------
    \sigma_L = arg Gamma(L+1+i\eta)
    \eta = z1 z2 / k
    -------
    """
    return cmath.phase(sympy.gamma(1+L+1.0j*eta))

