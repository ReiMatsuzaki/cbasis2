"""
wkb solution for coulomb function
2015/12/7
"""

import cmath
import numpy as np
from scipy.special import gamma


def coulomb_phase_shift(k, L):
    return cmath.phase(gamma(1.0 + L - 1.0j/k))

def wkb_local_momentum(k, L, r):
    c = 1.0 * L * (L + 1)
    return np.sqrt(k**2 + 2.0/r - c/(r**2))

def wkb_zeta(k, L, r):
    p = wkb_local_momentum(k, L, r)
    c = 1.0*L*(L+1)
    term1 = 1.0/(4.0*p**3)*(3.0*c/(r**4)-2.0/(r**3))
    term2 = 5.0/(8.0*p**5)*(1.0/(r*r)-c/(r**3))**2
    return p + term1 + term2

def wkb_zeta_prime(k, L, r):
    L = 1.0 * L
    return -(3.0*L*(L + 1)/r**5 - 3.0/(2*r**4))/(-L*(L + 1)/r**2 + k**2 + 2/r)**(3.0/2.0) + 5.0*(-12*L*(L + 1)/r**4 + 8/r**3)*(2*L*(L + 1)/r**3 - 2/r**2)/(32*(-L*(L + 1)/r**2 + k**2 + 2/r)**(5.0/2.0)) - (-3*L*(L + 1)/(4*r**4) + 1/(2*r**3))*(-3*L*(L + 1)/r**3 + 3/r**2)/(-L*(L + 1)/r**2 + k**2 + 2/r)**(5.0/2.0) + 5*(-5*L*(L + 1)/r**3 + 5/r**2)*(2*L*(L + 1)/r**3 - 2/r**2)**2/(32*(-L*(L + 1)/r**2 + k**2 + 2/r)**(7.0/2.0)) + (L*(L + 1)/r**3 - 1/r**2)/np.sqrt(-L*(L + 1)/r**2 + k**2 + 2/r)

def wkb_phi(k, L, r):
    p = wkb_local_momentum(k, L, r)
    if L == 0:
        phi0 = p*r + 1.0/k*np.log(2+2*k*k*r+2*p*k*r)
        term1 = (9+6*k**2*r+2*k**4*r**2)/(24*p**3*r**2)
        return phi0 + term1
    else:
        c = 1.0 * L * (L + 1)
        phi0 = p*r - np.sqrt(c)*np.arctan((r-c)/(np.sqrt(c)*r*p)) + 1.0/k*np.log(2+2*k*k*r+2*p*k*r)
        term1 = - (c + 2*c**2*k**2 + 3*c*k**4*r**2 + 3*r + 6*k**2*r**2 + k**4*r**3) / (24.0*r**3*p**3*(1+c*k*k))
        term2 = - 1.0/(8*np.sqrt(c)) * np.arctan((r-c)/(np.sqrt(c)*p*r))
        return phi0 + term1 + term2

def wkb_delta_and_A(k, L, r, psi, d_psi):
    zeta = wkb_zeta(k, L, r)
    zeta_prime = wkb_zeta_prime(k, L, r)
    arg_arctan = zeta*psi/(d_psi + psi*zeta_prime/(2.0*zeta))
    delta = np.arctan(arg_arctan)
    A = np.sqrt(zeta) * psi / np.sin(delta)
    return (delta, A)

def wkb_func(k, L, r_match, psi, d_psi):
    (delta, A) = wkb_delta_and_A(k, L, r_match, psi, d_psi)
    phi0 = wkb_phi(k, L, r_match)
    def __func__(r):
        phi = wkb_phi(k, L, r) - phi0
        zeta = wkb_zeta(k, L, r)
        return A * zeta**(-0.5) * np.sin(phi+delta)
    return __func__

def wkb_phase_shift(k, L, r, psi, d_psi, r_asym):
    (delta, A) = wkb_delta_and_A(k, L, r, psi, d_psi)
    phi = wkb_phi(k, L, r_asym) - wkb_phi(k, L, r)
    phase_shift = phi + delta - k*r_asym - (1.0/k)*np.log(2.0*k*r_asym) + L*np.pi/2.0

    if A > 0.0:
        return phase_shift % (2*np.pi)
    else:
        return (phase_shift + np.pi) % (2*np.pi)
    
def phase_shift_by_coulomb_func(psi_r, r, r_match, k):
    """
    psi_r: caluculated wave function written in sympy symbol r
    r:     sympy symbol used in psi_r
    r_match: matching point 
    k:       wave number
    """
    
