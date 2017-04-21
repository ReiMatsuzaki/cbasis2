"""
derivative basis for STO/GTO
"""
from linspace import *
from set_l2func import *

def nprime_over_n(basis_type, order, n, z):
    """ compute
    """
    if(basis_type == STO and order==1):
        return (n+0.5)/z
    if(basis_type == GTO and order==1):
        return (n * 0.5 + 0.25)/z
    if(basis_type == STO and order==2):
        return (4 * n * n - 1) / (4.0 * z * z);
    if(basis_type == GTO and order==2):
        return (-3.0/4.0+n/2.0) * (1.0/4.0 +n/2.0) / (z*z)
    msg = "basis_type:{0}, order:{1}".format(str(basis_type), order)
    raise Exception(msg)

def get_m(basis_type):
    if(basis_type == STO):
        return 1
    if(basis_type == GTO):
        return 2
    

def d0_basis(bt, n, z):
    n_term = 1.0/cnorm(bt(1.0, n, z))
    return bt(n_term, n, z)

def d1_basis(bt, n, z):
    m = get_m(bt)
    n_term = 1.0/cnorm(bt(1.0, n, z))
    cp = nprime_over_n(bt, 1, n, z)
    # print "aaaAA: ", n, m, n_term, cp
    return bt(n_term*cp, n, z) + bt(-n_term, n+m, z)

def d2_basis(bt, n, z):
    m = get_m(bt)
    n_term = 1.0/cnorm(bt(1.0, n, z))
    cp = nprime_over_n(bt, 1, n, z)
    cpp= nprime_over_n(bt, 2, n, z)
    return bt(n_term*cpp, n, z) + bt(-2*n_term*cp, n+m, z) + bt(n_term, n+2*m, z)

        
