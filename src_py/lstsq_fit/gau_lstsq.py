import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
from sympy.mpmath import polyroots


def lstsq(m, n, w, pn, f, 
          maxsteps=50, extraprec=10):
    """compute Gauss functions to fit function f
    
    f(x) = sum_j c_j x^pn exp(a_jx^2)

    Inputs
    ------
    m : Int
        Number of GTOs to expand function f
    n : Int 
        Number of points to use fitting. 
	require :  n>2m
    w : Scalar
        Paramaeter. Points used for fitting are
        x_k = Sqrt(kw) k = 1,...,n
    pn : Int
        Principle number of GTOs
    f  : lambda 
       Target function

    Returns
    -------
    (cs, zs, {"pn": pn, "f": f})
    cs : [Scalar]
        List of coefficient of GTOs
    zs : [Scalar]
        List of orbital exponents of GTOs
    """

    # grid points
    xn_s  = [ np.sqrt((k+1)*w) for k in range(n)]

    # discreted function values
    fn_s  = [ f(x)/(x**pn) for x in xn_s]

    # equation for determine S
    f_mat = [ fn_s[j:j+m] for j in range(n-m)]
    f_vec = [ -fn_s[j+m] for j in range(n-m)]
    (sm_s, resi, d1, d2) = np.linalg.lstsq(f_mat, f_vec)

    # solve polynomial of S
    coef_for_v = np.hstack([[1], sm_s[::-1]])
    vm_s = polyroots(coef_for_v, maxsteps, extraprec)

    # calculate orbital exponents
    am_s = np.array([np.log(np.real_if_close(complex(v)))/w for v in vm_s])

    # calculate coefficient
    v_mat = np.array([[ v**(k+1) for v in vm_s] for k in range(n)])
    (cm_s, resi_c, d1, d2) = np.linalg.lstsq(v_mat, fn_s)
    return (cm_s, am_s, {"pn": pn, "f": f, "xn_s":xn_s})


def lstsq_from_region(m, xmin, xmax, pn, f):
    w = xmin*xmin
    num = int(xmax*xmax/w)
    return lstsq(m, num, w, pn, f)

def to_func(cs_as_n):
    """ from the results of gau_least_square function, compute lambda.
    """
    (cm_s, am_s, other) = cs_as_n
    pn = other["pn"]
    def func(x):
        return sum([c*x**pn*np.exp(a*x*x) for (c, a) in zip(cm_s, am_s)])
    return func

def abs_err(m, n, pn, f, x1):
    """ compute distance of GTO linear combination and original function in L2 norm.
    this function is used for scipy.optimize.minimize function.

    input exp_ws is assumed [log(w)] 
    """
    def __func__(log_ws):
        w = np.exp(log_ws[0])
        res = lstsq(m, n, w, pn, f)
        f0 = res[-1]["f"]
        f1 = gau_func(res)
        def diff(x):
            return abs(f0(x) - f1(x))**2
        (val, err) = quad(diff, 0.0, x1)
        return val
    return __func__

def lstsq_opt(m, n, pn, f, w0, x1):
    err = abs_err(m, n, pn, f, x1)
    minimize_res =  minimize(err, np.log(w0), method='Nelder-Mead')
    w_opt = np.exp(minimize_res['x'][0])
    (cs, zs, other) = gau_least_square(m, n, w_opt, pn, f)
    other['minimize'] = minimize_res
    return (cs, zs, other)

