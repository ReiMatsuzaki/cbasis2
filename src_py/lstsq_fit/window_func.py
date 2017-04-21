""" define switching function for window function.
2015/12/5 R.Matsuzaki
"""

import numpy as np

def exp_switch(f, g, a, x0):
    """gives switching function which behave f if x<<x0 
    and behave g if x>>x0. a is behavior change parameter

    Define function S:
    S(x) = exp(a(x-x0))

    exp_switch function is defined as follow
    (f + gS)/(1+S)

    Inputs
    ------
    f,g : lambda (Real->Scalar)
    a   : Positive Real number
    x0  : Real

    Returns
    -----
    out : lambda (Real->Scalar)
    """

    def __func__(x):
        exp_term = np.exp(a*(x-x0))
        return (f(x) + g(x)*exp_term) / (1.0 + exp_term)
    return __func__

def exp_switch_zero(f, g, a, x0):
    """gives switching function which behave f if x<<x0 
    and behave g if x>>x0. a is behavior change parameter

    Define function S:
    S(x) = exp(a(x-x0))

    exp_switch function is defined as follow
     f*(1 + gS)
    -----------
        1+fS

    Inputs
    ------
    f,g : lambda (Real->Scalar)
    a   : Positive Real number
    x0  : Real

    Returns
    -----
    out : lambda (Real->Scalar)
    """
    def __func__(x):
        exp_term = np.exp(a*(x-x0))
        return f(x) * (1.0 + g(x) * exp_term) / (1.0 + f(x)*exp_term)
    return __func__

def double_switch(f, a, x0, g, b, x1):
    """gives swhitch function which behave 
    f if x << x0
    1 if x0 < x < x1
    g if x >> x1
    """
    one = lambda x:1.0
    left_term = exp_switch_zero(f, one, a, x0)
    right_term = exp_switch(one, g, b, x1)
    res = lambda x: left_term(x) * right_term(x)
    return res
