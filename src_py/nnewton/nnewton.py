import numpy as np
import scipy.linalg as la
from utils import *


# ==== Numerical partial derivative ====
# ---- impl ----
def hs_for_num_pd(h, num, k):
    """ compute finite difference vector.
    
    Inputs
    ------
    h : real, finite difference
    num : int, dimendion
    k   : int, index 

    Returns
    -------
    out : 
    """
    def one(i):
	if i==k:
	    return 1
	else:
	    return 0
    hs = np.array([ one(i)*h for i in range(num)])
    return hs

def num_pd_r1(g, x, h, k):
    """partial derivative with first order real"""
    num = len(x)
    hh = h(abs(x[k])) if callable(h) else h
    hs = hs_for_num_pd(hh, num, k)
    return (g(x+hs)-g(x-hs))/(2*hh)
    
def num_pd_c1(g, x, h, k):    
    """ partial derivative with one first order complex"""
    num = len(x)
    hh = h(abs(x[k])) if callable(h) else h
    hs = hs_for_num_pd(hh, num, k)
    ig = lambda y: 1.0j * g(y)
    ih = 1.0j*hs
    return (g(x+hs)-g(x-hs)-ig(x+ih)+ig(x-ih))/(4*hh)

def num_pd_r3(g, x, h, k):
    """ partial derivative with one third order real"""
    num = len(x)
    hh = h(abs(x[k])) if callable(h) else h
    hs = hs_for_num_pd(hh, num, k)
    m1 = (g(x+hs)   - g(x-hs)  )/(2*hh)
    m2 = (g(x+2*hs) - g(x-2*hs))/(4*hh)
    m3 = (g(x+3*hs) - g(x-3*hs))/(6*hh)
    return 3.0/2.0*m1-3.0/5.0*m2+1.0/10.0*m3


# ---- numerical partial derivative ----
def num_pd(g, x, h, k, method=0):
    """numerical partial derivative
    
    Inputs 
    ------
    g : [Scalar]->Scalar
    x : [Scalar]
    h : Real || Real->Real (compute finite difference from absolute value of x)
    k : int 
    method : "r1" || "c1" || "r3"

    Outpus
    ------
    result : Scalar
    """
    if(method == 0 or method == "r1"):
        return num_pd_r1(g, x, h, k)
    if(method == "c1"):
	return num_pd_c1(g, x, h, k)
    if(method == "r3"):
	return num_pd_r3(g, x, h, k)

def num_pd2(g, x, h, k1, k2, method=0):
    """numerical second partial derivative
    
    Inputs 
    ------
    g : [Scalar]->Scalar
    x : [Scalar]
    h: Real || 
       Real->Real (compute finite difference from absolute value of x)
    
    k1,k2 : int 
    method : "r1" || "c1" || "r3"

    Outpus
    ------
    result : Scalar
    """
    dg1 = lambda y: num_pd(g, y, h, k1, method)
    return num_pd(dg1, x, h, k2, method)


# ---- grad/hess ----
def num_grad_hess_r1(g, xs, h):
    num = len(xs)
    if(num != 1):
        raise Exception("len(xs) != 1 is not supported")
    x = xs[0]
    gp = g([x+h])
    gm = g([x-h])
    g0 = g([x])
    grad = np.array([(gp-gm)/(2*h)])
    hess = np.array([[(gp+gm-2*g0)/(h*h)]])
    return (grad, hess)

def num_grad_hess_c1(g, xs, h):
    num = len(xs)
    if(num != 1):
        raise Exception("len(xs)!=1 is not supported")
    x = xs[0]
    ii = 1.0j
    ih= ii*h
    gp = g([x+h])
    gm = g([x-h])
    gip= g([x+ih])
    gim= g([x-ih])
    grad = (gp - gm - ii*gip + ii*gim) / (4.0*h)
    hess = (gp + gm - gip - gim) / (2.0*h*h)
    return (np.array([grad]), np.array([[hess]]))

def num_grad_hess(g, x, h, method=0):
    if(method == "r1" or method==0):
        return num_grad_hess_r1(g, x, h)
    if(method == "c1"):
        return num_grad_hess_c1(g, x, h)
    else:
        raise Exception("mehtod without c1 and r1 is not implemented")

def ngrad(g, x, h, method=0):
    if(h == None):
        raise(Exception('h is None'))    
    return np.array([num_pd(g, x, h, k, method)
                     for k in range(len(x))])

def nhess(g, x, h, method=0):
    return np.array([[num_pd2(g, x, h, k1, k2, method)
                      for k1 in range(len(x))]
                     for k2 in range(len(x))])

# ==== Check function ====
def is_reasonable_val_grad_hess(f, xs0, eps = 10.0**(+8.0),
                                method="r1", h=0.001):
    (val, grad, hess) = f(xs0)
    target = lambda xs: f(xs)[0]
    num = len(xs0)
    ngrad = np.array([num_pd(target, xs0, h, k, method=method)
                      for k in range(num)])
    nhess = np.array([[num_pd2(target, xs0, h, k1, k2, method=method)
             for k1 in range(num)] for k2 in range(num)])
    return (abs(grad-ngrad)/num < eps and
            abs(flatten((hess-nhess)))/(num*num) < eps)
    

# ==== Numerical Newton-Raphson ====
def nnewton(f, x0, iter_max=100, eps=0.0000001,
            h=0.001, show_lvl=0, method=0, func="v"):
    """Newton method with numerical derivative.

    Inputs
    ------
    f : lambda [scalar]->scalar, target function
    x0: [scalar], initial guess
    iter_max: int, max of iteration
    eps: real, convergence torrelance
    show_lvl: int, 0=>no print, 1=>print value, 2=>print grad, 3=>print hess
    method: int, 0=>first order complex, 1=>third order real
    func: character: "v"=> f return value, "vfh" => f return (value, grad, hess)

    Outputs
    -------
    result: bool, convergence or not
    xs: [scalar], results
    f(xs): scalar, value 
    grad: [scalar], gradient 
    """
    num = len(x0)
    ns = range(num)

    # gradient/Hessian calculation
    if(func == "vgh"):
        def gh(xs):
            (val, grad, hess) = f(xs)
            return (grad, hess)

    elif(func == "v" and len(x0) == 1 and (method=="r1" or method=="c1")):
        def gh(xs):
            return num_grad_hess(f, xs, h, method)
    elif(func == "v"):
        def gh(xs):
            grad = ngrad(f, xs, h, method)
            hess = nhess(f, xs, h, method)
            return (grad, hess)
    else:
        raise Exception("func<-{vgh, v}")

    # variable in loop
    if(isinstance(x0, list)):
        xs = np.array(x0)
    else:
        xs = x0

    # optimization loop
    for i in range(iter_max):
        (grad, hess) = gh(xs)
	dx = -la.solve(hess, grad)

	# check convergence
	if(show_lvl==1):
	    print i, xs
	if(show_lvl==2):
	    print i, xs, grad
	if(show_lvl==3):
	    print i, xs, grad, hess

	if(max(abs(grad)) < eps):
            print "OPT", xs, f(xs), grad
	    return (True, xs, f(xs), grad)
	xs += dx

    return (False, xs, f(xs), grad)

# ==== Numerical Newton-Raphson with variables ====
def nnewton_params(param_vs_to_val, params, param0_guess0, **args):
    """
    Inputs
    ------
    param_vs_to_val : (Real, [Scalar]) -> Scalar
        Object function
    params: [Real]
        Parameter list
    param0_guess0: (Real, [Scalar])
        parameter and initial guess at the parameter
    **args: nnewton options
    """
    (param0, guess0) = param0_guess0
    params_down = sorted([p for p in params if p < param0],
                         reverse=True)
    params_up   = sorted([p for p in params if p >= param0])
    
    if("show_lvl" in args):
        show_lvl = args["show_lvl"]
    else:
        show_lvl = 0

    calc_res = [] 
    for params_part in [params_up, params_down]:
        guess = guess0
        for p in params_part:
            if(show_lvl > 0):
                print "param: {0}".format(p)
            target = lambda vs: param_vs_to_val(p, vs)
            (convq, vs, val, grad) = nnewton(target, guess, **args)
            if(not convq):
                print "optimization failed. at param = {0}".format(p)
                break
            calc_res.append((p, vs, val))
            guess = list(vs)		  

    calc_res = sorted(calc_res, key=lambda x:x[0])
    return calc_res	
