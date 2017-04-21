"""
<exp[ikr] | cart GTO(W)> = exp[ikw] <exp[ikr] | GTO(0)> 
.                        = exp[ikw] <r^nx exp[ixk_x-ax^2]>
.                                   <r^ny exp[iyk_y-ay^2]>
.                                   <r^ny exp[izk_z-az^2]>
-axx+ikx = -a{xx - ik/a x + (ik/2a)^2 - (ik/2a)^2}
.        = -a(x-ik/2a)^2 + a(ik/2a)^2
<r^n exp[ikx-axx]> = exp[-kk/4a] <(x+ik/2a)^n exp[-ax^2]>
.                  = exp[-kk/4a] sum_m nCm (ik/2a)^(n-m) <x^m exp[-ax^2]>
"""

from sympy import *

x = Symbol('x')
k = Symbol('k')
a = Symbol('a')
n = Symbol('n')

inte = simplify(integrate(x**n * exp(I*k*x -a*x*x), (x,0,oo), conds='none'))
print inte
print N(inte.subs({k:1.1, a:1.2, n:3}))
inte = simplify(integrate(x**n * exp(-a*x*x), (x,0,oo), conds='none'))
print inte
print N(inte.subs({a:1.2, n:3}))

"""
a**(-n/2 - 1/2)*meijerg(((-n/2 + 1/2,), ()), ((0, 1/2), ()), exp_polar(-I*pi)*polar_lift(k)**2/(4*a))/(2*sqrt(pi))
0.0744991235617568 + 0.299682718346912*I
a**(-n/2 - 1/2)*gamma(n/2 + 1/2)/2
0.347222222222222
"""

