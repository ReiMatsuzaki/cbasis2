from sympy import *

r = Symbol('r')
s = 2.2 * r**3 * exp(-1.1*r)
g = 1.3 * r**2 * exp(-1.2*r*r)
print "s2s_ref=", integrate(s*r*r*s, (r,0,oo)).evalf()
print "g2g_ref=", integrate(g*r*r*g, (r,0,oo)).evalf()
print "s2g_ref=", integrate(s*r*r*g, (r,0,oo)).evalf()
print "sDs_ref=", integrate(s*diff(s,r,2), (r,0,oo)).evalf()
print "gDg_ref=", integrate(g*diff(g,r,2), (r,0,oo)).evalf()
print "sDg_ref=", integrate(s*diff(g,r,2), (r,0,oo)).evalf()
print "------------------------"

"""output
s2s_ref= 161.644807242673
g2g_ref= 0.131127436620057
s2g_ref= 0.663645309086432
sDs_ref= -3.38091660405710
gDg_ref= -0.352470549634713
sDg_ref= 0.208872645967760
"""

s = 2.0 * r**3 * exp(-(1.2-0.3j)*r)
g = 1.3 * r**2 * exp(-(1.1-0.1j)*r*r)
print "s2s_ref=", integrate(s*r*r*s, (r,0,oo)).evalf()
print "g2g_ref=", integrate(g*r*r*g, (r,0,oo)).evalf()
print "s2g_ref=", integrate(s*r*r*g, (r,0,oo)).evalf()
print "sDs_ref=", integrate(s*diff(s,r,2), (r,0,oo)).evalf()
print "gDg_ref=", integrate(g*diff(g,r,2), (r,0,oo)).evalf()
print "sDg_ref=", integrate(s*diff(g,r,2), (r,0,oo)).evalf()
print "------------------------"

a = Symbol('a')
b = Symbol('b')
prod = lambda f: simplify(integrate(f, (r,0,oo), conds='none'))
print 0, prod(exp(-a*r-b*r*r))
print 1, prod(r*exp(-a*r-b*r*r))
print 2, prod(r**2*exp(-a*r-b*r*r))
print 3, prod(r**3*exp(-a*r-b*r*r))

""" output
s2s_ref= -27.5296456055511 + 37.4411871137165*I
g2g_ref= 0.16651663387627 + 0.0546850960763247*I
s2g_ref= 56.9050748227212 - 9.17658487944802*I
sDs_ref= -0.526917955926652 - 1.46206690245985*I
gDg_ref= -0.395454606004856 - 0.0541117842324456*I

0 sqrt(pi)*(-erf(a/(2*sqrt(b))) + 1)*exp(a**2/(4*b))/(2*sqrt(b))
1 (pi*a*exp(a**2/(4*b))*erf(a/(2*sqrt(b))) - pi*a*exp(a**2/(4*b)) + 2*sqrt(pi)*sqrt(b))/(4*sqrt(pi)*b**(3/2))
2 (pi*(a**2 + 2*b)*exp(a**2/(4*b)) - sqrt(pi)*(2*a*sqrt(b) + sqrt(pi)*(a**2 + 2*b)*exp(a**2/(4*b))*erf(a/(2*sqrt(b)))))/(8*sqrt(pi)*b**(5/2))
3 (-pi*a*(a**2 + 6*b)*exp(a**2/(4*b)) + sqrt(pi)*(-4*b**(3/2) + (a**2 + 6*b)*(sqrt(pi)*a*exp(a**2/(4*b))*erf(a/(2*sqrt(b))) + 2*sqrt(b))))/(16*sqrt(pi)*b**(7/2))
"""
