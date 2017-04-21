from sympy import *

r = Symbol("r")
z = Symbol("z")
n = Symbol("n")
N = Function('N')(z)
print "==== STO ===="
nfunc = N * (r**n) * exp(-z*r)
print "STO'/STO : ", simplify(diff(nfunc, z)/nfunc)
print "STO''/STO : ", simplify(diff(nfunc, z, 2)/nfunc)
func = (r**n)*exp(-z*r)
nterm = 1/sqrt(integrate(func*func, (r,0,oo), conds='none'))
print "N'/N = ", simplify(diff(nterm, z)/nterm)
print "N''/N = ", simplify(diff(nterm, z, 2)/nterm)
print ""
print "============="
print ""
print "==== GTO ===="

nfunc = N * (r**n) * exp(-z*r*r)
print "GTO'/GTO : ", simplify(diff(nfunc, z)/nfunc)
print "GTO''/GTO : ", simplify(diff(nfunc, z, 2)/nfunc)

func = (r**n)*exp(-z*r*r)
nterm = 1/sqrt(integrate(func*func, (r,0,oo), conds='none'))
#print n, "GTO:", nterm, func
print "N'/N = ", simplify(diff(nterm, z)/nterm)
print "N''/N = ", simplify(diff(nterm, z, 2)/nterm)

for n0 in [1,2,3,4]:
    func = (r**n0)*exp(-z*r*r)
    nterm = 1/sqrt(integrate(func*func, (r,0,oo), conds='none'))
    print n0, "GTO:", nterm, func

""" output
==== STO ====
STO'/STO :  -r + Derivative(N(z), z)/N(z)
STO''/STO :  (r**2*N(z) - 2*r*Derivative(N(z), z) + Derivative(N(z), z, z))/N(z)
N'/N =  (n + 1/2)/z
N''/N =  (n**2 - 1/4)/z**2

=============

==== GTO ====
GTO'/GTO :  -r**2 + Derivative(N(z), z)/N(z)
GTO''/GTO :  (r**4*N(z) - 2*r**2*Derivative(N(z), z) + Derivative(N(z), z, z))/N(z)
N'/N =  (2*n + 1)/(4*z)
N''/N =  (4*n**2 - 4*n - 3)/(16*z**2)
1 GTO: 2*2**(3/4)/(pi**(1/4)*sqrt(z**(-3/2))) r*exp(-r**2*z)
2 GTO: 4*2**(3/4)*sqrt(3)/(3*pi**(1/4)*sqrt(z**(-5/2))) r**2*exp(-r**2*z)
3 GTO: 8*sqrt(15)*2**(3/4)/(15*pi**(1/4)*sqrt(z**(-7/2))) r**3*exp(-r**2*z)
4 GTO: 16*sqrt(105)*2**(3/4)/(105*pi**(1/4)*sqrt(z**(-9/2))) r**4*exp(-r**2*z)
"""    


    
