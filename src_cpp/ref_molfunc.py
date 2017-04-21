from scipy.special import hyp1f1
from numpy import exp
def inc_gamma(m, z):
    return 1.0/(2*m+1) * hyp1f1(m+0.5, m+1.5, -z)

def exp_inc_gamma(m, z):
    return 1.0/(2*m+1) * hyp1f1(1, m+1.5, -z)

z0 = 40.0+40.0j
print exp(-z0) * inc_gamma(0, -z0), exp_inc_gamma(0, z0)
print exp_inc_gamma(1, z0);
print exp_inc_gamma(2, z0);

print ""
z1 = 5.0-4.3j
print exp_inc_gamma(0, z1)
print exp_inc_gamma(1, z1)
print exp_inc_gamma(2, z1)

