from linspace import *
from set_l2func import *

import numpy as np

class HAtom:
    def __init__(self, z):
        self.z = z

    def eigenfunc(self, n, l):
        if(n == 1 and l == 0):
            return STO(2.0, 1, 1.0)
        if(n == 2 and l == 0):
            return STO(1.0/np.sqrt(2.0), 1, 0.5) + STO(-1.0/(2.0*np.sqrt(2.0)), 2, 0.5)
        if(n == 2 and l == 1):
            return STO(1.0/(np.sqrt(6.0)*2.0), 2, 0.5)
        if(n == 3 and l == 2):
            return STO(4.0/(81.0*np.sqrt(30.0)), 3, 1.0/3.0)
        else:
            raise Exception("not implemented yet: n={0}, l={1}".format(n, l))

    def eigenenergy(self, n):
        return -1.0/(2.0*n*n)

    def length(self, n, l0, l1):
        if(abs(l0-l1) != 1):
            return FuncZero()
        f0 = self.eigenfunc(n, l0)
        return op_apply(rm(1), f0)

    def velocity(self, n, l0, l1):
        f0 = self.eigenfunc(n, l0)
        if(l0 + 1 == l1):
            o = dr(1) + (-1.0-l0)*rm(-1)
            return op_apply(o, f0)
        if(l0 - 1 == l1):
            o = dr(1) + (l0)*rm(-1)
            return op_apply(o, f0)
        else:
            return FuncZero()

    def hop(self, l):
        h0 = -0.5*D2() + (-1.0)*Rm(-1)
        if(l == 0):
            return h0
        else:
            return h0 + l*(l+1)*0.5*Rm(-2)

    def h_minus_ene_op(self, l, ene):
        return self.hop(l) + (-ene) * Rm(0)
