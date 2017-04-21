from l2func_bind import *
from linspace import cip_dict, cip_op_dict, set_as_func, set_as_op, OpId, op_apply_dict

# ==== Set Func ====
# ---- general ----
map(set_as_func, [GTO, STO, CutSTO])

# ---- repr/str ----
def STO_repr(self):
    return "STO({0},{1},{2})".format(self.c.__repr__(), self.n, self.z.__repr__())

def GTO_repr(self):
    return "GTO({0},{1},{2})".format(self.c.__repr__(), self.n, self.z.__repr__())

def CutSTO_repr(self):
    return "CutSTO({0},{1},{2},{3})".format(self.c.__repr__(), self.n, self.z.__repr__(),
                                            self.r0)

STO.__repr__ = STO_repr
GTO.__repr__ = GTO_repr
CutSTO.__repr__ = CutSTO_repr

# ==== Set Operator ====
# ---- general ----
map(set_as_op, [D1, D2, Rm, Cut])

# ---- repr ----
def d1_repr(self):
    return "D1"

def d2_repr(self):
    return "D2"

def rm_repr(self):
    return "Rm({0})".format(self.m)

D1.__repr__ = d1_repr
D2.__repr__ = d2_repr
Rm.__repr__ = rm_repr

# ==== inner product ====
cip_dict[(STO, STO)] = cip_ss
cip_dict[(GTO, STO)] = cip_gs
cip_dict[(STO, GTO)] = cip_sg
cip_dict[(GTO, GTO)] = cip_gg
cip_dict[(CutSTO, CutSTO)] = cip_cut_ss

cip_op_dict[(STO, D2, STO)] = cip_s_d2_s
cip_op_dict[(STO, D2, GTO)] = cip_s_d2_g   
cip_op_dict[(GTO, D2, STO)] = cip_g_d2_s
cip_op_dict[(GTO, D2, GTO)] = cip_g_d2_g
cip_op_dict[(CutSTO, D2, CutSTO)] = cip_cut_s_d2_s
                               
cip_op_dict[(STO, D1, STO)] = cip_s_d1_s
cip_op_dict[(STO, D1, GTO)] = cip_s_d1_g
cip_op_dict[(GTO, D1, STO)] = cip_g_d1_s
cip_op_dict[(GTO, D1, GTO)] = cip_g_d1_g
cip_op_dict[(CutSTO, D1, CutSTO)] = cip_cut_s_d1_s
                               
cip_op_dict[(STO, Rm, STO)] = cip_s_rm_s
cip_op_dict[(STO, Rm, GTO)] = cip_s_rm_g
cip_op_dict[(GTO, Rm, STO)] = cip_g_rm_s
cip_op_dict[(GTO, Rm, GTO)] = cip_g_rm_g
cip_op_dict[(CutSTO, Rm, CutSTO)] = cip_cut_s_rm_s

def cip_sss(a, o, b):
    ao = STO(a.c * o.c, a.n + o.n, a.z + o.z)
    return cip_ss(ao, b)

cip_op_dict[(STO, STO, STO)] = cip_sss


# ==== op func ====
def d1_sto(dum, sto):
    c = sto.c
    n = sto.n
    z = sto.z
    s1 = STO(c*n,  n-1, z)
    s2 = STO(-z*c, n,   z)
    return s1+s2

def d1_cutsto(dum, sto):
    c = sto.c
    n = sto.n
    z = sto.z
    r0 = sto.r0
    s1 = CutSTO(c*n,  n-1, z, r0)
    s2 = CutSTO(-z*c, n,   z, r0)
    return s1+s2

def d1_gto(dum, gto):
    c = gto.c
    n = gto.n
    z = gto.z
    s1 = GTO(c*n,      n-1, z)
    s2 = GTO(-2.0*z*c, n+1,   z)
    return s1+s2

op_apply_dict[(D1, STO)] = d1_sto
op_apply_dict[(D1, CutSTO)] = d1_cutsto
op_apply_dict[(D1, GTO)] = d1_gto


def rm_sto(rm, s):
    f = STO(s.c, s.n+rm.m, s.z)
    return f

def rm_cutsto(rm, s):
    f = CutSTO(s.c, s.n+rm.m, s.z, s.r0)
    return f

def rm_gto(rm, s):
    f = GTO(s.c, s.n+rm.m, s.z)
    return f
    
op_apply_dict[(Rm, STO)] = rm_sto
op_apply_dict[(Rm, CutSTO)] = rm_cutsto
op_apply_dict[(Rm, GTO)] = rm_gto

# ==== convenient ====
def dr(m):
    if(m==1):
        return D1()
    if(m==2):
        return D2()
    else:
        raise Exception("m=1 or 2")

def rm(m):
    return Rm(m)


def rhankel(k, L):
    if(L != 0):
        raise Exception("Not impl")
    return STO(1.0, 0, -1.0j * k)

def rbessel(k, L):
    if(L != 0):
        raise Exception("Not impl")
    ii = 1.0j
    return STO(-0.5j, 0, -ii*k) + STO(0.5j, 0, ii*k)

def __cut_exp_func(type_func):
    def __func__(cut, exp_func):
        f = exp_func
        return type_func(f.c, f.n, f.z, cut.r0)
    return __func__

op_apply_dict[(Cut, STO)] = __cut_exp_func(CutSTO)


# ==== Basis set 3D ====
"""
def lambda_s_mat(self, other = None):
    if(other == None):
        other = self
    mat = calc_s_mat(self, other)
    return mat.reshape((self.size, other.size))

def lambda_t_mat(self, other = None):
    if(other == None):
        other = self
    mat = calc_t_mat(self, other)
    return mat.reshape((self.size, other.size))

def lambda_v_mat(self, q, xyz, other = None):
    if(other == None):
        other = self
    (x, y, z) = xyz
    mat = calc_v_mat(self, q, x, y, z, other)
    return mat.reshape((self.size, other.size))

def lambda_xyz_mat(self, lmn, other = None):
    if(other == None):
        other = self
    (nx, ny, nz) = lmn
    mat = calc_xyz_mat(self, nx, ny, nz, other)
    return mat.reshape((self.size, other.size))

def lambda_add_one_basis(self, L, M, xyz, zeta):
    (x, y, z) = xyz
    self.add_one_basis_cpp(L, M, x, y, z, zeta)

def lambda_add_basis(self, L, xyz, zeta):
    (x, y, z) = xyz
    self.add_basis_cpp(L, x, y, z, zeta)
    
SphericalGTOSet.s_mat = lambda_s_mat
SphericalGTOSet.t_mat = lambda_t_mat
SphericalGTOSet.v_mat = lambda_v_mat
SphericalGTOSet.xyz_mat = lambda_xyz_mat
SphericalGTOSet.add_one_basis = lambda_add_one_basis
SphericalGTOSet.add_basis = lambda_add_basis
"""
