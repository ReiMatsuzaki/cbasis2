from l2func_molint_bind import *
import scipy.linalg as la
import numpy as np

def sort_eigs(eigvals, eigvecs):

    val_vec_list = sorted(zip(eigvals, eigvecs.T), key=lambda x: x[0].real)
    vals_out = np.array([val for (val, vecs) in val_vec_list])
    vecs_out = np.array([vecs for (val, vecs) in val_vec_list]).T
    return (vals_out, vecs_out)

def cnorm(x):
    return np.sqrt(np.sum(x*x))

def matrix_inv_sqrt(s):
    (eigs, vecs) = la.eig(s)

    # -- normalize the eigen vector --
    vecs = np.array([vt/cnorm(vt) for vt in vecs.T]).T

    # -- compute lambda_ij = delta_ij / sqrt(eig_i)
    lambda_mat = np.diag(np.array([1.0/np.sqrt(eig) for eig in eigs]))

    # -- compute S^(-1/2) = D lambda D^T --
    s_inv_sqrt = np.dot(vecs, np.dot(lambda_mat, vecs.transpose()))
    return s_inv_sqrt

def sym_eig(h, s):

    s2inv = matrix_inv_sqrt(s)
    hp = np.dot(s2inv, np.dot(h, s2inv))
    
    (eigs, vecs) = la.eig(hp)
    vecs = np.dot(s2inv, vecs)
    vecs = np.array([vt/np.sqrt(np.dot(vt, np.dot(s, vt))) for vt in vecs.T]).T
    return sort_eigs(eigs, vecs)

def c2_eig(h, s):
    (eigs, vecs) = la.eig(h, s)
    vecs = np.array([vt/np.sqrt(np.dot(vt, np.dot(s, vt))) for vt in vecs.T]).T
    return sort_eigs(eigs, vecs)


# ==== GTOs ====
def GTOs_at_r_ylm(self, l, m, rs, cs):
    rs_val = np.array(rs, dtype=complex)
    cs_val = np.array(cs, dtype=complex)
    vs = self.at_r_ylm_cpp(l, m, rs_val, cs_val)
    return vs

def GTOs_add_atom(self, q, xyz):
    (x, y, z) = xyz
    self.add_atom_cpp(q, x, y, z)

def GTOs_add_gtos(self, L, xyz, zeta):
    (x, y, z) = xyz
    self.add_gtos_cpp(L, x, y, z, zeta)

def GTOs_add_one_gto(self, L, M, xyz, zeta):
    (x, y, z) = xyz
    self.add_one_gto_cpp(L, M, x, y, z, zeta)    

GTOs.at_r_ylm = GTOs_at_r_ylm
GTOs.add_atom = GTOs_add_atom
GTOs.add_one_gto = GTOs_add_one_gto
GTOs.add_gtos = GTOs_add_gtos
