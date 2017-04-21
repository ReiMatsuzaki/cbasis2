"""
remove linear dependence of basis function set.
2015/12/6 R.Matsuzaki
"""
import scipy.linalg as la
from linspace import cnormalize, cnorm
from basis_set import build_matrix


def eig_overlap(fs):
    normalized_fs = [cnormalize(f) for f in fs]
    s = build_matrix(normalized_fs)
    return sorted(map(abs, la.eig(s)[0]))

def with_index(xs):
    return zip(range(len(xs)), xs)

def clength(f, g):
    return abs(cnorm(cnormalize(f) - cnormalize(g)))

def clen_list(fs):
    """return length list in the meaning of c-Norm"""
    return [(clength(fi, fj), fi, i, fj, j)
            for (i, fi) in with_index(fs) 
            for (j, fj) in with_index(fs) 
            if i<j]

def close_pair(fs):
    return sorted(clen_list(fs), key=lambda x:x[0])[0]

def contract_one(fs):
    (d, fi, i, fj, j) = close_pair(fs)
    others = [fk for (k, fk) in with_index(fs) if k!=i and k!=j]
    return [fi+fj] + others

def remove_lindep(fs, eps):
    if len(fs) == 1:
	return fs
    eig = eig_overlap(fs)[0]
    if eig > eps:
	return fs
    fs_prime = contract_one(fs)
    return remove_lindep(fs_prime, eps)    
