import numpy as np
from symmolint_bind import *
def SymGTOs_at_r_ylm(gtos, L, M, irrep, cs, rs):
    (ys, ds) = gtos.at_r_ylm_cpp(L, M, irrep, cs, rs)
    return (np.array(ys), np.array(ds))
    #return np.array(gtos.at_r_ylm_cpp(L, M, irrep, cs, rs))

ClassSymGTOs.at_r_ylm = SymGTOs_at_r_ylm

def calc_v_mat(g1, g2, xyz, q):

    g1p = g1.clone()
    g2p = g2.clone()
    mole = Molecule(g1.sym_group, [Atom("X", q, [xyz])])

    g1p.molecule = mole; g1p.setup()
    g2p.molecule = mole; g2p.setup()

    mat = calc_mat(g1p, g2p, True)
    return mat["v"]


def Atom(name, q, xyz_list=[]):
    at = new_Atom(name, q)
    for xyz in xyz_list:
        at.add(xyz)
    return at

def Molecule(sym, atoms=[]):
    mole = new_Molecule(sym)
    for atom in atoms:
        mole.add(atom)
    return mole

def SymGTOs(mole, subs=[]):
    gtos = new_SymGTOs(mole)
    for (atom_name, ns_list, rds_list, zeta) in subs:
        atom = mole.atom(atom_name)
        sub = SubSymGTOs(mole.sym_group(), atom)
        sub.zeta(zeta)
        for ns in ns_list:
            sub.ns(ns)
        for rds in rds_list:
            sub.rds(rds)
        gtos.sub(sub)
    return gtos
    
        
    
