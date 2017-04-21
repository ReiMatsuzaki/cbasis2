import os
import numpy as np
import scipy.linalg as la
import pandas as pd

import minieigen as me
from local import *

sym = D2h()
R0 = 2.0
mole = Molecule(sym,
                [Atom("H", 1.0, [(0,0,+R0/2),
                                 (0,0,-R0/2)]),
                 Atom("Cen", 0.0, [(0, 0, 0)])])

def get_gi_RM(reduce_num=0):
    """ see T.Rescigno and C.McCurdy(1985) PRA 31, 624
    """
    mult = lambda x: 3.0*x
    zeta_s_h = map(mult,
		   [0.0285649, 0.0812406, 0.190537, 0.463925, 1.202518,
		    3.379649, 10.60720, 38.65163, 173.5822, 1170.498])
    
    ## see lin_dip/
    if(reduce_num > 3):
        raise(Exception("not implemented"))

    zeta_p_h = map(mult,
		   [0.015442, 0.035652, 0.085676,
		    0.227763, 0.710128, 3.009711])
    zeta_d_h = [15.0, 5.0, 5.0/3.0, 5.0/9.0]
    zeta_s_cen = [14.4/(3.0**n) for n in range(12)]
    zeta_d_cen = [10.0/(3.0**n) for n in range(9)]

    if(reduce_num==1):
        del zeta_s_cen[3]
    if(reduce_num==2):
        del zeta_s_cen[3]
        del zeta_s_cen[5]
    if(reduce_num==3):
        del zeta_s_cen[3]
        del zeta_s_cen[4]
        del zeta_s_cen[5]
        
    g_i = (SymGTOs(mole)
	   .sub(SubSymGTOs(sym, mole.atom("H"))
		.ns( (0,0,0))
		.rds(Reduction(sym.irrep_s, [[1], [1]]))
		.zeta(zeta_s_h))
	   .sub(SubSymGTOs(sym, mole.atom("H"))
		.ns( (0,0,1))
		.rds((Reduction(sym.irrep_s, [[1], [-1]])))
		.zeta(zeta_p_h))
	   .sub(SubSymGTOs(sym, mole.atom("H"))
		.ns( (0,0,2))
		.rds(Reduction(sym.irrep_s, [[1], [1]]))
		.zeta(zeta_d_h))
	   .sub(SubSymGTOs(sym, mole.atom("Cen"))
		.ns( (0,0,0))
		.rds(Reduction(sym.irrep_s, [[1]]))
		.zeta(zeta_s_cen))
           .sub(SubSymGTOs(sym, mole.atom("Cen"))
		.ns( (0,0,2))
		.rds(Reduction(sym.irrep_s, [[1]]))
		.zeta(zeta_d_cen))
	   .setup())
    return g_i

def et_zeta(num, z0, z1):
    r = (z1/z0)**(1.0/(num-1))
    return [z0*r**n for n in range(num)]


## below basis is from /Users/rei/src/git/opt_cbf/lstsq_fit/calc/win_irr_coulomb_atan/5basis
def fitted_zeta(L, ene):
    zeta_p_05 = [0.0812755955262-0.0613237786222j,
	         0.00497147387796-0.0113737972763j,
                 0.0323712622673-0.0451300037076j,
                 0.00317417887792-0.022582254987j,
                 0.0118391719646-0.0327847352576j]
    zeta_p_10 = [0.126638175526-0.105030118631j,
	         0.00829648213915-0.0206661218928j,
	         (0.0518178348978-0.0778901537559j),
	         (0.00655155739796-0.0427765414833j),
	         (0.0204086224266-0.0573315334604j)]
    zeta_f_5 = [0.114846838205-0.0528053215279j,
                0.0502495991466-0.0443217040212j,
                0.00291922532665-0.0190617777251j,
                0.0223369214731-0.0358454508219j,
                0.0083110955684-0.0268941896615j]
    zeta_f_10 = [0.176407650439-0.092875535931j,
	         0.0775939426597-0.0773920329099j,
	         0.0072361335156-0.0306232949503j,
	         0.0344401450273-0.0624797301291j,
	         0.0118345934315-0.0471053022308j]
    zeta_map = {}
    zeta_map[1, 0.5] = zeta_p_05
    zeta_map[1, 1.0] = zeta_p_10
    zeta_map[3, 0.5] = zeta_f_5
    zeta_map[3, 1.0] = zeta_f_10
    zeta_map[5, 0.5] = zeta_p_05
    zeta_map[5, 1.0] = zeta_p_10
    return zeta_map[L, ene]


def get_gpsi1_simple(ls, zs, add_sub_list=[]):
    g_psi1 = SymGTOs(mole)

    for L in ls:
        g_psi1.sub(sub_solid_sh(sym, mole.atom("Cen"), L, [-1,0,1], zs))

    for sub in add_sub_list:
        g_psi1.sub(sub)
    g_psi1.setup()
    return g_psi1

def get_gpsi1(ls, add_zs, add_sub_list=[], ene=0.5):
    g_psi1 = SymGTOs(mole)

    for L in ls:
        zs = fitted_zeta(L, ene) + add_zs
        g_psi1.sub(sub_solid_sh(sym, mole.atom("Cen"), L, [-1,0,1], zs))

    for sub in add_sub_list:
        g_psi1.sub(sub)
    g_psi1.setup()
    return g_psi1


def get_gpsi0_L_simple(ls, zs):
    gpsi0_L = {}

    for L in ls:
        gpsi0_L[L] = SymGTOs(mole)
        gpsi0_L[L].sub(sub_solid_sh(sym, mole.atom("H"), L, [-1,0,1], zs))
        gpsi0_L[L].setup()
    
    return gpsi0_L
    
    
def get_gpsi0_L(ls, add_zs, ene=0.5):
    gpsi0_L = {}

    for L in ls:
        zs = fitted_zeta(L, ene) + add_zs
        gpsi0_L[L] = SymGTOs(mole)
        gpsi0_L[L].sub(sub_solid_sh(sym, mole.atom("Cen"), L, [-1,0,1], zs))
        gpsi0_L[L].setup()
    
    return gpsi0_L


def get_gchi_L(ls, zeta_chi):
    res = {}
    for L in ls:
        gtos = (SymGTOs(mole)
                .sub(sub_solid_sh(sym, mole.atom("Cen"), L, [-1,0,1], [zeta_chi]))
                .setup())
        res[L] = gtos
    return res

