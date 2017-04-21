import numpy as np
import os
import os.path
import sys

def append_on_git(dirname):
    dn = os.environ['HOME'] + "/src/" + dirname
    sys.path.append(dn)

append_on_git("cbasis/src_py/coulomb")
import coulomb

append_on_git("cbasis/src_py/lstsq_fit")
import gau_lstsq

def lstsq_from_region(num_gto, num_grid, w, pn, target):
    try:
        return gau_lstsq.lstsq(m = num_gto,
                               n = num_grid,
                               w = w,
                               pn = pn,
                               f = target)
    except:
        return None
