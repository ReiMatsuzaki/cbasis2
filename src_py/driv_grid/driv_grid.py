import os
import sys
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg

"""
Solve radial driven equation with simple grid method.
"""

# ==== utility ====
def flatten(xss):
    """
    Inputs
    -------
    xss : [[x]]
    .     list of list of something
    
    Returns
    --------
    xs   : [x]
    .    flattened list
    """
    return reduce(lambda a,b:a+b, xss)

def laplacian_mat(n):
    """
    gives laplacian matrix with n grid points
    
    Inputs
    ------
    n : int
    .   number of grid points

    Returns
    -------
    mat : sparse matrix
    .     laplacian matrix
    """
    data = [1, -2, 1]*n
    i = flatten([[k,k,k] for k in range(n)])
    j = flatten([[k-1, k, k+1] for k in range(n)])
    return scipy.sparse.coo_matrix((data[1:-1], (i[1:-1], j[1:-1])))

def bc_outgoing_mat(n, h, k):
    """
    gives matrix representing outgoing boundary condition
    
    Inputs
    ------
    n : int
    .   number of grid
    h : double
    .   width of grid
    k : double
    .   wave number
    """
    
    d = [1.0, 2.0j*k*h]
    i = [n-1, n-1]
    j = [n-2, n-1]
    return scipy.sparse.coo_matrix((d, (i, j)))

# ==== solver ====
def solve_driv(v, ene, s, n, h):
    """
    solve (T+V-E)f=s with outgoing boundary condition
    
    Inputs
    ------
    v: lambda scalar->scalar
      potential
    ene: scalar
      E in (T+V-E)
    s: lambda scalar->scalar
      driven term
    n: int
      number of grid points
    h: double
      grid width

    Return
    xs: [double]
    .   grid points
    ys  [complex]
    .   function values on grid points
    -------
    
    """

    xs = np.array([(k+1)*h for k in range(n)])
    h2 = h*h
    k = np.sqrt(2.0*ene)
    
    vs = [v(x)-ene for x in xs]

    mat = laplacian_mat(n) -2.0 * h2 * scipy.sparse.diags(vs, 0) + bc_outgoing_mat(n, h, k)
    vec = np.array([-2.0*h*h*s(x) for x in xs])

    ys = scipy.sparse.linalg.spsolve(mat, vec)
    return (xs, ys)

def solve_h_length(ene, channel, n, h):

    raise(Exception("Not supported. Use solve_h instead"))

    (n0, l0, l1) = channel
    v = lambda x: -1.0/x + l1*(l1+1)*0.5/(x*x)
    driv_term = HAtom(1.0).length(n0, l0, l1)
    s = lambda x: driv_term.at(x)
    return solve_driv(v, ene, s, n, h)

def solve_h(ene, channel, dipole, n, h):
    """
    solve (H-E)f=s with outgoing boundary condition.

    Inputs
    ------
    ene : scalar
    .   E in operator (H-E)
    dipole : string
    .   dipole type ("length" or "velocity")
    channel : integer tuple 
    .   (n0, l0, l1) where n0 is initial quantum number, l0 and l1 are 
    .   initial and final angular quantum number 
    n : number of grids
    h : grid width

    Returns
    -------
    xs : [double]
    .     grid points
    ys : [complex]
    .     function values on each grid points
    """

    (n0, l0, l1) = channel
    v = lambda x: -1.0/x + l1*(l1+1)*0.5/(x*x)
    
    if(dipole == "length"):
        driv_term = HAtom(1.0).length(n0, l0, l1)
    elif(dipole == "velocity"):
        driv_term = HAtom(1.0).velocity(n0, l0, l1)
    else:
        raise(Exception("dipole must be length or velocity"))
    
    s = lambda x: driv_term.at(x)
    return solve_driv(v, ene, s, n, h)
