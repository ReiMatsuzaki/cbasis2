from linspace import *
from set_l2func import *

import numpy as np

def build_matrix(fs, op=None):
    if(op == None):
        op = OpId()
    return np.array([[cip(a, op, b) for a in fs] for b in fs])

def build_vector(fs, b):
    return np.array([cip(f, b) for f in fs])
        

class BasisSet:
    def __init__(self, fs):
        self.fs = fs

    def matrix(self, op=None):
        if(op == None):
            op = OpId()
        return np.array([[cip(a, op, b) for a in self.fs] for b in self.fs])

    def vector(self, b):
        return np.array([cip(a, b) for a in self.fs])

    def func(self, cs):
        return linear_combination(cs, self.fs)
