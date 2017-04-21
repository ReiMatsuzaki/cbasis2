import time
import sys
sys.path.append('../')
from l2func import *
import numpy as np

def create_basis(num) :
    return [ normalized_sto(2, 0.1 * (i + 1)) for i in range(num)]

us = create_basis(5)
hatom = HLikeAtom(1, 1.0, 0)
op = hatom.h_minus_energy_sto(0.5)

t0 = time.clock()
acc = 0.0
for i in us:
    l_i = op.operate(i)
    for j in us:
        acc += sym_ip_ss(i, l_i)
t1 = time.clock()
print "op_5STO : ", t1 - t0

us = create_basis(50)
t0 = time.clock()
acc = 0.0
for i in us:
    l_i = op.operate(i)
    for j in us:
        acc += sym_ip_ss(j, l_i)
t1 = time.clock()
print "op_50STO : ", t1 - t0

us = np.array(us)
t0 = time.clock()
acc = 0.0
for i in us:
    l_i = op.operate(i)
    for j in us:
        acc += sym_ip_ss(j, l_i)

t1 = time.clock()
print "op_50_np_STO : ", t1 - t0
    



