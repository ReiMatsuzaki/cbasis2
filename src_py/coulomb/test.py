import numpy as np
import coulomb

(f, g, fp, gp, err)= coulomb.coulomb(0.2, 0.3, 0.0, 0)
print (f, g)
print (fp, gp)
print err

y = coulomb.outgoing(-1.0, 100.0, 0.0, 0)
arg = coulomb.asym_arg(-1.0, 100.0, 0)
print y, np.exp(1.0j*arg)
