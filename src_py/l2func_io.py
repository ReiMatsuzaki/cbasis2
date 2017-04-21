"""
input/output in l2func
2015/12/6 R.Matsuzaki
"""

from set_l2func import *


def write_func_list(fs, fn):
    """
    Inputs
    ------
    fs : [function]
          function list to write

    fn : string 
         file name
    """

    with open(fn, "w") as fp:
        for f in fs:
            if(type(f) == STO):
                basis = "STO"
            if(type(f) == GTO):                
                basis = "GTO"
            fp.write("{0} {1} {2} {3}\n".format(basis, f.c, f.n, f.z))


def read_func_list(fn):
    def one_line(line):
        [str_type, str_c, str_n, str_z] = line.split(" ")
        c = complex(str_c)
        n = int(str_n)
        z = complex(str_z)
        if(str_type == "STO"):
            return STO(c, n, z)
        if(str_type == "GTO"):
            return GTO(c, n, z)

    with open(fn) as fp:
        return map(one_line, fp.readlines())
        
