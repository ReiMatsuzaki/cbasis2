from symmolint import *
import minieigen as me

## partial wave expansion using cartesian GTOs


## ==== Spherical Harmonic conversion ====
# s : Solid Spherical Harmonics
# y : Complex Spherical Harmonics

def s2y(s_map, L):
    """
    Inputs
    -------
    s_map : map<int, complex> : s_map[M] gives <0|S_LM> or <0|S_LM|1>.
    
    Returns
    -------
    y_map : map<int, complex> : y_map[M] gives <0|Y_LM> or <0|S_LM|1>.

    M > 0 case
    <0 | Y_LM> = (-)^M / sqrt(2) * (<0|S_LM + 1.0j S_L-M> )
    .          = (-)^M / sqrt(2) * (<0|S_LM> + 1.0j <0|S_L-M> )
    """
    
    y_map = {}
    y_map[0] = s_map[0]
    for M in range(1, L+1):
        y_map[M] = ((-1)**M) * (s_map[M] + 1.0j * s_map[-M]) / np.sqrt(2.0)
        y_map[-M]= 1.0       * (s_map[M] - 1.0j * s_map[-M]) / np.sqrt(2.0)

    return y_map

def s2y_bra(s_map, L):
    """
    <Y_LM | 0> = <0 | Y_LM>^* = <0^* | Y_LM^* > 
    .          = (-)^m <0^* | Y_L-M>
    .          = (-)^m /sqrt(2) <0^* | S_LM - 1.0j *S_L-M>
    .          = (-)^m /sqrt(2) (<0^* | S_LM> -1.0j <0^*|S_L-M>)
    .          = (-)^m /sqrt(2) (<S_LM|0> -1.0j <S_L-M|0>)
    """
    y_map = {}
    y_map[0] = s_map[0]
    for M in range(1, L+1):
        y_map[+M] = ((-1)**M) * (s_map[M] - 1.0j * s_map[-M]) / np.sqrt(2.0)
        y_map[-M] =             (s_map[M] + 1.0j * s_map[-M]) / np.sqrt(2.0)
    return y_map

def ss0_2_yy0(ss0, L, Lp):
    """
    { <S_LM | S_L'M'|const> | MM'} 
    to
    { <Y_LM | Y_L'M'|const> | MM'}
    
    S : Normalized Solid Spherical Harmonics
    Y : Normalized Complex Spherical Harmonics

    Basic relations between two type of spehrical harmonic functions are:
    .   Y(L0) = S(L0)
    .   Y(LM) = (-)^m / sqrt(2) * [S(LM) + iS(L-M)]
    .   Y(L-M)= 1     / sqrt(2) * [S(LM) - iS(L-M)]
    where M > 0. 

    Inputs
    ------
    ss0 : map<pair<int, int>, complex> : ss0[M,M'] gives <aLM|bL'M'|i> 

    Algorithm:
    <Y_LM | Y_L'M' | 0> = <0 | Y_L'M'^* | Y_LM^*> 
    """

    Mps = range(-Lp, Lp+1)
    Ms  = range(-L, L+1)

    for Mp in Mps:
        for M in Ms:
            if (M, Mp) not in ss0.keys():
                raise(Exception("({0}, {1}) is not found".format(M, Mp)))

    # result
    sy0 = {}
    for M in range(-L, L+1):
        tmp_s = {Mp:ss0[M, Mp] for Mp in Mps}
        tmp_y = s2y(tmp_s, Lp)
        sy0.update( {(M,Mp) : tmp_y[Mp] for Mp in Mps})
            
    yy0 = {}
    for Mp in range(-Lp, Lp+1):
        tmp_s = {M: sy0[M, Mp] for M in Ms}
        tmp_y = s2y_bra(tmp_s, L)
        yy0.update( { (M,Mp) : tmp_y[M] for M in Ms } )

    return yy0

def ss0_2_yy0_old(ss0, L, Lp):
    # result
    yy0 = {}
    
    for M in range(0,L+1):
        for Mp in range(0, Lp+1):
            if(M == 0 and Mp == 0):
                yy0[0, 0] = ss0[M, Mp]
            elif(M == 0):
                yy0[0,+Mp] = (((-1)**Mp)/np.sqrt(2.0)
                              *(ss0[0,Mp]+1.0j*ss0[0,-Mp]))
                yy0[0,-Mp] = (1.0/np.sqrt(2.0)
                              *(ss0[0,Mp]-1.0j*ss0[0,-Mp]))
            elif(Mp == 0):
                yy0[+M, 0] = (((-1)**M)/np.sqrt(2.0)
                              *(ss0[M,0]-1.0j*ss0[-M,0]))
                yy0[-M, 0] = (1.0/np.sqrt(2.0)
                              *(ss0[M,0]+1.0j*ss0[-M,0]))
            else:
                yy0[+M,+Mp] = (((-1)**(M+Mp))/2.0 * 
                               (ss0[M,Mp] -1.0j*ss0[-M,Mp] +1.0j*ss0[M,-Mp] +ss0[-M,-Mp])) # 
                yy0[-M,+Mp] = ((-1**(Mp))/2.0 * 
                               (ss0[M,Mp] +1.0j*ss0[-M,Mp] + 1.0j*ss0[M,-Mp]  -ss0[-M,-Mp]))
                yy0[+M,-Mp] = (((-1)**(M))/2.0 * 
                               (ss0[M,Mp] -1.0j*ss0[-M,Mp] - 1.0j*ss0[M,-Mp]  -ss0[-M,-Mp]))
                yy0[-M,-Mp] = (1.0/2.0 * 
                               (ss0[M,Mp] +1.0j*ss0[-M,Mp] - 1.0j*ss0[M,-Mp]  +ss0[-M,-Mp]))
    return yy0


## ==== partial wave (Solid Spherical Harmonics) ====
def sub_solid_sh_Ms_old(L, Ms, zeta, sym):

    for M in Ms:
        if abs(M) > 1:
            raise(Exception("|M|>1 case is not supported "))
        
    if L not in [0,1,2,3,5]:
        raise(Exception("implemented  for L=0,1,2,3 only"))

    if(sym.irrep_x == sym.irrep_y or
       sym.irrep_y == sym.irrep_z or
       sym.irrep_x == sym.irrep_z):
        raise(Exception("this symmetry can not be separate x,y,z"))

    irrep_map = {-1: sym.irrep_y,
                 0:  sym.irrep_z,
                 1:  sym.irrep_x}    
    sub = SubSymGTOs().xyz((0, 0, 0)).zeta(zeta)
    
    if L == 0:
        sub.ns((0, 0, 0))
        sub.rds(Reduction(sym.irrep_x, [[1]]))
         
    elif L == 1:
        sub.ns((1, 0, 0)).ns((0, 1, 0)).ns((0, 0, 1))
        if 1 in Ms:
            sub.rds(Reduction(sym.irrep_x, [[1, 0, 0]]))
        if -1 in Ms:
            sub.rds(Reduction(sym.irrep_y, [[0, 1, 0]]))
        if 0 in Ms:
            sub.rds(Reduction(sym.irrep_z, [[0, 0, 1]]))

    elif L == 2:
        sub.ns((2,0,0)).ns((0,2,0)).ns((0,0,2))
        sub.ns((1,1,0)).ns((0,1,1)).ns((1,0,1))
        if 1 in Ms:
            sub.rds(Reduction(sym.irrep_x, [[0,0,0,
                                             0,0,1]]))
        if 0 in Ms:
            sub.rds(Reduction(sym.irrep_z, [[-1,-1,2,
                                             0,0,0]]))
        if -1 in Ms:
            sub.rds(Reduction(sym.irrep_y, [[0,0,0,
                                             0,1,0]]))
    elif L == 3:
        sub.ns((1,0,2)).ns((3,0,0)).ns((1,2,0))
        sub.ns((0,0,3)).ns((2,0,1)).ns((0,2,1))
        sub.ns((0,1,2)).ns((2,1,0)).ns((0,3,0))
        if 1 in Ms:
            sub.rds(Reduction(sym.irrep_x, [[4,-1,-1,
                                             0,0,0,
                                             0,0,0]]))
        if 0 in Ms:
            sub.rds(Reduction(sym.irrep_z, [[0,0,0,
                                             2,-3,-3,
                                             0,0,0]]))
        if -1 in Ms:
            sub.rds(Reduction(sym.irrep_y, [[0,0,0,
                                             0,0,0,
                                             4,-1,-1]]))
    elif L == 5:
        sub.ns((5,0,0)).ns((3,2,0)).ns((1,4,0))
        sub.ns((3,0,2)).ns((1,2,2)).ns((1,0,4))

        sub.ns((4,0,1)).ns((2,2,1)).ns((0,4,1))
        sub.ns((2,0,3)).ns((0,2,3)).ns((0,0,5))

        sub.ns((4,1,0)).ns((2,3,0)).ns((0,5,0))
        sub.ns((2,1,2)).ns((0,3,2)).ns((0,1,4))

        if 1 in Ms:
            sub.rds(Reduction(sym.irrep_x, [[1, 2, 1,-12,-12,8,
                                             0, 0, 0, 0, 0,0,
                                             0, 0, 0, 0, 0,0]]))
        if 0 in Ms:
            sub.rds(Reduction(sym.irrep_z, [[0, 0, 0, 0, 0,0,
                                             15,30,15,-40,-40,8,
                                             0, 0, 0, 0, 0,0]]))
        if -1 in Ms:
            sub.rds(Reduction(sym.irrep_y, [[0, 0, 0, 0, 0,0,
                                             0, 0, 0, 0, 0,0,
                                             1,2,1,-12,-12,8]]))
        
    else:
        raise(Exception("not implemented yet"))
    return sub
    
def sub_solid_sh_M_old(L, M, zeta, irrep):

    sub = SubSymGTOs().xyz((0,0,0)).zeta(zeta)
    
    if(L == 0 and M == 0):
        sub.rds(Reduction(irrep, [[1]]))
        sub.ns((0,0,0))
    elif(L == 1):
        sub.rds(Reduction(irrep, [[1]]))
        if(M == 1):
            sub.ns((1,0,0))   # S11 = x
        elif(M == 0):
            sub.ns((0,0,1))   # S10 = z
        elif(M == -1):
            sub.ns((0,1,0))   # S1-1 = y

    elif(L == 2):
        if(M == 2):
            # S22 = x^2 - y^2
            sub.ns((2,0,0)).ns((0,2,0))
            sub.rds(Reduction(irrep, [[1, -1]]))
        elif(M == 1):
            # S21 = xz
            sub.ns((1,0,1))
            sub.rds(Reduction(irrep, [[1]]))
        elif(M == 0):
            # S20 = -x^2 - y^2 + 2z^2 = 3z^2 - r^2
            sub.ns((2,0,0)).ns((0,2,0)).ns((0,0,2))
            sub.rds(Reduction(irrep, [[-1, -1, 2]]))
        elif(M == -1):
            # S2-1 = yz
            sub.ns((0,1,1))
            sub.rds(Reduction(irrep, [[1]]))
        elif(M == -2):
            # S2-2 = xy
            sub.ns((1,1,0))
            sub.rds(Reduction(irrep, [[1]]))
    elif(L == 3):
        if(M == 1):
            # S31 = 4xz^2 - x^3 - xy^2
            sub.ns((1,0,2)).ns((3,0,0)).ns((1,2,0))
            sub.rds(Reduction(irrep, [[4, -1, -1]]))
        elif(M == 0):
            # S30 = 2z^3 -3 x^2z -3 y^2z
            sub.ns((0,0,3)).ns((2,0,1)).ns((0,2,1))
            sub.rds(Reduction(irrep, [[2, -3, -3]]))
        elif(M == -1):
            # S3-1 = 4yz^2 -x^2y - y^3
            sub.ns((0,1,2)).ns((2,1,0)).ns((0,3,0))
            sub.rds(Reduction(irrep, [[4, -1, -1]]))
        else:
            raise(Exception("not implemented yet"))
    elif(L == 5):
        if(M == 1):
            # S51 = x^4y + 2x^2y^3 + y^5 -12x^2yz^2 -12y^3z^2 + 8yz^4
            sub.ns((5,0,0)).ns((3,2,0)).ns((1,4,0)).ns((3,0,2)).ns((1,2,2)).ns((1,0,4))
            sub.rds(Reduction(irrep, [[1,2,1,-12,-12,8]]))
        elif(M == 0):
            # S50 = 15x^4z + 30x^2y^2z + 15y^4z - 40x^2z^3 - 40y^2z^3 + 8z^5
            sub.ns((4,0,1)).ns((2,2,1)).ns((0,4,1)).ns((2,0,3)).ns((0,2,3)).ns((0,0,5))
            sub.rds(Reduction(irrep, [[15,30,15,-40,-40,8]]))
        elif(M == -1):
            # S5m1 = x^4y + 2x^2y^3 + y^5 -12x^2yz^2 -12y^3z^2 + 8yz^4
            sub.ns((4,1,0)).ns((2,3,0)).ns((0,5,0)).ns((2,1,2)).ns((0,3,2)).ns((0,1,4))
            sub.rds(Reduction(irrep, [[1,2,1,-12,-12,8]]))
        else:
            raise(Exception("not implemented yet"))
    else:
        raise(Exception("not implemented yet"))

    return sub
    
def sub_solid_sh_old(L, M_or_Ms, zeta, irrep_or_sym):
    """ gives solid spherical harmonics
    see http://www.f-denshi.com/000TokiwaJPN/14bibnh/436cub.html
    S_l0 = Y_l0
    S_lm = ((-)^m Y_lm + Y_l-m)/sqrt(2) = sqrt(2)(-)^m Re[Y_lm]            (m>0)
    S_lm = i^-1 {(-1)^|m| - Y_{l,-|m|}}/sqrt(2) = sqrt(2) (-)^m Im[Y_l|m|] (m<0)
    """
    
    if(type(M_or_Ms) == list):
        Ms = M_or_Ms
        sym = irrep_or_sym
        return sub_solid_sh_Ms_old(L, Ms, zeta, sym)
    else:
        M = M_or_Ms
        irrep = irrep_or_sym
        return sub_solid_sh_M_old(L, M, zeta, irrep)

def solid_part_waves(sym, L, Ms, zeta):

    gtos = SymGTOs.create().sym(sym)
    gtos.sub(sub_solid_sh(sym, L, Ms, (0,0,0), zeta))
    """
    for M in Ms:
        print "M:", M
    gtos.sub(sub_solid_sh(sym, L, M, (0,0,0), zeta))
    """
    gtos.setup()
    
    return gtos
        
            
