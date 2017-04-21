from symmolint import *
from part_wave import *
from coulomb import *
import scipy.linalg as la
from numpy.linalg import inv


## ==== constant ====
c_light = 137.035999258
au2mb = pow(5.291772, 2)
au2ev = 27.2114

## ==== calculation ====
def calc_coef_Azeta0(w, Alm, zeta, ls):
    """ see
    J. C. Tully, R. S. Berry, and B. J. Dalton, 
    Physical Review 176, 95 (1968)
    """

    c0 = 4.0*np.pi*np.pi/(c_light*w)
    M_g = 0
    M_i = 0
    cumsum = 0.0
    for L1 in ls:
	for L2 in ls:
	    c1 = np.sqrt((2*L1+1)*(2*L2+1)*1.0)/(4.0*np.pi*(2*zeta+1))
	    for M1 in range(-L1,L1+1):
		for M2 in range(-L2,L2+1):
		    for mzeta in range(-zeta, zeta+1):
			c2 = pow(-1, M_g+M_i)
			c3 = (cg_coef(1, 1, M_i-M1-M_g, M2+M_g-M_i, zeta, mzeta)*
			      cg_coef(1, 1, 0,          0,          zeta, 0)*
			      cg_coef(L2, L1, M2,       -M1,        zeta, mzeta)*
			      cg_coef(L2, L1, 0,        0,          zeta, 0))
			c4 = Alm[L1,M1] * Alm[L2,M2].conjugate()
#                        if zeta==2:
#                            print "c:" , c2, c3, c4
		        cumsum += c0 * c1 * c2 * c3 * c4
    return cumsum

def calc_coulomb_cck_mat(S00, H00, ene, c0_chi):
    L0 = S00 * ene - H00
    c0_psi0 = la.solve(L0, c0_chi)
    
    # see the comment in the below function "one_e_pi"
    im_psi0_chi = np.dot(c0_psi0, c0_chi).imag
    
    k = np.sqrt(2.0 * ene)
    sign = im_psi0_chi / abs(im_psi0_chi)    
    c0_npsi =  sign * np.sqrt(k/2.0) / np.sqrt(sign * im_psi0_chi) * c0_psi0
    return c0_npsi

def calc_coulomb_cck(g_psi0, g_chi, ir, Z, ene):
    """
    gives the coefficient of Im[~Psi0] to representing coulomb function
    with unit amplitude.
    
    Inputs
    -------
    g_psi0 : SymGTOs : basis for psi0
    g_chi  : SymGTOs : basis for driven term
    ir     : Irrep   
    Z      : double  : charge of nucleus
    ene    : double  : energy
    """

    mat_00 = calc_mat_complex(g_psi0, False)
    v0 = calc_v_mat(g_psi0, g_psi0, [0,0,0], Z)

    c0_chi = calc_mat(g_psi0, g_chi, False)["s"][ir, ir].col(0)

    res = calc_coulomb_cck_mat(mat_00["s"][ir, ir],
                               mat_00["t"][ir, ir] + v0[ir, ir],
                               ene, c0_chi)
    return res
    
    """
    L0 = (mat_00["s"][irrep,irrep]*ene
    -mat_00["t"][irrep,irrep]
    -v0[irrep, irrep])

    # compute coefficient of psi0
    c0_psi0 = la.solve(L0, c0_chi)

    # see the comment in the below function "one_e_pi"
    im_psi0_chi = np.dot(c0_psi0, c0_chi).imag

    k = np.sqrt(ene*2.0)
    sign = im_psi0_chi / abs(im_psi0_chi)    
    c0_npsi =  sign * np.sqrt(k/2.0) / np.sqrt(sign * im_psi0_chi) * c0_psi0
    #c0_npsi =  -np.sqrt(k/2.0) / np.sqrt(-im_psi0_chi) * c0_psi0
    return c0_npsi
    """

def solve_init(g_i):
    mat = calc_mat_complex(g_i, True)
    h = mat.get_matrix("t", 0, 0) + mat.get_matrix("v", 0, 0)
    s = mat.get_matrix("s", 0, 0)
    tmp = ceig(h,s)
    ei = tmp[0][0].real
    ci = tmp[1].col(0)
    return (ei, ci)

def solve_init_canonical(g_i, num0):
    mat = calc_mat_complex(g_i, True)
    h = mat["t"][0, 0] + mat["v"][0, 0]
    s = mat["s"][0, 0]
    tmp = ceig_canonical(h, s, num0)
    ei = tmp[0][0].real
    ci = tmp[1].col(0)
    return (ei, ci)

def dip_list(sym, dipole):
    if dipole == "length":
        xyz = ["y", "z", "x"]
    elif dipole == "velocity":
        xyz = ["dy", "dz", "dx"]
    else:
        raise(Exception("dipole must be length or velocity"))
    return zip([-1,0,1],
               xyz,
               [sym.irrep_y, sym.irrep_z, sym.irrep_x])        

## ==== main ====
def t2k(t_mat):
    (numi, numj) = t_mat.shape
    if numi != numj:
        raise(Exception("t_mat is not diagonal"))
    num = numi
    i_mat = np.identity(num)
    s_mat = t_mat + i_mat
    k_mat = 1.0j * (i_mat-s_mat) * inv(i_mat+s_mat)
    return k_mat
    
class PhotoIonizationCoulombPlus():
    """
    Solve photoionization problem by separating potential 
    to coulomb and other potential.
    H = H0 + V
    H0 : kinetic energy + coulomb potential
    V  : other potential

    TDM = <phi_0 | mu phi_i + V psi_1>
    """
    
    def __init__(self, Z, ne, ls, dipole="length", Mi=0, ir_i=0):

        if(type(ls) != list):
            raise Exception("ls must be list of integer")
        if(len(ls) == 0):
            raise Exception("ls is null")
        if(type(ls[0]) != int ):
            raise Exception("ls must be list of integer")
        if(dipole not in ["length", "velocity"]):
            raise Exception("dipole must be length or velocity")

        self.Z = Z    ## charge for Hamiltonian H0
        self.ne = ne  ## number of electron
        self.ls = ls  ## L values for partial wave expansions
        self.dipole = dipole # dipole operator("length" or "velocity")
        self.Mi = Mi  ## M for initial electron orbital
        self.ir_i = ir_i ## irreducible representation for init
        self.ei = 0.0 ## minus ionization energy
        self.use_2 = True ## True => use 2 term for TDM

        self.Sii = None
        self.Hii = None
        
        self.S11 = {}
        self.H11 = {}
        self.mu_phi_1 = {}

        self.S00 = {}
        self.H00 = {}
        self.c0_chi = {}        
        self.mu_phi_0H = {}
        self.mu_phi_0C = {}
        self.v_H = {}
        self.v_C = {}        

    def precalc_common(self, g_i, g_chi_L, g_psi1, g_psi0_L):
        
        """
        g_i     : SymGTOs : basis for initial orbitals
        g_chi_L :{SymGTOs}: g_chi_L[L] gives for L th partial wave basis
        g_psi1  : SymGTOs : basis for psi1
        g_psi0_L:{SymGTOs}: g_psi0_L[L] gives for L partial wave basis
        
        Assumption:
        g_chi_L[L] and g_psi0_L[L] are assumed to distingish each 
        partial wave by M
        """

        ls = self.ls

        ## for psi1/initial
        for (mpp, xyz, ir) in dip_list(g_i.sym_group, self.dipole):
            m = mpp + self.Mi
            xyz_1i = calc_mat(g_psi1, g_i, False)[xyz][ir, self.ir_i]
            self.mu_phi_1[m] = np.dot(xyz_1i, self.ci)

        ## for psi0        
        for L in ls:
            g_psi0 = g_psi0_L[L]
            g_chi  = g_chi_L[L]
            mat_00 = calc_mat_complex(g_psi0, False)
            v_00 = calc_v_mat(g_psi0, g_psi0, [0,0,0], self.Z)
            
            mat_H_0i = calc_mat(g_psi0.conj(), g_i, False)
            mat_C_0i = calc_mat(g_psi0,        g_i, False)
            
            mat_0chi = calc_mat(g_psi0, g_chi, False)
            
            for (mpp, xyz, ir) in dip_list(g_i.sym_group, self.dipole):
                m = mpp + self.Mi
                xyz_H_0i = mat_H_0i[xyz][ir, self.ir_i]
                xyz_C_0i = mat_C_0i[xyz][ir, self.ir_i]

                self.S00[L, m] = mat_00["s"][ir, ir]
                self.H00[L, m] = mat_00["t"][ir, ir] + v_00[ir, ir]
                
                self.mu_phi_0H[L, m] = np.dot(xyz_H_0i, self.ci)
                self.mu_phi_0C[L, m] = np.dot(xyz_C_0i, self.ci)
                
                self.c0_chi[L, m] = mat_0chi["s"][ir, ir].col(0)

    def precalc_1ele(self, g_i, g_chi_L, g_psi1, g_psi0_L, num0=None):

        ls = self.ls

        if self.ne != 1:
            raise(Exception("only ne = 1 is supported"))

        if num0==None:
            (self.ei, self.ci) = solve_init(g_i)
        else:
            (self.ei, self.ci) = solve_init_canonical(g_i, num0)

        ## psi1
        mat_11 = calc_mat_complex(g_psi1, True)        
        for (mpp, xyz, ir) in dip_list(g_i.sym_group, self.dipole):
            m = mpp + self.Mi
            self.S11[m] = mat_11["s"][ir, ir]
            self.H11[m] = mat_11["t"][ir, ir] + mat_11["v"][ir, ir]

        ## V between psi1 and psi0
        for L in ls:
            g_psi0 = g_psi0_L[L];
            Z = self.Z
            v_H_01 = calc_mat(g_psi0.conj(), g_psi1, True)["v"]
            v_C_01 = calc_mat(g_psi0,        g_psi1, True)["v"]
            v0_H_01= calc_v_mat(g_psi0.conj(),g_psi1, [0,0,0], Z )
            v0_C_01= calc_v_mat(g_psi0, g_psi1, [0,0,0], Z )
            
            
            for (mpp, xyz, ir) in dip_list(g_i.sym_group, self.dipole):
                m = mpp + self.Mi
                self.v_H[L, m] = v_H_01[ir, ir] - v0_H_01[ir,ir]
                self.v_C[L, m] = v_C_01[ir, ir] - v0_C_01[ir,ir]
                print "v_C",  L,m, self.v_C[L,m][1,2]
                print "v_H",  L,m, self.v_H[L,m][1,2]
                """
                print "v_C",  L,m,v_C_01[ir,ir][1,2]
                print "v0_C", L,m,v0_C_01[ir,ir][1,2]
                print "v_H",  L,m,v_H_01[ir,ir][1,2]
                print "v0_H", L,m,v0_H_01[ir,ir][1,2]
                """
        
        self.precalc_common( g_i, g_chi_L, g_psi1, g_psi0_L)

    def precalc_stex(self, g_i, g_chi_L, g_psi1, g_psi0_L, 
                     eri_method, ei_in = None, mo_in = None,
                     coef_J=1.0, coef_K=1.0):

        ls = self.ls
        
        if self.ne != 2:
            raise(Exception("only ne = 2 is supported"))

        ## set ei and ci
        if mo_in == None:
            sym = g_i.sym_group
            mat = calc_mat_complex(g_i, True)
            eri = calc_ERI_complex(g_i, eri_method)
            mo  = calc_RHF(sym, mat, eri, self.ne, 20, 0.0001, 0)
        else:
            mo = mo_in
            
        if ei_in == None:
            self.ei = mo.eigs[0][0]
        else:
            self.ei = ei_in
        self.ci = mo.C[self.ir_i, self.ir_i].col(0)
        
        ## psi1
        mat_11 = calc_mat_complex(g_psi1, True)
        eri_J_11 = calc_ERI(g_psi1, g_psi1, g_i, g_i, eri_method)
        J_11 = calc_J(eri_J_11, self.ci, 0, g_psi1, g_psi1, coef_J)
        eri_K_11 = calc_ERI(g_psi1, g_i, g_i, g_psi1, eri_method)
        K_11 = calc_K(eri_K_11, self.ci, 0, g_psi1, g_psi1, coef_K)
        
        for (mpp, xyz, ir) in dip_list(g_i.sym_group, self.dipole):
            m = mpp + self.Mi
            self.S11[m] = mat_11["s"][ir, ir]
            self.H11[m] = (mat_11["t"][ir, ir] + mat_11["v"][ir, ir]
                           + J_11[ir, ir] + K_11[ir, ir])
            
        ## V between psi1 and psi0
        for L in ls:
            g_psi0 = g_psi0_L[L]
            g_psi0_H = g_psi0.conj()
            Z = self.Z
            
            eri_JC = calc_ERI(g_psi0, g_psi1, g_i, g_i, eri_method)
            JC = calc_J(eri_JC, self.ci, 0, g_psi0, g_psi1, coef_J)
            eri_KC = calc_ERI(g_psi0, g_i, g_i, g_psi1, eri_method)
            KC = calc_K(eri_KC, self.ci, 0, g_psi0, g_psi1, coef_K)
            v1_C = calc_mat(g_psi0, g_psi1, True)["v"]
            v0_C = calc_v_mat(g_psi0, g_psi1, [0,0,0], Z)

            eri_JH = calc_ERI(g_psi0_H, g_psi1, g_i, g_i, eri_method)
            JH = calc_J(eri_JH, self.ci, 0, g_psi0, g_psi1, coef_J)
            eri_KH = calc_ERI(g_psi0_H, g_i, g_i, g_psi1, eri_method)
            KH = calc_K(eri_KH, self.ci, 0, g_psi0, g_psi1, coef_K)
            v1_H = calc_mat(g_psi0_H,   g_psi1, True)["v"]
            v0_H = calc_v_mat(g_psi0_H, g_psi1, [0,0,0], self.Z)
            
            for (mpp, xyz, ir) in dip_list(g_i.sym_group, self.dipole):
                m = mpp + self.Mi
                self.v_H[L,m] = v1_H[ir, ir]+coef_J*JH[ir, ir]+coef_K*KH[ir, ir] -v0_H[ir,ir]
                self.v_C[L,m] = v1_C[ir, ir]+coef_J*JC[ir, ir]+coef_K*KC[ir, ir] -v0_C[ir,ir]

        self.precalc_common(g_i, g_chi_L, g_psi1, g_psi0_L)

    def calc_one(self, w, use_2, sigma_kind):
        ene = w + self.ei; k = np.sqrt(2.0 * ene)

        ls = self.ls

        ## compute psi1 if necessary
        if self.use_2:
            self.c1_psi1_m = {}
            c1_psi1_m = self.c1_psi1_m
            for m in [self.Mi-1, self.Mi, self.Mi+1]:
                L1 = self.S11[m] * ene - self.H11[m]
                c1_psi1_m[m] = la.solve(L1, self.mu_phi_1[m])

        ## compute Dl(TDM)
        Alm = {}
        for L in ls:
            
            ## init
            dl_ss = {(m,mp):0.0 for m in range(-L,L+1) for mp in [-1,0,1]}
            mc_list = range(-L,L+1) if sigma_kind=="total" else [0]
            
            for m in mc_list:
                mppc = m-self.Mi
                if abs(mppc) > 1:
                    continue
                
                ## compute
                c_npsi0 = calc_coulomb_cck_mat(self.S00[L,m], self.H00[L,m], ene,
                                               self.c0_chi[L,m])
                psi0H_muphi = np.dot(c_npsi0.conj(), self.mu_phi_0H[L, m])
                psi0C_muphi = np.dot(c_npsi0       , self.mu_phi_0C[L, m])
                impsi0_muphi = (psi0H_muphi-psi0C_muphi)/(-2.0j)
                
                if use_2:
                    c1_psi1 = c1_psi1_m[m]
                    v_H = self.v_H[L, m]; v_C = self.v_C[L, m]
                    psi0H_v_psi1 = np.dot(c_npsi0.conj(), np.dot(v_H, c1_psi1))
                    psi0C_v_psi1 = np.dot(c_npsi0,        np.dot(v_C, c1_psi1))
                    impsi0_v_psi1 = (psi0H_v_psi1-psi0C_v_psi1)/(-2.0j)
                    impsi0_other = impsi0_muphi + impsi0_v_psi1
                else:
                    impsi0_other = impsi0_muphi

                print "matele: ", L, m, impsi0_other
                dl = np.sqrt(1.0*self.ne) * 2j/np.sqrt(k) * impsi0_other
                dl_ss[m, mppc] = np.sqrt(3.0/(4.0*np.pi)) * dl
                
            etal = coulomb_phase(L, -self.Z/k)
            if self.dipole == "length":
                coef = w * 1.0j * np.sqrt(2.0/3.0) * (1.0j**(-L)) * np.exp(1.0j*etal)
            if self.dipole == "velocity":
                coef = 1.0j * np.sqrt(2.0/3.0) * (1.0j**(-L)) * np.exp(1.0j*etal)
            dl_yy = ss0_2_yy0(dl_ss, L, 1)
            for m in range(-L, L+1):
                if abs(m-self.Mi) > 1:
                    Alm[L, m] = 0.0
                else:
                    Alm[L, m] = coef * dl_yy[m, m]
                    print "Alm: ",L,m, Alm[L,m]
                    
        A00 = calc_coef_Azeta0(w, Alm, 0, ls)
	A20 = calc_coef_Azeta0(w, Alm, 2, ls)
	cs = (4.0 * np.pi * A00 * au2mb).real
	beta = (A20/A00).real
        return (cs, beta)

    def calc_cs_imalpha(self, w, sigma_kind):
        ene = w + self.ei; k = np.sqrt(2.0 * ene)

        mpp_list = [-1,0,1] if sigma_kind=="total" else [0]
        
        alpha = 0.0
        for mpp in mpp_list:
            m = mpp + self.Mi
            L1 = self.S11[m] * ene - self.H11[m]
            c1_psi1_m = la.solve(L1, self.mu_phi_1[m])
            alpha_m = np.dot(c1_psi1_m,  self.mu_phi_1[m]) / 3.0
            alpha += alpha_m
        cs = 0.0
        if(self.dipole == "length"):
            cs = 4.0 * np.pi * w / c_light * alpha.imag * au2mb
        elif(self.dipole == "velocity"):
            cs = 4.0 * np.pi / (w * c_light) * alpha.imag * au2mb
        return cs
       
    def phase_shift_easy(self, w, L, irrep_Dooh):
        """compute phase shift
        
        Inputs
        ------
        w : double : photon energy
        L : int    : quantum number L
        irrep_Dooh : 'sigu' or 'piu'
        """

        ene = w + self.ei; k = np.sqrt(2.0 * ene)
        if(irrep_Dooh == "sigu"):
            mpp = 0
        elif(irrep_Dooh == "piu"):
            mpp = 1
        else:
            raise("invalid irrep_Dooh")        
        m = mpp + self.Mi
        
        ## psi1
        L1 = self.S11[m] * ene - self.H11[m]
        c1_psi1 = la.solve(L1, self.mu_phi_1[m])        

        c_npsi0 = calc_coulomb_cck_mat(self.S00[L,m], self.H00[L,m],
                                       ene, self.c0_chi[L,m])
        psi0H_muphi = np.dot(c_npsi0.conj(), self.mu_phi_0H[L, m])
        psi0C_muphi = np.dot(c_npsi0       , self.mu_phi_0C[L, m])
        impsi0_muphi = (psi0H_muphi-psi0C_muphi)/(-2.0j)

        v_H = self.v_H[L, m]; v_C = self.v_C[L, m]
        psi0H_v_psi1 = np.dot(c_npsi0.conj(), np.dot(v_H, c1_psi1))
        psi0C_v_psi1 = np.dot(c_npsi0,        np.dot(v_C, c1_psi1))
        impsi0_v_psi1 = (psi0H_v_psi1-psi0C_v_psi1)/(-2.0j)
        
        dd = impsi0_muphi + impsi0_v_psi1
        return cmath.phase(dd)
    

    def precalc_phase_shift_1ele(self, g_chi_L, g_psi0_L, ls):

        self.v_0H_0 = {}
        self.v_0_0 = {}
        for Li in ls:
            for Lj in ls:
                gi = g_psi0_L[Li]
                gj = g_psi0_L[Li]
                giH = gi.conj()
                Z = self.Z
                v1_H = calc_mat(giH,   gj, True)["v"]
                v1_C = calc_mat(gi,    gj, True)["v"]
                v0_H = calc_v_mat(giH, gj, [0,0,0], Z)
                v0_C = calc_v_mat(gi,  gj, [0,0,0], Z)

                for (mpp, xyz, ir) in dip_list(gi.sym_group):
                    m = mpp + self.Mi
                    self.v_0H_0[Li, Lj, m] = v1_H[ir, ir] - v0_H[ir,ir]
                    self.v_0_0[ Li, Lj, m] = v1_C[ir, ir] - v0_C[ir,ir]

    def t_matrix(self, ene, Ls, irrep_Dooh):
        """ Compute phase shift from calculated T-matrix
        
        T_ij = <phi_i | V + VGV | phi_j>
        
        <im y | V | im y> = <(y-y*)/2i | V | (y-y*)/2i>
        .                 = (<y|V|y>+<y*|V|y*>-<y*|V|y>-<y|V|y*>)/4
        .                 = ((y*|V|y)+(y|V|y*)-(y|V|y)-(y*|V|y*))/4
        .                 = ((y*|V|y)+(y*|V|y)*-(y|V|y)-(y|V|y)*)/4
        .                 = (Re[(y*|V|y)] -Re[(y|V|y)])/2
        <im y |VGV| im y> = ((y*|VGV|y)+(y|VGV|y*)-(y|VGV|y)-(y*|VGV|y*))/4
        """
        
        k = np.sqrt(2.0 * ene)
        if irrep_Dooh == "sigu":
            mpp = 0
        elif irrep_Dooh == "piu":
            mpp = 1

        m = mpp + self.Mi    
        ## psi1
        L11 = self.S11[m] * ene - self.H11[m]            

        num = len(Ls)
        res = np.zeros((num, num), dtype=complex)

        ci_L = {}
        for L in Ls:
            ci_0 = calc_coulomb_cck_mat(self.S00[L,m], self.H00[L,m],
                                        ene, self.c0_chi[L,m])
            ci_L[L] = ci_0
        
        for (Li, i) in zip(Ls, range(num)):
            for (Lj, j) in zip(Ls, range(num)):
                
                ci_0  = me.VectorXc(ci_L[Li])
                ci_0H = me.VectorXc(ci_L[Li].conj())
                cj_0  = me.VectorXc(ci_L[Lj])
                cj_0H = me.VectorXc(ci_L[Lj].conj())
                
                v_0H_0 = self.v_0H_0[Li, Lj, m];
                v_0_0 = self.v_0_0[Li, Lj, m]
                V_0H_1 = self.v_H[Li, m]
                V_1_0H = self.v_H[Lj, m].transpose()
                V_0_1 = self.v_C[Li, m]
                V_1_0 = self.v_C[Lj, m].transpose()

                v1 = (np.dot(ci_0H, np.dot(v_0H_0, cj_0)).real-
                      np.dot(ci_0,  np.dot(v_0_0,  cj_0)).real) / 2.0
                
                t1 = np.dot(ci_0H, np.dot(V_0H_1, la.solve(L11, np.dot(V_1_0,  cj_0))))
                t2 = np.dot(ci_0,  np.dot(V_0_1,  la.solve(L11, np.dot(V_1_0H, cj_0H))))
                t3 = np.dot(ci_0,  np.dot(V_0_1,  la.solve(L11, np.dot(V_1_0,  cj_0))))
                t4 = np.dot(ci_0H, np.dot(V_0H_1, la.solve(L11, np.dot(V_1_0H, cj_0H))))
                                                   
                v2 = (t1+t2-t3-t4)/4.0
                res[i, j] = (v1 + v2) / np.sqrt(k)
        return res
        
    def phase_shift_t_matrix(self, ene, Ls, irrep_Dooh):
        t_mat = self.t_matrix(ene, Ls, irrep_Dooh)
        k_mat = t2k(t_mat)
        
     
    def calc_psi1(self, w, mpp):
        L = self.S11[mpp] * (self.ei + w) - self.H11[mpp]
        c1 = la.solve(L, self.mu_phi_1[mpp])
        return c1
