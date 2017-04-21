# Solve driven Schrodinger equation with basis function expansion.

from r1gtoint import *
import scipy.linalg as la

class TwoPotDrivEq():
    """
    Solve two potential driven Schrodinger equation with complex basis functions.
    .    (E-T-V0-V1)\psi = s
    where
    .    E : energy
    .    T : kinetic operator
    .    V0 : Long range potential(only pure coulomb with charge Z is supported)
    .    V1 : Short range potential(only STO type potential is supported)
    .    s : source term or driven term (only STO type function is supported)
    """
    def __init__(self, gtos, energy, z, v1, s):
        self.gtos = gtos
        self.gtos.normalize()
        self.gtos_cc = R1GTOs()
        self.gtos_cc.set_conj(self.gtos)
        self.energy = energy
        self.z = z
        self.v1 = v1
        self.s = s

    def amp_int_one(self):
        """
        Compute transition amplitude of wave function by one particle operator 
        integration form. 
        .    A = 2k^(-1) <F, [E-H0]\psi>
        """
        v1 = self.gtos.calc_mat_sto(self.v1)
        (s, t, v0)    = self.gtos.calc_mat_stv(1)
        (sH, tH, v0H) = self.gtos_cc.calc_mat_stv(self.gtos, 1)
        m = self.gtos.calc_vec_sto(self.s)
        
        # solve driven eq
        c0 = la.solve(self.energy*s-self.z*v0-t,    m)
        c1 = la.solve(self.energy*s-self.z*v0-v1-t, m)

        # compute amplitude 
        k = np.sqrt(2.0*self.energy)
        plmx_y0p = np.dot(m, c0)
        j_plmx = np.sqrt(-k*plmx_y0p.imag)

        # compute amplitude (E-H0)
        y0p_esmh0_yp = np.dot(c0.conj(), np.dot(sH*self.energy-tH-v0H, c1))
        y0m_esmh0_yp = np.dot(c0       , np.dot(s *self.energy-t -v0 , c1))
        amp2 = (y0p_esmh0_yp - y0m_esmh0_yp)/(2.0j*j_plmx)
        return 2.0/k*amp2        

    def amp_int_two(self):
        """
        Compute transition amplitude of wave function by two particle operator form.
        .    A = 2k^(-1)<F, s + v1psi> 
        In this function, F is evaluated by TIWP method and CBF method.
        TIWP is based on F is proportional to the solution of H0 driven equation.
        We put driven term as s.
        .    (E-H0)psi0 = s
        .    F = Im[psi0]/Sqrt(k Im<s, psi0>)
        """
        
        v1 = self.gtos.calc_mat_sto(self.v1)
        v1H = self.gtos_cc.calc_mat_sto(self.gtos, self.v1)
        (s, t, v0) = self.gtos.calc_mat_stv(1)
        v0 = self.z * v0
        m = self.gtos.calc_vec_sto(self.s)

        # solve driven eq
        c0 = la.solve(self.energy*s-v0-t,    m)
        c1 = la.solve(self.energy*s-v0-v1-t, m)

        # compute amplitude (V2psi + muphi)
        # plmx = -m (wave packet)
        k = np.sqrt(2.0*self.energy)
        plmx_y0p = np.dot(m, c0)
        y0p_muphi = plmx_y0p.conj()
        y0p_v_yp = np.dot(c0.conj(), np.dot(v1H,c1))
        y0m_v_yp = np.dot(c0,        np.dot(v1, c1))
        imy0p_v_yp = (y0p_v_yp - y0m_v_yp)/(2.0j)
        j_plmx = np.sqrt(-np.pi*plmx_y0p.imag)
        amp = (y0p_muphi.imag + imy0p_v_yp)/j_plmx
        return amp        

    def amp_c_matching(self, r0):
        """
        Compute transition amplitude by matching with Coulomb function.
        A = psi(r0)/H^+(r0)
        """

        # compute matrix/vector
        v1 = self.gtos.calc_mat_sto(self.v1)
        (s, t, v0) = self.gtos.calc_mat_stv(1)
        v0 = self.z * v0
        m = self.gtos.calc_vec_sto(self.s)

        # solve driven eq
        c1 = la.solve(self.energy*s-v0-v1-t, m)

        # matching 
        psi0 = gtos.at_r()
        k = np.sqrt(2.0*self.energy)
        amp = psi0/coulomb.outgoing(-1.0/k, k*r0, 1)

        # return
        return amp

    """
    def amp_im_matching(r0):
        # compute matrix/vector
        v1 = self.gtos.calc_mat_sto(self.v1)
        (s, t, v0) = self.gtos.calc_mat_stv(1)
        m = self.gtos.calc_vec_sto(self.s)

        # solve driven eq
        c1 = la.solve(self.energy*s-v0-v1-t, m)

        # matching
        dlog_impsi = (self.gtos.deriv_at_r(c1, [r0])[0].imag /
                      self.gtos.at_r(      c1, [r0])[0].imag)
        k0 = np.sqrt(2.0*ene)
        (f, g, fp, gp, err) = coulomb.coulomb(-1.0/k0, k0*r0, 1)
        phase = np.arctan(-(f*dlog_impsi - fp)/(g* dlog_impsi-gp))
    """
    def solve_full(self):
        v1 = self.gtos.calc_mat_sto(self.v1)
        (s, t, v0) = self.gtos.calc_mat_stv(1)
        v0 = self.z * v0
        m = self.gtos.calc_vec_sto(self.s)

        # solve driven eq
        c1 = la.solve(self.energy*s-v0-v1-t, m)
        return c1


