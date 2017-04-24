#include <gtest/gtest.h>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include "../utils/eigen_plus.hpp"
#include "../utils/gtest_plus.hpp"
#include "../utils/typedef.hpp"
#include "one_int.hpp"
#include "two_int.hpp"
#include "mo.hpp"
#include "ivo.hpp"

using namespace std;
using namespace cbasis;
using namespace Eigen;
using boost::format;
using namespace boost::assign;


class TestPracticeTwoPot : public ::testing::Test {
public:
  // basis
  SymGTOs us_g, us_0, us_1;
  ERIMethod method;

  // ground state
  MO mo_g;
  const static Irrep irr0=0;
  BMat U_g;
  VectorXcd c0_g;
  BVec E_g;

  // first order wave function
  BMat U_HF_1;
  CalcRPA rpa; // contain Y,Z,wi

  // 0th order wave function
  BVec eig0;
  BMat U_0;
    
  void PrepareBasis() {
    
    // -- symmetry --
    SymmetryGroup sym = SymmetryGroup_D2h();

    // -- target molecule --
    Molecule mole = NewMolecule(sym);
    mole->Add(NewAtom("He", 2.0)->Add(0,0,0));
    Molecule mole0 = NewMolecule(sym);
    mole0->Add(NewAtom("dum", 1.0)->Add(0,0,0));
       
    // -- basis for ground state --
    us_g = NewSymGTOs(mole);
    VectorXcd zeta_s(10);
    //    zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723, 488.8941, 3293.694;
    zeta_s << 0.107951, 0.552610, 3.522261, 30.17990, 488.8941;  
    us_g->NewSub("He").SolidSH_M(0,0,zeta_s);
    us_g->SetUp();
    
    // -- basis for excited state --
    VectorXcd zs(8);
    //zs << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053, 30.17990, 108.7723;
    zs << 0.240920, 1.352436, 9.789053, 108.7723;
    us_1 = NewSymGTOs(mole);
    us_1->NewSub("He").SolidSH_M(1, 0, zs);
    us_1->NewSub("He").SolidSH_M(1, -1, zs);
    us_1->NewSub("He").SolidSH_M(1, +1, zs);
    us_1->SetUp();

    // -- basis for zeroth --
    us_0 = NewSymGTOs(mole0);
    us_0->NewSub("dum").SolidSH_M(1, 0, zs);
    us_0->NewSub("dum").SolidSH_M(1, -1, zs);
    us_0->NewSub("dum").SolidSH_M(1, +1, zs);
    us_0->SetUp();

    // -- other --
    method.set_symmetry(1); method.set_coef_R_memo(1);

  }
  void Solve_ground() {
    // ==== HF for ground state ====
    SCFOptions opts; opts.max_iter = 10;
    
    bool conv;
    mo_g = CalcRHF(us_g, 2, method, opts, &conv);
    if(not conv) {
      cout << "failed to conversion" << endl;
    }
    U_g = mo_g->C;
    c0_g = U_g(irr0, irr0).col(0);
    E_g = mo_g->eigs;
  }
  void Solve_first () {
    // ==== Eigen values for first psi ====    
    Molecule mole = us_1->molecule();
    // -- AO --
    BMat S1, T1, V1;
    CalcSTVMat(us_1, us_1, &S1, &T1, &V1);
    B2EInt eri_J_11 = CalcERI(us_1, us_1, us_0, us_0, method);
    B2EInt eri_K_11 = CalcERI(us_1, us_0, us_0, us_1, method);
    BMat K1; Copy(S1, K1); K1.SetZero();
    AddK(eri_K_11, c0_g, 0, 1.0, K1);
    // -- HF --
    BMat H_HF_1;
    CalcFock(T1, V1, 0, c0_g, eri_J_11, eri_K_11, &H_HF_1);
    BVec eig_HF_1;
    BMatEigenSolve(H_HF_1, S1, &U_HF_1, &eig_HF_1);
    // -- RPA eigen system --
    BMat H_IVO_1;
    CalcSTEXHamiltonian(T1, V1, irr0, c0_g, eri_J_11, eri_K_11, &H_IVO_1);
    CalcRPA rpa;
    rpa.CalcH_HF(H_IVO_1, E_g(irr0)(0), K1, eig_HF_1, U_HF_1);
    rpa.SolveEigen();
    // -- dipole moment
    SymmetryGroup sym = us_1->sym_group();    
    Irrep x = sym->irrep_x();
    Irrep y = sym->irrep_y();
    Irrep z = sym->irrep_z();
    BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i, s1_ao, sD1_ao;
    CalcDipMat(us_1, us_g, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
    s1_ao(x,0) = X1i(x,0); sD1_ao(x,0) = DX1i(x,0);
    s1_ao(y,0) = Y1i(y,0); sD1_ao(y,0) = DY1i(y,0);
    s1_ao(z,0) = Z1i(z,0); sD1_ao(z,0) = DZ1i(z,0);
    BMat r1g_HF, d1g_HF;
    BMatCtAD(U_HF_1, U_g, s1_ao, &r1g_HF);
    BMatCtAD(U_HF_1, U_g, sD1_ao, &d1g_HF);
    BVec r1_RPA, d1_RPA;
    rpa.CalcOneInt(r1g_HF,  &r1_RPA, true);
    rpa.CalcOneInt(d1g_HF, &d1_RPA, false);
    

    
  }
  void Solve_zero() {
    // -- AO --
    BMat S, T, H;
    CalcSTVMat(us_0, us_0, &S, &T, &H);
    H.Add(1.0, T);
    BMatEigenSolve(H, S, &U_0, &eig0);
  }
  void Solve() {
    this->Solve_ground();
    this->Solve_first();
  }
  TestPracticeTwoPot() {
    this->PrepareBasis();
    this->Solve();
  }
};
TEST_F(TestPracticeTwoPot, first) {

}

class TestIVO : public ::testing::Test {
public:  
  SymGTOs us_0;
  SymGTOs us_1;

  MO mo_0;
  VectorXcd c0;
  dcomplex E0;
  CalcRPA RPA_HF;
  CalcRPA RPA_IVO;
  BVec de_HF;
  BMat U_HF;
  BVec de_IVO;
  BMat U_IVO;
  
  void PrepareBasis() {
    U_HF.set_name("U_HF");
    U_IVO.set_name("U_IVO");
    
    // -- symmetry --
    SymmetryGroup sym = SymmetryGroup_D2h();

    // -- target molecule --
    Molecule mole = NewMolecule(sym);
    mole->Add(NewAtom("He", 2.0)->Add(0,0,0));
       
    // -- basis for ground state --
    us_0 = NewSymGTOs(mole);
    VectorXcd zeta_s(10);
    zeta_s << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053,
      30.17990, 108.7723, 488.8941, 3293.694;  
    us_0->NewSub("He").SolidSH_M(0,0,zeta_s);
    us_0->SetUp();

    // -- HF for ground state --
    SCFOptions opts; opts.max_iter = 10;
    ERIMethod method; method.set_symmetry(1); method.set_coef_R_memo(1);
    bool conv;
    mo_0 = CalcRHF(us_0, 2, method, opts, &conv);
    if(not conv) {
      cout << "failed to conversion" << endl;
    }
    c0 = mo_0->C(0, 0).col(0);
    E0 = mo_0->eigs(0)(0);
    
    // -- basis for excited state --
    VectorXcd zs(8);
    zs << 0.107951, 0.240920, 0.552610, 1.352436, 3.522261, 9.789053,
      30.17990, 108.7723;
    us_1 = NewSymGTOs(mole);
    us_1->NewSub("He").SolidSH_M(1, 0, zs);
    us_1->NewSub("He").SolidSH_M(1, -1, zs);
    us_1->NewSub("He").SolidSH_M(1, +1, zs);
    us_1->SetUp();

  }
  void PrepareRPA() {
    // -- computa AO matrix --
    BMat S, T, V; 
    CalcSTVMat(us_1, us_1, &S, &T, &V);
    
    ERIMethod m;
    B2EInt eri_J = CalcERI(us_1, us_1, us_0, us_0, m);
    B2EInt eri_K = CalcERI(us_1, us_0, us_0, us_1, m);
  
    // -- extract initial orbital --
    VectorXcd c0 = mo_0->C(0, 0).col(0);
    dcomplex eig0_HF = mo_0->eigs(0)(0);
    BMat K; Copy(S, K); K.SetZero();
    AddK(eri_K, c0, 0, 1.0, K);
    
    // -- Hartree Fock --
    BMat H_HF;
    CalcFock(T, V, 0, c0, eri_J, eri_K, &H_HF);
    BVec eig_HF("eig_HF");
    BMatEigenSolve(H_HF, S, &U_HF, &eig_HF);
    de_HF = eig_HF; de_HF.set_name("de_HF");
    de_HF.Shift(-eig0_HF);
    
    // -- IVO --
    BMat H_IVO_AO;
    CalcSTEXHamiltonian(T, V, 0, c0, eri_J, eri_K, &H_IVO_AO);
    BVec eig_IVO("eig_IVO");  
    BMatEigenSolve(H_IVO_AO, S, &U_IVO, &eig_IVO);  
    de_IVO = eig_IVO; de_IVO.set_name("de_IVO");
    de_IVO.Shift(-eig0_HF);
    
    // -- RPA on HF --    
    RPA_HF.CalcH_HF(H_IVO_AO, eig0_HF, K, eig_HF, U_HF);
    RPA_HF.SolveEigen();
    
    // -- RPA on IVO --
    RPA_IVO.CalcH_IVO(eig0_HF, eig_IVO, U_IVO, K);
    RPA_IVO.SolveEigen();
  }
  TestIVO(): RPA_HF("RPA_HF"), RPA_IVO("RPA_IVO") {
    this->PrepareBasis();
    this->PrepareRPA();
  }
};
TEST_F(TestIVO, Delta_E) {
  Irrep irr1 = us_1->sym_group()->irrep_z();
  int n = us_1->size_basis_isym(irr1);  
  cout << "E0 = " << mo_0->eigs(0)(0).real() << endl;
  cout << "Delta E" << endl;
  for(int I = 0; I < n; I++) {
    cout << format("%10.5f, %10.5f, %10.5f, %10.5f\n")
      % de_HF(irr1)(I).real()
      % de_IVO(irr1)(I).real()
      % RPA_HF.w(irr1)(I).real()
      % RPA_IVO.w(irr1)(I).real();
  }
}
TEST_F(TestIVO, Norm) {
 
  // -- check RPA norm --
  SymmetryGroup sym = us_1->sym_group();
  Irrep irr1 = sym->irrep_z();
  int n = us_1->size_basis_isym(irr1);
  for(int i = 0; i < n; i++) {
    EXPECT_C_EQ(1.0, RPA_HF.CalcYZNorm(irr1, i));
    EXPECT_C_EQ(1.0, RPA_IVO.CalcYZNorm(irr1, i));
  }
  
}
TEST_F(TestIVO, Moment) {

  // -- AO basis dipole matrix element --
  BMat x,y,z("z"),dx,dy,dz;
  CalcDipMat(us_1, us_0, &x, &y, &z, &dx, &dy, &dz);
  BVec zz;
  Irrep irr1 = us_1->sym_group()->irrep_z();
  zz(irr1) = z(irr1,0) * c0;  

  // -- many electron --
  BVec z_HF("z_HF");
  BVecCtx(U_HF, zz, &z_HF);
  BVec z_IVO("z_IVO");
  BVecCtx(U_IVO, zz, &z_IVO);
  BVec z_RPA_HF;
  RPA_HF.CalcOneInt(z_HF, &z_RPA_HF, true);
  BVec z_RPA_IVO;
  RPA_IVO.CalcOneInt(z_IVO, &z_RPA_IVO, true);

  // -- print out --  
  int n = us_1->size_basis_isym(irr1);
  for(int I = 0; I < n; I++) {
    cout << format("%10.5f, %10.5f, %10.5f, %10.5f\n")
      % z_HF(irr1)(I).real()
      % z_IVO(irr1)(I).real()
      % z_RPA_HF(irr1)(I).real()
      % z_RPA_IVO(irr1)(I).real();
  }

  // -- check sum rule --
  dcomplex acc_HF(0), acc_IVO(0), acc_RPA_HF(0), acc_RPA_IVO(0);
  for(int I = 0; I < n; I++) {
    dcomplex ele_z = z_HF(irr1)(I);
    acc_HF += 2.0/3.0 * de_HF(irr1)(I) * ele_z * ele_z;

    ele_z = z_IVO(irr1)(I);
    acc_IVO += 2.0/3.0 * de_IVO(irr1)(I) * ele_z * ele_z;
    
    ele_z = z_RPA_HF(irr1)(I);
    acc_RPA_HF += 2.0/3.0 * RPA_HF.w(irr1)(I) * ele_z * ele_z;

    ele_z = z_RPA_IVO(irr1)(I);
    acc_RPA_IVO += 2.0/3.0 * RPA_IVO.w(irr1)(I) * ele_z * ele_z;
    
  }
  cout << "sum of oscillator strength" << endl;
  cout << format("%10.5f, %10.5f, %10.5f, %10.5f\n")
    % acc_HF.real()
    % acc_IVO.real()
    % acc_RPA_HF.real()
    % acc_RPA_IVO.real();

}
TEST_F(TestIVO, Moment2) {
  
  // -- AO basis dipole matrix element --
  BMat x,y,z("z"),dx,dy,dz;
  CalcDipMat(us_1, us_0, &x, &y, &z, &dx, &dy, &dz);
  BVec vecl, vecv;
  Irrep irrx = us_1->sym_group()->irrep_x();
  Irrep irry = us_1->sym_group()->irrep_y();
  Irrep irrz = us_1->sym_group()->irrep_z();
  vecl(irrx) = x(irrx,0) * c0;
  vecl(irry) = y(irry,0) * c0;
  vecl(irrz) = z(irrz,0) * c0;
  vecv(irrx) = dx(irrx,0) * c0;
  vecv(irry) = dy(irry,0) * c0;
  vecv(irrz) = dz(irrz,0) * c0;  

  // -- many electron --
  BVec vecl_HF, vecv_HF;
  BVecCtx(U_HF, vecl, &vecl_HF);
  BVecCtx(U_HF, vecv, &vecv_HF);
  BVec vecl_IVO, vecv_IVO;
  BVecCtx(U_IVO, vecl, &vecl_IVO);
  BVecCtx(U_IVO, vecv, &vecv_IVO);
  BVec vecl_RPA_HF, vecv_RPA_HF;
  RPA_HF.CalcOneInt(vecl_HF, &vecl_RPA_HF, true);
  RPA_HF.CalcOneInt(vecv_HF, &vecv_RPA_HF, false);
  BVec vecl_RPA_IVO, vecv_RPA_IVO;
  RPA_IVO.CalcOneInt(vecl_IVO, &vecl_RPA_IVO, true);
  RPA_IVO.CalcOneInt(vecv_IVO, &vecv_RPA_IVO, false);

  // -- check sum rule --
  vector<Irrep> irreps; irreps += irrx, irry, irrz;
  dcomplex acc_HF(0), acc_IVO(0), acc_RPA_HF(0), acc_RPA_IVO(0);
  dcomplex accv_HF(0), accv_IVO(0), accv_RPA_HF(0), accv_RPA_IVO(0);
  BOOST_FOREACH(Irrep irr, irreps) {
    int n = us_1->size_basis_isym(irr);
    dcomplex d, w;
    for(int I = 0; I < n; I++) {
      d = vecl_HF(irr)(I); w = de_HF(irr)(I);      
      acc_HF += 2.0/3.0*w*d*d;
      d = vecl_IVO(irr)(I); w = de_IVO(irr)(I);      
      acc_IVO += 2.0/3.0*w*d*d;
      d = vecl_RPA_HF(irr)(I); w = RPA_HF.w(irr)(I);      
      acc_RPA_HF += 2.0/3.0*w*d*d;
      d = vecl_RPA_IVO(irr)(I); w = RPA_IVO.w(irr)(I);      
      acc_RPA_IVO += 2.0/3.0*w*d*d;

      d = vecv_HF(irr)(I); w = de_HF(irr)(I);      
      accv_HF += 2.0/3.0/w*d*d;
      d = vecv_IVO(irr)(I); w = de_IVO(irr)(I);      
      accv_IVO += 2.0/3.0/w*d*d;
      d = vecv_RPA_HF(irr)(I); w = RPA_HF.w(irr)(I);      
      accv_RPA_HF += 2.0/3.0/w*d*d;
      d = vecv_RPA_IVO(irr)(I); w = RPA_IVO.w(irr)(I);      
      accv_RPA_IVO += 2.0/3.0/w*d*d;
    }
  }

  // -- print out --
  cout << "sum of oscillator strength" << endl;
  cout << format("%10.5f, %10.5f, %10.5f, %10.5f\n")
    % acc_HF.real()
    % acc_IVO.real()
    % acc_RPA_HF.real()
    % acc_RPA_IVO.real();
  cout << format("%10.5f, %10.5f, %10.5f, %10.5f\n")
    % accv_HF.real()
    % accv_IVO.real()
    % accv_RPA_HF.real()
    % accv_RPA_IVO.real();
}


int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}


