#include <gtest/gtest.h>
#include "../utils/eigen_plus.hpp"
#include "../utils/gtest_plus.hpp"
#include "../utils/typedef.hpp"
#include "one_int.hpp"
#include "two_int.hpp"

#include "ivo.hpp"

using namespace std;
using namespace cbasis;
using namespace Eigen;
class TestIVO : public ::testing::Test {
public:  
  SymGTOs us_0;
  SymGTOs us_1;

  MO mo_0;
  CalcRPA RPA_HF;
  CalcRPA RPA_IVO;
  BVec de_HF;
  BMat U_HF;
  BVec de_IVO;
  BMat U_IVO;
  
  void PrepareBasis() {
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

    
    // -- basis for excited state --
    VectorXcd zs = zeta_s;
    us_1 = NewSymGTOs(mole);
    us_1->NewSub("He").SolidSH_M(1, 0, zs);
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
 
  cout << de_HF << endl;
  cout << de_IVO << endl;
  cout << RPA_HF.w << endl;
  cout << RPA_IVO.w << endl;

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

  /**
     AO basis dipole matrix element
     <u1_i | z | u0_j>
   */
  BMat x,y,z,dx,dy,dz;
  CalcDipMat(us_1, us_0, &x, &y, &z, &dx, &dy, &dz);

  BMat z_HF, z_IVO, z_RPA_HF, z_RPA_IVO;
  BMatCtAD(U_HF,      mo_0->C, z,     &z_HF);
  BMatCtAD(U_IVO,     mo_0->C, z,     &z_IVO);
  BMatCtAD(RPA_HF.U,  U_HF,    z_HF,  &z_RPA_HF);
  BMatCtAD(RPA_IVO.U, U_IVO,   z_IVO, &z_RPA_IVO);

  Irrep irr1 = us_1->sym_group()->irrep_z();
  for(int I = 0; I < 10; I++) {
    cout << z_HF(irr1, 0)(I, 0) << ", "
	 << z_IVO(irr1, 0)(I, 0) << ", "
	 << z_RPA_HF(irr1, 0)(I, 0) << ", "
	 << z_RPA_IVO(irr1, 0)(I, 0) 
	 << endl;
  }
}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}


