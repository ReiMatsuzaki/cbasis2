#include "ivo.hpp"
#include "../utils/eigen_plus.hpp"
#include "mo.hpp"

using namespace Eigen;

namespace cbasis {
  
  void CalcFock(const BMat& T, const BMat& V,
		Irrep irr0, const Eigen::VectorXcd& c0,
		B2EInt eri_J, B2EInt eri_K, BMat *H) {
    if(not T.is_same_structure(V)) {
      THROW_ERROR("size mismatch (T and V)");
    }
    if(not H->is_same_structure(T)) {
      Copy(T, *H);
    }

    H->SetZero();
    H->Add(1.0, T);
    H->Add(1.0, V);
    AddJ(eri_J, c0, irr0, +2.0, *H);
    AddK(eri_K, c0, irr0, -1.0, *H);
  }
  void CalcSTEXHamiltonian(const BMat& T, const BMat& V,
			   Irrep irrep0, const Eigen::VectorXcd& c0,
			   B2EInt eri_J, B2EInt eri_K, BMat *H) {

    if(not T.is_same_structure(V)) {
      THROW_ERROR("size mismatch (T and V)");
    }
    if(not H->is_same_structure(T)) {
      Copy(T, *H);      
    }

    H->SetZero();
    H->Add(1.0, T);
    H->Add(1.0, V);
    AddJ(eri_J, c0, irrep0, 1.0, *H);
    AddK(eri_K, c0, irrep0, 1.0, *H);

  }
  CalcRPA::CalcRPA(std::string name) {
    A.set_name(name+"_A");
    B.set_name(name+"_B");
    ApB.set_name(name+"_A-B");
    AmB.set_name(name+"_A+B");
    sqrt_AmB.set_name(name+"_(A-B)^(1/2)");
    inv_sqrt_AmB.set_name(name+"_(A-B)^(-1/2)");
    H.set_name(name+"_H");
    U.set_name(name+"_U");
    w2.set_name(name+"_w2");
    w.set_name(name+"_w");
    Y.set_name(name+"_Y");
    Z.set_name(name+"_Z");
  }
  void CalcRPA::CalcH() {
    Copy(A, AmB); AmB.Add(-1.0, B);
    Copy(A, ApB); ApB.Add(+1.0, B);
    BMatSqrt(AmB, sqrt_AmB);
    BMatInvSqrt(AmB, inv_sqrt_AmB);
    Multi3(sqrt_AmB, ApB, sqrt_AmB, H);
  }
  void CalcRPA::CalcH_HF(const BMat& H_IVO_AO, dcomplex eig0_HF,
			 const BMat& K,
			 const BVec& eig_HF, const BMat& U_HF) {
    // -- A Mat --
    BMatCtAC(U_HF, H_IVO_AO, &A); // change to HF basis
    A.Shift(-eig0_HF); 
    
    // -- B Mat --
    BMatCtAC(U_HF, K, &B);  // B is now HF basis exchange matrix
    
    // -- H --
    this->CalcH();
  }
  void CalcRPA::CalcH_IVO(dcomplex eig0_HF,
			  const BVec& eig_IVO, const BMat& U_IVO, 
			  const BMat& K) {
    // -- A Mat --
    BMatDiag(eig_IVO, &A);
    A.Shift(-eig0_HF); 
    
    // -- B Mat --
    BMatCtAC(U_IVO, K, &B);
    
    // -- Hamiltonian --
    this->CalcH();
  }
  void CalcRPA::SolveEigen() {
    
    // -- solve eigen value problem of RPA --
    BMatEigenSolve(H, &U, &w2);
    BVecSqrt(w2, &w);
    
    // -- compute Y+Z and Y-Z --
    BMat YpZ, YmZ;
    Multi(sqrt_AmB,     U, YpZ);
    Multi(inv_sqrt_AmB, U, YmZ);
    for(BMat::iterator it = U.begin(); it != U.end(); ++it) {
      BMat::Key k = it->first;
      Irrep irr = it->first.first;
      Eigen::MatrixXcd& m(it->second);
      int n = m.cols();
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < n; j++) {
	  dcomplex sqrt_w = sqrt(sqrt(w2(irr)(j)));
	  YpZ[k](i, j) = 1.0/sqrt_w * YpZ[k](i, j);
	  YmZ[k](i, j) = sqrt_w * YmZ[k](i, j);
	  }
      }
    }
    
    // -- compute Y and Z --
    Copy(YpZ, Y); Y.Add(+1.0, YmZ); Y.Scale(0.5);
    Copy(YpZ, Z); Z.Add(-1.0, YmZ); Z.Scale(0.5);
    
  }
  dcomplex CalcRPA::CalcYZNorm(Irrep irr, int I) {
    dcomplex acc_y2(0), acc_z2(0);
    for(BMat::const_iterator it = Y.begin(); it != Y.end(); ++it) {
      BMat::Key k = it->first;
      if(k.first != irr or k.second != irr)
	continue;
      const Eigen::MatrixXcd& m(it->second);
      acc_y2 += TDot(m.col(I), m.col(I));
    }
    for(BMat::const_iterator it = Z.begin(); it != Z.end(); ++it) {
      BMat::Key k = it->first;
      if(k.first != irr or k.second != irr)
	continue;
      const Eigen::MatrixXcd& m(it->second);
      acc_z2 += TDot(m.col(I), m.col(I));
    }
    return acc_y2 - acc_z2;
  }
  void CalcRPA::CalcOneInt(const BVec& z_HF, BVec *z_RPA, bool symmetric) {
    /**
       z_HF  : HF basis one particle operator matrix.
       z_HF(i) = <i|z|0> in HF basis
       <I_RPA|o|0> = sqrt(2) sum_i (Y_iI <i|z|0> + Z_iI <0|o|i>)
    */
    double sign = (symmetric ? 1.0 : -1.0);
    typedef BVec::const_iterator It;
    for(It it = z_HF.begin(); it != z_HF.end(); ++it) {
      Irrep irr = it->first;
      const VectorXcd& z = it->second;      
      const MatrixXcd& ZZ(Z(irr, irr));
      const MatrixXcd& YY(Y(irr, irr));
      int nI = ZZ.cols();
      int ni = ZZ.rows();
      (*z_RPA)(irr) = VectorXcd::Zero(nI);
      VectorXcd& zz = (*z_RPA)(irr);
      for(int I = 0; I < nI; I++) {
	for(int i = 0; i < ni; i++) {
	  zz(I) += YY(i,I)*z(i) + sign*ZZ(i, I)*z(i);
	}
      }
    }      
  }  
}
