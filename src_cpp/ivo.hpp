#ifndef IVO_H
#define IVO_H

#include "symgroup.hpp"
#include "mo.hpp"
#include "bmatset.hpp"
#include "b2eint.hpp"

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
  class CalcRPA {
  public:
    // -- intermediate --
    BMat A, B;
    BMat ApB, AmB, sqrt_AmB, inv_sqrt_AmB;
    BMat H, U;
    BVec w2;
    // -- final results --
    BVec w;
    BMat Y, Z;
    // -- routines --
  public:
    CalcRPA(std::string name) {
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
    void CalcH() {
      Copy(A, AmB); AmB.Add(-1.0, B);
      Copy(A, ApB); ApB.Add(+1.0, B);
      BMatSqrt(AmB, sqrt_AmB);
      BMatInvSqrt(AmB, inv_sqrt_AmB);
      Multi3(sqrt_AmB, ApB, sqrt_AmB, H);
    }
    void CalcH_HF(const BMat& H_IVO_AO, dcomplex eig0_HF,		    
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
    void CalcH_IVO(dcomplex eig0_HF,
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
    void SolveEigen() {

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
    dcomplex CalcYZNorm(Irrep irr, int I) {
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
  };
  /*
  void CalcRPAHamiltonian_AO(const BMat& H_IVO_AO, const BMat& S, const BMat& K,
			     dcomplex eig0_HF, BMat *H) {
    BMat A("A");
    Copy(H_IVO_AO, A);
    A.Add(-eig0_HF, S);

    CalcRPAHamiltonian(A, K, H);
  }
  */
}












#endif
