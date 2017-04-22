#ifndef IVO_H
#define IVO_H

#include "symgroup.hpp"
#include "bmatset.hpp"
#include "b2eint.hpp"

namespace cbasis {
  
  void CalcFock(const BMat& T, const BMat& V,
		Irrep irr0, const Eigen::VectorXcd& c0,
		B2EInt eri_J, B2EInt eri_K, BMat *H);
  void CalcSTEXHamiltonian(const BMat& T, const BMat& V,
			   Irrep irrep0, const Eigen::VectorXcd& c0,
			   B2EInt eri_J, B2EInt eri_K, BMat *H);
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
    CalcRPA(std::string name);
    void CalcH();
    void CalcH_HF(const BMat& H_IVO_AO, dcomplex eig0_HF,	    
		  const BMat& K,
		  const BVec& eig_HF, const BMat& U_HF);
    void CalcH_IVO(dcomplex eig0_HF,
		   const BVec& eig_IVO, const BMat& U_IVO, 
		   const BMat& K);
    void SolveEigen();
    dcomplex CalcYZNorm(Irrep irr, int I);
    void CalcOneInt(const BVec& z_HF, BVec *z_RPA, bool symmetric);
  };
}

#endif
