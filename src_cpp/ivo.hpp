#ifndef IVO_H
#define IVO_H

#include "symmolint.hpp"
#include "symgroup.hpp"
#include "bmatset.hpp"
#include "b2eint.hpp"

namespace cbasis {

  void CalcCoreHamiltonian(SymGTOs basis_a, Molecule mole, SymGTOs basis_b, BMat *H);
  void CalcFock(const BMat& T, const BMat& V,
		Irrep irr0, const Eigen::VectorXcd& c0,
		B2EInt eri_J, B2EInt eri_K, BMat *H);
  void CalcSTEXHamiltonian(const BMat& T, const BMat& V,
			   Irrep irrep0, const Eigen::VectorXcd& c0,
			   B2EInt eri_J, B2EInt eri_K, BMat *H);
  void CalcSTEXHamiltonian(const BMat& T, const BMat& V,
			   Irrep irrep0, const BMat& C0, int num_oo,
			   B2EInt eri_J, B2EInt eri_K, BMat *H);  
  void CalcSTEXHamiltonian(SymGTOs basis_a, SymGTOs basis_b , SymGTOs basis_0,
			   Molecule mole,
			   Irrep irrep0, const Eigen::VectorXcd& c0,
			   BMat *H);
  void CalcSTEXHamiltonian(SymGTOs basis_a, SymGTOs basis_b , SymGTOs basis_0,
			   Molecule mole,
			   Irrep irrep0, const BMat& C, int num_oo, ERIMethod m,
			   BMat *H);    
  class CalcRPA {
  public:
    // -- intermediate --
    BMat A, B, S;
    BMat ApB, AmB, sqrt_AmB, inv_sqrt_AmB;
    BMat H, U;
    BVec w2;
    // -- final results --
    bool is_ortho;
    BVec w;
    BMat Y, Z;
    // -- routines --
  public:
    CalcRPA();
    CalcRPA(std::string name);
    void CalcH();
    void SetS(const BMat& in_S);
    void CalcH_AO(const BMat& H_IVO_AO, const BMat& S,
		  dcomplex eig0_HF, const BMat& K);
    void CalcH_HF(const BMat& H_IVO_AO, dcomplex eig0_HF,	    
		  const BMat& K,
		  const BVec& eig_HF, const BMat& U_HF);
    void CalcH_IVO(dcomplex eig0_HF,
		   const BVec& eig_IVO, const BMat& U_IVO, 
		   const BMat& K);
    void SolveEigen();
    dcomplex CalcYZNorm(Irrep irr, int I);
    void CalcOneInt(const BVec& z_HF, BVec *z_RPA, bool symmetric);
    void CalcOneInt(const BMat& z_HF, BVec *z_RPA, bool symmetric);
    void CalcOneIntRight(const BMat& x, BMat *y, bool symmetric,
			 dcomplex cy=1.0, dcomplex cz=1.0) const;
  };
}

#endif
