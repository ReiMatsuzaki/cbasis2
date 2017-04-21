#ifndef ONE_INT_H
#define ONE_INT_H

#include <Eigen/Core>
#include "../utils/typedef.hpp"
#include "mol_func.hpp"
#include "symmolint.hpp"

namespace cbasis {

  // ==== Definition ====
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<Reduction>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<Reduction>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 1> A1dc;
  typedef MultArray<dcomplex, 2> A2dc;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;


  // ==== Slow routines ====
  dcomplex SMatEle(CartGTO& a, CartGTO& b);
  dcomplex TMatEle(CartGTO& a, CartGTO& b);
  dcomplex VMatEle(CartGTO& a, Eigen::Vector3cd at, CartGTO& b);
  dcomplex XMatEle(CartGTO& a, CartGTO& b);
  dcomplex YMatEle(CartGTO& a, CartGTO& b);
  dcomplex ZMatEle(CartGTO& a, CartGTO& b);  
  dcomplex DXMatEle(CartGTO& a, CartGTO& b);
  dcomplex DYMatEle(CartGTO& a, CartGTO& b);
  dcomplex DZMatEle(CartGTO& a, CartGTO& b);
  dcomplex PWVecEle(const Eigen::Vector3cd& k, CartGTO& a);
  dcomplex PWXVecEle(const Eigen::Vector3cd& k, CartGTO& a);
  dcomplex PWYVecEle(const Eigen::Vector3cd& k, CartGTO& a);
  dcomplex PWZVecEle(const Eigen::Vector3cd& k, CartGTO& a);

  // ==== SymGTOs ====
  BMatSet CalcMat(SymGTOs a, SymGTOs b, bool calc_coulomb);
  BMatSet CalcMat_Complex(SymGTOs g, bool calc_coulomb);
  BMatSet CalcMat_Hermite(SymGTOs g, bool calc_coulomb);

  // ==== SymGTOs(new) ====
  void InitBVec(SymGTOs a, BVec *ptr_bvec);
  void InitBVec(SymGTOs a, Irrep irrep, BVec *ptr_bvec);
  void InitBVec(SymGTOs a, const std::vector<Irrep>& irrep_list, BVec *ptr_bvec);
  void InitBMat(SymGTOs a, Irrep krrep, SymGTOs b, BMat *ptr_mat);
  void CalcSTMat(SymGTOs a, SymGTOs b, BMat *S, BMat *T);
  void CalcSTVMat(SymGTOs a, SymGTOs b, BMat *S, BMat *T, BMat *V);
  void CalcSMat(SymGTOs a, SymGTOs b, BMat *S);
  void CalcVMat(SymGTOs a, Molecule mole, SymGTOs b, BMat *V);
  void CalcDipMat(SymGTOs a, SymGTOs b,
		  BMat* X, BMat* Y, BMat* Z, BMat* DX, BMat* DY, BMat* DZ);
  void CalcDipMat(SymGTOs a, SymGTOs b, BMat* L, BMat *V);
  void CalcPWVec(SymGTOs a, const Eigen::Vector3cd& k,
		 BVec *S, BVec *X, BVec *Y, BVec *Z);
  
  
}

#endif
