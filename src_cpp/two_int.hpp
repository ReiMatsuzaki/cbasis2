#ifndef TWO_INT_H
#define TWO_INT_H

#include "mol_func.hpp"
#include "symmolint.hpp"

namespace cbasis {

  // ==== Slow routines ====
  dcomplex ERIEle(CartGTO& a, CartGTO& b, CartGTO& c, CartGTO& d);

  // ==== aa ====
  void coef_R_eri_switch(dcomplex zarg,
			 dcomplex wPx, dcomplex wPy, dcomplex wPz,
			 dcomplex wPpx,dcomplex wPpy,dcomplex wPpz,
			 int max_n, dcomplex *Fjs, dcomplex mult_coef,
			 MultArray<dcomplex, 3>& res, ERIMethod method);
			 

  // ==== SymGTOs ====
  B2EInt CalcERI_Complex(SymGTOs i, ERIMethod m);
  B2EInt CalcERI_Hermite(SymGTOs i, ERIMethod m);
  B2EInt CalcERI(SymGTOs i, SymGTOs j, SymGTOs k, SymGTOs l, ERIMethod method);  
	       
}

#endif
