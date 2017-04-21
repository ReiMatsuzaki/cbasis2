#ifndef INT_EXP_H
#define INT_EXP_H

#include <vector>
#include "../utils/typedef.hpp"

namespace cbasis {

  // -- slow --
  dcomplex STOInt_Rplus(int n, dcomplex a);
  dcomplex GTOInt_Rplus(int n, dcomplex a);  
  dcomplex STO_GTOInt_Rplus(int n, dcomplex a, dcomplex b);
  dcomplex STO_GTOInt_Rplus_2(int n, dcomplex a, dcomplex b);
  
  dcomplex GTOInt_R(int n, dcomplex a);
  dcomplex STO_GTOInt_R(int n, dcomplex a, dcomplex b);  

  // -- rapid --
  template<class Scalar>
  void STOInt_Rplus_array(int n, Scalar a, std::vector<Scalar> *res);
  template<class Scalar>
  void GTOInt_Rplus_array(int n, Scalar a, std::vector<Scalar> *res);
  template<class Scalar>
  void STO_GTOInt_Rplus_array(int n, Scalar a, Scalar b, std::vector<Scalar> *res);
   
}

#endif
