#include <iostream>
#include <stdexcept>
#include "../utils/macros.hpp"
#include "../utils/fact.hpp"
#include "int_exp.hpp"
#include "erfc.hpp"

using namespace std;
using namespace erfc_mori;

namespace cbasis {

  dcomplex STOInt_Rplus(int n, dcomplex a) {

    if(n < 0) {
      string msg; SUB_LOCATION(msg);
      msg += "\nn must be bigger than 0";
      throw runtime_error(msg);
    }

    if(n == 0)
      return 1.0/a;
    
    return STOInt_Rplus(n-1, a) * (1.0*n) / a;
    
    //    return DFactorial(n)/pow(a, n+1)

  }
  dcomplex GTOInt_Rplus(int n, dcomplex a) {

    if(n < 0) {
      std::string msg; SUB_LOCATION(msg);
      msg += "\nn must be bigger than 0";
      throw std::runtime_error(msg);
    }

    if(n == 0)
      return sqrt(M_PI)/(2.0*sqrt(a));
    if(n == 1)
      return 0.5/a;

    return dcomplex(n-1)/(2.0*a) * GTOInt_Rplus(n-2, a);

  }
  dcomplex STO_GTOInt_Rplus(int n, dcomplex a, dcomplex b) {

    if(n < 0) {
      std::string msg; SUB_LOCATION(msg);
      msg += "\nn must be bigger than 0";
      throw std::runtime_error(msg);
    }

    if(n == 0) {
      return sqrt(M_PI) * erfcx(a/(2.0*sqrt(b))) /(2.0*sqrt(b));
    }
    if(n == 1) {
      return (-sqrt(M_PI)*a*erfcx(a/(2.0*sqrt(b)))+2.0*sqrt(b)) / (4.0*pow(b,1.5));
	      
    }

    return ( (n-1.0)/(2.0*b) * STO_GTOInt_Rplus(n-2,a,b)
	     - a/(2.0*b) * STO_GTOInt_Rplus(n-1,a,b)
	     );	     
  }
  dcomplex GTOInt_R(int n, dcomplex a) {

    if(n % 2 == 1) {
      return 0.0;
    }

    return GTOInt_Rplus(n, a) * 2.0;
  }
  dcomplex STO_GTOInt_R(int n, dcomplex a, dcomplex b) {

    /*
      Int(r^n Exp[-ar -brr], (-oo,+oo)) = 
      Int(r^n Exp[-b(r+a/2b)^2 + a^2/4b], (-oo,+oo)) =
      Exp[a^2/4b] Int((r-a/2b)^n Exp[-br^2], (-oo,+oo)) =
      Exp[a^2/4b] sum_m nCm (a/2b)^(n-m) Int(r^m Exp[-br^2], (-oo, +oo)) = 
    */

    dcomplex res = 0.0;
    for(int m = 0; m <= n; m++) {
      res += dcomplex(Combination(n, m)) * pow(a/(2.0*b), n-m) * GTOInt_R(m, b);
    }
    res *= exp(a*a/(4.0*b));
    return res;
  }

  template<class Scalar>  
  void STOInt_Rplus_array(int maxn, Scalar a, vector<Scalar> *res) {
    /**
       compute { Int_0^oo r^n exp[-ar] dr   | n=0,...,maxn }

       Inputs
       ------
       maxn   : maximum of integer n  (maxn>0)
       a      : orbital exponent

       Returns
       -------
       res[n] : n th result 
     */
    
    if(maxn < 0) {
      string msg; SUB_LOCATION(msg);
      msg += "\n" + msg + "maxn must be positive";
      throw runtime_error(msg);
    }
    if((int)res->size() != maxn+1) {
      res->resize(maxn+1);
    }

    vector<Scalar>& vec = *res;
    vec[0] = 1.0/a;

    for(int n = 1; n <= maxn; n++) {
      vec[n] = vec[n-1] * (Scalar(n)/a);
    }
    
  }  
  template<class Scalar>  
  void GTOInt_Rplus_array(int maxn, Scalar a, vector<Scalar> *res) {
    /**
       compute { Int_0^oo r^n exp[-ar^2] dr   | n=0,...,maxn }

       Inputs
       ------
       maxn   : maximum of integer n  (maxn>0)
       a      : orbital exponent

       Returns
       -------
       res[n] : n th result 
     */
    
    if(maxn < 0) {
      string msg; SUB_LOCATION(msg);
      msg += "\n" + msg + "maxn must be positive";
      throw runtime_error(msg);
    }
    if((int)res->size() != maxn+1) {
      res->resize(maxn+1);
    }

    vector<Scalar>& vec = *res;
    
    vec[0] = sqrt(M_PI)/(2.0*sqrt(a));

    if(maxn >= 1)
      vec[1] = 0.5/a;

    for(int n = 2; n <= maxn; n++) {
      vec[n] = Scalar(n-1)/(2.0*a) * vec[n-2];
    }    
  }
  template<class Scalar>
  void STO_GTOInt_Rplus_array(int maxn, Scalar a, Scalar b, vector<Scalar> *res) {
    
    if(maxn < 0) {
      string msg; SUB_LOCATION(msg);
      msg += "\n" + msg + "maxn must be positive";
      throw runtime_error(msg);
    }
    if((int)res->size() != maxn+1) {
      res->resize(maxn+1);
    }
    
    Scalar val_erfcx =  erfcx(a/(2.0*sqrt(b)));
    vector<Scalar>& vec = *res;
    vec[0] = sqrt(M_PI) * val_erfcx/(2.0*sqrt(b));

    if(maxn >= 1) {
      vec[1] = (-sqrt(M_PI)*a*val_erfcx +2.0*sqrt(b)) / (4.0*pow(b,1.5));
		
    }
    for(int n = 2; n <= maxn; n++) {
      vec[n] = (n-1.0)/(2.0*b) * vec[n-2] - a/(2.0*b) * vec[n-1];
    }
  }
  template void STOInt_Rplus_array<dcomplex>(int, dcomplex, vector<dcomplex>*);
  template void GTOInt_Rplus_array<dcomplex>(int, dcomplex, vector<dcomplex>*);
  template void STO_GTOInt_Rplus_array<dcomplex>(int, dcomplex, dcomplex, vector<dcomplex>*);
}
