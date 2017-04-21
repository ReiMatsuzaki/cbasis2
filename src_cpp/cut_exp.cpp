#include "cut_exp.hpp"
#include <limits>

namespace cbasis {

  // ==== Constructors ====
  template<class F, int m> 
  CutExpFunc<F,m>::CutExpFunc(): func(), 
				 r0_(0) {}
  
  template<class F, int m> 
  CutExpFunc<F,m>::CutExpFunc(F _c, int _n, F _z, double r0): func(_c, _n, _z),
							      r0_(r0) {}

  template<class F, int m> 
  CutExpFunc<F,m>::CutExpFunc(const CutExpFunc<F, m>& o): func(o.func),
							  r0_(o.r0()) {}


  // ==== Method ====
  template<class F, int m> 
  F CutExpFunc<F,m>::at(F x) const {
    
    /*
      std::numeric_limits<F>::epsilon() returns machine epsilon(the difference between 1 and the least value greater than 1 that is representable)
     */

    if(std::abs(x) <= this->r0_ * std::abs(F(1)+F(10)*std::numeric_limits<F>::epsilon()))
      return this->func.at(x);
    else
      return F(0);
  }

  template<class F, int m>
  std::string CutExpFunc<F,m>::str() const {
    std::stringstream ss;
    if(m == 1)
      ss << "CutSTO[";
    else
      ss << "CutGTO[";
    ss << c() << "," << n() << "," << z() <<  "," << r0() << "]";
    return ss.str();    
  }

  template<class F, int m> 
  void CutExpFunc<F,m>::SetComplexConjugate() {
    
    this->func.SetComplexConjugate();

  }

  template<class F, int m> 
  void CutExpFunc<F,m>::SetDerivParam() {

    this->func.SetDerivParam();

  }

  template<class F, int m> 
  void CutExpFunc<F,m>::SetScalarProd(F c) {

    this->func.SetScalarProd(c);

  }

  template<class F, int m> 
  void CutExpFunc<F,m>::SetRmProd(int n) {

    this->SetRmProd(n);

  }

  template<class F, int m>
  typename CutExpFunc<F,m>::FuncDerivOne CutExpFunc<F,m>::DerivParamOne() const {
    FuncDerivOne s(*this);
    s.SetDerivParam();
    return s;
  }
  template<class F, int m>
  typename CutExpFunc<F,m>::FuncDerivTwo CutExpFunc<F,m>::DerivParamTwo() const {
    FuncDerivOne s(*this);
    s.SetDerivParam();
    s.SetDerivParam();
    return s;
  }

  // ==== Externals ====
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, const CutExpFunc<F,m>& a) {

    os << "Cut(" << a.r0() << ")" << a.get_func();
    return os;

  }

  
  // ==== Expllicit Declation ====
  template class CutExpFunc<double, 1>;
  template class CutExpFunc<std::complex<double>, 1>;
  template class CutExpFunc<double, 2>;
  template class CutExpFunc<std::complex<double>, 2>;

  template std::ostream& operator <<<std::complex<double>,1> 
  (std::ostream& os, const CutCSTO& a);
  template std::ostream& operator <<<double,1> 
  (std::ostream& os, const CutRSTO& a);
  template std::ostream& operator <<<std::complex<double>,2> 
  (std::ostream& os, const CutCGTO& a);
  template std::ostream& operator <<<double,2>
  (std::ostream& os, const CutRGTO& a);
  

}
