#ifndef D_EXP_FUNC_TEMPLATE_H
#define D_EXP_FUNC_TEMPLATE_H

/**
   derivative of orbital exponents of STO/GTO
 */

#include "op.hpp"

#include "exp_func.hpp"

namespace cbasis {
  
  // N th derivative of normalization constant over normalization constant
  template<class FuncT, int N>
  struct CalcNPrimOverN {
    typedef typename FuncT::Field Field;
    Field operator()(int n,  Field z);
  };

  template<class F>
  struct CalcNPrimOverN<ExpFunc<F, 1>, 1> {
    F operator()(int n, F z) {
      return (F(n) + F(1)/F(2)) / z;
    }
  };

  template<class F>
  struct CalcNPrimOverN<ExpFunc<F, 1>, 2> {
    F operator()(int n, F z) {
      return F(4 * n * n - 1) / (F(4) * z * z);
    }
  };

  template<class F>
  struct CalcNPrimOverN<ExpFunc<F, 2>, 1> {
    F operator()(int n, F z) {
      return (n * F(1)/F(2) + F(1)/F(4)) / z;
    }
  };

  template<class F>
  struct CalcNPrimOverN<ExpFunc<F, 2>, 2> {
    F operator()(int n, F z) {
      return (-F(3)/F(4) + n / F(2)) * (F(1)/F(4) + n/F(2)) / (z*z);
    }
  };

  /*
  template<class F>
  F CalcNOnePrimeOverN(int n, F z, sto_tag) {
    return (F(n) + F(1)/F(2)) / z;
  }
  template<class F>
  F CalcNOnePrimeOverN(int n, F z, gto_tag) {
    return (n * F(1)/F(2) + F(1)/F(4)) / z;
  }
  template<class F>
  F CalcNTwoPrimeOverN(int n, F z, sto_tag) {
    return F(4 * n * n - 1) / (F(4) * z * z);
  }
  template<class F>
  F CalcNTwoPrimeOverN(int n, F z, gto_tag) {
    return (-F(3)/F(4) + n / F(2)) * (F(1)/F(4) + n/F(2)) / (z*z);
  }
  */

  template<class F, int m>
  class NormalExpFunc :public ExpFunc<F,m> {

  public:
    typedef ExpFunc<F,m> Func;
    typedef FuncAdd<Func, Func> FuncDerivOne;
    typedef FuncAdd<Func, FuncAdd<Func, Func> > FuncDerivTwo;

  public:
    void setNormalize() {
      F cc = CNorm(*this);
      this->SetScalarProd(F(1)/cc);
    }
    NormalExpFunc() : ExpFunc<F,m>(1, 1, 1) { setNormalize(); }
    NormalExpFunc(int n, F z) : ExpFunc<F,m>(1, n, z) { setNormalize(); }
    NormalExpFunc(const NormalExpFunc<F,m>& o) : ExpFunc<F,m>(o) {}
    ~NormalExpFunc() {}

    void set_z(F z) { Func::set_z(z); setNormalize(); }
    void set_n(int n) { Func::set_n(n); setNormalize(); }

    FuncDerivOne DerivParamOne() const {

      /*
	N r^n e^(-zr^m) -> N'/N Nr^ne^(-zr^m)  -   N r^(n+m)e^(-zr^m)
       */
      
      // F  cp = CalcNOnePrimeOverN(this->n(), this->z(), 
      //			 typename func_traits<Func>::func_tag());
      F cp = CalcNPrimOverN<Func, 1>()(this->n(), this->z());

      Func a(*this); a.SetScalarProd(cp);
      Func b(*this); b.SetScalarProd(-1); b.SetRmProd(m);

      return func_add_func(a, b);
    }

    FuncDerivTwo DerivParamTwo() const {

      /*
	N r^n e^(-zr^m) -> N'/N Nr^ne^(-zr^m)  -   N r^(n+m)e^(-zr^m)
	                -> N''/N u - 2N' r^m u + r^(2m)u
       */

      //F  cpp = CalcNTwoPrimeOverN(this->n(), this->z(), 
      //typename func_traits<Func>::func_tag());
      //F  cp  = CalcNOnePrimeOverN(this->n(), this->z(), 
      //			  typename func_traits<Func>::func_tag());
      F cp = CalcNPrimOverN<Func, 1>()(this->n(), this->z());
      F cpp = CalcNPrimOverN<Func, 2>()(this->n(), this->z());

      Func u1(*this); u1.SetScalarProd(cpp);
      Func u2(*this); u2.SetScalarProd(-F(2)*cp); u2.SetRmProd(m);
      Func u3(*this); u3.SetRmProd(2*m);

      return func_add_func(u1, func_add_func(u2, u3));
    }

  };

  // ==== STO ====
  typedef NormalExpFunc<double, 1> NormalRSTO;
  typedef NormalExpFunc<std::complex<double>, 1> NormalCSTO;

  // ==== GTO ====
  typedef NormalExpFunc<double, 2> NormalRGTO;
  typedef NormalExpFunc<std::complex<double>, 2> NormalCGTO;
  
  
}



#endif
