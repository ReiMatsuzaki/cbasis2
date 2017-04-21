#ifndef EXP_FUNC_TEMPLATE_H
#define EXP_FUNC_TEMPLATE_H

#include <iostream>
#include <complex>
#include "linspace.hpp"
#include "op.hpp"

namespace cbasis {

  // ==== STO or GTO ====
  template<class F, int m>
  class ExpFunc :public Func<F, double> {
  public:
    // ---- type ----
    typedef F Field;
    enum EExpPower { exp_power=m };
    typedef ExpFunc<F,m> FuncDerivOne;
    typedef ExpFunc<F,m> FuncDerivTwo;

  private:    
    // ---- Member Field ----
    F c_;
    int n_;
    F z_;

  public:
    // ---- Constructors ----
    ExpFunc();
    ExpFunc(F c, int n, F z);
    ExpFunc(const ExpFunc<F, m>& o);
    template<class F2>
    ExpFunc(const ExpFunc<F2,m>& o): c_(o.c()), n_(o.n()), z_(o.z()) {}
    ~ExpFunc();

    // ---- Accessor ----
    F c() const { return this->c_; }
    int n() const {return this->n_; }
    F z() const { return this->z_; }
    
    void set_c(F c) { this->c_ = c; }
    void set_n(int n) { this->n_ = n;}
    void set_z(F z) { this->z_ = z;}
    
    // ---- Method ----
    F at(F x) const;
    std::string str() const;
    void SetComplexConjugate();
    void SetDerivParam();
    void SetScalarProd(F c);
    void SetRmProd(int n);

    FuncDerivOne DerivParamOne();
    FuncDerivTwo DerivParamTwo();
  };

  /*
  // ==== Operator ====
  // ---- Rm ----
  template<class F, int m>
  struct op_func_res<OpRm<F, double>, ExpFunc<F, m> > {
    typedef ExpFunc<F, m>  type;
  };

  template<class F, int m>
  typename op_func_res<OpRm<F, double>, ExpFunc<F, m> >::type 
  OP_prim(const OpRm<F, double>& o, const ExpFunc<F, m>& a) {
    ExpFunc<F, m> res(a);
    res.n += o.m();
    return res;
  }

  // ---- D1 ----
  template<class F, int m>
  struct op_func_res<OpD<F, double, 1>, ExpFunc<F, m> > {
    typedef typename NTermFunc<2, ExpFunc<F, m> >::type  type;
  };

  template<class F, int m>
  typename op_func_res<OpD<F, double, 1>, ExpFunc<F, m> >::type 
  OP_prim(const OpD<F, double, 1>&, const ExpFunc<F, m>& a) {
    typedef ExpFunc<F, m> Func;
    Func f(a); f.set_n(f.n()-1);
    Func g(a); g.set_c(-a.z()*F(m)*f.c()); g.set_n(g.n()-1);
    return func_add_func(f, g);
  }
	  

  // ---- D2 ----
  template<class F, int m>
  struct op_func_res<OpD<F, double, 2>, ExpFunc<F, m> > {
    typedef typename NTermFunc<3, ExpFunc<F, m> >::type  type;
  };

  template<class F, int m>
  NTermFunc<3, ExpFunc<F, m> > OP_prim(const OpD<F, double, 2>&,
				       const ExpFunc<F, m>& a) {
    int n = a.n();
    F   c = a.c();
    F   z = a.z();
    ExpFunc<F, m> s(a); s.set_c(c* F(n*n-n));            s.set_n(n-2);
    ExpFunc<F, m> t(a); t.set_c(c* (-z*F(2*n*m+m*m-m))); t.set_n(n+m-2);
    ExpFunc<F, m> u(a); u.set_c(c* (z*z*F(m*m)));        u.set_n(n+2*m-2);
    return func_add_func(s, func_add_func(t, u));
  }

*/
  // ==== External ====
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, const ExpFunc<F,m>& a);

  typedef ExpFunc<double, 1> RSTO;
  typedef ExpFunc<std::complex<double>, 1> CSTO;
  
  typedef ExpFunc<double, 2> RGTO;
  typedef ExpFunc<std::complex<double>, 2> CGTO;
}

#endif
