#ifndef CUT_EXP_TEMPLATE_H
#define CUT_EXP_TEMPLATE_H

#include "exp_func.hpp"

/**
   Cutted STO/GTO
*/

namespace cbasis {

  template<class Field, int m>
  class CutExpFunc :public Func<Field, double> {
    
  public:
    // ---- type ----
    typedef Field F;
    enum EExpPower {exp_power=m};
    typedef CutExpFunc<Field,  m> type;
    typedef type FuncDerivOne;
    typedef type FuncDerivTwo;

  private:
    // ---- Member Field ----
    ExpFunc<Field, m> func;
    double r0_;          // R-matrix radius

  public:
    // ---- Constructors ----
    CutExpFunc();
    CutExpFunc(F _c, int _n, F _z, double r0);
    CutExpFunc(const type& o);
    template<class U> CutExpFunc(const CutExpFunc<U,m>& o):
      func(o.func), r0_(o.r0_)  {}

    // ---- Accessor ----
    F c() const { return func.c(); }
    int n() const { return func.n(); }
    F z() const { return func.z(); }
    double r0() const { return this->r0_; }
    
    void set_c(F c) { func.set_c(c); }
    void set_n(int n) { func.set_n(n);}
    void set_z(F z) { func.set_z(z);}
    void set_r0(double r0) { this->r0_ = r0;}

    const ExpFunc<F,m>& get_func() const { return func; }
    
    // ---- Method ----
    F at(F x) const;
    std::string str() const;
    void SetComplexConjugate();
    void SetDerivParam();
    void SetScalarProd(F c);
    void SetRmProd(int n);

    FuncDerivOne DerivParamOne() const;
    FuncDerivTwo DerivParamTwo() const;
  };


  // ==== External ====
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, const CutExpFunc<F,m>& a);

  // ==== Operator ====
  // ---- ExpFunc->CutExpFunc ----
  template<class Field, class Coord>
  class Cut:public Op<Field, Coord> {
  private:
    Coord r0_;
  public:
    Cut(Coord _r0): r0_(_r0){}
    Coord r0() const { return r0_; }
  };

  // ==== STO ====
  typedef CutExpFunc<double, 1> CutRSTO;
  typedef CutExpFunc<std::complex<double>, 1> CutCSTO;

  // ==== GTO ====
  typedef CutExpFunc<double, 2> CutRGTO;
  typedef CutExpFunc<std::complex<double>, 2> CutCGTO;
}
#endif
