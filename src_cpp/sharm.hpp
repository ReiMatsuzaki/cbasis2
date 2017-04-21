#ifndef SHARM_TEMPLATE_H
#define SHARM_TEMPLATE_H

#include "linspace.hpp"

namespace cbasis {

  // ==== Sheraical Harmonics ====
  template<class RadFunc>
  class SHarm :public Func<typename RadFunc::Field, R3d> {
  public:
    // ---- type ----
    typedef Func<typename RadFunc::Field, R3d> base;
    typedef typename base::Field F;

  private:
    // ---- member field ----
    int L_;
    int M_;
    RadFunc rad_func_;
    
  public:
    // ---- Constructors ----
    SHarm(int _L, int _M, const RadFunc& f): L_(_L), M_(_M), rad_func_(f) {}
    
    // ---- Accessor ----
    int L() const { return L_; }
    int M() const { return M_; }
    const RadFunc rad_func() const { return rad_func_; }

    // ---- Method ----
    //F at(R3d x) const;
    
  };

  // ==== Spherical tensor operator ====
  template<class RadOp>
  class TensorOp: public Op<typename RadOp::Field, R3d> {
  public:
    typedef Func<typename RadOp::Field, R3d> base;
    typedef typename base::Field F;
  public:
    // ---- member field ----
    int K_;
    int Q_;
    RadOp op_;
    
    // ---- Construcors ----
    TensorOp(int _K, int _Q, const RadOp& _op): K_(_K), Q_(_Q), op_(_op) {}
  };

  // ==== CIP ====
  template<class RadFunc1, class RadFunc2>
  typename RadFunc1::Field CIP_impl_prim(const SHarm<RadFunc1>& a,
					 const One&,
					 const SHarm<RadFunc2>& b) {
    typedef typename RadFunc1::Field F;
    if(a.L() != b.L() || a.M() != b.M())
      return F(0);

    return CIP(a, One(), b);
  }
}
#endif
