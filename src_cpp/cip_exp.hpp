#ifndef CIP_EXP_TEMPLATE_H
#define CIP_EXP_TEMPLATE_H

#include "exp_func.hpp"
#include "cut_exp.hpp"
#include "normal_exp.hpp"
#include "exp_int.hpp"

namespace cbasis {
  
  // ==== static dispatch ====
  // ---- ExpFunc ----
  template<class F, int m1, int m2>
  struct ExpInt {
    F operator()(F za, F zb, int n);
  };
  template<class F>
  struct ExpInt<F, 1, 1> {
    F operator()(F za, F zb, int n) {
      return STO_Int(za+zb, n);
    }
  };
  template<class F>
  struct ExpInt<F, 1, 2> {
    F operator()(F za, F zb, int n) {
      return STO_GTO_Int(za, zb, n);
    }
  };
  template<class F>
  struct ExpInt<F, 2, 1> {
    F operator()(F za, F zb, int n) {
      return STO_GTO_Int(zb, za, n);
    }
  };
  template<class F>
  struct ExpInt<F, 2, 2> {
    F operator()(F za, F zb, int n) {
      return GTO_Int(zb+za, n);
    }
  };

  // ---- CutExp ----
  template<class F, int m1, int m2>
  struct CutExpInt;
  template<class F>
  struct CutExpInt<F, 1, 1> {
    F operator()(F za, F zb, int n, double r0) {
      return CutSTO_Int(za+zb, n, r0);
    }
  };


  // ==== set CIP_impl_prim ====
  // ---- ExpFunc, ExpFunc ----
  template<class F, int m1, int m2>
  F CIP_impl_prim(const ExpFunc<F, m1>& a, 
		  const One&, 
		  const ExpFunc<F, m2>& b) {
    return CIP_impl_prim(a, OpRm<F, double>(0), b);
  }
  
  template<class F, int m1, int m2>
  F CIP_impl_prim(const ExpFunc<F, m1>& a, 
		  const OpRm<F, double>& o, 
		  const ExpFunc<F, m2>& b) {
    return ExpInt<F, m1, m2>()(a.z(), b.z(), o.m()+a.n()+b.n()) * a.c() * b.c();
  }

  template<class F, int m1, int m2>
  F CIP_impl_prim(const ExpFunc<F, m1>& a, 
		  const OpD<F, double, 1>&, 
		  const ExpFunc<F, m2>& b) {
    int n = b.n();
    int m = b.exp_power;
    F z = b.z();
    typedef OpRm<F, double> Rm;
    return F(n) * CIP_impl_prim(a, Rm(-1), b) - z*F(m)*CIP_impl_prim(a, Rm(m-1), b);
  }

  template<class F, int m1, int m2>
  F CIP_impl_prim(const ExpFunc<F, m1>& a, 
		  const OpD<F, double, 2>&, 
		  const ExpFunc<F, m2>& b) {
    int n = b.n();
    int m = b.exp_power;
    F z = b.z();
    typedef OpRm<F, double> Rm;
    F cumsum(0);
    
    if(n!=1)
      cumsum += F(n*n-n)        * CIP_impl_prim(a, Rm(-2), b);
    cumsum += -z*F(2*n*m+m*m-m) * CIP_impl_prim(a, Rm(m-2), b);
    cumsum += F(m*m)*z*z        * CIP_impl_prim(a, Rm(2*m-2), b);
    return cumsum;
  }


  // ---- CutExp, CutExp ----
  template<class F, int m1, int m2>
  F CIP_impl_prim(const CutExpFunc<F, m1>& a, 
		  const One&, 
		  const CutExpFunc<F, m2>& b) {
    return CIP_impl_prim(a, OpRm<F, double>(0), b);
  }  
  template<class F, int m1, int m2>
  F CIP_impl_prim(const CutExpFunc<F, m1>& a,
		  const OpRm<F, double>& o,
		  const CutExpFunc<F, m2>& b) {
    double r0 = a.r0() < b.r0() ? a.r0() : b.r0();
    return CutExpInt<F, 1, 1>()(a.z(), 
				b.z(), 
				a.n() + b.n() + o.m(), 
				r0) * a.c() * b.c();
  }

  template<class F, int m1, int m2>
  F CIP_impl_prim(const CutExpFunc<F, m1>& a, 
		  const OpD<F, double, 1>&, 
		  const CutExpFunc<F, m2>& b) {
    int n = b.n();
    int m = b.exp_power;
    F z = b.z();
    typedef OpRm<F, double> Rm;
    return F(n) * CIP_impl_prim(a, Rm(-1), b) - z*F(m)*CIP_impl_prim(a, Rm(m-1), b);
  }

  template<class F, int m1, int m2>
  F CIP_impl_prim(const CutExpFunc<F, m1>& a, 
		  const OpD<F, double, 2>&, 
		  const CutExpFunc<F, m2>& b) {
    int n = b.n();
    int m = b.exp_power;
    F z = b.z();
    typedef OpRm<F, double> Rm;
    F cumsum(0);
    
    if(n!=1)
      cumsum += F(n*n-n) * CIP_impl_prim(a, Rm(-2), b);
    cumsum += -z*F(2*n*m+m*m-m) * CIP_impl_prim(a, Rm(m-2), b);
    cumsum += F(m*m)*z*z        * CIP_impl_prim(a, Rm(2*m-2), b);
    return cumsum;
  }

  // ---- CutExp, ExpFunc ----
  template<class F, int m1, int m2, class OpT>
  F CIP_impl_prim(const CutExpFunc<F, m1>& a,
		  const OpT& o, 
		  const ExpFunc<F, m2>& b) {
    CutExpFunc<F, m2>  cut_b(b.c(), b.n(), b.z(), a.r0());
    return CIP_impl_prim(a, o, cut_b);
  }

  template<class F, int m1, int m2, class OpT>
  F CIP_impl_prim(const ExpFunc<F, m1>& a,
		  const OpT& o, 
		  const CutExpFunc<F, m2>& b) {
    CutExpFunc<F, m1>  cut_a(a.c(), a.n(), a.z(), b.r0());
    return CIP_impl_prim(cut_a, o, b);
  }

  // ---- NormalExp ----
  template<class F, int m1, int m2, class OpT>
  F CIP_impl_prim(const NormalExpFunc<F, m1>& a,
		  const OpT& o,
		  const NormalExpFunc<F, m2>& b) {
    ExpFunc<F, m1> aa = a;
    ExpFunc<F, m2> bb = b;
    return CIP_impl_prim(aa, o, bb);
  }
  /*
  template<class F, int m1, class OpT, class FuncA>
  F CIP_impl_prim(const FuncA& a,
		  const OpT& o,
		  const NormalExpFunc<F, m1>& b) {
    return CIP_impl_prim<FuncA, OpT, ExpFunc<F, m1> >(a, o, b);
  }
  */
}
#endif
