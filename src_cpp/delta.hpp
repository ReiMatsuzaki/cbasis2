#ifndef DELTA_TEMPLATE_H
#define DELTA_TEMPLATE_H

#include "linspace.hpp"
#include <complex>

/**
   Represent Dirac's Delta function.
   δ(r-r_0)
   whose width is very small value h.
*/

namespace cbasis {

  // ==== Delta function ====
  template<class Field, class Coord>
  class DiracDelta :public Func<Field, Coord>{
    
  public:
  private:
    // ---- Field Member ----
    double r0_; // location of delta function

  public:
    // ---- Constructors ----
    DiracDelta();
    DiracDelta(Coord _r0);
    DiracDelta(const DiracDelta<Field, Coord>& o);

    // ---- Accessor ----
    double r0() const { return this->r0_; }
    void set_r0(double r0) { this->r0_ = r0; }

    // ---- Method ----
    //    Field at(Coord x) const;
  };

  // ==== External ====
  template<class Field, class Coord>
  std::ostream& operator << (std::ostream& os, const DiracDelta<Field, Coord>& a);

  // ==== typedef ====
  typedef DiracDelta<double,double> RDelta;
  typedef DiracDelta<std::complex<double>,double > CDelta;

  // ==== CIP ====
  template<class F, class FuncA>
  F CIP_impl_prim(const DiracDelta<F, double>& d,
		  const One&,
		  const FuncA& a) {
    return a.at(d.r0());
  }
  template<class F, class FuncA>
  F CIP_impl_prim(const FuncA& a,
		  const One&,
		  const DiracDelta<F, double>& d
		  ) {
    return a.at(d.r0());
  }

}

#endif
