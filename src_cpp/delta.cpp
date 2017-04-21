#include "delta.hpp"

namespace cbasis {

  // ==== Constructors ====
  template<class F, class C>
  DiracDelta<F, C>::DiracDelta() : r0_(0.0) {}
  template<class F, class C> 
  DiracDelta<F, C>::DiracDelta(C _r0) :r0_(_r0) {}
  template<class F, class C> 
  DiracDelta<F, C>::DiracDelta(const DiracDelta& o) : r0_(o.r0()) {}
    

  // ==== Method ====
  /*
  template<class F, class C> F DiracDelta<F, C>::at(C x) const {
    if( std::abs(x-this->r0()) < this->h_) 
      return F(1);
    else
      return F(0);
  }
  */

  // ==== Externals ====
  template<class F, class C>
  std::ostream& operator << (std::ostream& os, const DiracDelta<F, C>& a) {
    os << "Delta(r0=" << a.r0() << ")";
    return os;
  }

  // ==== Explicit Decralation ====
  template class DiracDelta<double, double>;
  template class DiracDelta<std::complex<double>, double>;
}
