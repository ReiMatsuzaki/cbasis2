#include "math_utils.hpp"
#include<complex>

typedef std::complex<double> CD;

namespace cbasis {

  template<class F> F ConjugateIfPossible(F x) { return x; }
  template<> double ConjugateIfPossible(double x) { return x; }
  template<> CD ConjugateIfPossible<CD>(CD x) { return std::conj(x); }
}
