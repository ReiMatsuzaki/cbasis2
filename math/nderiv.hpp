#include <boost/function.hpp>

namespace cbasis {
  typedef boost::function<dcomplex (dcomplex)> CFunc;
  dcomplex NDerivOne_R1(CFunc f, dcomplex x, dcomplex h);
  dcomplex NDerivOne_C1(CFunc f, dcomplex x, dcomplex h);
  dcomplex NDerivOne_R3(CFunc f, dcomplex x, dcomplex h);

  dcomplex NDerivTwo_R1(CFunc f, dcomplex x, dcomplex h);
  dcomplex NDerivTwo_C1(CFunc f, dcomplex x, dcomplex h);
  dcomplex NDerivTwo_R3(CFunc f, dcomplex x, dcomplex h);
}
