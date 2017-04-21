#include <iostream>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include "../utils/typedef.hpp"
#include "nderiv.hpp"

using namespace boost::lambda;

namespace cbasis {
  dcomplex NDerivOne_R1(CFunc f, dcomplex x, dcomplex h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
  }
  dcomplex NDerivOne_C1(CFunc f, dcomplex x, dcomplex h) {
    dcomplex i(0,1);
    return (f(x+h) - f(x-h) - i*f(x+i*h) + i*f(x-i*h)) / (4.0*h);
  }
  dcomplex NDerivOne_R3(CFunc f, dcomplex x, dcomplex h) {
    dcomplex m1 = (f(x+h) - f(x-h)) / (2.0*h);
    dcomplex m2 = (f(x+2.0*h) - f(x-2.0*h)) / (4.0*h);
    dcomplex m3 = (f(x+3.0*h) - f(x-3.0*h)) / (6.0*h);
    return 3.0/2.0*m1 -3.0/5.0*m2 + 1.0/10.0*m3;
  }

  dcomplex NDerivTwo_R1(CFunc f, dcomplex x, dcomplex h) {
    //CFunc ff = bind(NDerivOne_R1, f, _1, h);
    //    return NDerivOne_R1(ff , x, h);
    return NDerivOne_R1(bind(NDerivOne_R1, f, _1, h) , x, h);
  }
  dcomplex NDerivTwo_C1(CFunc f, dcomplex x, dcomplex h) {
    return NDerivOne_R1(bind(NDerivOne_C1, f, _1, h) , x, h);
    /*
      f(x+h) = f(x) + h f'(x)  + hh/2 f''(x) + hhh/6 f'''(x) + 
      f(x+h) = f(x) - h f'(x)  + hh/2 f''(x) - hhh/6 f'''(x) + 
      f(x+ih)= f(x) + ih f'(x) - hh/2 f''(x) - ihhh/6 f'''(x) + 
      f(x-ih)= f(x) - ih f'(x) - hh/2 f''(x) + ihhh/6 f'''(x) + 
     */
  }
  dcomplex NDerivTwo_R3(CFunc f, dcomplex x, dcomplex h) {
    return NDerivOne_R3(bind(NDerivOne_C1, f, _1, h) , x, h);
  }
}
