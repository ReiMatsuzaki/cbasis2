#include "cfunc.hpp"

dcomplex casin(dcomplex z) {
  dcomplex ii(0, -1);
  return -ii * log(z + sqrt(z*z-1.0));
}
dcomplex cacos(dcomplex z) {
  dcomplex ii(0, -1);
  return -ii * log(z + ii * sqrt(1.0 - z*z));
}
dcomplex catan(dcomplex z) {
  dcomplex ii(0, -1);
  return ii * 0.5 * log((ii + z) / (ii - z));
}
