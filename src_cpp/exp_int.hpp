#ifndef EXP_INT_H
#define EXP_INT_H

namespace cbasis {

  template<class F> F STO_Int(F z, int n);
  template<class F> F GTO_Int(F z, int n);
  template<class F> F STO_GTO_Int(F as, F ag, int n);
  template<class F> F CutSTO_Int(F z, int n, double r0);

}

#endif
