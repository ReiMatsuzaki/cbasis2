#ifndef OP_TEMPLATE_H
#define OP_TEMPLATE_H

#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>

#include "linspace.hpp"

namespace cbasis {

  template<int N, class OpT>
  struct NTermOp {
    typedef OpAdd<OpT, typename NTermOp<N-1, OpT>::type> type;
  };
  template<class OpT>
  struct NTermOp<1, OpT> {
    typedef OpT type;
  };

  template<int N, class OpT>
  struct NLinOp {
    typedef OpAdd<ScalarOpMult<OpT>,
		  typename NLinOp<N-1, OpT>::type> type;
  };
  template<class OpT>
  struct NLinOp<1, OpT> {
    typedef ScalarOpMult<OpT> type;
  };
  
  // ==== deriv operator ====
  template<class Field, class Coord, int N>
  struct OpD :public Op<Field, Coord>{};
  template<class Field, class Coord, int N>
  std::ostream& operator << (std::ostream& os, const OpD<Field, Coord, N>&) {
    os << (char*)"OpD" << N;
    return os;
  }

  typedef OpD<double, double, 1> RD1;
  typedef OpD<std::complex<double>, double, 1> CD1;
  typedef OpD<double, double, 2> RD2;
  typedef OpD<std::complex<double>, double, 2> CD2;

  

  // ==== Rm operator ====
  template<class Field, class Coord>  
  class OpRm :public Op<Field, Coord> {
  private:
    int m_;
  public:
    OpRm(int _m):m_(_m) {}
    int m() const {return m_; }
  };
  template<class Field, class Coord>  
  std::ostream& operator << (std::ostream& os, const OpRm<Field, Coord>& a) {
  os << "OpRm(" << a.m() << ")";
  return os;
}

  typedef OpRm<double, double> RRm;
  typedef OpRm<std::complex<double>, double> CRm;
}
#endif
