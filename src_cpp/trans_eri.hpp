#ifndef TRANS_ERI_H
#define TRANS_ERI_H

#include "mo.hpp"
#include "b2eint.hpp"

namespace cbasis {

  void TransformERI_Slow(IB2EInt* ao, MO mo, IB2EInt* res);
  void TransformERI(IB2EInt* ao, MO mo, IB2EInt* res);
  void CalcJK_MO(IB2EInt* eri_ao, BMat& C, int A0, int a0, int B0, int b0,
		 dcomplex cJ, dcomplex cK, BMat* bmat);
}

#endif
