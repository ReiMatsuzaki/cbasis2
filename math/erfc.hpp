#ifndef ERFC_H
#define ERFC_H

namespace {
}

namespace erfc_mori {

  // ================ Utils =========================
  // compilation raise error
  template<typename T> struct raise_error;
  
  // norm is less than
  template<class F> bool NormIsLessThan(F x, double y);
  
  // eps
  double machine_eps();
  
  // erfc start num term
  template<class F> int erfc_start_num_term();

  // add Eh term or not
  template<class F>  bool erfc_add_Eh_q(F x, F h);
  
  // n term erfc 
  template<class F>
  int ErfcAtN(F x, unsigned int n, F& y, F& yn );

  // erfc calculation data
  struct ErfcCalcData {
    int  num_term;
    bool convergence;
  };

  // erfc
  template<class F>
  int Erfc(F x, F & y, ErfcCalcData& data);

  template<class F>
  F erfc(F x);

  // Exp(z^2) Erfc(z)  with n term
  template<class F>
  int Exp2ErfcAtN(F x, unsigned int n, F& y, F& yn);

  // Exp(z^2) Erfc(z) or erfcx (scaling complement error function)
  template<class F>
  int Exp2Erfc(F x, F& y, ErfcCalcData& data);

  template<class F>
  F erfcx(F x);
}  
#endif
