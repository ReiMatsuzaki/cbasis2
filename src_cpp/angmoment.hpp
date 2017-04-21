#ifndef ANGMOMENT_HPP
#define ANGMOMENT_HPP

#include <exception>
#include "../utils/typedef.hpp"

namespace cbasis {

  // ==== Exception class ====
  // -- to be removed
  /*
  class ExceptionBadYlm : public std::exception {
  public:
    std::string msg_;
    ExceptionBadYlm(int L, int M, std::string msg);
    ~ExceptionBadYlm() throw();
    virtual const char* what() const throw();
  };
  */

  // ==== Utilities ====
  double cg_coef(int j1, int j2, int m1, int m2, int j3, int m3);
  int lm_index(int L, int M);
  int num_lm_pair(int max_l);
  bool is_lm_pair(int L, int M);

  /*
    -- to removed --
  class Array5Dim {
  private:
    double* xs_;
    int n0_;
    int n1_;
    int n2_;
    int n3_;
    int n4_;    
  public:
    Array5Dim(int n0, int n1, int n2, int n3, int n4) {
      n0_ = n0;
      n1_ = n1;
      n2_ = n2;
      n3_ = n3;
      n4_ = n4;
      xs_ = new double[n0*n1*n2*n3*n4];
    }
    ~Array5Dim() {
      delete xs_;
    }
    double get(int i0, int i1, int i2, int i3, int i4) {
      return xs_[i0 + i1*n0_ + i2*n0_*n1_ + i3*n0_*n1_*n2_ + i4*n0_*n1_*n2_*n3_];
    }
    void set(int i0, int i1, int i2, int i3, int i4, double x) {
      xs_[i0 + i1*n0_ + i2*n0_*n1_ + i3*n0_*n1_*n2_ + i4*n0_*n1_*n2_*n3_] = x;
    }
  };
  */

  // ==== Spherical GTO single expansion ====
  // -- not used now
  double GTOExpansionCoef(int l, int m, int lp, int lppp,
			  int Jp, int Mp, int Jpp, int Mpp);
  // ---- Modified spherical bessel function ----
  // values and derivative of Modified Spherical Bessel function i_n(x)
  //
  // Inputs
  // -------
  // max_n >= 0
  // y_dy_array : new dcomplex[2*max_n+2]
  //
  // Returns
  // --------
  // y_dy_array[n]         : i_n(x)  (n=0,...,max_n)
  // y_dy_array[n+max_n+1] : i'_n(x) (n=0,...,max_n)
  //
  //
  // 
  // Implementation.
  // for small x, I used power series.
  // i_n^(1)(x) = x^n \sum_{k=0}^oo {(0.5x^2)^k} / {k! (2n+2k+1)!!}
  void ModSphericalBessel(dcomplex x, int max_n, dcomplex* y_dy_array);
  // ---- Asscociate Legendre function ----
  // compute Associated Legendre functions {P_{LM}(x)}_{LM}
  // 
  // Inputs
  // -------
  // max_l >=0 
  // res  : new dcomplex[num_lm_pair(max_l)]
  // 
  // Returns 
  // --------
  // res[lm_index(L, M)] gives function value P_{LM}(x).
  void AssociatedLegendre(dcomplex x, int max_l, dcomplex* res);
  void RealSphericalHarmonics(dcomplex theta, dcomplex phi, int max_l, dcomplex* res);  
  
  // -- <Y(Jpp,Mpp) | GTO(zeta; (x,y,z); (l, m)>
  void gto_00_r(dcomplex x, dcomplex y, dcomplex z,		    
		     int Jpp, int Mpp,
		     dcomplex zeta,
		     dcomplex* r, int num_r, dcomplex* work, dcomplex* res);
  void gto_lm_r_center(int l, int m, int Jpp, int Mpp,
		       dcomplex* r, int num_r, dcomplex zeta, dcomplex* res);
  void gto_lm_r(int l, int m,       
		dcomplex x, dcomplex y, dcomplex z,      
		int Jpp, int Mpp,   
		dcomplex zeta,
		int lppp_max,
		dcomplex* r, int num_r,
		dcomplex* work,
		dcomplex* res);

  /*
  class GTOCenterExpansion {
  private:
    int lmax_;
    double* coef;
  public:
    void Init(int lmax) {
      lmax_ = lmax;
    }
  };
  */
  
}

#endif
