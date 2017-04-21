#include <math.h>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_sf_coupling.h>

#include "../utils/typedef.hpp"
#include "../utils/macros.hpp"
#include "../utils/fact.hpp"

#include "cfunc.hpp"
#include "angmoment.hpp"

using namespace std;

namespace cbasis {

  // ==== Exception class ====
  // -- to be removed
  /*
  ExceptionBadYlm::ExceptionBadYlm(int L, int M, std::string msg) :std::exception() {
      std::stringstream ss;
      ss << "\nUnphysical (L, M) pair. (L, M) = (" << L << ", " << M << ")";
      msg_ = msg + ss.str();
  }
  ExceptionBadYlm::~ExceptionBadYlm() throw() {}
  const char* ExceptionBadYlm::what() const throw() {
    return msg_.c_str();
  }
  */

  // ==== Utilities ====
  double cg_coef(int j1, int j2, int m1, int m2, int j3, int m3) {
    return ipow(-1, j1-j2+m3) * sqrt(2*j3+1.0) *
      gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, -2*m3);
  }
  int lm_index(int L, int M) {
    return L * L + (M+L);
  }
  int num_lm_pair(int max_l) {

    return lm_index(max_l+1, -max_l-1);

  }
  bool is_lm_pair(int L, int M) {
    return 0<=L && iabs(M) <= L;
  }
  
  // -- not used now
  double GTOExpansionCoef(int l, int m, int lp, int lppp,
			  int Jp, int Mp, int Jpp, int Mpp) {

    // See
    // APPLICATION OF THE SCHWINGER VARIATIONAL PRINCIPLE TO ELECTRON-MOLECULE COLLISIONS AND MOLECULAR PHOTOIONIZATION
    // R.Lucchese, K.Takatsuka, V.McKoy
    // Phys. Rep. 131, (1986), 147

    double cumsum(0);

    for(int mp = -lp; mp <= lp; mp++)
      for(int mpp = -lp; mpp <= lp; mpp++)
	for(int J = l-lp; J <= l+lp; J++)
	  for(int M = -J; M <= J; M++)
	    for(int mppp = -lppp; mppp <= lppp; mppp++) {
	      double t1 = ipow(-1, l+m-(lp+mp)) * (2.0*lppp+1) / DFactorial(l-lp);
	      double t2 = sqrt((4.0 * M_PI * (2*l+1)*DFactorial(l-mp)*DFactorial(l+mp)) /
			       ((2.0*Jp+1)*(2.0*Jpp+1)*DFactorial(lp-mp)*DFactorial(lp+mp)));
	      double t3 = (cg_coef(lp,l,mpp,-m,J,M) *
			   cg_coef(lp,l,mp,-mp,J,0) * 
			   cg_coef(J,lppp,0,0,Jp,0) *
			   cg_coef(J,lppp,M,mppp,Jp,Mp) * 
			   cg_coef(lp,lppp,0,0,Jpp,0) *
			   cg_coef(lp,lppp,mpp,mppp,Jpp,Mpp));
		
	      cumsum += t1*t2*t3;
	    }
    return cumsum;
  }

  // ==== Special functions ====
  void ModSphericalBessel(dcomplex x, int max_n, dcomplex* y_dy_array) {

    if(max_n < 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": max_n must be 0 or positive.";
      throw runtime_error(msg);
    }


    dcomplex* ys = y_dy_array;
    dcomplex* dys= &y_dy_array[max_n + 1];

    double eps(0.02);

    // -- compute { ys[n] | n=0,...,max_n+1 }
    if(abs(x) < eps) {

      // >>> small x
      dcomplex xx(0.5*x*x);
      dcomplex xn(1);
      
      for(int n = 0; n <= max_n+1; n++) {

	int kf(20);
	dcomplex c(1.0/dcomplex(2*n+1));
	dcomplex cumsum(0.0);
	dcomplex xxk(1.0);
	bool convq(false);
	double machine_eps(pow(10.0, -15.0));
	for(int k = 0; k < kf; k++) {
	  if(abs(c) < machine_eps) {
	    convq = true;
	    break;
	  }
	  c = xxk / dcomplex(Factorial(k) * DoubleFactorial(2*n+2*k+1));
	  xxk *= xx;
	  cumsum += c;	  
	}
	/*
	c = 1.0 / dcomplex(2*n+1);
	for(int k = 0; k < kf; k++) {
	  if(abs(c) < machine_eps) {
	    convq = true;
	    break;
	  }
	  cumsum += c;
	  c *= xx / dcomplex((k+1) * (2*n+2*k+3));
	}
	*/
	if(!convq) {
	  string msg; SUB_LOCATION(msg);
	  ostringstream oss;
	  oss << msg << ": " << endl;
	  oss << "x=" << x << endl;
	  oss << "n=" << n << endl;
	  oss << "c=" << c << endl;
	  oss << ": not converged";
	  throw runtime_error(oss.str());
	}
	ys[n] = xn * cumsum;
	xn *= x;

      }
      // <<< small x

    } else {

      // >>> large x
      ys[0] = sinh(x) / x;
      ys[1] = (x*cosh(x)-sinh(x)) / (x*x);
      if(max_n > 1)
	for(int n = 1; n <= max_n; n++)  {
	  ys[n+1] = ys[n-1] - (2.0*n+1.0)/x*ys[n];
	}
      // <<< large x

    }

    // -- compute {dys[n] | n=0,...,max_n}
    for(int n = max_n; n >= 0; n--) {
      dys[n] = dcomplex(n+1)*ys[n+1];
      if(n != 0)
	dys[n] += dcomplex(n) * ys[n-1];
      dys[n] *= 1.0 / dcomplex(2*n+1);
    }

  }
  void AssociatedLegendre(dcomplex x, int max_l, dcomplex* res) {

    if(max_l < 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": max_l must be 0 or positive.";
      throw runtime_error(msg);
    }

    double eps(0.000001);
    if(abs(x-1.0) < eps|| abs(x+1.0) < eps) {
      // -- (L-M+1) P(L+1,M) = (2L+1) x P(LM) - (L+M) P(L-1,M)

      res[0] = 1.0;

      for(int L = 0; L < max_l; L++) {
	for(int M = -L-1; M <= L+1; M++) {
	  if(M == 0) {
	    dcomplex c(1.0/double(L+1));
	    res[lm_index(L+1, 0)] = c * dcomplex(2*L+1) * x * res[lm_index(L,0)];
	    if(L != 0)
	      res[lm_index(L+1, 0)] -= c*dcomplex(L+M) * res[lm_index(L-1,0)];
	  } else {
	    res[lm_index(L+1, M)] = 0.0;
	  }
	}
      }

    } else {
    
      res[0] = 1.0;
      for(int L = 0; L < max_l; L++) {
	dcomplex val = -(2.0*L+1.0) * sqrt(1.0-x*x) * res[lm_index(L, L)];
	res[lm_index(L+1, L+1)] = val;
	res[lm_index(L+1, -L-1)] = pow(-1.0, L+1) / DFactorial(2*(L+1)) * val;
      }

      for(int L = 1; L <= max_l; L++) {
	for(int M = -L; M <= L-2; M++) {
	  dcomplex term(0);
	  if(L-M != 0) {
	    term += (L-M)*1.0*x*res[lm_index(L, M)];
	  }
	  if(L+M != 0) {
	    term -= (L+M)*1.0*res[lm_index(L-1, M)];
	  }
	  res[lm_index(L, M+1)] =  term / sqrt(1.0-x*x);
	}
      }
    }
  }
  void RealSphericalHarmonics(dcomplex theta, dcomplex phi, int max_l, dcomplex* res) {

    // WARNING
    // The expression used for Real valued Spherical Harmonics is different from
    // one in Wiki.
    // 
    // res must be allocate before this function:
    // // res = new dcomplex[num_lm_pair(max_l)];
    
    AssociatedLegendre(cos(theta), max_l, res);
    for(int l = 0; l <= max_l; l++) {
      
      res[lm_index(l, 0)] = sqrt((2*l+1)/(4.0*M_PI)) * res[lm_index(l, 0)];
      for(int m = 1; m <= l; m++) {
	int am = iabs(m);
	dcomplex plm = res[lm_index(l, am)];
	dcomplex t0 = ipow(-1, m)*sqrt(2.0)*
	  sqrt((2*l+1) / (4.0*M_PI)*DFactorial(l-am)*1.0/DFactorial(l+am));
	res[lm_index(l, -m)] = t0 * plm * sin(am*1.0*phi);
	
	dcomplex t1 = (ipow(-1, m) * sqrt(2.0) *
		       sqrt((2.0*l+1)/(4.0*M_PI) * DFactorial(l-m)/DFactorial(l+m)));
	res[lm_index(l,m)] = t1 * plm * cos(m*1.0*phi);
      }
    }
  }
  
  void gto_00_r(dcomplex x, dcomplex y, dcomplex z, int Jpp, int Mpp,dcomplex zeta,
		dcomplex* rs, int num_r,  dcomplex* work, dcomplex* res) {
    /*
      gives one center expansion values of (Jpp,Mpp) wave for s-GTO at (x,y,z).

      Inputs
      ------
      x,y,z : location of s-GTO to expand
      Jpp,Mpp: partial wave to obtain 
      zeta   : orbital exponents
      *rs    : location array for evaluation. (new dcomplex[num_r])
      num_r  : number of location rs.
      *work  : working space (new dcomplex[num_lm(Jpp) + 2*Jpp +2])
     */

    dcomplex a2= x*x+y*y+z*z;
    dcomplex a = sqrt(a2);
    dcomplex theta = cacos(z / a);
    dcomplex phi   = cacos(x/sqrt(x*x+y*y));

    dcomplex* ylm = &work[0];
    dcomplex* il =  &work[num_lm_pair(Jpp)];

    RealSphericalHarmonics(theta, phi, Jpp, ylm);

    for(int i = 0; i < num_r; i++) {
      dcomplex r(rs[i]);
      ModSphericalBessel(2.0*zeta*a*r, Jpp, il);
      res[i] = (4.0*M_PI * exp(-zeta*(r*r+a2)) * il[Jpp] *
		pow(-1.0, Mpp) * ylm[lm_index(Jpp, -Mpp)]);
    }
  }
  void gto_lm_r_center(int l, int m, int Jpp, int Mpp, dcomplex zeta,
			    dcomplex* rs, int num_r, dcomplex* res) {

    for(int i = 0; i < num_r; i++) {
      dcomplex r(rs[i]);
      if(l != Jpp || m != Mpp) {
	res[i] = dcomplex(0.0, 0.0);
      } else {
	res[i] = pow(r, l) * exp(-zeta*r*r);
      }
    }
  }

  void gto_lm_r_general(int l, int m,
			     dcomplex x, dcomplex y, dcomplex z,
			     int Jpp, int Mpp,
			     dcomplex zeta,
			     int lppp_max,
			     dcomplex* rs, int num_r,
			     dcomplex* work,
			     dcomplex* res) {
    std::string msg;
    SUB_LOCATION(msg);
    msg += "not implemented for general case";
    throw std::runtime_error(msg);


    dcomplex *ylm = &work[0];
    dcomplex *il  = &work[num_r];
    dcomplex theta = cacos(z / sqrt(x*x+y*y+z*z));
    dcomplex phi = cacos(x / sqrt(x*x+y*y));
    dcomplex A = sqrt(x*x + y*y + z*z);
    RealSphericalHarmonics(theta, phi, 2*l+lppp_max, ylm);

    for(int i = 0; i < num_r; i++) {
      dcomplex r(rs[i]);
      ModSphericalBessel(2.0*A * r, lppp_max, il);
      dcomplex e_term = exp(-zeta*(r*r+A*A));

      dcomplex cumsum(0);
      for(int lp = 0; lp <= l; lp++)
	for(int lppp = 0; lppp <= lppp_max; lppp++)
	  for(int Jp = iabs(l-lp-lppp); Jp <= l+lp+lppp; Jp++) 
	    for(int Mp = -Jp; Mp <= Jp; Mp++)
	      cumsum += GTOExpansionCoef(l, m, lp, lppp, Jp, Mp, Jpp, Mpp) *
		pow(A, l) * pow(r/A, lp) * e_term * il[lppp] *
		ylm[lm_index(Jp, Mp)];
      res[i] = cumsum;
    }

  }
  void gto_lm_r(int l, int m,
		     dcomplex x, dcomplex y, dcomplex z,
		     int Jpp, int Mpp,
		     dcomplex zeta,
		     int lppp_max,
		     dcomplex* rs, int num_r,
		     dcomplex* work, dcomplex* res) {

    double eps = pow(10.0, -10.0);
    if(l == 0 && m == 0) {
      gto_00_r(x, y, z, Jpp, Mpp, zeta, rs, num_r, work, res);
    } else if(abs(x) < eps && abs(y) < eps && abs(z) < eps) {
      gto_lm_r_center(l, m, Jpp, Mpp, zeta, rs, num_r, res);
    } else {
      gto_lm_r_general(l, m, x, y, z, Jpp, Mpp, zeta, lppp_max,
		       rs, num_r, work, res);
    }
  }


}
