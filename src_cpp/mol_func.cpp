#include <string>
#include <stdexcept>
#include <sstream>
#include "mol_func.hpp"
#include "../utils/fact.hpp"

using namespace std;

namespace cbasis {

  typedef MultArray<dcomplex, 3> A3dc;

  string CartGTO::str() const {
    ostringstream oss;
    oss << nx << ny << nz << ": " << x << y << z << ": " << zeta;
    return oss.str();
  }
  dcomplex dist2(dcomplex dx, dcomplex dy, dcomplex dz) {
    return dx*dx+dy*dy+dz*dz;
  }


  void IncompleteGamma_F1(int max_m, dcomplex z, dcomplex* res_list) {

    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));

    if(x < +eps || y < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "Re[z] and Im[z] must positive in F1 algorithms";
      throw std::runtime_error(msg);
    }

    double pi(M_PI);
    double z2 = x*x+y*y;
    double az = abs(z);
    double c = sqrt(pi*(az+x)/(8.0*z2));
    double s = sqrt(pi*(az-x)/(8.0*z2)) + eps;

    int NF_max(50);
    //    int NF(0);
    int NF_0(10);
    double* anR = new double[NF_max+1];
    double* anI = new double[NF_max+1];
    
    double phiR = -x/(2.0*z2)*exp(-x)*cos(y) + y/(2.0*z2)*exp(-x)*sin(y);
    double phiI = +y/(2.0*z2)*exp(-x)*cos(y) + x/(2.0*z2)*exp(-x)*sin(y);
    anR[0] = phiR;
    anI[0] = phiI;
    double sum_anR(anR[0]);
    double sum_anI(anI[0]);
    double delta(pow(10.0, -15.0));
    bool convq(false);
    for(int n = 1; n <= NF_max; n++) {
      anR[n] = -(2.0*n-1) * (x/(2.0*z2)*anR[n-1] + y/(2.0*z2)*anI[n-1]);
      anI[n] = +(2.0*n-1) * (y/(2.0*z2)*anR[n-1] - x/(2.0*z2)*anI[n-1]);      
      sum_anR += anR[n];
      sum_anI += anI[n];
      if(n > NF_0)
	if(c * delta > std::abs(anR[n]) && s*delta > std::abs(anI[n])) {
	  convq = true;
	  //	  NF = n;
	  break;
	}
    }
    
    if(!convq) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss; 
      oss << msg << ": not converged" << endl;
      oss << "z = " << z << endl;
      oss << "c = " << c << endl;
      oss << "s = " << s << endl;
      oss << "delta = " << delta << endl;
      oss << "anR[NF_max] = " << anR[NF_max] << endl;
      
      oss << "anI[NF_max] = " << anI[NF_max] << endl;
      
      throw runtime_error(oss.str());
    }

    double* fmR = new double[max_m+1];
    double* fmI = new double[max_m+1];
    fmR[0] = +c + sum_anR;
    fmI[0] = -s + sum_anI;
    for(int m = 1; m <= max_m; m++) {
      fmR[m] = (2*m-1.0) * (x/(2.0*z2)*fmR[m-1] + y/(2.0*z2)*fmI[m-1]);
      fmI[m] = (2*m-1.0) * (x/(2.0*z2)*fmI[m-1] - y/(2.0*z2)*fmR[m-1]);
    }

    double* bmR = new double[max_m+1];
    double* bmI = new double[max_m+1];
    bmR[0] = 0.0;
    bmI[0] = 0.0;
    //    bmR[1] = phiR;
    //    bmI[1] = phiI;
    for(int m = 1; m <= max_m; m++) {
      bmR[m] = phiR + (2*m-1.0) * (x/(2.0*z2)*bmR[m-1] + y/(2.0*z2)*bmI[m-1]);
      bmI[m] = phiI + (2*m-1.0) * (x/(2.0*z2)*bmI[m-1] - y/(2.0*z2)*bmR[m-1]);
    }
    //    std::cout << "bmR[2]: " << bmR[2] << std::endl;    

    for(int m = 0; m <= max_m; m++) 
      res_list[m] = dcomplex(fmR[m]+bmR[m], fmI[m]+bmI[m]);

    delete[] anR; delete[] anI;
    delete[] fmR; delete[] fmI;
    delete[] bmR; delete[] bmI;
  }
  void IncompleteGamma_F2(int max_m, dcomplex z, dcomplex* res) {

    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));
    
    if(x < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "Re[z] must positive in F2 algorithms";
      throw std::runtime_error(msg);
    }

    if(y < -eps) {
      IncompleteGamma_F2(max_m, dcomplex(x, -y), res);
      for(int m = 0; m <= max_m; m++)
	res[m] = conj(res[m]);
    }

    static const int NR(47);
    dcomplex Bn[NR+1];
    dcomplex An[NR+1];

    Bn[0] = 1.0;
    Bn[1] = 1.0 + 0.5*z;
    for(int n = 2; n <= NR; n++) 
      Bn[n] = Bn[n-1] + z*z/(4*(2*n-1)*(2*n-3)*1.0)*Bn[n-2];

    for(int m = 0; m <= max_m; m++) {

      An[0] = 1.0;
      double t1 = (2.0*m+1) / (2.0*m+3);
      double t2 = (2.0*m+1) / ((2*m+3)*(2*m+5));
      double t3 = (double((2*m+1)*((2*m+1)*(2*m+1)+44)) /
		   double(60*(2*m+3)*(2*m+5)*(2*m+7)));
	
      An[1] = Bn[1] - t1*z;
      An[2] = Bn[2] - t1*z - t2*z*z;
      An[3] = Bn[3] - t1*z - t2*z*z - t3*z*z*z;
      for(int n = 4; n <= NR; n++) {
	double F1 = double(2*n-2*m-5)/double(2*(2*n-3)*(2*n+2*m+1));
	double F2 = double(1) / double(4*(2*n-1)*(2*n-3));
	double F3 = double(-F1) / double(4*(2*n-3)*(2*n-5));
	double E = -F1;
	An[n] = (1.0+F1*z)*An[n-1] + (E + F2*z)*z*An[n-2] + F3*z*z*z*An[n-3];
      }
      res[m] = 1.0 / (2*m+1) * An[NR]/Bn[NR];
    }
    
  }
    
  void IncompleteGamma(int max_m, dcomplex z, dcomplex* res_list) {
    
    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -15.0));

    if(x < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += ": negative Re[z] is not supported. Use ExpIncompleteGamma with -z";
      cout << "z: " << z << endl;
      throw std::runtime_error(msg);
    } 
    
    if(x > -eps && y > -eps) {
      if(x < 21.0 && x+y < 37.0) 
	IncompleteGamma_F2(max_m, z, res_list);
      else
	IncompleteGamma_F1(max_m, z, res_list);
    } else {
      IncompleteGamma(max_m, dcomplex(x, -y), res_list);
      for(int m = 0; m <= max_m; m++)
	res_list[m] = conj(res_list[m]);
    }

  }

  void ExpIncompleteGamma_G1(int max_m, dcomplex z, dcomplex *res_list) {
    cout << "G1" << endl;
    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));
    double delta(pow(10.0, -15.0));
    static const int MAX_M(100);
    if(max_m >= MAX_M) {
      string msg; SUB_LOCATION(msg); msg += ": max_m must be lesser than 100";
      throw runtime_error(msg);
    }

    if(x < +eps || y < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "Re[z] and Im[z] must positive in F1 algorithms";
      throw std::runtime_error(msg);
    }
    double xxyy(x*x + y*y);
    double az(sqrt(xxyy));
    double two_xxyy = 2.0 * xxyy;

    double c(sqrt((M_PI * (az+x) / (8.0 * xxyy))));
    double s(sqrt((M_PI * (az-x) / (8.0 * xxyy))));
    double an_R(0), an_I(0), anm1_R(0), anm1_I(0);
    bool conv(false);
    int NF_max(50);
    anm1_R = x /two_xxyy; anm1_I = -y / two_xxyy;
    double Sn_R((c*exp(-x)*sin(y) + s*exp(-x)*cos(y)) + anm1_R);
    double Sn_I((c*exp(-x)*cos(y) - s*exp(-x)*sin(y)) + anm1_I);
    for(int n = 1; n < NF_max; n++) {      
      an_R = (2 * n - 1) * ( x / two_xxyy * anm1_R + y / two_xxyy * anm1_I);
      an_I = (2 * n - 1) * ( x / two_xxyy * anm1_I - y / two_xxyy * anm1_R);
      Sn_R += an_R;
      Sn_I += an_I;
      anm1_R = an_R;
      anm1_I = an_I;
      if(delta > abs(an_R/Sn_R) && delta > abs(an_I/Sn_I)) {
	conv = true;
	break;
      }
    }
    
    if(!conv) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss; 
      oss << msg << ": not converged" << endl;
      throw runtime_error(oss.str());
    }

    double g_R[MAX_M];
    double g_I[MAX_M];
    g_R[0] = Sn_R;
    g_I[0] = Sn_I;
    for(int m = 1; m <= max_m; m++) {
      g_R[m] = -(2*m-1) * (x/two_xxyy * g_R[m-1] + y/two_xxyy * g_I[m-1]);
      g_I[m] = +(2*m-1) * (y/two_xxyy * g_R[m-1] - x/two_xxyy * g_I[m-1]);
    }

    double b_R[MAX_M];
    double b_I[MAX_M];
    b_R[0] = 0.0; // b_R[1] = +x/two_xxyy;
    b_I[0] = 0.0; // b_I[1] = -y/two_xxyy;
    for(int m = 1; m <= max_m; m++) {
      b_R[m] = +x/two_xxyy - (2*m-1)*(x/two_xxyy * b_R[m-1] + y/two_xxyy*b_I[m-1]);
      b_I[m] = -y/two_xxyy + (2*m-1)*(y/two_xxyy * b_R[m-1] - x/two_xxyy*b_I[m-1]);
    }
    
    for(int m = 0; m <= max_m; m++) {
      res_list[m] = dcomplex(g_R[m] + b_R[m], g_I[m] + b_I[m]);
    }

  }
  void ExpIncompleteGamma_G2(int max_m, dcomplex z, dcomplex *res_list) {

    //    cout << "G2" << endl;

    double x = real(z);
    double y = imag(z);
    double delta(pow(10.0, -15.0));
    /*
    static const int MAX_M(100);
    if(max_m >= MAX_M) {
      string msg; SUB_LOCATION(msg); msg += ": max_m must be lesser than 100";
      throw runtime_error(msg);
    }
    */
    if(x < -delta) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "Re[z] must positive in G2 algorithms";
      throw std::runtime_error(msg);
    }

    for(int m = 0; m <= max_m; m++) {

      bool conv(false);
      int max_n(40);
      double An_R(0), An_I(0), Bn_R(0), Bn_I(0);
      double Anm1_R(0), Anm1_I(0), Bnm1_R(0), Bnm1_I(0);
      double Anm2_R(0), Anm2_I(0), Bnm2_R(0), Bnm2_I(0);
      Bnm2_R = 1;                      Bnm2_I = 0.0;
      Bnm1_R = 1 + 2.0/(2.0*m+5.0)*x;  Bnm1_I = 2.0/(2.0*m+5.0)*y;
      Anm2_R = 1;                      Anm2_I = 0; 
      Anm1_R = Bnm1_R - 2.0*x/(2*m+3); Anm1_I = Bnm1_I - 2.0*y/(2*m+3); 
      dcomplex Bn(0), An(0);

      for(int n = 2; n < max_n; n++) {
	double F1 = (double(2*(2*m+1)) /
		     double((4*n+2*m+1)*(4*n+2*m-3)));
	double F2 = (double(8*(n-1)*(2*n+2*m-1)) /
		     double((4*n+2*m-1)*(4*n+2*m-3)*(4*n+2*m-3)*(4*n+2*m-5)));
	An_R = (1+F1*x)*Anm1_R - F1*y*Anm1_I + (x*x-y*y)*F2*Anm2_R - 2*x*y*F2*Anm2_I;
	An_I = (1+F1*x)*Anm1_I + F1*y*Anm1_R + (x*x-y*y)*F2*Anm2_I + 2*x*y*F2*Anm2_R;
	Bn_R = (1+F1*x)*Bnm1_R - F1*y*Bnm1_I + (x*x-y*y)*F2*Bnm2_R - 2*x*y*F2*Bnm2_I;
	Bn_I = (1+F1*x)*Bnm1_I + F1*y*Bnm1_R + (x*x-y*y)*F2*Bnm2_I + 2*x*y*F2*Bnm2_R;

	Bn = 0.0;
	for(int j = 0; j <= n; j++)
	  Bn += (double(DoubleFactorial(2*n+2*m+2*j+1)) / 
		 double(DoubleFactorial(4*n+2*m+1)) *
		 Combination(n, j) *
		 pow(2.0*z, n-j));
	
	dcomplex Gn(dcomplex(An_R, An_I) / dcomplex(Bn_R, Bn_I));
	dcomplex Gnm1(dcomplex(Anm1_R, Anm1_I) / dcomplex(Bnm1_R, Bnm1_I));
	if(abs((Gn-Gnm1)/Gn) < delta) {
	  conv = true; break;
	}
	Anm2_R = Anm1_R; Anm2_I = Anm1_I; Anm1_R = An_R; Anm1_I = An_I; 
	Bnm2_R = Bnm1_R; Bnm2_I = Bnm1_I; Bnm1_R = Bn_R; Bnm1_I = Bn_I; 
      }

      if(!conv) {
	string msg; SUB_LOCATION(msg); msg += ": not converged.";
	cout <<
	  ((abs((Anm1_R-Anm2_R)/An_R) < delta) ? "Yes" : "No") << endl << 
	  ((abs((Anm1_I-Anm2_I)/An_I) < delta) ? "Yes" : "No") << endl << 
	  ((abs((Bnm1_R-Bnm2_R)/Bn_R) < delta) ? "Yes" : "No") << endl << 
	  ((abs((Bnm1_I-Bnm2_I)/Bn_I) < delta) ? "Yes" : "No") << endl;
	cout <<
	  (Anm1_R-Anm2_R)/An_R << endl << 
	  (Anm1_I-Anm2_I)/An_I << endl << 
	  (Bnm1_R-Bnm2_R)/Bn_R << endl << 
	  (Bnm1_I-Bnm2_I)/Bn_I << endl;
	throw runtime_error(msg);
      }

      res_list[m] = 1.0/(2*m+1) * dcomplex(An_R, An_I) / dcomplex(Bn_R, Bn_I);
    }
  }
  
  void ExpIncompleteGamma(int max_m, dcomplex z, dcomplex* res_list) {
    //    cout << "ExpIncGamma: " << z << endl;
    double x = real(z);
    double y = imag(z);
    double eps(pow(10.0, -10.0));
    

    if(x < -eps) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += ": negative Re[z] is not supported. use IncompleteGamma with -z.";
      cout << "z: " << z << endl;
      throw std::runtime_error(msg);
    }
 
    if(x > -eps && y > -eps) {
      if( (x > 36.0 && y > 36) || x + y > 51) {
	ExpIncompleteGamma_G1(max_m, z, res_list);
      } else {
	ExpIncompleteGamma_G2(max_m, z, res_list);
      }
    } else {
      ExpIncompleteGamma(max_m, dcomplex(x, -y), res_list);
      for(int m = 0; m <= max_m; m++)
	res_list[m] = conj(res_list[m]);
    }
  }

  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk) {
    
    //    std::cout << nAk << nBk << Nk <<std::endl;
		  
    if(nAk == 0 && nBk == 0 && Nk == 0)
      return 1.0;

    if(Nk < 0 || Nk > nAk + nBk) 
      return 0.0;

    if(nAk > 0) 
      return 1.0/(2.0*zetap) * coef_d(zetap, wPk, wAk, wBk, nAk-1, nBk, Nk-1) +
	(wPk - wAk)          * coef_d(zetap, wPk, wAk, wBk, nAk-1, nBk, Nk)   +
	(Nk + 1.0)           * coef_d(zetap, wPk, wAk, wBk, nAk-1, nBk, Nk+1);
    else 
      return 1.0/(2.0*zetap) * coef_d(zetap, wPk, wAk, wBk, nAk, nBk-1, Nk-1) +
	(wPk - wBk)          * coef_d(zetap, wPk, wAk, wBk, nAk, nBk-1, Nk)   +
	(Nk + 1.0)           * coef_d(zetap, wPk, wAk, wBk, nAk, nBk-1, Nk+1);
    
  }

  void calc_d_coef(int max_ni, int max_nj, int max_n,
		   dcomplex zetaP, dcomplex wPx,
		   dcomplex xi, dcomplex xj,
		   MultArray<dcomplex, 3>& res) {
    
    //
    //  compute coef d by recursion formula.
    //
    //  d(nAk+1, nBk,   Nk) = (2zetaP)^-1 d(nAk,nBk,Nk-1) + (WPk-WAk)d(nAk,nBk,Nk) + (Nk+1)d(nAk,nBk,Nk+1)
    //  d(nAk,   nBk+1, Nk) = (2zetaP)^-1 d(nAk,nBk,Nk-1) + (WPk-WBk)d(nAk,nBk,Nk) + (Nk+1)d(nAk,nBk,Nk+1)
    //

    res.SetRange(0, max_ni,
		 0, max_nj,
		 -max_ni-max_nj, max_ni+max_nj+max_n);

    res(0, 0, 0) = 1.0;
    for(int n = -max_ni -max_nj; n < 0; n++) 
      res(0, 0, n) =  0.0;
    for(int n = 1; n <= max_ni+max_nj+max_n; n++) 
      res(0, 0, n) = 0.0;

    for(int ni_nj = 1; ni_nj <= max_ni + max_nj; ni_nj++) {
      for(int ni = 0; ni <= std::min(max_ni, ni_nj); ni++) {
	int nj = ni_nj - ni;
	if(nj <= max_nj) {
	  for(int n = -(max_ni + max_nj - ni_nj);
	      n <= max_ni + max_nj + max_n - ni_nj; n++) {
	    if(ni > 0) {
	      dcomplex v = (1.0/(2.0*zetaP) * res(ni-1, nj, n-1) +
			    (wPx - xi)      * res(ni-1, nj, n) + 
			    (n + 1.0)      * res(ni-1, nj, n+1));
	      res(ni, nj, n) = v;
	    } else if(nj > 0) {
	      dcomplex v = (1.0/(2.0*zetaP) * res(ni, nj-1, n-1) +
			    (wPx - xj)      * res(ni, nj-1, n) + 
			    (n + 1.0)      * res(ni, nj-1, n+1));
	      res(ni, nj, n) = v;
	    }
	  }
	}
      }
    }
  }

  dcomplex coef_R(dcomplex zetaP,
		  dcomplex wPx, dcomplex wPy, dcomplex wPz,
		  dcomplex cx,  dcomplex cy,  dcomplex cz,
		  int mx, int my, int mz, int j, dcomplex* Fjs) {
    /* Compute function R_{mx,my,mz,j}(zetaP, Wp, C) where
       wP = (wPx, wPy, wPz)
       C  = (Cx,  Cy,  Cz).
       by recursion formula. The last term can be obtained by 
       R_{0,0,0,j} = (-2zetaP)^j * Fj(zetaP(wP-C)^2).
     */

    if(mx == 0 && my == 0 && mz == 0) {
      //      dcomplex arg = zetaP * dist2(wPx-cx, wPy-cy, wPz-cz);
      //      dcomplex* Fjs = IncompleteGamma(j, arg);
      dcomplex Fj = Fjs[j];
      return pow(-2.0*zetaP, j) * Fj;
    }
    
    if(mx > 0) {
      dcomplex res(0.0);
      if(mx > 1) 
	res += (mx-1.0) * coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx-2, my, mz, j+1, Fjs);
      res += (wPx-cx) *   coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx-1, my, mz, j+1, Fjs);
      return res;
    }
    if(my > 0) {
      dcomplex res(0.0);
      if(my > 1) 
	res += (my-1.0) * coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my-2, mz, j+1, Fjs);
      res += (wPy-cy) *   coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my-1, mz, j+1, Fjs);
      return res;
    }
    if(mz > 0) {
      dcomplex res(0.0);
      if(mz > 1) 
	res += (mz-1.0) * coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my, mz-2, j+1, Fjs);
      res += (wPz-cz) *   coef_R(zetaP, wPx, wPy, wPz, cx, cy, cz, mx, my, mz-1, j+1, Fjs);
      return res;
    }
    std::string msg;
    SUB_LOCATION(msg);
    msg += " one of mx, my, mz is negative integer.";
    throw std::runtime_error(msg);    
    return 0.0;    
  }  

  dcomplex GTOOverlap(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB) {

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx = (zetaA*wAx+zetaB*wBx)/zetaP;
    dcomplex wPy = (zetaA*wAy+zetaB*wBy)/zetaP;
    dcomplex wPz = (zetaA*wAz+zetaB*wBz)/zetaP;
    A3dc dxmap(100);
    A3dc dymap(100); 
    A3dc dzmap(100); 
    calc_d_coef(nAx,nBx,0,zetaP,wPx,wAx,wBx,dxmap);
    calc_d_coef(nAy,nBy,0,zetaP,wPy,wAy,wBy,dymap);
    calc_d_coef(nAz,nBz,0,zetaP,wPz,wAz,wBz,dzmap);
    
    dcomplex d2 = pow(wAx-wBx,2) + pow(wAy-wBy,2) + pow(wAz-wBz,2);	
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*d2);
    dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);

    dcomplex res = (ce *
		    dxmap(nAx, nBx, 0) * 
		    dymap(nAy, nBy, 0) * 
		    dzmap(nAz, nBz, 0));
    return res;

  }
  dcomplex GTONormConst(int nAx, int nAy, int nAz, dcomplex zetaA) {

    dcomplex norm2 = GTOOverlap(nAx, nAy, nAz, 0.0, 0.0, 0.0, zetaA,
				nAx, nAy, nAz, 0.0, 0.0, 0.0, zetaA);
    return 1.0 / sqrt(norm2);
    
  }
  dcomplex GTODipZ(int nAx, int nAy, int nAz, 
		   dcomplex wAx, dcomplex wAy, dcomplex wAz,
		   dcomplex zetaA,
		   int nBx, int nBy, int nBz, 
		   dcomplex wBx, dcomplex wBy, dcomplex wBz,
		   dcomplex zetaB) {

    return (GTOOverlap(nAx, nAy, nAz,   wAx, wAy, wAz, zetaA,
		       nBx, nBy, nBz+1, wBx, wBy, wBz, zetaB) +
	    wBz * 
	    GTOOverlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
		       nBx, nBy, nBz, wBx, wBy, wBz, zetaB));
	    

  }
  dcomplex GTOKinetic(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB) {

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx = (zetaA*wAx+zetaB*wBx)/zetaP;
    dcomplex wPy = (zetaA*wAy+zetaB*wBy)/zetaP;
    dcomplex wPz = (zetaA*wAz+zetaB*wBz)/zetaP;
    A3dc dxmap(100); calc_d_coef(nAx,nBx+2,0,zetaP,wPx,wAx,wBx,dxmap);
    A3dc dymap(100); calc_d_coef(nAy,nBy+2,0,zetaP,wPy,wAy,wBy,dymap);
    A3dc dzmap(100); calc_d_coef(nAz,nBz+2,0,zetaP,wPz,wAz,wBz,dzmap);
    dcomplex d2 = pow(wAx-wBx,2) + pow(wAy-wBy,2) + pow(wAz-wBz,2);	
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*d2);
    dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);

    dcomplex dx00 = dxmap(nAx, nBx,   0);
    dcomplex dx02 = dxmap(nAx, nBx+2, 0);
    dcomplex dy00 = dymap(nAy, nBy,   0);
    dcomplex dy02 = dymap(nAy, nBy+2, 0);
    dcomplex dz00 = dzmap(nAz, nBz,   0);
    dcomplex dz02 = dzmap(nAz, nBz+2, 0);

    dcomplex res = 4.0 * zetaB * zetaB *
      (dx02*dy00*dz00 + dx00*dy02*dz00 + dx00*dy00*dz02);
    res += -2.0*(2*nBx+2*nBy+2*nBz+3)*zetaB*dx00*dy00*dz00;
    
    if(nBx > 1)
      res += 1.0*nBx*(nBx-1) * dxmap(nAx, nBx-2, 0) * dy00 * dz00;
    if(nBy > 1)
      res += 1.0*nBy*(nBy-1) * dymap(nAy, nBy-2, 0) * dx00 * dz00;
    if(nBz > 1)
      res += 1.0*nBz*(nBz-1) * dzmap(nAz, nBz-2, 0) * dx00 * dy00;

    res *= -0.5 * ce;
    return res;    
  }
  dcomplex GTONuclearAttraction(int nAx, int nAy, int nAz, 
				dcomplex wAx, dcomplex wAy, dcomplex wAz,
				dcomplex zetaA,
				int nBx, int nBy, int nBz, 
				dcomplex wBx, dcomplex wBy, dcomplex wBz,
				dcomplex zetaB,
				dcomplex wCx, dcomplex wCy, dcomplex wCz) {

    dcomplex* Fjs = new dcomplex[100];

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx = (zetaA*wAx+zetaB*wBx)/zetaP;
    dcomplex wPy = (zetaA*wAy+zetaB*wBy)/zetaP;
    dcomplex wPz = (zetaA*wAz+zetaB*wBz)/zetaP;

    dcomplex dx = wPx-wCx;
    dcomplex dy = wPy-wCy;
    dcomplex dz = wPz-wCz;
    dcomplex d2p = dx*dx + dy*dy + dz*dz;

    A3dc dxmap(100); calc_d_coef(nAx,nBx,nAx+nBx,zetaP,wPx,wAx,wBx,dxmap);
    A3dc dymap(100); calc_d_coef(nAy,nBy,nAy+nBy,zetaP,wPy,wAy,wBy,dymap);
    A3dc dzmap(100); calc_d_coef(nAz,nBz,nAz+nBz,zetaP,wPz,wAz,wBz,dzmap);
    IncompleteGamma(nAx+nBx+nAy+nBy+nAz+nBz, zetaP * d2p, Fjs);

    dcomplex d2 = pow(wAx-wBx,2) + pow(wAy-wBy,2) + pow(wAz-wBz,2);	
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*d2);

    dcomplex res(0);
    for(int nx = 0; nx <= nAx + nBx; nx++)
      for(int ny = 0; ny <= nAy + nBy; ny++)
	for(int nz = 0; nz <= nAz + nBz; nz++)
	  res += (dxmap(nAx, nBx, nx) *
		  dymap(nAy, nBy, ny) *
		  dzmap(nAz, nBz, nz) *
		  coef_R(zetaP, wPx, wPy, wPz,
			 wCx, wCy, wCz,
			 nx, ny, nz, 0, Fjs));

    res *= -2.0 * M_PI/zetaP * eAB;
    delete[] Fjs;
    return res;    
  }

}
