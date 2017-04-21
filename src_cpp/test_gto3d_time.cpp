//#include "utils.hpp"
//#include "cints.hpp"
//#include "angmoment.hpp"
//#include "gto3d.hpp"
#include "gto3dset.hpp"
#include "timer.hpp"
#include "molint.hpp"

using namespace cbasis;

void gto_set() {
  /*
    2016/3/2
    ./test_gto3d_time
    label : time / s
    S : 0.108516
    T : 0.568004
    V : 1.93674
    basis : 0.004714
   */

  Timer timer;
  
  timer.Start("basis");
  SphericalGTOSet gtos;
  for(int L = 0; L <= 2; L++)
    for(int n = -5; n < 5; n++)
      gtos.AddBasis(L, 0.1, 0.2, 0.3, pow(2.0, n));
  timer.End("basis");

  timer.Start("S");
  dcomplex* smat = gtos.SMat();
  timer.End("S");
  timer.Start("V");
  dcomplex* vmat = gtos.VMat(1.0, 0.0, 0.0, 0.0);
  timer.End("V");
  timer.Start("T");
  dcomplex* tmat = gtos.TMat();
  timer.End("T");

  timer.Display();

  delete smat;
  delete vmat;
  delete tmat;

}
void molint() {

  // 2016/3/4  : 2.70645
  //           : 2.66855 (reduce new)
  //           : 1.25852 (add transpose)
  //           : 0.40530 (reduce number of calling IncompleteGamma)
  // 2016/3/4  : 0.152456 (pre calculate coef d)
  //             0.089503 (pre calculate for all d)
  //                in molint
  //                    27 %: coef_d_coef
  //                    24 %: 
  //                    19 %: Incomplete Gamma
  //                    16 %: coef_R
  //           : 0.054927 (pre calculation of Fjs)
  //           :: 0.047015 (get_safe -> get)
  //           : 0.118636 (forgot H atom)
  //           : 0.297726 (in debug)
  //             in molint
  //                   27 : x
  //                   22 : calc_d_coef 
  //                   20 : coef_R
  //                   20 : IncompleteGamma

  Timer timer;
  GTOs gtos;
  for(int L = 0; L <= 2; L++)
    for(int n = -5; n < 5; n++)
      gtos.AddSphericalGTOs(L, 0.1, 0.2, 0.3, pow(2.0, n));
  gtos.AddAtom(1.0, 0.0, 0.0, 0.0);
  dcomplex* smat;
  dcomplex* tmat;
  dcomplex* zmat;
  dcomplex* vmat;

  timer.Start("molint");
  std::cout << 1;
  gtos.CalcMat(&smat, &tmat, &zmat, &vmat);
  timer.End("molint");
  timer.Display();

  delete[] smat;
  delete[] tmat;
  delete[] zmat;
  delete[] vmat;
}
void molint2() {

  Timer timer;
  GTOs gtos;
  for(int nx = 0; nx < 1; nx++) {
    dcomplex x = 0.1*nx;
    gtos.AddAtom(1.0, x, 0.0, 0.0);
    for(int L = 0; L <= 2; L++)
      for(int n = -3; n < 3; n++)
	gtos.AddSphericalGTOs(L, x, 0.2, 0.3, pow(2.0, n));
  }

  dcomplex* smat;
  dcomplex* tmat;
  dcomplex* zmat;
  dcomplex* vmat;

  timer.Start("molint");
  std::cout << 1;
  gtos.CalcMat(&smat, &tmat, &zmat, &vmat);
  timer.End("molint");
  timer.Display();

  delete[] smat;
  delete[] tmat;
  delete[] zmat;
  delete[] vmat;
}

int main(){

  gto_set();  
  molint();
  //  molint2();
  
}

