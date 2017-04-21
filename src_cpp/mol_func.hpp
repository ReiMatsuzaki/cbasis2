#ifndef MOL_FUNC_H
#define MOL_FUNC_H

#include <Eigen/Core>
#include "../utils/typedef.hpp"
#include "mult_array.hpp"

namespace cbasis {

  class CartGTO {
  public:
    int nx, ny, nz;
    dcomplex x, y, z;
    dcomplex zeta;
    CartGTO(){}
    CartGTO(int _nx, int _ny, int _nz, dcomplex _x, dcomplex _y, dcomplex _z, dcomplex _zeta):
      nx(_nx), ny(_ny), nz(_nz), x(_x), y(_y), z(_z), zeta(_zeta) {}
    std::string str() const;
  };

  dcomplex dist2(dcomplex dx, dcomplex dy, dcomplex dz);

  // ==== Incomplete Gamma ====
  // K.Ishida J.Comput.Chem. 25, (2004), 739
  // F1 and F2 algorithms 
  // Warning:
  // Eq. (24)(25)(26) are incorrect.
  void IncompleteGamma(int max_m, dcomplex z, dcomplex* res);

  // K.Ishida J.Comput.Chem. 25, (2004), 739
  // G1 and G2 algorithms 
  // Evaluate G(z) = Exp(-z)F(-z) for Re[z] > 0
  void ExpIncompleteGamma(int max_m, dcomplex z, dcomplex* res_list);

  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk);

  void calc_d_coef(int max_ni, int max_nj, int max_n,
		   dcomplex zetaP, dcomplex wPx,
		   dcomplex xi, dcomplex xj, MultArray<dcomplex, 3>& res);


  dcomplex coef_d(dcomplex zetap,
		  dcomplex wPk, dcomplex wAk, dcomplex wBk,
		  int nAk, int nBk, int Nk);

  
  dcomplex coef_R(dcomplex zetaP,
		  dcomplex wPx, dcomplex wPy, dcomplex wPz,
		  dcomplex cx,  dcomplex cy,  dcomplex cz,
		  int mx, int my, int mz, int j, dcomplex* Fjs);

  // ==== molecular integral(To be removed) ====
  /*
  dcomplex GTOOverlap(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB);
  dcomplex GTONormConst(int nAx, int nAy, int nAz, dcomplex zetaA);
  dcomplex GTODipZ(int nAx, int nAy, int nAz, 
		   dcomplex wAx, dcomplex wAy, dcomplex wAz,
		   dcomplex zetaA,
		   int nBx, int nBy, int nBz, 
		   dcomplex wBx, dcomplex wBy, dcomplex wBz,
		   dcomplex zetaB);
  dcomplex GTOKinetic(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB);
  dcomplex GTONuclearAttraction(int nAx, int nAy, int nAz, 
				dcomplex wAx, dcomplex wAy, dcomplex wAz,
				dcomplex zetaA,
				int nBx, int nBy, int nBz, 
				dcomplex wBx, dcomplex wBy, dcomplex wBz,
				dcomplex zetaB,
				dcomplex wCx, dcomplex wCy, dcomplex wCz);
  */
}

#endif
