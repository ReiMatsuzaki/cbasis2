#ifndef GTO3DSET_H
#define GTO3DSET_H

#include <vector>
#include "gto3d.hpp"

namespace cbasis {

  typedef SphericalGTO<dcomplex, dcomplex> SGTO;

  class SphericalGTOSet {
  private:
    
    std::vector<SGTO*> basis_list_;
  public:
    typedef vector<SGTO*>::const_iterator const_iterator;
    typedef vector<SGTO*>::iterator iterator;

  public:
    SphericalGTOSet();
    ~SphericalGTOSet();
    int size() const;
    const SGTO& basis(int i) const;
    void AddOneBasis(int L, int M, dcomplex x, dcomplex y, dcomplex z, dcomplex zeta);
    void AddBasis(int L, dcomplex x, dcomplex y, dcomplex z, dcomplex zeta);
    dcomplex SMatEle(int i, int j) const;
    dcomplex TMatEle(int i, int j) const;
    dcomplex VMatEle(int i, int j, dcomplex q, dcomplex x, dcomplex y, dcomplex z) const;
    dcomplex* SMat() const;
    dcomplex* TMat() const;
    dcomplex* VMat(dcomplex q, dcomplex x, dcomplex y, dcomplex z) const;
    dcomplex* XyzMat(int nx, int ny, int nz) const;
		   
    // np::ndarray SMatWtihOther(const SphericalGTOSet& o);
  };

}

#endif
