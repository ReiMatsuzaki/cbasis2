#include "gto3dset.hpp"

namespace cbasis {

  // typedef std::complex<double> dcmplx;
  // typedef SphericalGTO<dcmplx, array3<dcmplx> > Basis;

  SphericalGTOSet::SphericalGTOSet() {}
  SphericalGTOSet::~SphericalGTOSet() {
    for(iterator it = basis_list_.begin(); 
	it != basis_list_.end(); ++it) {
      delete(*it);
    }
  }
  int SphericalGTOSet::size() const { return basis_list_.size(); }
  const SGTO& SphericalGTOSet::basis(int i) const {
    return *basis_list_[i];
  }
  void SphericalGTOSet::AddOneBasis(int L, int M, dcomplex x, dcomplex y, dcomplex z, dcomplex zeta) {
    SGTO* ptr;
    try {
      ptr = new SGTO(L, M, c3(x, y, z), zeta);
    } catch(const ExceptionBadYlm& e) {
      throw(e);
    }
    this->basis_list_.push_back(ptr);
  }
  void SphericalGTOSet::AddBasis(int L, dcomplex x, dcomplex y, dcomplex z, dcomplex zeta) {

    for (int M = -L; M <= L; M++) {
      this->AddOneBasis(L, M, x, y, z, zeta);
    }

  }


  dcomplex SphericalGTOSet::VMatEle(int i, int j, dcomplex q, dcomplex x, dcomplex y, dcomplex z) const {
    OpNA<dcomplex, c3> op(q, c3(x, y, z));
    return CIP(this->basis(i), op, this->basis(j));
  }

/*
  template<class TOp>
  dcomplex* CalcMat(const SphericalGTOSet& a, const TOp& op, const SphericalGTOSet& b) {
    int numi = a.size();
    int numj = b.size();
    dcomplex* vs = new dcomplex[numi * numj];

    for(int i = 0; i < numi; i++) {      
      for(int j = 0; j < numj; j++) {
	vs[i*numj + j] = CIP(a.basis(i), op, a.basis(j));;
      }
    }
    return vs;
  }
*/
  dcomplex* SphericalGTOSet::SMat() const {

    int num = this->basis_list_.size();
    dcomplex* vs = new dcomplex[num * num];
    int i(0);
    for(const_iterator it_i = basis_list_.begin();
	it_i != basis_list_.end(); ++it_i) {      
      int j(0);
          for(const_iterator it_j = basis_list_.begin();
	      it_j != basis_list_.end(); ++it_j) {
	    vs[i*num + j] = CIP(**it_i, **it_j);
	    j++;
	  }
	  i++;
    }
    return vs;

    /*
      int num = this->size();
      dcomplex* vs = new dcomplex[num * num];
    for(int i = 0; i < num; i++)
      for(int j = 0; j < num; j++)
	vs[i*num+j] = CIP(this->basis(i), this->basis(j));
    return vs;
    */
  }
  dcomplex* SphericalGTOSet::TMat() const {
    int num = this->size();
    dcomplex* vs = new dcomplex[num * num];
    for(int i = 0; i < num; i++)
      for(int j = 0; j < num; j++)
	vs[i*num+j] = CIP(this->basis(i), OpKE<dcomplex, c3>(), this->basis(j));
    return vs;
  }
  
  dcomplex* SphericalGTOSet::VMat(dcomplex q, dcomplex x, dcomplex y, dcomplex z) const {
    OpNA<dcomplex, c3> op(q, c3(x, y, z));
    int num = this->size();
    dcomplex* vs = new dcomplex[num * num];    
    for(int i = 0; i < num; i++)
      for(int j = 0; j < num; j++)
	vs[i*num+j] = CIP(this->basis(i), op, this->basis(j));
    return vs;
  }
  dcomplex* SphericalGTOSet::XyzMat(int nx, int ny, int nz) const {
    
    OpXyz<dcomplex, c3> op(nx, ny, nz);
    int num = this->size();
    dcomplex* vs = new dcomplex[num * num];    
    for(int i = 0; i < num; i++)
      for(int j = 0; j < num; j++)
	vs[i*num+j] = CIP(this->basis(i), op, this->basis(j));
    return vs;
  }

}
