#ifndef MOLINT_HPP
#define MOLINT_HPP

#include <vector>
#include <boost/array.hpp>
#include "math_utils.hpp"
#include "typedef.hpp"
#include <map>

namespace cbasis {

  using std::vector;

  // -- old 
  dcomplex gto_overlap(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB);

  // -- old 
  dcomplex gto_kinetic(int nAx, int nAy, int nAz,
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB);

  // -- old 
  dcomplex gto_moment_z(int nAx, int nAy, int nAz,
			dcomplex wAx, dcomplex wAy, dcomplex wAz,
			dcomplex zetaA,
			int nBx, int nBy, int nBz, 
			dcomplex wBx, dcomplex wBy, dcomplex wBz,
			dcomplex zetaB);

  // -- old 
  dcomplex gto_nuclear_attraction(int nAx, int nAy, int nAz, 
				  dcomplex wAx, dcomplex wAy, dcomplex wAz,
				  dcomplex zetaA,
				  int nBx, int nBy, int nBz, 
				  dcomplex wBx, dcomplex wBy, dcomplex wBz,
				  dcomplex zetaB,
				  dcomplex wCx, dcomplex wCy, dcomplex wCz);

  class MatrixSet {
  public:
    int nbasis_i_;
    int nbasis_j_;
    std::map<std::string, dcomplex*> mat_map_;
  public:
    MatrixSet(int nbi, int nbj);
    ~MatrixSet();
    int size_basis_i() const;
    int size_basis_j() const;
    void set(std::string label, dcomplex* ptr);    
    dcomplex* get(std::string label);
    std::map<std::string, dcomplex*>& get_map() { return mat_map_; }
    //np::ndarray get_py(std::string label);
  };  
  class GTOs {
  private:
    void Normalize();
    void Add(dcomplex _zeta,
	     dcomplex x, dcomplex y, dcomplex z,
	     vector<int> _nx, vector<int> _ny, vector<int> _nz, 
	     vector<vector<dcomplex> > _coef, int L, vector<int> M);

  public:
    // ish   : index of shell
    // iprim : index of primitive GTO in ish shell.
    // icont : index of contracted basis in ish shell.
    //
    vector<dcomplex>    zeta_ish; 
    vector<dcomplex>    x_ish;
    vector<dcomplex>    y_ish;
    vector<dcomplex>    z_ish;
    vector<vector<int> > nx_ish_iprim;
    vector<vector<int> > ny_ish_iprim;
    vector<vector<int> > nz_ish_iprim;
    vector<int> l_ish;
    vector<vector<int> > m_ish_ibasis;
    
    vector<int> offset_ish; // offset_ish[ish] + ibasis gives global index    
    vector<dcomplex>                   coef_ylm_ish;
    vector<vector<vector<dcomplex> > > coef_ish_ibasis_iprim;

    vector<dcomplex> x_iat;
    vector<dcomplex> y_iat;
    vector<dcomplex> z_iat;
    vector<dcomplex> q_iat;

  public:
    GTOs();
    int size_basis() const {
      int cumsum(0);
      for(int ish = 0; ish < this->size_sh(); ish++) {
	cumsum += this->size_basis_ish(ish);      
      }
      return cumsum;
    }
    int size_sh()  const {
      return zeta_ish.size();    
    }
    int size_basis_ish(int ish) const {
      return coef_ish_ibasis_iprim[ish].size();    
    }
    int size_prim_ish(int ish) const {
      return nx_ish_iprim[ish].size();
    }
    int size_prim() const {
      int cumsum(0);
      for(int ish = 0; ish < this->size_sh(); ish++) {
	cumsum += this->size_prim_ish(ish);
      }
      return cumsum;
    }
    int size_atom() const {
      return x_iat.size();
    }

    void AddSphericalGTOs(int L, dcomplex x, dcomplex y, dcomplex z, dcomplex _zeta);
    void AddOneSphericalGTO(int L, int M,
			    dcomplex x, dcomplex y, dcomplex z, dcomplex _zeta);
    void AddAtom(dcomplex q, dcomplex x, dcomplex y, dcomplex z);

    MatrixSet Calc();
    void CalcMat(dcomplex** s, dcomplex** t, dcomplex** dz, dcomplex** v);
    MatrixSet CalcZMatOther(const GTOs&);
    
    void AtR_Ylm(int l, int m, dcomplex* rs, int num_r, dcomplex* cs, dcomplex* res);

    void Show() const;
    
    // void VMat(dcomplex** v);
  };

}


#endif
