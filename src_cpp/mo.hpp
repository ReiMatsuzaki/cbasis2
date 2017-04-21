#ifndef MO_H
#define MO_H

#include <vector>
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include "../utils/typedef.hpp"
#include "symgroup.hpp"
#include "symmolint.hpp"

namespace cbasis {

  struct SCFOptions {
    bool use_real;
    int max_iter;
    double tol;
    int debug_lvl;
    SCFOptions():
      use_real(false), max_iter(50), tol(0.0000001), debug_lvl(0) {}
  };
  
  class _MO {
  public:
    _MO();
    SymmetryGroup sym;
    BMat H, S, J, K, F, C, P;
    BVec eigs;
    std::vector<int> num_occ_irrep;
    std::vector<Irrep> irrep_list;
    dcomplex energy;
    
    const BMat& HMat() { return H; }
    BMat& SMat() { return S; }
    BMat& JMat() { return J; }
    BMat& KMat() { return K; }
    BMat& FMat() { return F; }
    BMat& CMat() { return C; }
    BMat& PMat() { return P; }
    BVec& Eigs() { return eigs; }
  };
  typedef boost::shared_ptr<_MO> MO;

  // ---- RHF ----
  MO NewMO(SymmetryGroup sym, BMat& _H, BMat& _C, BVec& _eigs, int num_ele);
  vector<int> CalcOccNum(const BVec& eigs, int num_sym, int num_orb);
  MO CalcOneEle(SymmetryGroup sym, BMatSet mat_set, int debug_lvl = 0);
  void AddJK(B2EInt eri,  BMat& C, int I0, int i0,
	     dcomplex coef_J, dcomplex coef_K, BMat& JK, bool use_real=false);
  void AddJK_Slow(B2EInt eri, BMat& C, int I0, int i0,
		  dcomplex coef_J, dcomplex coef_K, BMat& H);
  void AddJ(B2EInt eri, const Eigen::VectorXcd& Ca, Irrep ir_a, dcomplex coef, BMat& J, bool use_real=false);
  void AddK(B2EInt eri, const Eigen::VectorXcd& Ca, Irrep ir_a, dcomplex coef, BMat& K, bool use_real=false);
  void CalcJ(B2EInt eri, const Eigen::VectorXcd& Ca, Irrep ir_a, BMat *J);
  void CalcK(B2EInt eri, const Eigen::VectorXcd& Ca, Irrep ir_a, BMat *K);
  MO CalcRHF(SymGTOs gtos, int nele, ERIMethod method, const SCFOptions& opts, 
	     bool *is_conv);
  MO CalcRHF(SymmetryGroup sym, BMatSet mat_set, B2EInt eri, int nele, 
	     const BVec& E0, const BMat& C0,const SCFOptions& opts,  bool *is_conv);
  void CalcSEHamiltonian(MO mo, B2EInt eri, Irrep I0, int i0, BMat* hmat,
			 int method = 0);
  dcomplex CalcAlpha(MO mo, BMatSet mat_set, Irrep I0, int i0, BMat& h_stex, double w, Coord coord, int method=0);
  double PITotalCrossSection(dcomplex alpha, double w, int num_occ_ele);  
}

#endif

