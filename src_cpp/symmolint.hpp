#ifndef SYMMOLINT_H
#define SYMMOLINT_H

#include <string>
#include <vector>
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include "../utils/typedef.hpp"
#include "bmatset.hpp"
#include "symgroup.hpp"
#include "molecule.hpp"
#include "b2eint.hpp"

namespace cbasis {
  
  // ==== ERI method ====
  class ERIMethod {
  public:
    int symmetry;
    int coef_R_memo;
    int perm;
    ERIMethod();
    void set_symmetry(int s);
    void set_coef_R_memo(int s);
    void set_perm(int s);
  };

  // ==== AO Reduction ====
  struct Reduction {
    Irrep irrep;
    Eigen::MatrixXcd coef_iat_ipn;
    bool is_solid_sh; // true=> this Reduction represent solid spherical harmonics
    dcomplex coef_sh; // coefficient for Solid Spherical Harmonics
    int L;    
    int M;
    
    // ---- for calculation  ----
    Eigen::VectorXcd coef_iz; // normalization constant for each iz
    int offset;

    // ---- constructor ----
    Reduction(int _irrep, Eigen::MatrixXcd _coef):
      irrep(_irrep), coef_iat_ipn(_coef), is_solid_sh(false),
      coef_sh(0), L(-1), M(0),offset(0) {}

    void SetLM(int _L, int _M, dcomplex _coef_sh);
    std::string str() const;
    void Display() const;
    
    const Eigen::MatrixXcd& get_coef_iat_ipn() const { return coef_iat_ipn; }
    inline int size_at() const { return coef_iat_ipn.rows(); }
    inline int size_pn() const { return coef_iat_ipn.cols(); }
    void set_zs_size(int num_zs) {
      coef_iz = Eigen::VectorXcd::Zero(num_zs);
    }
  };

  // ==== Sub sets of SymGTOs ====
  struct SubSymGTOs {
    // ---- type ----
    typedef std::vector<Reduction>::const_iterator cRdsIt;
    typedef std::vector<Reduction>::iterator RdsIt;

    // ---- Calculation data ----
    SymmetryGroup sym_group_;
    Atom atom_;
    std::vector<int> nx_ipn;
    std::vector<int> ny_ipn;
    std::vector<int> nz_ipn;
    Eigen::VectorXcd zeta_iz;
    std::vector<Reduction> rds;

    // ---- for calculation  ----
    Eigen::MatrixXi ip_iat_ipn;
    Eigen::MatrixXi ip_jg_kp;       // GTO[j] = G[i][GTO[k]]
    Eigen::MatrixXi sign_ip_jg_kp; // sign of above relation
    
    bool setupq;
    int maxn;
    //    int maxnx;

    // ---- Constructors ----
    SubSymGTOs(SymmetryGroup _sym, Atom _atom);

    // ---- Accessors ----
    std::string str() const;
    inline SymmetryGroup sym_group() { return sym_group_; }
    inline Atom atom() { return atom_; }
    inline int nx(int ipn) const { return nx_ipn[ipn]; }
    inline int ny(int ipn) const { return ny_ipn[ipn]; }
    inline int nz(int ipn) const { return nz_ipn[ipn]; }
    inline dcomplex x(int iat) const { return atom_->xyz_list()[iat][0]; }
    inline dcomplex y(int iat) const { return atom_->xyz_list()[iat][1]; }
    inline dcomplex z(int iat) const { return atom_->xyz_list()[iat][2]; }    
    inline dcomplex zeta(int iz) const { return zeta_iz[iz]; }
    inline cRdsIt begin_rds() const { return rds.begin(); }
    inline cRdsIt end_rds() const { return rds.end(); }
    inline RdsIt begin_rds() { return rds.begin(); }
    inline RdsIt end_rds() { return rds.end(); }
    //    inline void SetSym(pSymmetryGroup _sym_group) {
    //      sym_group = _sym_group;
    //    }
   
    SubSymGTOs& AddNs(int nx, int ny, int nz);
    SubSymGTOs& AddNs(Eigen::Vector3i ns);
    SubSymGTOs& AddZeta(const Eigen::VectorXcd& zs);
    SubSymGTOs& AddRds(const Reduction& rds);

    void Mono(Irrep irrep, Eigen::Vector3i ns, Eigen::VectorXcd zs);
    void SolidSH_Ms(int L, Eigen::VectorXi _Ms, Eigen::VectorXcd zs);
    void SolidSH_M(int L, int M, Eigen::VectorXcd zs);
    
    int size_at() const;// { return this->atom_->size(); }
    inline int size_pn() const { return nx_ipn.size(); }
    inline int size_rds() const { return rds.size(); }
    inline int size_prim() const { return this->size_at() * this->size_pn(); }
    inline int size_zeta() const { return zeta_iz.rows(); }

    // ---- SetUp ----   
    // -- calculate inner information and check values.
    void SetUp();

  };

  // ==== SymGTOs ====
  class _SymGTOs;
  typedef boost::shared_ptr<_SymGTOs> SymGTOs;
  class _SymGTOs {
  public:
    SymmetryGroup sym_group_;
    std::vector<SubSymGTOs> subs_;
    Molecule molecule_;
    bool setupq;
  public:  
    // ---- Constructors ----
    _SymGTOs(Molecule mole);
    //_SymGTOs(pSymmetryGroup _sym_group);

    // ---- Accessors ----    
    int size_atom() const;
    int size_basis() const;
    int size_basis_isym(Irrep isym) const;
    int size_subs() const {return subs_.size(); }
    std::string str() const;
    std::string show();
    int max_num_prim() const;
    SymmetryGroup sym_group() { return sym_group_; }
    SymmetryGroup sym_group() const { return sym_group_; }
    void set_sym_group(SymmetryGroup sym) { sym_group_ = sym; }
    Molecule molecule() { return molecule_; }
    Molecule molecule() const { return molecule_; }
    void set_molecule(Molecule mole) { molecule_ = mole; }
    std::vector<SubSymGTOs>& subs() { return subs_; };
    const std::vector<SubSymGTOs>& subs() const { return subs_; };
    SubSymGTOs& sub(int i) { return subs_[i]; };
    
    // ---- Other basis ----
    SymGTOs Clone() const;
    SymGTOs Conj() const;

    // ---- SetUp ----
    _SymGTOs* AddSub(SubSymGTOs sub);
    SubSymGTOs& NewSub(Atom atom);
    SubSymGTOs& NewSub(std::string atom_name);
    _SymGTOs* SetUp();
    
  private:
    void SetOffset();
    void Normalize();

  public:
    // ---- Calculation ----
    void InitBVec(BVec *vec);
    void InitBMat(Irrep irrep, BMat *vec);
    // -- to be removed --
    void loop();
    // -- not uesd now --
    int max_n() const;
    // -- Radial wave function --
    void AtR_Ylm(int L, int M,  int irrep,
		 const Eigen::VectorXcd& cs_ibasis,
		 const Eigen::VectorXcd& rs,
		 Eigen::VectorXcd* vs,
		 Eigen::VectorXcd* dvs);
    void AtR_YlmOne(int L, int M,  int irrep,
		    const Eigen::VectorXcd& cs_ibasis,
		    dcomplex r, dcomplex *v, dcomplex *dv);
    // -- Correction of wave function sign --
    void CorrectSign(int L, int M, int irrep, Eigen::VectorXcd& cs);
  };
  SymGTOs NewSymGTOs(Molecule mole);
  
  // ==== Add Sub ====
  //  SubSymGTOs Sub_Mono(SymmetryGroup sym, Atom atom, Irrep irrep, 
  //			 Eigen::Vector3i ns, Eigen::VectorXcd zs);
  //  SubSymGTOs Sub_SolidSH_Ms(SymmetryGroup sym, Atom atom, int L,
  //			       Eigen::VectorXi _Ms, Eigen::VectorXcd zs);
  //  SubSymGTOs Sub_SolidSH_M(SymmetryGroup sym, Atom atom, int L, int M, 
  //			      Eigen::VectorXcd zs);
  
}

#endif
