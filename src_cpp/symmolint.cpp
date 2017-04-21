#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <numeric>
#include <boost/foreach.hpp>
#include "../utils/eigen_plus.hpp"
#include "cfunc.hpp"
#include "angmoment.hpp"
#include "mol_func.hpp"
#include "one_int.hpp"
#include "two_int.hpp"
#include "symmolint.hpp"

namespace cbasis {

  using namespace std;
  using namespace Eigen;
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<Reduction>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<Reduction>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 1> A1dc;
  typedef MultArray<dcomplex, 2> A2dc;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;


  // ==== ERI method ====
  ERIMethod::ERIMethod(): symmetry(0), coef_R_memo(0), perm(0) {}
  void ERIMethod::set_symmetry(int s) {symmetry = s; }
  void ERIMethod::set_coef_R_memo(int s) {coef_R_memo = s; }
  void ERIMethod::set_perm(int s) {perm = s; }

  // ==== Reduction ====
  void Reduction::SetLM(int _L, int _M, dcomplex _coef_sh) {
     L=_L;
     M=_M;
     is_solid_sh=true;
     coef_sh = _coef_sh;     
  }
  string Reduction::str() const {
    ostringstream oss;
    oss << "==== Redcution ====" << endl;
    oss << "irrep : " << irrep << endl;
    oss << "coef_iat_ipn: " << endl;
    oss << coef_iat_ipn << endl;
    oss << "coef_iz: ";
    for(int iz = 0; iz < coef_iz.size(); iz++)
      oss << coef_iz[iz] << " ";
    oss << endl;
    oss << "offset:  " << offset << endl;
    oss << "===================" << endl;
    return oss.str();
  }
  void Reduction::Display() const {
    cout << this->str() ;
  }
  
  // ==== Sub ====
  SubSymGTOs::SubSymGTOs(SymmetryGroup _sym, Atom _atom) {
    sym_group_ = _sym;
    atom_ = _atom;
    zeta_iz = VectorXcd::Zero(0);
    setupq = false;
  }
  void SubSymGTOs::SetUp() {

    // ---- check values ----
    if(sym_group().get() == NULL) {
      string msg; SUB_LOCATION(msg);
      msg += ": sym_group is not set.";
      throw runtime_error(msg);
    }

    if(this->size_at() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": size_at=0";
      throw runtime_error(msg);
    }

    if( this->size_pn() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": size_pn==0";
      throw runtime_error(msg);
    }

    if(this->zeta_iz.size() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": zeta is not set.";
      throw runtime_error(msg);
    }

    if(this->rds.size() == 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": ReductionSet is not set.";
      throw runtime_error(msg);
    }

    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      if(it->size_at() != this->size_at() ||
	 it->size_pn() != this->size_pn()) {
	string msg; SUB_LOCATION(msg);
	ostringstream oss; oss << msg;
	oss << ": size mismatch."
	    << endl
	    << "size_at (Reduction, xyz) = "
	    << "(" << it->size_at() << ", " << this->size_at() << ")"
	    << endl
	    << "size_pn (Reduction, ns) = "
	    << "(" << it->size_pn() << "," << this->size_pn() << ")"
	    << endl << endl;
	oss << "object:" << this->str() << endl;
	throw runtime_error(oss.str());
      }
    }

    // ---- compute internal values ----
    ip_iat_ipn = MatrixXi::Zero(this->size_at(), this->size_pn());
    int ip(0);
    for(int iat = 0; iat < this->size_at(); ++iat)
      for(int ipn = 0; ipn <  this->size_pn(); ++ipn) {
	ip_iat_ipn(iat, ipn) = ip;
	ip++;
      }
    
    maxn = 0;
    for(int ipn = 0; ipn < size_pn(); ipn++) {
      if(maxn < nx(ipn) + ny(ipn) + nz(ipn))
	maxn = nx(ipn) + ny(ipn) + nz(ipn);
    }
    for(RdsIt it = rds.begin(); it != rds.end(); ++it) {
      it->set_zs_size(zeta_iz.size());
    }

    // -- build PrimGTO set.
    int nat(this->size_at());
    int npn(this->size_pn());
    vector<PrimGTO> gtos(nat * npn);
    for(int iat = 0; iat < nat; iat++) 
      for(int ipn = 0; ipn < npn; ipn++)  {
	int ip = this->ip_iat_ipn(iat, ipn);
	gtos[ip] = PrimGTO(this->nx(ipn),
			   this->ny(ipn),
			   this->nz(ipn),
			   this->x(iat),
			   this->y(iat),
			   this->z(iat));
      }
    // -- compute symmetry operation and primitive GTO relations --
    sym_group()->CalcSymMatrix(gtos, this->ip_jg_kp, this->sign_ip_jg_kp);

    // -- set flag --
    setupq = true;
  }
  SubSymGTOs& SubSymGTOs::AddNs(int nx, int ny, int nz) {

    setupq = false;
    this->nx_ipn.push_back(nx);
    this->ny_ipn.push_back(ny);
    this->nz_ipn.push_back(nz);
    return *this;
  }
  SubSymGTOs& SubSymGTOs::AddNs(Vector3i ns) {
    
    this->AddNs(ns[0], ns[1], ns[2]);
    return *this;
    
  }
  SubSymGTOs& SubSymGTOs::AddZeta(const VectorXcd& zs) {
    setupq = false;
    VectorXcd res(zeta_iz.size() + zs.size());
    for(int i = 0; i < zeta_iz.size(); i++)
      res(i) = zeta_iz(i);
    for(int i = 0; i < zs.size(); i++)
      res(i+zeta_iz.size()) = zs(i);
    zeta_iz.swap(res);
    return *this;
  }
  SubSymGTOs& SubSymGTOs::AddRds(const Reduction& _rds) {
    setupq = false;
    rds.push_back(_rds);
    return *this;
  }
  int SubSymGTOs::size_at() const {
    if(this->atom_)
      return this->atom_->size();
    else
      return 0;
  }
  void SubSymGTOs::Mono(Irrep irrep, Vector3i ns, VectorXcd zs) {

    if(this->atom_->size() != 1) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + ": atom have to have only one xyz";
      throw runtime_error(msg);
    }

    this->AddNs(ns);
    this->AddZeta(zs);
    this->AddRds(Reduction(irrep, MatrixXcd::Ones(1, 1)));

  }
  void SubSymGTOs::SolidSH_Ms(int L, VectorXi _Ms, VectorXcd zs) {
    // -- memo --
    // from wiki
    // real form and complex form:
    // Ylm = sqrt(2)(-)^m Im[Yl^|m|]   (m<0)
    // Yl0 = Yl^0   (m=0)
    // Ylm = sqrt(2)(-)^m Re[Yl^m]     (m>0)
    // or do google "spherical harmonics table"
    
    // Us(r) = N exp[-ar^2] = N exp[-ar^2] Y00 sqrt(4pi)
    // upz(r) = N z exp[-ar^2] = N exp[-ar^2] Y10 sqrt(4pi/3)
    // ud_z2(r) = N (2z^2-x^2-y^2) exp[-ar^2] = ...sqrt(16pi/5)
    // ud_zx(r) = N zx exp[-ar^2] = ...sqrt(4pi/15)

    if(this->atom_->size() != 1) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + ": atom have to have only one xyz";
      throw runtime_error(msg);
    }
    
    this->AddZeta(zs);
    SymmetryGroup sym = this->sym_group();
    
    vector<int> Ms;
    for(int i = 0; i < _Ms.size(); i++) {
      int m = _Ms[i];
      if(abs(m) > 1) {
	string msg; SUB_LOCATION(msg); msg += ": now |m| > 1 is not implemented";
	throw runtime_error(msg);
      }
      Ms.push_back(m);
    }

    if(L < 0 || L > 3) {
      string msg; SUB_LOCATION(msg); msg += ": only L=0,1,2,3 is supported";
      throw runtime_error(msg);
    }

    bool find_0 = find(Ms.begin(), Ms.end(),   0) != Ms.end();
    bool find_p1 = find(Ms.begin(), Ms.end(), +1) != Ms.end();
    bool find_m1 = find(Ms.begin(), Ms.end(), -1) != Ms.end();

    if(L == 0) {
      this->AddNs(Vector3i::Zero());
      if(find_0) {
	Reduction rds(sym->irrep_s(), MatrixXcd::Ones(1, 1));
	rds.SetLM(0, 0, sqrt(4.0*M_PI)); this->AddRds(rds);
      }
    }
    if(L == 1) {
      this->AddNs(1,0,0); this->AddNs(0,1,0); this->AddNs(0,0,1);
      if(find_p1) {
	Reduction rds(sym->irrep_x(), m13cd(1,0,0));
	rds.SetLM(1, 1, sqrt(4.0 * M_PI / 3.0)); this->AddRds(rds);
      }      
      if(find_0) {
	Reduction rds(sym->irrep_z(), m13cd(0,0,1));
	rds.SetLM(1, 0, sqrt(4.0 * M_PI / 3.0)); this->AddRds(rds);
      }
      if(find_m1) {
	Reduction rds(sym->irrep_y(), m13cd(0,1,0));
	rds.SetLM(1, -1, sqrt(4.0 * M_PI / 3.0)); this->AddRds(rds);
      }
    }
    if(L == 2) {
      this->AddNs(2,0,0); this->AddNs(0,2,0); this->AddNs(0,0,2);
      this->AddNs(1,1,0); this->AddNs(0,1,1); this->AddNs(1,0,1);
      if(find_p1) {
	MatrixXcd m(1, 6); m << 0,0,0, 0,0,1;
	Reduction rds(sym->irrep_x(), m);
	rds.SetLM(2, 1, sqrt(4.0*M_PI/15.0)); this->AddRds(rds);
      }      
      if(find_0) {
	MatrixXcd m(1, 6); m << -1,-1,2, 0,0,0;
	Reduction rds(sym->irrep_z(), m);
	rds.SetLM(2, 0, sqrt(16.0*M_PI/5.0)); this->AddRds(rds);
      }
      if(find_m1) {
	MatrixXcd m(1, 6); m << 0,0,0, 0,1,0;
	Reduction rds(sym->irrep_y(), m);
	rds.SetLM(2, -1, sqrt(4.0*M_PI/15.0)); this->AddRds(rds);
      }
    }
    if(L == 3) {
      this->AddNs(1,0,2); this->AddNs(3,0,0); this->AddNs(1,2,0);
      this->AddNs(0,0,3); this->AddNs(2,0,1); this->AddNs(0,2,1);
      this->AddNs(0,1,2); this->AddNs(2,1,0); this->AddNs(0,3,0);
      if(find_p1) {
	MatrixXcd m(1,9); m << 4,-1,-1,  0,0,0,  0,0,0;
	Reduction rds(sym->irrep_x(), m);
	rds.SetLM(3, +1, sqrt(32.0*M_PI/21.0)); this->AddRds(rds);
      }
      if(find_0) {
	MatrixXcd m(1,9); m << 0,0,0,  2,-3,-3,  0,0,0;
	Reduction rds(sym->irrep_z(), m);
	rds.SetLM(3, 0, sqrt(16.0*M_PI/7.0)); this->AddRds(rds);
      }
      if(find_m1) {
	MatrixXcd m(1,9); m << 0,0,0,  0,0,0,  4,-1,-1;
	Reduction rds(sym->irrep_y(), m);
	rds.SetLM(3, -1, sqrt(32.0*M_PI/21.0)); this->AddRds(rds);
      }
    }
  }    
  void SubSymGTOs::SolidSH_M(int L, int M, Eigen::VectorXcd zs) {
    VectorXi Ms(1);
    Ms(0) = M;
    this->SolidSH_Ms(L, Ms, zs);
  }
  string SubSymGTOs::str() const {
    ostringstream oss;
    oss << "==== SubSymGTOs ====" << endl;
    //oss << "sym : " << sym_group->str() << endl;    
    for(int iat = 0; iat < size_at(); iat++)
      oss << "xyz" << iat << ": " << x(iat) << y(iat) << z(iat) << endl;
    for(int ipn = 0; ipn < size_pn(); ipn++) 
      oss << "ns" << ipn << ": " << nx(ipn) << ny(ipn) << nz(ipn) << endl;
    oss << "zeta: ";
    for(int iz = 0; iz < zeta_iz.size(); iz++)
      oss << zeta_iz[iz] << " ";
    oss << endl;
    oss << "maxn: " << maxn << endl;
    oss << "ip_iat_ipn: " << endl << ip_iat_ipn << endl;
    oss << "ip_jg_kp: " << endl << ip_jg_kp << endl;
    oss << "sign_ip_jg_kp: " << endl << sign_ip_jg_kp << endl;    
    for(cRdsIt it = rds.begin(); it != rds.end(); ++it)
      oss << it->str();
    oss << "====================" << endl;
    return oss.str();
  }

  // ==== SymGTOs ====
  // ---- Constructors ----
  _SymGTOs::_SymGTOs(Molecule mole):
    sym_group_(mole->sym_group()), molecule_(mole), setupq(false) {}

  // ---- Accessors ----
  int _SymGTOs::size_basis() const {
    int cumsum(0);
    for(cSubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {
      cumsum += isub->size_zeta() * isub->rds.size();
    }
    return cumsum;
  }
  int _SymGTOs::size_basis_isym(Irrep isym) const {

    int cumsum(0);
    for(cSubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {

      for(cRdsIt irds = isub->rds.begin(); irds != isub->rds.end(); 
	  ++irds) {

	if(irds->irrep == isym) 
	  cumsum += isub->size_zeta();
      }
    }
    return cumsum;

  }
  string _SymGTOs::str() const {
    ostringstream oss;
    oss << "==== SymGTOs ====" << endl;
    oss << "Set Up?" << (setupq ? "Yes" : "No") << endl;
    oss << "sym:";
    if(this->sym_group()) {
      oss << endl << this->sym_group()->str();
    } else {
      oss << "No" << endl;
    }
    
    for(cSubIt it = subs_.begin(); it != subs_.end(); ++it) {
      oss << "sub" << distance(subs_.begin(), it) << ": " << endl;
      oss << it->str();
    }
    oss << "molecule:";
    if(this->molecule()) {
      oss << endl << this->molecule()->str();
    } else {
      oss << "No";
    }
    oss << endl;
    oss << "=================" << endl;
    return oss.str();
  }
  string _SymGTOs::show() {
    
    string sep = " | ";
    string line = "-|-";
    ostringstream oss;

    ostringstream lines;
    lines << setfill('-')
	  << line << setw(4) << right << ""
	  << line << setw(4) << right << ""
	  << line << setw(5) << right << ""
	  << line << setw(3) << right << ""
	  << line << setw(20) << right << ""
	  << line << endl;
    oss << lines.str();
    
    oss << setfill(' ')
	<< sep << setw(4) << right << "atom"
	<< sep << setw(4) << right << "pn"
	<< sep << setw(5) << right << "irrep"
	<< sep << setw(3) << right << "idx"
	<< sep << setw(20)<< right << "zeta"
	<< sep << endl;
    oss << lines.str();
    
    BOOST_FOREACH(SubSymGTOs& sub, this->subs_) {
      

      for(RdsIt irds = sub.rds.begin(); irds != sub.rds.end(); ++irds) {

	ostringstream pn_oss;
	if(sub.size_pn() == 1) {
	  pn_oss << sub.nx(0)<<sub.ny(0)<<sub.nz(0);
	} else {
	  pn_oss << sub.maxn;
	}
	oss << sep << setw(4) << right << sub.atom()->name()
	    << sep << setw(4) << right << pn_oss.str()
	    << sep << setw(5) << right << sym_group_->GetIrrepName(irds->irrep)
	    << sep << setw(3) << right << ""
	    << sep << setw(20) << right << ""
	    << sep << endl;

	for(int iz = 0; iz < sub.size_zeta(); iz++) {
	  oss << sep << setw(4) << right << ""
	      << sep << setw(4) << right << ""
	      << sep << setw(5) << right << ""
	      << sep << setw(3) << right << irds->offset+iz
	      << sep << setw(20) << right <<  sub.zeta(iz)
	      << sep << endl;
	}
      }
    }
    oss << lines.str();
    return oss.str();
  }
  
  // ---- Other ----
  SymGTOs _SymGTOs::Clone() const {

    SymGTOs gtos = NewSymGTOs(this->molecule());
    gtos->subs_ = this->subs_;
    gtos->SetUp();
    return gtos;

  }
  SymGTOs _SymGTOs::Conj() const {

    SymGTOs gtos = this->Clone();
    
    for(SubIt isub = gtos->subs_.begin(); isub != gtos->subs_.end(); ++isub) {
      isub->zeta_iz = isub->zeta_iz.conjugate();
      for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
	irds->coef_iat_ipn = irds->coef_iat_ipn.conjugate();
	irds->coef_iz = irds->coef_iz.conjugate();
      }
    }
    gtos->SetUp();
    return gtos;
  }

  // ---- SetUp ----
  _SymGTOs* _SymGTOs::AddSub(SubSymGTOs sub) {
    subs_.push_back(sub);
    return this;
  }
  SubSymGTOs& _SymGTOs::NewSub(Atom atom) {
    SubSymGTOs new_sub(this->sym_group(), atom);
    subs_.push_back(new_sub);
    return subs_[subs_.size()-1];
  }
  SubSymGTOs& _SymGTOs::NewSub(std::string atom_name) {
    return this->NewSub(this->molecule()->atom(atom_name));
  }
  _SymGTOs* _SymGTOs::SetUp() {
    for(SubIt it = subs_.begin(); it != subs_.end(); ++it) {
      
      // -- setup sub --
      it->sym_group_ = this->sym_group();
      if(it->setupq == false) {
	
	try {
	  it->SetUp();
	} catch(exception& e) {
	  ostringstream oss;
	  string loc; SUB_LOCATION(loc);
	  oss << endl << loc << ": error in SetUp of SubSymGTOs" << endl;
	  oss << "error on " << distance(subs().begin(), it) << " th sub" << endl;
	  oss << e.what();
	  throw runtime_error(oss.str());
	}
      }

      // -- check each sub --
      if(it->ip_jg_kp.rows() != sym_group_->order()) {
	string msg; SUB_LOCATION(msg); 
	msg += ": size of symmetry transformation matrix and num_class of symmetry group is not matched.";
	cout << "ip_jg_kp:" << endl;
	cout << it->ip_jg_kp << endl;
	throw runtime_error(msg);
      }

      // -- init offset --
      for(RdsIt irds = it->begin_rds(); irds != it->end_rds(); ++irds) {
	sym_group_->CheckIrrep(irds->irrep);
	irds->offset = 0;
      }
    }
    try {
      this->SetOffset();
    } catch(exception& e) {
      string loc; SUB_LOCATION(loc);
      string msg = "\n" + loc + ": error in SetOffset\n" + e.what();
      throw runtime_error(msg);
    }
    try {
      this->Normalize();
    } catch(exception& e) {
      string loc; SUB_LOCATION(loc);
      string msg = "\n" + loc + ": error in Normalize\n" + e.what();
      throw runtime_error(msg);      
    }
    setupq = true;
    return this;
  }
  void _SymGTOs::SetOffset() {
    
    vector<int> num_irrep(10, 0);
    
    for(SubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(), end = isub->rds.end();
	  irds != end; ++irds) {
	irds->offset = num_irrep[irds->irrep];
	num_irrep[irds->irrep] += isub->size_zeta();
      }
    }

  }
  void _SymGTOs::Normalize() {

    /*
  A3dc dxmap(100), dymap(100), dzmap(100);
    for(SubIt isub = subs.begin(), end = subs.end(); isub != end; ++isub) {
      for(int iz = 0; iz < isub->size_zeta(); iz++) {

	dcomplex zetai = isub->zeta_iz[iz];
	dcomplex zetaP = zetai + zetai;

	int mi = isub->maxn;
	calc_d_coef(mi,mi,0, zetaP, 0.0, 0.0, 0.0, dxmap);
	calc_d_coef(mi,mi,0, zetaP, 0.0, 0.0, 0.0, dymap);
	calc_d_coef(mi,mi,0, zetaP, 0.0, 0.0, 0.0, dzmap);

	int nipn(isub->size_pn());
	for(int ipn = 0; ipn < nipn; ipn++) {
	  int nxi, nyi, nzi;
	  nxi = isub->nx(ipn);
	  nyi = isub->ny(ipn);
	  nzi = isub->nz(ipn);
	  dcomplex dx00, dy00, dz00;
	  dx00 = dxmap(nxi, nxi ,0);
	  dy00 = dymap(nyi, nyi ,0);
	  dz00 = dzmap(nzi, nzi ,0);
	  dcomplex s_ele = dx00 * dy00 * dz00;
	  dcomplex ce = pow(M_PI/zetaP, 1.5);	  
	  for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
	    irds->coef_iz(iz) = 1.0/sqrt(ce*s_ele);
	  }
	}
      }
    }
    */

    A3dc dxmap(100), dymap(100), dzmap(100);
    // >>> Irrep Adapted GTOs >>>
    for(SubIt isub = subs_.begin(), end = subs_.end(); isub != end; ++isub) {
      for(int iz = 0; iz < isub->size_zeta(); iz++) {
	vector<Reduction>& rds(isub->rds);
	for(RdsIt irds = rds.begin(); irds != rds.end(); ++irds) {
	  dcomplex zetai = isub->zeta_iz[iz];
	  dcomplex zetaP = zetai + zetai;
	  int niat(isub->size_at()); int nipn(isub->size_pn());
	  dcomplex norm2(0.0);
	  

	  // >>> Primitive GTOs >>>
	  for(int iat = 0; iat < niat; iat++) {
	    for(int jat = 0; jat < niat; jat++) {
	      dcomplex xi, xj, yi, yj, zi, zj, wPx, wPy, wPz;
	      xi = isub->x(iat); xj = isub->x(jat);
	      yi = isub->y(iat); yj = isub->y(jat);
	      zi = isub->z(iat); zj = isub->z(jat);
	      wPx = (zetai*xi+zetai*xj)/zetaP;
	      wPy = (zetai*yi+zetai*yj)/zetaP;
	      wPz = (zetai*zi+zetai*zj)/zetaP;
	      dcomplex d2 = pow(xi-xj,2) + pow(yi-yj,2) + pow(zi-zj,2);	
	      dcomplex eAB = exp(-zetai*zetai/zetaP*d2);
	      dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);
	      int mi = isub->maxn;
	      calc_d_coef(mi,mi,0, zetaP,wPx,xi,xj,dxmap);
	      calc_d_coef(mi,mi,0, zetaP,wPy,yi,yj,dymap);
	      calc_d_coef(mi,mi,0, zetaP,wPz,zi,zj,dzmap);
	      
	      for(int ipn = 0; ipn < nipn; ipn++) {
		for(int jpn = 0; jpn < nipn; jpn++) {
		  int nxi, nxj, nyi, nyj, nzi, nzj;
		  nxi = isub->nx(ipn); nxj = isub->nx(jpn);
		  nyi = isub->ny(ipn); nyj = isub->ny(jpn);
		  nzi = isub->nz(ipn); nzj = isub->nz(jpn); 
		  dcomplex dx00, dy00, dz00;
		  dx00 = dxmap(nxi, nxj ,0);
		  dy00 = dymap(nyi, nyj ,0);
		  dz00 = dzmap(nzi, nzj ,0);
		  dcomplex s_ele = dx00 * dy00 * dz00;
		  norm2 += ce * s_ele *
		    irds->coef_iat_ipn(iat, ipn) *
		    irds->coef_iat_ipn(jat, jpn);
		}
	      }
	    }
	  }

	  if(abs(norm2) < pow(10.0, -14.0)) {
	    string msg; SUB_LOCATION(msg);
	    ostringstream oss; oss << msg << ": " << endl;
	    oss << "norm is too small" << endl;
	    oss << "isub:"  << distance(subs_.begin(), isub) << endl;
	    oss << "iz  : " << iz << endl;
	    oss << "zeta: " << zetai << endl;
	    oss << "irds:" << distance(isub->rds.begin(), irds) << endl;
	    oss << "nipn: " << nipn << endl;
	    oss << "niat: " << niat << endl;
	    oss << "norm2: " << norm2 << endl;	    
	    throw runtime_error(oss.str());
	  }

	  // <<< Primitive GTOs <<<
	  irds->coef_iz(iz) = 1.0/sqrt(norm2);
	}
      }
    }
    // <<< Irrep Adapted GTOs <<<

  }

  // ---- utils ----
  // -- not used now --
  int _SymGTOs::max_n() const {
    int max_n(0);
    for(cSubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {
      for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
	int nx(isub->nx(ipn));
	int ny(isub->ny(ipn));
	int nz(isub->nz(ipn));
	if(max_n < nx)
	  max_n = nx;
	if(max_n < ny)
	  max_n = ny;
	if(max_n < nz)
	  max_n = nz;
      }
    }    
    return max_n;
  }
  
  // ---- CalcMat ---- 
  void _SymGTOs::loop() {

    for(SubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {
      for(SubIt jsub = subs_.begin(); jsub != subs_.end(); ++jsub) {

	// -- loop over each zeta --
	for(int iz = 0; iz < isub->size_zeta(); iz++) {
	  for(int jz = 0; jz < jsub->size_zeta(); jz++) {
	    dcomplex zetai, zetaj;
	    zetai = isub->zeta_iz[iz]; zetaj = jsub->zeta_iz[jz];
	    int niat(isub->size_at()); int njat(jsub->size_at());
	    int nipn(isub->size_pn()); int njpn(jsub->size_pn());

	    // -- primitive basis --
	    for(int iat = 0; iat < niat; iat++) {
	      for(int jat = 0; jat < njat; jat++) { 
		for(int ipn = 0; ipn < nipn; ipn++) {
		  for(int jpn = 0; jpn < njpn; jpn++) {

		  }}}}

	    // -- contractions --
	    for(RdsIt irds = isub->rds.begin();irds != isub->rds.end();++irds) {
	      for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end();++jrds) {
	      }}

	  }}
      }}
    }
  int _SymGTOs::max_num_prim() const {
    
    int max_num_prim(0);
    for(cSubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {
      max_num_prim = std::max(max_num_prim, isub->size_prim());
    }
    return max_num_prim;

  }

  // ---- AtR ----
  bool IsCenter(SubIt isub, double eps) {

    if(isub->size_at() != 1)
      return false;
    
    dcomplex x  = isub->x(0);
    dcomplex y  = isub->y(0);
    dcomplex z  = isub->z(0);
    dcomplex a2 = x*x+y*y+z*z;
    dcomplex a = sqrt(a2);
    return (abs(a) < eps);
  }
  void AtR_Ylm_cen(RdsIt irds, dcomplex z,
		   dcomplex r, int L, int M, dcomplex* v, dcomplex* dv) {

    if(!irds->is_solid_sh) {
      string msg; SUB_LOCATION(msg);
      msg += ": it is not Solid spherical harmonics";
      throw runtime_error(msg);
    }

    if(L > 3 || L < 0  || abs(M) > L) {
      string msg; SUB_LOCATION(msg);
      msg += ": not supported (L,M)";
      throw runtime_error(msg);
    }
    
    if(irds->L == L && irds->M == M) {
      *v = irds->coef_sh  * pow(r, L+1) * exp(-z*r*r);
      *dv = irds->coef_sh * ((L+1.0)*pow(r, L) -2.0*z*pow(r,L+2)) * exp(-z*r*r);
    } else {
      *v = 0.0;
      *dv=0.0;
    }
    
  }
  void AtR_Ylm_noncen(SubIt isub, int iz, int iat, int ipn,
		      dcomplex r, int L, int M, dcomplex* v, dcomplex *dv) {

    dcomplex* il_ipl  = new dcomplex[2*L+2];
    dcomplex* il = &il_ipl[0];
    dcomplex* ipl= &il_ipl[L+1];
    dcomplex* ylm = new dcomplex[num_lm_pair(L)];

    dcomplex res;

    int nn = (isub->nx(ipn) +
	      isub->ny(ipn) +
	      isub->nz(ipn));
    dcomplex zeta = isub->zeta_iz(iz);
    dcomplex x  = isub->x(iat);
    dcomplex y  = isub->y(iat);
    dcomplex z  = isub->z(iat);
    dcomplex xxyy = x*x+y*y;
    dcomplex a2 = xxyy+z*z;
    dcomplex a = sqrt(a2);

    if(abs(a) < 0.000001) {
      string msg; SUB_LOCATION(msg);
      msg += ": (x,y,z) is centered on origin."; 
      throw runtime_error(msg);
    }

    dcomplex expz = exp(-zeta*(r*r+a2));
    dcomplex theta = cacos(z / a);
    dcomplex phi   = (abs(xxyy) < 0.00001) ? 0.0 : cacos(x / sqrt(xxyy));
    ModSphericalBessel(2.0*zeta*a*r, L, il_ipl);
    RealSphericalHarmonics(theta, phi, L, ylm);
    if(nn == 0) {
      
      dcomplex c((4.0*M_PI) * pow(-1.0, M) * ylm[lm_index(L, -M)]);
      *v = c * r * il[L] * expz;
      *dv = (c*   il[L]*expz +
	     c*r*ipl[L]*expz*(+2.0*zeta*a) +
	     c*r* il[L]*expz*(-2.0*zeta*r));
    } else {
      string msg; SUB_LOCATION(msg);
      msg += ": not implemented yet for p or higher orbital";
      throw runtime_error(msg);
    }		  

    delete[] il_ipl;
    delete[] ylm;
    
  }      
  void _SymGTOs::AtR_Ylm(int L, int M, int irrep, const VectorXcd& cs_ibasis,
			const VectorXcd& rs,
			VectorXcd* res_vs, VectorXcd* res_dvs ) {

    if(not setupq) 
      this->SetUp();

    if(!is_lm_pair(L, M)) {
      string msg; SUB_LOCATION(msg);
      msg += ": invalid L,M pair";
      throw runtime_error(msg);
    }

    try{
      sym_group_->CheckIrrep(irrep);
    } catch(const runtime_error& e) {
      string msg; SUB_LOCATION(msg);
      msg += ": Invalid irreq.\n";
      msg += e.what(); throw runtime_error(msg);
    }

    if(cs_ibasis.size() != this->size_basis_isym(irrep)) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << endl << ": size of cs must be equal to basis size" << endl;
      oss << "(L, M, irrep)     = " << L << M << irrep << endl;
      oss << "csibasis.size()   = " << cs_ibasis.size() << endl;
      oss << "size_basis_isym() = " << this->size_basis_isym(irrep) << endl;
      oss << endl;
      oss << "print SymGTOs:" << endl;
      oss << this->str();
      
      throw runtime_error(oss.str());
    }

    VectorXcd vs  = VectorXcd::Zero(rs.size());   // copy
    VectorXcd dvs = VectorXcd::Zero(rs.size());   // copy
    dcomplex* ylm = new dcomplex[num_lm_pair(L)];
    dcomplex* il  = new dcomplex[2*L+2];
    double eps(0.0000001);

    // Y00   = 1/sqrt(4pi)
    // r Y10 = sqrt(3/4pi) z
    for(SubIt isub = subs_.begin(); isub != subs_.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {
	if(irds->irrep != irrep) 
	  continue;
	if(IsCenter(isub, eps) && irds->is_solid_sh) {
	  if(irds->L == L && irds->M == M) {
	    for(int ir = 0; ir < rs.size(); ir++) {
	      for(int iz = 0; iz < isub->size_zeta(); iz++) {
		int ibasis = irds->offset + iz;
		dcomplex r = rs[ir];
		dcomplex z = isub->zeta_iz(iz);
		dcomplex expo=exp(-z*r*r);
		dcomplex v, dv, c;
		c = cs_ibasis(ibasis) * irds->coef_iz(iz) * irds->coef_sh;
		vs[ir] += c * pow(r, L+1) * expo;
		dvs[ir]+= c * ((L+1.0)*pow(r, L) -2.0*z*pow(r,L+2)) * expo;
	      }
	    }
	  }
	} else {
	  for(int iat = 0; iat < isub->size_at(); iat++) {
	    for(int ipn = 0; ipn < isub->size_pn(); ipn++) {
	      for(int ir = 0; ir < rs.size(); ir++) {
		for(int iz = 0; iz < isub->size_zeta(); iz++) {
		  int ibasis = irds->offset + iz;
		  dcomplex c = (cs_ibasis(ibasis) *
				irds->coef_iat_ipn(iat, ipn) *
				irds->coef_iz(iz));
		  dcomplex v, dv;
		  AtR_Ylm_noncen(isub, iz, iat, ipn, rs[ir], L, M, &v, &dv);
		  vs[ir] += c * v;
		  dvs[ir]+= c * dv;
		}
	      }
	    }
	  }
	}
      }
    }
    delete[] ylm;
    delete[] il;
    res_vs->swap(vs);
    res_dvs->swap(dvs);
  }
  void _SymGTOs::AtR_YlmOne(int L, int M,  int irrep,
			    const Eigen::VectorXcd& cs_ibasis,
			    dcomplex r, dcomplex *v, dcomplex *dv) {
    VectorXcd rs(1); rs << r;
    VectorXcd vs, dvs;
    this->AtR_Ylm(L, M, irrep, cs_ibasis, rs, &vs, &dvs);
    *v = vs[0];
    *dv = dvs[0];
  }
  // ---- Correct Sign ----
  void _SymGTOs::CorrectSign(int L, int M, int irrep, Eigen::VectorXcd& cs) {

    if(not setupq) {
      string msg; SUB_LOCATION(msg);
      msg += ": call SetUp() before calculation.";
      throw runtime_error(msg);
    }

    if(!is_lm_pair(L, M)) {
      string msg; SUB_LOCATION(msg);
      msg += ": invalid L,M pair";
      throw runtime_error(msg);
    }

    try{
      sym_group_->CheckIrrep(irrep);
    } catch(const runtime_error& e) {
      string msg; SUB_LOCATION(msg);
      msg += ": Invalid irreq.\n";
      msg += e.what(); throw runtime_error(msg);
    }

    if(cs.size() != this->size_basis_isym(irrep)) {
      string msg; SUB_LOCATION(msg);
      msg += ": size of cs must be equal to basis size";
      throw runtime_error(msg);
    }

    double r = 0.3;
    dcomplex v, dv;
    this->AtR_YlmOne(L, M, irrep, cs, r, &v, &dv);

    if(v.real() < 0) {
      cs = -cs;
    }
  }
  
  // ---- Calculation ----
  void _SymGTOs::InitBVec(BVec *ptr_bvec) {
    throw(runtime_error("dont use _SymGTOs::InitBVec"));
    BVec& bvec = *ptr_bvec;
    for(SubIt isub = this->subs_.begin(); isub != this->subs_.end(); ++isub) {
      for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {

	bool need_update = false;
	Irrep irrep = irds->irrep;
	if(not bvec.has_block(irrep))
	  need_update = true;
	else if(bvec[irrep].size() != this->size_basis_isym(irrep)) 
	  need_update = true;

	if(need_update)
	  bvec[irrep] = VectorXcd::Zero(this->size_basis_isym(irrep));
      }
    }
  }
  void _SymGTOs::InitBMat(Irrep krrep, BMat *ptr_mat) {
    throw(runtime_error("dont use _SymGTOs::InitBMat"));
    BMat& mat = *ptr_mat;
    for(SubIt isub = this->subs_.begin(); isub != this->subs_.end(); ++isub) {
      for(SubIt jsub = this->subs_.begin(); jsub != this->subs_.end(); ++jsub) {
	for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds) {	
	  for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end(); ++jrds) {
	    Irrep irrep = irds->irrep;
	    Irrep jrrep = jrds->irrep;
	    int ni = this->size_basis_isym(irrep);
	    int nj = this->size_basis_isym(jrrep);
	    
	    if(this->sym_group_->Non0_3(irds->irrep, krrep, jrds->irrep)) {
	      pair<Irrep, Irrep> ijrrep(irrep, jrrep);
	      if(mat.has_block(ijrrep))
		mat[ijrrep] = MatrixXcd::Zero(ni, nj);
	      if(mat[ijrrep].rows() != ni || mat[ijrrep].cols() != nj) 
		mat[ijrrep] = MatrixXcd::Zero(ni, nj);
	    }
	  }
	}
      }
    }
  }

  SymGTOs NewSymGTOs(Molecule mole) { return SymGTOs(new _SymGTOs(mole)); }
  
}

