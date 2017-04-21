#include <iostream>
#include <algorithm>
#include <Eigen/Dense>
#include <iomanip>

#include "macros.hpp"
#include "op_driv.hpp"
#include "r1gtoint.hpp"
#include "l_algebra.hpp"
#include "eigen_plus.hpp"

#include "opt_alpha.hpp"


namespace cbasis {
  
  using namespace Eigen;
  using namespace std;

  // ==== calculation for driv ===
  dcomplex CalcAlpha(IDriv* driv, IOp* op, const R1GTOs& gs) {
    VectorXcd m;
    MatrixXcd L;
    driv->Calc(gs, m);
    op->Calc(gs, L);

    return (m.array() * L.fullPivLu().solve(m).array()).sum();
  }
  void SolveAlpha(IDriv* driv, IOp* op, const R1GTOs& gs, Eigen::VectorXcd& c) {
    VectorXcd m;
    MatrixXcd L;
    driv->Calc(gs, m);
    op->Calc(gs, L);
    //    int n = gs.size_basis();
    c = L.fullPivLu().solve(m);
    
  }
  
  // ==== optimization result ====
  OptResult::OptResult(int n):
    conv_q(false),zs(n), val(0.0), dz(n), grad(n), hess(n, n) {}  
  ostream& operator<<(ostream& out, const OptResult& a) {
    out << "conv_q: " << (a.conv_q ? "Yes" : "No") << endl;
    if(a.zs.size() < 3) {
      out << "zs    : ";
      for(int i = 0; i < a.zs.size(); i++)
	out << a.zs[i] << " ";
      out << endl;
    } else {
      out << "zs   : " << a.zs[0] << " " << a.zs[1] << " ... " << a.zs[a.zs.size()-1];
      out << endl;
    }
    return out;
  }

  // ==== optimizer ====
  // ---- Interface ----
  IOptimizer::IOptimizer(int _max_iter, double _eps, IOptTarget* _target, int lvl):
    max_iter(_max_iter), eps(_eps), target(_target), debug_lvl(lvl) {

    if(max_iter <= 0) {
      string msg; SUB_LOCATION(msg);
      msg += ": max_iter must be positive integer";
      throw runtime_error(msg);            
    }

    if(eps < 0.0) {
      string msg; SUB_LOCATION(msg);
      msg += ": eps must be positive";
      throw runtime_error(msg);      
    }

    if(target == NULL) {
      string msg; SUB_LOCATION(msg);
      msg += ": target is null";
      throw runtime_error(msg);
    }

  }
  IOptimizer::~IOptimizer() {}
  OptResult* IOptimizer::Optimize(const Eigen::VectorXcd& z0s) {

    // -- Prepare results object --
    OptResult* p_res = new OptResult(target->dim());
    p_res->zs = z0s;

    for(int i = 0; i < max_iter; i++) {
      target->ValGradHess(p_res->zs, &p_res->val, p_res->grad, p_res->hess);
      this->OneStep(p_res->grad, p_res->hess, p_res->dz);
      if(debug_lvl > 0) {
	cout << "i: " << i << endl;
	cout << "zs: " << p_res->zs << endl;
	cout << "dz: " << p_res->dz << endl;
	cout << "grad: " << p_res->grad << endl;
	cout << "hess: " << p_res->hess << endl;
      }
      if(p_res->dz.array().abs().matrix().maxCoeff() < eps) {
	p_res->conv_q = true;
	break;
      }
      p_res->zs += p_res->dz;
    }
    return p_res;
  }
  OptResult* IOptimizer::Optimize(dcomplex z0) {
    
    if(target->dim() != 1) {
      string msg; SUB_LOCATION(msg);
      msg += ": Optimizae(dcomplex z0) can be used for dimension 1. ";
      throw runtime_error(msg);
    }
    VectorXcd zs(1);
    zs << z0;
    OptResult* p_res = this->Optimize(zs);
    return p_res;
  }

  // ---- External ----
  bool CheckOptTarget(IOptTarget* target, const VectorXcd& zs, double h, double eps) {
    
    int n(zs.size());
    dcomplex ii(0.0, 1.0);
    VectorXcd hvec = VectorXcd::Zero(n);
    VectorXcd ihvec = VectorXcd::Zero(n);
    hvec[0] = h;
    ihvec[0] = ii*h;

    dcomplex val;
    VectorXcd grad;
    MatrixXcd hess;
    target->ValGradHess(zs, &val, grad, hess);

    dcomplex ap0, am0, a0p, a0m;
    ap0 = target->Val(zs+hvec);
    am0 = target->Val(zs-hvec);
    a0p = target->Val(zs+ihvec);
    a0m = target->Val(zs-ihvec);

    dcomplex num_grad = (ap0 - am0 + ii*a0m - ii*a0p) / (4.0 * h);
    dcomplex num_hess = (ap0 + am0 - a0p - a0m) / (2.0*h*h);

    if(abs(num_grad - grad[0]) > eps) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss <<  msg;
      oss << ": test failed\n";
      oss << "grad[0]:  " << grad[0] << endl;
      oss << "num_grad: " << num_grad << endl;
      oss << "diff    : " << abs(num_grad-grad[0]) << endl;
      throw runtime_error(oss.str());
    }
    if(abs(num_hess - hess(0, 0)) > eps) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss <<  msg;
      oss << ": test failed\n";
      oss << "hess[0]:  " << hess(0, 0) << endl;
      oss << "num_hess: " << num_hess << endl;
      oss << "diff    : " << abs(hess(0, 0)-num_hess) << endl;
      throw runtime_error(oss.str());
    }
    return true;

  }
  bool CheckOptTarget(IOptTarget* target, dcomplex z, double h, double eps) {
    if(target->dim() != 1) {
      string msg; SUB_LOCATION(msg);
      msg += ": only one dim is supported.";
      throw runtime_error(msg);
    }
    VectorXcd zs(1); zs << z;
    try {
      CheckOptTarget(target, zs, h, eps);
    } catch(const exception& e) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << "\n";
      oss << e.what();
      throw runtime_error(oss.str());
    }
    return true;
  }

  // ---- Newton method ----
  OptNewton::OptNewton(int _max_iter, double _eps, IOptTarget* _target, int _lvl):
    IOptimizer(_max_iter, _eps, _target, _lvl) {}
  OptNewton::~OptNewton() {}
  void OptNewton::OneStep(const VectorXcd& _grad, const MatrixXcd& _hess,
			  VectorXcd& _dz) {
	       
    _dz = -_hess.fullPivLu().solve(_grad);

  }

  // ==== Optimize target ====
  // ---- Interface ----
  //IOptTarget::IOptTarget() { cout << "In IOptTarget" << endl;}
  IOptTarget::~IOptTarget() {}
  dcomplex IOptTarget::Val(const VectorXcd& zs) {
    dcomplex val;
    VectorXcd g;
    MatrixXcd h;
    this->ValGradHess(zs, &val, g, h);
    return val;
  }

  // ---- optimize for alpha ----
  OptAlpha::OptAlpha(IDriv* _driv, IOp* _op, R1GTOs& _g0): ref_g0s(_g0){

    if(_driv == NULL) {
      string msg; SUB_LOCATION(msg);
      msg += ": driv is null";
      throw runtime_error(msg);
    }
    if(_op == NULL) {
      string msg; SUB_LOCATION(msg);
      msg += ": op is null";
      throw runtime_error(msg);
    }

    driv = _driv;
    op = _op;
    ref_g0s.Normalize();
    g1s = NULL;
    g2s = NULL;

  }
  OptAlpha::~OptAlpha() {
    if(g1s != NULL)      
      delete g1s;
    if(g2s != NULL)
      delete g2s;
  }
  dcomplex OptAlpha::Val(const VectorXcd& zs) {

    this->Update(zs);

    VectorXcd m;
    MatrixXcd L;
    this->op->Calc(ref_g0s, L);
    this->driv->Calc(ref_g0s, m);

    return (m.array() * L.fullPivLu().solve(m).array()).sum();
  }
  int  OptAlpha::dim() { return ref_g0s.size_basis(); }
  void OptAlpha::Solve(VectorXcd& res) {
    
    if(!ref_g0s.setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": not set up.";
      throw runtime_error(msg);
    }

    static VectorXcd m;
    static MatrixXcd L;

    this->op->Calc(ref_g0s, L);
    this->driv->Calc(ref_g0s, m);
    res = L.fullPivLu().solve(m);

  }
  void OptAlpha::ValGradHess(const VectorXcd& zs, dcomplex* val,
			     VectorXcd& g, MatrixXcd& h) {
    
    this->Update(zs);

    int n(ref_g0s.size_basis());
    MatrixXcd D00(n,n), D10(n,n), D20(n,n), D11(n,n);
    VectorXcd m0(n), m1(n), m2(n);

    R1GTOs* g0s = &ref_g0s;
    
    op->Calc(*g0s, D00);
    op->Calc(*g1s, *g0s, D10);
    op->Calc(*g2s, *g0s, D20);
    op->Calc(*g1s, D11);

    driv->Calc(*g0s, m0);
    driv->Calc(*g1s, m1);
    driv->Calc(*g2s, m2);

    MatrixXcd Dinv = D00.inverse();
    VectorXcd c = D00.fullPivLu().solve(m0);

    // alpha
    *val = (m0.array() * c.array()).sum();

    // grad
    g = (-Calc_a_Aj_b(c, D10, c) +
	 2.0*Calc_ai_b(m1, c));

    // hess
    MatrixXcd tmp1 = (2*m2.array()*c.array()).matrix().asDiagonal();
    MatrixXcd tmp2 = 2 * (m1 * m1.transpose()).array() * Dinv.array();
    MatrixXcd tmp3 = -Calc_a_Aij_a(c, D20, D11);
    MatrixXcd tmp4 = -2.0*Calc_ai_A_Bj_b(m1, Dinv, D10, c);
    MatrixXcd tmp5 = tmp4.transpose();
    MatrixXcd tmp6 = 2.0*Calc_a_Ai_B_Aj_b(c, D10, Dinv, c);
    h = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6;
    
  }
  void OptAlpha::Update(const VectorXcd& zs) {
    try {
      ref_g0s.Set(zs);
    } catch(const exception& e) {
      string msg; SUB_LOCATION(msg);
      msg += ": update orbital exponents is failed.\n";
      msg += e.what();
      throw runtime_error(msg);
    }

    if(g1s == NULL) {
      g1s = new R1GTOs();
    }
    if(g2s == NULL) {
      g2s = new R1GTOs();
    }

    ref_g0s.Normalize();
    g1s->SetOneDeriv(ref_g0s);
    g2s->SetTwoDeriv(ref_g0s); 
  }

  // ---- opt alpha by shift ----
  /*
    .  z_i = x_i + z (if i is in A)
    .      = x_i     (otherwise)
    df          dz_i  df
    -- = SUM(i) ---- ---- = SUM(i in A)df_i
    dz          dz    dz_i
    d2f/dz2 = SUM(i in A,j) dzj/dz df_i/dzj = SUM((i,j) in A) df_ij
   */
  OptAlphaShift::OptAlphaShift(IDriv* _driv, IOp* _op, R1GTOs& _g0,
			       const VectorXi& _index):
    opt_alpha(new OptAlpha(_driv, _op, _g0)), index(_index) {
    
    int n(_g0.size_basis());
    if(index.minCoeff() < 0 || n <= index.maxCoeff()) {
      string msg; SUB_LOCATION(msg);
      msg += ": index out of range\n";
      ostringstream oss;      
      oss << msg << endl;
      oss << "index.minCoeff(): " << index.minCoeff() << endl;
      oss << "index.maxCoeff(): " << index.maxCoeff() << endl;
      oss << "# of basis      : " << n << endl;
    }
    z0s = VectorXcd::Zero(n);
    for(int i = 0; i < n; i++) {
      z0s[i] = _g0.prim(i).z;
    }
  }
  OptAlphaShift::~OptAlphaShift() {
    delete opt_alpha;
  }
  int OptAlphaShift::dim() { return 1; }
  //  dcomplex OptAlphaShift::Val(const VectorXcd& zs) { return 1.0;}
  void OptAlphaShift::ValGradHess(const VectorXcd& zs, dcomplex* a,
				  VectorXcd& g, MatrixXcd& h) {
    if(zs.size() != this->dim()) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << endl;
      oss << "zs.size(): " << zs.size()   << endl;
      oss << "dim      : " << this->dim() << endl;
      throw runtime_error(oss.str());
    }

    // convert shift to full
    VectorXcd zs_full(z0s);
    for(int i = 0; i < index.size(); i++) {
      zs_full[index(i)] += zs[0];
    }

    // compute grad/hess in full space
    VectorXcd g_full;
    MatrixXcd h_full;
    opt_alpha->ValGradHess(zs_full, a, g_full, h_full);

    // convert grad/hess to restricted space
    dcomplex gval(0);
    dcomplex hval(0);
    for(int i = 0; i < index.size(); i++) {
      gval += g_full(index[i]);
      for(int j = 0; j < index.size(); j++) {
	hval += h_full(index[i], index[j]);
      }
    }
    
    // set results
    g = VectorXcd::Zero(1);
    h = MatrixXcd::Zero(1, 1);
    g(0) = gval;
    h(0, 0) = hval;
  }  
  
  // ---- alpha by partial ----
  /*
    optimize variable xi (i<-Range(n)).
    xi = z_I(i)     (i<-Range(n))
   */
  OptAlphaPartial::OptAlphaPartial(IDriv* _driv, IOp* _op, R1GTOs& _g0,
				   const VectorXi& _index):
    opt_alpha(new OptAlpha(_driv, _op, _g0)), index(_index) {

    int n(_g0.size_basis());
    z0s = VectorXcd::Zero(n);
    for(int i = 0; i < n; i++) {
      z0s[i] = _g0.prim(i).z;
    }
  }
  OptAlphaPartial::~OptAlphaPartial() {
    delete opt_alpha;
  }
  int OptAlphaPartial::dim() { return index.size(); }
  //  dcomplex OptAlphaPartial::Val(const VectorXcd& zs) { return 1.0;}
  void OptAlphaPartial::ValGradHess(const VectorXcd& zs, dcomplex* a,
				    VectorXcd& g, MatrixXcd& h) {
    int n_part(this->dim());
    if(zs.size() != n_part) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << endl;
      oss << "zs.size(): " << zs.size()   << endl;
      oss << "dim      : " << this->dim() << endl;
      throw runtime_error(oss.str());
    }

    // convert restricted to full
    VectorXcd zs_full(z0s);
    for(int i = 0; i < n_part; i++) {
      zs_full[index(i)] = zs[i];
    }

    // compute grad/hess in full space
    VectorXcd g_full;
    MatrixXcd h_full;
    opt_alpha->ValGradHess(zs_full, a, g_full, h_full);

    // convert grad/hess to restricted space
    g = VectorXcd::Zero(n_part);
    h = MatrixXcd::Zero(n_part, n_part);    

    for(int i = 0; i < n_part; i++) {
      g(i) = g_full(index(i));
      for(int j = 0; j < n_part; j++) {
	h(i, j) = h_full(index(i), index(j));
      }
    }
  }

  void VectorShift(VectorXcd& zs_in, VectorXcd& zs_out,
		   const VectorXi& opt_idx, dcomplex z_shift) {

    int num(zs_out.size());
    for(int i = 0; i < num; i++)
      zs_out[i] = zs_in[i];
    for(int i = 0; i < opt_idx.size(); i++) {
      zs_out[opt_idx[i]] += z_shift;
    }
  }

  void AlphaGrad(const VectorXcd& d, const VectorXcd& c,
   		 const MatrixXcd& D10,
		 const VectorXcd& r1,
		 const VectorXcd& s0, const VectorXcd& s1,
		 dcomplex *a, VectorXcd *g) {

    // alpha
    *a = (d.array() * s0.array()).sum();

    // alpha
    *a = (d.array() * s0.array()).sum();

    // grad
    VectorXcd g1 = Calc_a_Aj_b(d, D10, c);
    VectorXcd g2 = Calc_ai_b(r1, c);
    VectorXcd g3 = Calc_ai_b(s1, d);
      
    *g = g2 + g3 - g1;
  }
  
  void AlphaGradHess(const VectorXcd& c, const VectorXcd& d, const MatrixXcd& U,
		     const MatrixXcd& D00, const MatrixXcd& D10,
		     const MatrixXcd& D20, const MatrixXcd& D11,
		     const VectorXcd& r0,  const VectorXcd& r1,  const VectorXcd& r2,
		     const VectorXcd& s0,  const VectorXcd& s1,  const VectorXcd& s2,
		     dcomplex* a, VectorXcd* g, MatrixXcd* h) {

    // alpha
    *a = (d.array() * s0.array()).sum();

    // grad
    VectorXcd g1 = Calc_a_Aj_b(d, D10, c);
    VectorXcd g2 = Calc_ai_b(r1, c);
    VectorXcd g3 = Calc_ai_b(s1, d);
      
    *g = g2 + g3 - g1;

    // hess
    MatrixXcd rij_c       = (r2.array() * c.array()).matrix().asDiagonal();
    MatrixXcd d_sij       = (s2.array() * d.array()).matrix().asDiagonal();
    MatrixXcd d_Lij_c     = Calc_a_Aij_b(d, D20, D11, c);

    MatrixXcd tmp = Calc_ai_A_Bj_b(r1, U, D10, c);
    MatrixXcd ri_U_Lj_c   = tmp + tmp.transpose();    

    tmp = (r1 * s1.transpose()).array() * U.array();
    MatrixXcd ri_U_sj = tmp + tmp.transpose();

    tmp = Calc_a_Ai_B_Aj_b(d, D10, U, c);
    MatrixXcd d_Li_U_Lj_c = tmp + tmp.transpose();

    tmp = Calc_ai_A_Bj_b(s1, U, D10, d);
    MatrixXcd d_Li_U_sj   = tmp + tmp.transpose();

    *h = rij_c - d_Lij_c + d_sij
      - ri_U_Lj_c + d_Li_U_Lj_c - ri_U_Lj_c + ri_U_sj;

  }

  VectorXcd SolveAlpha(const R1STOs& driv,  R1GTOs& gtos, int L, dcomplex ene,
		       MatVecMap& d, double eps) {

    gtos.CalcMatSTV(L, d, "s", "t", "v");
    gtos.CalcVec(driv, d, "m");

    VectorXcd& m = d.vec("m");
    MatrixXcd& t = d.mat("t");
    MatrixXcd& v = d.mat("v");
    MatrixXcd& s = d.mat("s");

    if(eps > pow(10.0, -14.0)) {
      MatrixXcd X;
      CanonicalMatrix(s, eps, &X);

      MatrixXcd Lp = (X.transpose() *
		      (t+v-ene*s) *
		      X);
      VectorXcd mp = X.transpose() * m;
      return X*Lp.fullPivLu().solve(mp);
    } else {
      return (t+v-ene*s).fullPivLu().solve(m);
    }
  }

/*
  void SolveAlpha(const R1STOs& driv, R1GTOs& gtos,
		  int L, dcomplex ene,
		  VectorXcd& c, dcomplex* a) {

    static MatVecMap d;

    gtos.CalcMatSTV(L, d, "s", "t", "v");
    //    gtos.CalcMatSTO(sto_pot, d, "v2");
    gtos.CalcVec(driv, d, "m");

    VectorXcd& m = d.vec("m");
    MatrixXcd& t = d.mat("t");
    MatrixXcd& v = d.mat("v");
    //    MatrixXcd& v2= d.mat("v2");
    MatrixXcd& s = d.mat("s");
    return (m.array() * 
	    ((t+v-ene*s).fullPivLu().solve(m)).array()).sum();
  }

  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos,
		     int L, dcomplex ene, double eps) {

    
    static MatVecMap mat_vec;
    
    gtos.CalcMatSTV(L, mat_vec, "s", "t", "v");
    gtos.CalcVec(driv, mat_vec, "m");
    MatrixXcd X;
    CanonicalMatrix(mat_vec.mat("s"), eps, &X);
    MatrixXcd Lp = (X.transpose() *
		    (mat_vec.mat("t") + mat_vec.mat("v") - ene * mat_vec.mat("s")) *
		    X);
    VectorXcd m = X.transpose()*mat_vec.vec("m");
    dcomplex a = (m.array() *
		  Lp.fullPivLu().solve(m).array()).sum();
  
    return a;

  }  
  dcomplex CalcAlphaFull(const R1STOs& driv,R1GTOs& gtos, int L, dcomplex ene) {

    static MatVecMap mat_vec;
    
    gtos.CalcMatSTV(L, mat_vec, "s", "t", "v");
    gtos.CalcVec(driv, mat_vec, "m");

    VectorXcd& m = mat_vec.vec("m");
    MatrixXcd& t = mat_vec.mat("t");
    MatrixXcd& v = mat_vec.mat("v");
    MatrixXcd& s = mat_vec.mat("s");
    dcomplex a = (m.array() * 
		  ((t+v-ene*s).fullPivLu().solve(m)).array()).sum();
    return a;
  }
*/

/*
  void OptAlphaShiftFull(const R1STOs& driv,
			 const VectorXi& opt_idx,
			 int L,
			 dcomplex ene,
			 double h,
			 int max_iter,
			 double eps,
			 R1GTOs& gtos,
			 dcomplex& z_shift,
			 bool* convq,
			 dcomplex* alpha) {

    if(!gtos.setup_q()) {
      string msg; SUB_LOCATION(msg);
      msg += ": gtos is not SetUp.";
      throw runtime_error(msg);
    }

    if(opt_idx.array().minCoeff() < 0 ||
       opt_idx.array().maxCoeff() >= gtos.size_prim()) {
      string msg; SUB_LOCATION(msg);
      ostringstream oss;
      oss << msg << ": invalid opt_idx" << endl;
      oss << "min(opt_idx): " << opt_idx.array().minCoeff() << endl;
      oss << "max(opt_idx): " << opt_idx.array().maxCoeff() << endl;
      oss << "size_prim   : " << gtos.size_prim() << endl;
      throw runtime_error(oss.str());
    }

    VectorXcd zs0(gtos.size_prim());
    VectorXcd zs1(gtos.size_prim());
    dcomplex ih(0.0, h);
    dcomplex ii(0.0, 1.0);
    dcomplex a00;

    for(int i = 0; i < gtos.size_prim(); i++) {
      try {
	zs0[i] = gtos.prim(i).z;
      } catch(const exception& e) {
	cout << "failed to prim" << endl;
	cout << e.what() << endl;
      }
    }

    int i(0);
    for(i = 0; i < max_iter; i++) {
      dcomplex ap0, am0, a0p, a0m;
      VectorShift(zs0, zs1, opt_idx, z_shift + h);
      gtos.Set(2, zs1); gtos.Normalize();
      ap0 = CalcAlphaFull(driv, gtos, L, ene);

      VectorShift(zs0, zs1, opt_idx, z_shift - h);
      gtos.Set(2, zs1); gtos.Normalize();
      am0 = CalcAlphaFull(driv, gtos, L, ene);

      VectorShift(zs0, zs1, opt_idx, z_shift + ih);
      gtos.Set(2, zs1); gtos.Normalize();
      a0p = CalcAlphaFull(driv, gtos, L, ene);

      VectorShift(zs0, zs1, opt_idx, z_shift - ih);
      gtos.Set(2, zs1);       gtos.Normalize();
      a0m = CalcAlphaFull(driv, gtos, L, ene);
      
      dcomplex grad = (ap0 - am0 + ii*a0m - ii*a0p) / (4.0 * h);
      dcomplex hess = (ap0 + am0 - a0p - a0m) / (2.0*h*h);

      dcomplex dz(-grad/hess);
      z_shift += dz;
      if(abs(grad) < eps && abs(dz) < eps)
	break;
    }

    *convq = (i < max_iter-2);
    cout << i << endl;

    VectorShift(zs0, zs1, opt_idx, z_shift);
    gtos.Set(2, zs1); gtos.Normalize();
    a00 = CalcAlphaFull(driv, gtos, L, ene);

    *alpha = a00;

  }

  void OptAlphaShiftCanonical(const R1STOs& driv,
			      const VectorXi& opt_idx,
			      dcomplex ene,
			      double h,
			      int max_iter,
			      double eps,
			      double eps_canonical,
			      R1GTOs& gtos,
			      dcomplex& z_shift,
			      bool* convq,
			      dcomplex* alpha) {

    VectorXcd zs0(gtos.size_basis());
    VectorXcd zs1(gtos.size_basis());
    dcomplex ih(0.0, h);
    dcomplex ii(0.0, 1.0);
    
    MatrixXcd L;

    for(int i = 0; i < gtos.size_basis(); i++) {
      zs0[i] = gtos.prim(i).z;
    }

    int i(0);
    for(i = 0; i < max_iter; i++) {
      dcomplex ap0, am0, a0p, a0m;
      VectorShift(zs0, zs1, opt_idx, z_shift + h);
      gtos.Set(2, zs1);
      ap0 = CalcAlpha(driv, gtos, ene, eps_canonical);

      VectorShift(zs0, zs1, opt_idx, z_shift - h);
      gtos.Set(2, zs1);
      am0 = CalcAlpha(driv, gtos, ene, eps_canonical);

      VectorShift(zs0, zs1, opt_idx, z_shift + ih);
      gtos.Set(2, zs1);
      a0p = CalcAlpha(driv, gtos, ene, eps_canonical);	

      VectorShift(zs0, zs1, opt_idx, z_shift - ih);
      gtos.Set(2, zs1);      
      a0m = CalcAlpha(driv, gtos, ene, eps_canonical);	
      
      dcomplex grad = (ap0 - am0 + ii*a0m - ii*a0p) / (4.0 * h);
      dcomplex hess = (ap0 + am0 - a0p - a0m) / (2.0*h*h);

      dcomplex dz(-grad/hess);
      z_shift += dz;
      if(abs(grad) < eps && abs(dz) < eps)
	break;
    }

    cout << i << endl;
    *convq = (i < max_iter-2);
    VectorShift(zs0, zs1, opt_idx, z_shift);
    gtos.Set(2, zs1);
    *alpha = CalcAlpha(driv, gtos, ene, 0.0);;

  }

  void OptAlphaShift(const R1STOs& driv,
		     const VectorXi& opt_idx,
		     int L,
		     dcomplex ene,
		     double h,
		     int max_iter,
		     double eps,
		     double eps_canonical,
		     R1GTOs& gtos,
		     dcomplex& z_shift,
		     bool* convq,
		     dcomplex* alpha) {

    double eps_canonical_th(pow(10.0, -14.0));
    if(eps_canonical < eps_canonical_th)
      OptAlphaShiftFull(driv, opt_idx, L, ene, h, max_iter, eps,
			gtos, z_shift, convq, alpha);
    else {
      string msg; SUB_LOCATION(msg);
      msg += ": canonical calculation is not supported";
      throw runtime_error(msg);
      //OptAlphaShiftCanonical(driv, opt_idx, ene, h, max_iter, eps, eps_canonical,
      //gtos, z_shift, convq, alpha);
		    }
  }
*/      
}
