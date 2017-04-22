#include <sstream>
#include <stdexcept>
#include <iostream>
#include <boost/format.hpp>
#include "typedef.hpp"
#include "macros.hpp"

#include "eigen_plus.hpp"

using namespace std;
using namespace Eigen;
using boost::format;

MatrixXcd
m13cd(dcomplex m00, dcomplex m01, dcomplex m02) {
  MatrixXcd M(1, 3);
  M << m00, m01, m02;
  return M;
}

MatrixXcd
m33cd(dcomplex m00, dcomplex m01, dcomplex m02,
      dcomplex m10, dcomplex m11, dcomplex m12,
      dcomplex m20, dcomplex m21, dcomplex m22) {
  MatrixXcd M(3, 3);
  M << m00, m01, m02, m10, m11, m12, m20, m21, m22;
  return M;
}

VectorXi v1i(int i) {
  VectorXi v(1); v << i;
  return v;
}
VectorXi v3i(int i,int j,int k) {
  VectorXi v(3); v<<i,j,k;
  return v;
}
VectorXcd v1cd(dcomplex v) {
  VectorXcd vec(1); vec << v;
  return vec;
}
VectorXcd v3cd(dcomplex v0, dcomplex v1, dcomplex v2) {
  VectorXcd vec(3); vec << v0, v1, v2;
  return vec;
}

dcomplex TDot(const VectorXcd& xs, const VectorXcd& ys) {
  return (xs.array() * ys.array()).sum();
}
double TakeReal(dcomplex x) {
  return x.real();
}
double TakeAbs(dcomplex x) {
  return abs(x);
}
std::complex<double> cnorm(const CV& v) {
  return sqrt(TDot(v, v));
}
void complex_normalize(CV& v) {
  Eigen::ArrayXcd u = v.array();
  std::complex<double> cnorm = sqrt((u*u).sum());
  v = v / cnorm;
}
void col_cnormalize(CM& c) {

  int n = c.cols();
  for(int j = 0; j < n; j++) {

    std::complex<double> cn = cnorm(c.col(j));

    for(int i = 0; i < n; i++)
      c(i, j) /= cn;

  }
}
void matrix_inv_sqrt(const CM& s, CM* s_inv_sqrt) {
  // compute eigen value problem of S.
  /// S D = D V;
  Eigen::ComplexEigenSolver<CM> es;
  es.compute(s, true);
  CV v = es.eigenvalues();
  CM c = es.eigenvectors();

  // ensure ||c_i|| = 1;
  col_cnormalize(c);

  // lambda_ij = delta_ij / sqrt(v_i);
  CV tmp = v.array().inverse().sqrt();
  CM lambda_mat = tmp.asDiagonal();

  // S^(-1/2) = D diag{1/sqrt(v_i)} D^T
  CM c_tr = c;
  c_tr.transposeInPlace();
  *s_inv_sqrt = c * lambda_mat * c_tr;
}
void matrix_sqrt(const CM& s, CM& s_sqrt) {

  // -- structure check
  if(s.rows() != s.cols()) {
    THROW_ERROR("only square matrix");
  }

  // -- compute eigen values
  Eigen::ComplexEigenSolver<CM> es;
  es.compute(s, true);
  CV v = es.eigenvalues();
  CM c = es.eigenvectors();
  
  // -- ensure |c_i| =1
  col_cnormalize(c);

  // -- lamba_ij = delta_ij * sqrt(v_i)
  CV tmp = v.array().sqrt();
  CM lambda_mat = tmp.asDiagonal();

  // -- S^(1/2) = D diag(sqrt(v_i)) D^T
  CM c_tr = c; c_tr.transposeInPlace();  
  s_sqrt = c * lambda_mat * c_tr;
  
}
void SortEigs(Eigen::VectorXcd& eigs, Eigen::MatrixXcd& eigvecs,
	      double (*to_real)(dcomplex), bool reverse) {
  
  int n = eigs.size();
  if(n != eigvecs.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch";
    throw runtime_error(msg);
  }

  for(int i = 1; i < n; i++) {
    for(int j = i-1; j > -1; j--) {
      if(not reverse && to_real(eigs[j]) > to_real(eigs[j+1])) {
	dcomplex tmp = eigs[j];
	eigs[j] = eigs[j+1];
	eigs[j+1] = tmp;
	eigvecs.col(j).swap(eigvecs.col(j+1));
      }
      if(reverse && to_real(eigs[j]) < to_real(eigs[j+1])) {
	dcomplex tmp = eigs[j];
	eigs[j] = eigs[j+1];
	eigs[j+1] = tmp;
	eigvecs.col(j).swap(eigvecs.col(j+1));
      }
    }
  }
    
}
void generalizedComplexEigenSolve(const CM& f, const CM& s, CM* c, CV* eig){

  SymGenComplexEigenSolver solver(f, s);
  *c = solver.eigenvectors();
  *eig = solver.eigenvalues();
  
  /*
  THROW_ERROR("do not call it");
  if(f.rows() != s.rows() || f.rows() == 0) {
    string msg; SUB_LOCATION(msg);
    stringstream oss;
    oss << msg << ": invalid matrix size for f and s." << endl
	<< "s = " << s.rows() << s.cols() << endl
	<< "f = " << f.rows() << f.cols() << endl;
    throw runtime_error(oss.str());
  }

  // s2inv means S^(-1/2)
  CM s2inv;
  matrix_inv_sqrt(s, &s2inv);
  
  // fp means F' = S^(-1/2)FS^(-1/2)
  CM fp = s2inv * f * s2inv;

  // solve F'C' = C diag{e_i}
  Eigen::ComplexEigenSolver<CM> es;
  es.compute(fp, true);
  *c = s2inv *  es.eigenvectors();
  *eig = es.eigenvalues();
  SortEigs(*eig, *c, TakeReal);
  */
    
}
void CanonicalMatrix(const CM& S, double eps, CM* res) {

  int num(S.rows());
  if(num != S.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch."; 
    throw runtime_error(msg);
  }

  ComplexEigenSolver<CM> es;
  es.compute(S, true);
  CV s;
  CM U;
  s = es.eigenvalues();
  U = es.eigenvectors();
  SortEigs(s, U, TakeAbs, true);
  
  int num_non0(0);
  for(int j = 0; j < num; j++)
    if(abs(s[j]) > eps) {
      num_non0++;
    }
  
  MatrixXcd X(num, num_non0);
  int j_idx(0);
  for(int j = 0; j < num; j++) {
    if(abs(s[j]) > eps) {
      for(int i = 0; i < num; i++) {
	X(i, j_idx) = U(i, j) / sqrt(s(j));
      }
      j_idx++;
    }
  }
  
  *res = CM::Zero(1, 1);
  res->swap(X);

}
void CanonicalMatrixNum(const CM& S, int num0, CM* res) {

  int num_all(S.rows());
  
  if(num_all != S.cols()) {
    string msg; SUB_LOCATION(msg);
    msg += ": size mismatch."; 
    throw runtime_error(msg);
  }
  
  if(num_all < num0) {
    string msg; SUB_LOCATION(msg);
    msg += "num0 must be lesser than num_all";
    throw runtime_error(msg);
  }

  ComplexEigenSolver<CM> es;
  es.compute(S, true);
  CV s;
  CM U;
  s = es.eigenvalues();
  U = es.eigenvectors();
  SortEigs(s, U, TakeAbs, true);

  MatrixXcd X = MatrixXcd::Zero(num_all, num0);
  for(int j = 0; j < num0; j++) {
    for(int i = 0; i < num_all; i++ ) {
      X(i, j) = U(i, j) / sqrt(s(j));
    }
  }
  
  *res = CM::Zero(1, 1);
  res->swap(X);
  
}
void CEigenSolveCanonical(const CM& F, const CM& S, double eps, CM* c, CV* eig) {

  MatrixXcd X;
  CanonicalMatrix(S, eps, &X);
  CM Fp = X.transpose() * F * X;

  ComplexEigenSolver<CM> es;
  es.compute(Fp, true);
  *c = X * es.eigenvectors();
  *eig = es.eigenvalues();
  SortEigs(*eig, *c, TakeReal, false);

}
void CEigenSolveCanonicalNum(const CM& F, const CM& S, int num0,
			      CM* c, CV* eig) {

  MatrixXcd X;
  CanonicalMatrixNum(S, num0, &X);
  CM Fp = X.transpose() * F * X;

  ComplexEigenSolver<CM> es;
  es.compute(Fp, true);
  *c = X * es.eigenvectors();
  *eig = es.eigenvalues();
  SortEigs(*eig, *c, TakeReal, false);
  
}

void CtAD(const CM& C, const CM& D, const CM& A, CM *res) {
  /**
     gives C^T A C where T means translate
     (C^T A D)_{ij} = C_{ki} A_{kl} D(lj)
  */
 
  int ni = C.rows();
  int nj = D.cols();
  int nk = C.cols();
  int nl = D.rows();

  if(A.rows() != nk) {
    stringstream ss;
    ss << "size of A and C mismatch\n";
    ss << format("C: (%d, %d)\n") % C.rows() % C.cols();
    ss << format("A: (%d, %d)\n") % A.rows() % A.cols();
    ss << format("D: (%d, %d)\n") % D.rows() % D.cols();
    THROW_ERROR(ss.str());
  }
  if(A.cols() != nl) {
    stringstream ss;
    ss << "size of A and D mismatch\n";
    ss << format("C: (%d, %d)\n") % C.rows() % C.cols();
    ss << format("A: (%d, %d)\n") % A.rows() % A.cols();
    ss << format("D: (%d, %d)\n") % D.rows() % D.cols();
    THROW_ERROR(ss.str());
  }

  if(res->rows() != ni || res->cols() != nj) {
    *res = MatrixXcd::Zero(ni, nj);
  }

  for(int i = 0; i < ni; i++)
    for(int j = 0; j < nj; j++) {
      dcomplex acc(0);
      for(int k = 0; k < nk; k++)
	for(int l = 0; l < nl; l++) {
	  acc += C(k, i) * A(k, l) * D(l, j);
	}
      (*res)(i, j) = acc;
    }
}
void CtAD(const CM& C, const CM& D, CM& A) {
  int ni = C.rows();
  int nj = C.cols();
  static CM res(ni, nj);
  if(res.cols() != ni || res.rows() != nj) 
    res = CM::Zero(ni, nj);
  CtAD(C, D, A, &res);
  A.swap(res);
}
void CaAD(const CM& C, const CM& D, const CM& A, CM *res) {
 /**
     gives C^+ A C where T means translate
     (C^+ A D)_{ij} = C^*_{ki} A_{kl} D(lj)
  */

  int ni = C.rows();
  int nj = D.cols();
  int nk = C.cols();
  int nl = D.rows();

  if(A.rows() != nk) {
    THROW_ERROR("size of A and C mismatch");
  }
  if(A.cols() != nl) {
    THROW_ERROR("size of A and D miscmatch");
  }

  if(res->rows() != ni || res->cols() != nj) {
    *res = MatrixXcd::Zero(ni, nj);
  }

  for(int i = 0; i < ni; i++)
    for(int j = 0; j < nj; j++) {
      dcomplex acc(0);
      for(int k = 0; k < nk; k++)
	for(int l = 0; l < nl; l++) {
	  acc += conj(C(k, i)) * A(k, l) * D(l, j);
	}
      (*res)(i, j) = acc;
    }  
 
}
void CaAD(const CM& C, const CM& D,CM& A) {
  int ni = C.rows();
  int nj = D.cols();
  static CM res(ni, nj);
  if(res.rows() != ni || res.cols() != nj) 
    res = CM::Zero(ni, nj);
  CaAD(C, D, A, &res);
  A.swap(res);
}

void CtAC(const CM& C, const CM& A, CM *CtAC) {
  /**
     gives C^T A C where T means translate
     (C^T A C)_{ij} = C_{ki} A_{kl} C(lj)
   */
  CtAD(C, C, A, CtAC);
}
void CtAC(const CM& C, CM& A) {
  CtAD(C, C, A);
}
void CaAC(const CM& C, const CM& A, CM *CdAC) {
  CaAD(C, C, A, CdAC);
}
void CaAC(const CM& C, CM& A) {
  CaAD(C, C, A);
}

void Ctx(const CM& C, const CV& x, CV *Ctx) {
  /**
     gives C^T x where T means translate
     (C^T x)_{ij} = sum_k C_{ki}x_i
   */
 
  int n = C.rows();
  if(n != C.cols()) {
    THROW_ERROR("C is not square");
  }
  if(n != x.size()) {
    THROW_ERROR("size of x mismatch to A");
  }
  if(Ctx->size() != n) {
    *Ctx = VectorXcd::Zero(n);
  }

  for(int i = 0; i < n; i++) {
    dcomplex acc(0);
    for(int k = 0; k < n; k++) {
      acc += C(k, i) * x(k);
    }
    (*Ctx)(i) = acc;
  }    
}
void Ctx(const CM& C, CV& x) {
  int n = x.size();
  static CV res(n);
  if(res.size() != n) {
    res = CV::Zero(n);
  }
  Ctx(C, x, &res);
  res.swap(x);
}
SymGenComplexEigenSolver::SymGenComplexEigenSolver() {}
SymGenComplexEigenSolver::SymGenComplexEigenSolver(int _num0) {
  num0_ = _num0;
}
SymGenComplexEigenSolver::SymGenComplexEigenSolver(const CM& f, const CM& s) {
  num0_ = f.rows();
  this->compute(f, s);
}
SymGenComplexEigenSolver::SymGenComplexEigenSolver(const CM& f, const CM& s, int _num0) {
  num0_ = _num0;
  CEigenSolveCanonicalNum(f, s, num0_,
			  &this->eigenvectors_, &this->eigenvalues_);
}
void SymGenComplexEigenSolver::compute(const CM& f, const CM& s) {
  
  if(f.rows() != s.rows() || f.rows() != f.cols() || f.rows() == 0) {
    string msg; SUB_LOCATION(msg);
    stringstream oss;
    oss << msg << ": invalid matrix size for f and s." << endl
	<< "s = " << s.rows() << s.cols() << endl
	<< "f = " << f.rows() << f.cols() << endl;
    throw runtime_error(oss.str());
  }

  int n(f.rows());

  if(n != num0_) {
    CEigenSolveCanonicalNum(f, s, num0_,
			    &this->eigenvectors_, &this->eigenvalues_);
    return;
  }
  
  if(s2inv_.rows() != n) {
    s2inv_ = MatrixXcd::Zero(n, n);
  }
  if(fp_.rows() != n) {
    fp_ = MatrixXcd::Zero(n, n);
  }

  // -- s2inv means S^(-1/2) --
  this->compute_inv_sqrt(s);
  
  // -- fp means F' = S^(-1/2)FS^(-1/2) --
  fp_ = s2inv_ * f * s2inv_;

  // -- solve F'C' = C diag{e_i} --
  es_.compute(fp_, true);

  // -- C = S^{-1/2}C'
  eigenvectors_ = s2inv_ *  es_.eigenvectors(); 
  eigenvalues_ = es_.eigenvalues();

  // -- sort --
  SortEigs(eigenvalues_, eigenvectors_, TakeReal);

  // -- normalize to unity --
  for(int i = 0; i < n; i++) {
    const MatrixXcd& Ci =  eigenvectors_.col(i);
    dcomplex cSc = TDot(Ci, s*Ci);
    eigenvectors_.col(i) = Ci / sqrt(cSc);
  }
}
void SymGenComplexEigenSolver::compute_inv_sqrt(const CM& s) {

  int n(s.rows());
  if(v_.size() != s.rows())
    v_ = VectorXcd::Zero(s.rows());
  if(c_.rows() != n || c_.cols() != n)
    c_ = MatrixXcd::Zero(s.rows(), s.cols());
  if(lambda_mat_.rows() != n || lambda_mat_.cols() != n)
    lambda_mat_ = MatrixXcd::Zero(n, n);
  
  // -- compute eigen value problem of S. --
  // -- S D = D V --
  es_.compute(s, true);
  c_ = es_.eigenvectors();

  // -- ensure (c_i, c_i) = 1 --
  col_cnormalize(c_);

  // -- lambda_ij = delta_ij / sqrt(v_i) --
  lambda_mat_ = es_.eigenvalues().array().inverse().sqrt().matrix().asDiagonal();

  // S^(-1/2) = D diag{1/sqrt(v_i)} D^T

  s2inv_ = c_ * lambda_mat_ * c_.transpose();
}
const CV& SymGenComplexEigenSolver::eigenvalues() const {
  return eigenvalues_;
}
const CM& SymGenComplexEigenSolver::eigenvectors() const {
  return eigenvectors_;
}

LinearSolver::LinearSolver(string _method) {
  if(_method == "householderQr") {
    method_ = 0;
  } else if(_method == "colPivHouseholderQr") {
    method_ = 1;
  } else if(_method == "fullPivHouseholderQr") {
    method_ = 2;
  } else {
    string loc; SUB_LOCATION(loc);
    string msg = loc + "\nunsupported method. choose (householderQr, colPivHouseholderQr, fullPivHouseholderQr)";
    throw runtime_error(msg);
  }
  
}

//void LinearSolver::SetMatrix(Eigen::MatrixXcd& _mat) {
//  if(method_ == method_householderQr) {
//    householder = _mat.householderQr();
//  }
//  if(method_ == method_colPivHouseholderQr) {
//    col_piv = _mat.colPivHouseholderQr();
//  }
//  if(method_ == method_fullPivHouseholderQr) {
//    full_piv = _mat.fullPivHouseholderQr();
//  }
  //}

void LinearSolver::Solve(MatrixXcd& mat, VectorXcd& vec, VectorXcd *sol) {
    
  if(method_ == method_householderQr) {
    *sol = mat.householderQr().solve(vec);
  }
  if(method_ == method_colPivHouseholderQr) {
    *sol = mat.colPivHouseholderQr().solve(vec);
  }
  if(method_ == method_fullPivHouseholderQr) {
    *sol = mat.fullPivHouseholderQr().solve(vec);
  }
}
void LinearSolver::Inv(MatrixXcd& m, MatrixXcd *sol) {

  //  THROW_ERROR("do not use this");
  int n(m.rows());
  if(m.cols() != n) {
    THROW_ERROR("m must be square");
  }

  if(sol == NULL) {
    *sol = MatrixXcd::Zero(n, n);
  }

  if(sol->rows() != n && sol->rows() != n) {
    *sol = MatrixXcd::Zero(n, n);
  }

  VectorXcd id(n);
  VectorXcd sol_vec(n);

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++)
      id(j) = 0.0;
    id(i) = 1.0;
    
    if(method_ == method_householderQr) {
      sol_vec = m.householderQr().solve(id);
    }
    if(method_ == method_colPivHouseholderQr) {
      sol_vec = m.colPivHouseholderQr().solve(id);
    }
    if(method_ == method_fullPivHouseholderQr) {
      sol_vec = m.fullPivHouseholderQr().solve(id);
    }
    for(int j = 0; j < n; j++) {
      (*sol)(j, i) = sol_vec(j);
    }
  }
  
}
string LinearSolver::str() {
  if(method_ == method_householderQr) {
    return "householderQr";
  }
  if(method_ == method_colPivHouseholderQr) {
    return "colPivHouseholderQr";
  }
  if(method_ == method_fullPivHouseholderQr) {
    return "fullPivHouseholderQr";
  } else {
    return "unsupported";
  }
}
