#ifndef EIGNN_PLUS_H
#define EIGNN_PLUS_H

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "typedef.hpp"

template<class F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>
m33(F m00, F m01,
    F m10, F m11) {
  Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> M(2, 2);
  M << m00, m01, m10, m11;
  return M;
}

template<class F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>
m13(F m00, F m01, F m02) {
  Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> M(1, 3);
  M << m00, m01, m02;
  return M;
}
Eigen::MatrixXcd
m13cd(dcomplex m00, dcomplex m01, dcomplex m02);
Eigen::MatrixXcd
m33cd(dcomplex m00, dcomplex m01, dcomplex m02,
      dcomplex m10, dcomplex m11, dcomplex m12,
      dcomplex m20, dcomplex m21, dcomplex m22);
Eigen::VectorXi v1i(int i);
Eigen::VectorXi v3i(int,int,int);
Eigen::VectorXcd v1cd(dcomplex v);
Eigen::VectorXcd v3cd(dcomplex v1, dcomplex v2, dcomplex v4);

template<class F>
Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>
m33(F m00, F m01, F m02,
    F m10, F m11, F m12,
    F m20, F m21, F m22) {
  Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic> M(3, 3);
  M << m00, m01, m02, m10, m11, m12, m20, m21, m22;
  return M;
}

namespace {
  typedef Eigen::MatrixXcd CM;
  typedef Eigen::VectorXcd CV;
}

dcomplex TDot(const Eigen::VectorXcd& xs, const Eigen::VectorXcd& ys);
double TakeReal(dcomplex x);
double TakeAbs(dcomplex x);
std::complex<double> cnorm(const CV& v);
void complex_normalize(CV& v);
void col_cnormalize(CM& c);
void matrix_inv_sqrt(const CM& s, CM* s_inv_sqrt);
void matrix_sqrt(const CM& s, CM& s_sqrt);
void SortEigs(Eigen::VectorXcd& eigs, Eigen::MatrixXcd& eigvecs,
	      double (*to_real)(dcomplex), bool reverse=false);
void generalizedComplexEigenSolve(const CM& f, const CM& s, CM* c, CV* eig);
void CanonicalMatrix(const CM& S, double eps, CM* res);
void CEigenSolveCanonical(const CM& f, const CM& s, double eps, CM* c, CV* eig);
void CEigenSolveCanonicalNum(const CM& F, const CM& S, int num0,
			     CM* c, CV* eig);

void CtAD(const CM& C, const CM& D, const CM& A, CM *res);
void CtAD(const CM& C, const CM& D,CM& A);
void CaAD(const CM& C, const CM& D, const CM& A, CM *res);
void CaAD(const CM& C, const CM& D,CM& A);

void CtAC(const CM& C, const CM& A, CM *CtAC);
void CtAC(const CM& C, CM& A);
void CaAC(const CM& C, const CM& A, CM *CdAC);
void CaAC(const CM& C, CM& A);
void Ctx(const CM& C, const CV& x, CV *v);
void Ctx(const CM& C, CV& x);
class SymGenComplexEigenSolver {
private:
  // -- solution of calculation --
  CM eigenvectors_;
  CV eigenvalues_;  

  // -- intermediate --
  Eigen::ComplexEigenSolver<CM> es_;
  CM s2inv_;
  CM fp_;
  CM c_;
  CV v_;
  CM lambda_mat_;
  int num0_;    // size of calculations
  
public:
  SymGenComplexEigenSolver();
  SymGenComplexEigenSolver(int _num0);
  SymGenComplexEigenSolver(const CM& f, const CM& s);
  SymGenComplexEigenSolver(const CM& f, const CM& s, int _num0);
  void compute(const CM& f, const CM& s);
  void compute_inv_sqrt(const CM& s);
  const CV& eigenvalues() const;
  const CM& eigenvectors() const;
};

class LinearSolver {
private:
  static const int method_householderQr = 0;
  static const int method_colPivHouseholderQr = 1;
  static const int method_fullPivHouseholderQr = 2;
  int method_;
  //  Eigen::MatrixXcd *ptr_mat_;
  //  bool set_matrix;
  //  Eigen::HouseholderQR<Eigen::MatrixXcd> householder;
  //  Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> col_piv;
  //  Eigen::FullPivHouseholderQR<Eigen::MatrixXcd> full_piv;
public:
  LinearSolver(): method_(0) {}
  LinearSolver(std::string method);
  //  void SetMatrix(Eigen::MatrixXcd& _mat);
  //  void Solve(Eigen::VectorXcd& x, Eigen::VectorXcd *y);
  void Solve(Eigen::MatrixXcd& mat, Eigen::VectorXcd& vec, Eigen::VectorXcd *sol);
  void Inv(Eigen::MatrixXcd& m, Eigen::MatrixXcd *sol);
  std::string str();
};


#endif
