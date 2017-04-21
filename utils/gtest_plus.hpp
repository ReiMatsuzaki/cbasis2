#ifndef GTEST_PLUS_H
#define GTEST_PLUS_H

#include <complex>
#include <gtest/gtest.h>
#include <Eigen/Core>

bool complex_is_near(std::complex<double> a, std::complex<double> b, double eps) {
  return (abs(a-b) < eps);
}
bool complex_is_near(std::complex<double> a, std::complex<double> b) {
  double eps(pow(10.0, -10.0));
  return complex_is_near(a, b, eps);
}

::testing::AssertionResult AssertComplexEq(const char* a_expr,
					   const char* b_expr,
					   std::complex<double> a,
					   std::complex<double> b) {

  /*
  typedef ::testing::internal::FloatingPoint<float> FP;
  const FP a_r(a.real()), a_i(a.imag()), b_r(b.real()), b_i(b.imag());

  if(a_r.AlmostEquals(b_r) && a_i.AlmostEquals(b_i)) 
  return ::testing::AssertionSuccess();
  */
  if(complex_is_near(a, b))
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
    << a_expr << " and " << b_expr << " are not near." << std::endl
    << a_expr << " : " << a << std::endl 
    << b_expr << " : " << b << std::endl
    << "|" << a_expr << "-" << b_expr<< "| : " << abs(a-b) << std::endl;

}

::testing::AssertionResult AssertComplexNear(const char* a_expr,
					     const char* b_expr,
					     const char* eps_expr,
					     std::complex<double> a,
					     std::complex<double> b,
					     double eps) {


  if(abs(a-b) < eps)
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
    << a_expr << " and " << b_expr << " are not near." << std::endl
    << a_expr << " : " << a << std::endl 
    << b_expr << " : " << b << std::endl
    << eps_expr << " : " << eps << std::endl
    << "|" << a_expr << "-" << b_expr<< "| : " << abs(a-b) << std::endl;

}

#define EXPECT_C_EQ(a, b) EXPECT_PRED_FORMAT2(AssertComplexEq, a, b)
#define EXPECT_C_NEAR(a, b, c) EXPECT_PRED_FORMAT3(AssertComplexNear, (a), (b), (c))

::testing::AssertionResult AssertMatrixXcdEq(const char* a_expr,
					     const char* b_expr,
					     const Eigen::MatrixXcd& a,
					     const Eigen::MatrixXcd& b) {
  
  if(a.cols() != b.cols() || a.rows() != b.rows()) {

    return ::testing::AssertionFailure()
      << a_expr << " and " << b_expr << " are not same size." << std::endl
      << a_expr << " : (" << a.rows() << ", " << a.cols() << ")" << std::endl 
      << b_expr << " : (" << b.rows() << ", " << b.cols() << ")" << std::endl;

  }
  
  for(int i = 0; i < a.rows(); i++)
    for(int j = 0; j < a.cols(); j++) 
      if(! complex_is_near(a.coeff(i, j), b.coeff(i, j)))
	return ::testing::AssertionFailure()
	  << a_expr << " and " << b_expr << " are not near." << std::endl
	  << "at " << i << ", " << j << std::endl
	  << a_expr << " : " << a.coeff(i, j) << std::endl
	  << b_expr << " : " << b.coeff(i, j) << std::endl;

  return ::testing::AssertionSuccess();
  
}
#define EXPECT_MATXCD_EQ(a, b) EXPECT_PRED_FORMAT2(AssertMatrixXcdEq, a, b)

::testing::AssertionResult AssertMatrixXiEq(const char* a_expr,
					    const char* b_expr,
					    const Eigen::MatrixXi& a,
					    const Eigen::MatrixXi& b) {
  
  if(a.cols() != b.cols() || a.rows() != b.rows()) {

    return ::testing::AssertionFailure()
      << a_expr << " and " << b_expr << " are not same size." << std::endl
      << a_expr << " : (" << a.rows() << ", " << a.cols() << ")" << std::endl 
      << b_expr << " : (" << b.rows() << ", " << b.cols() << ")" << std::endl;

  }
  
  for(int i = 0; i < a.rows(); i++)
    for(int j = 0; j < a.cols(); j++) 
      if(!(a.coeff(i, j) == b.coeff(i, j)))
	return ::testing::AssertionFailure()
	  << a_expr << " and " << b_expr << " are not near." << std::endl
	  << "at " << i << ", " << j << std::endl
	  << a_expr << " : " << a.coeff(i, j) << std::endl
	  << b_expr << " : " << b.coeff(i, j) << std::endl;

  return ::testing::AssertionSuccess();
}
#define EXPECT_MATXI_EQ(a, b) EXPECT_PRED_FORMAT2(AssertMatrixXiEq, a, b)

#endif
