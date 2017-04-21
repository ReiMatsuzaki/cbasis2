#include <iostream>
#include <gtest/gtest.h>
#include <Eigen/Core>

#include "../utils/gtest_plus.hpp"
#include "../utils/eigen_plus.hpp"
#include "nderiv.hpp"
#include "erfc.hpp"
#include "int_exp.hpp"

using namespace std;
using namespace Eigen;
using namespace cbasis;
dcomplex Func1(dcomplex x) {
  return sin(x);
}
TEST(NDeriv, OneTwo) {
  dcomplex h(0.0001);
  dcomplex x(0.2, 0.1);
  EXPECT_C_NEAR(cos(x),
		NDerivOne_R1(Func1, x, h),
		pow(10.0, -8.0));
  EXPECT_C_NEAR(cos(x),
		NDerivOne_C1(Func1, x, h),
		pow(10.0, -12.0));
  EXPECT_C_NEAR(cos(x),
		NDerivOne_R3(Func1, x, h),
		pow(10.0, -12.0));

  EXPECT_C_NEAR(-sin(x),
		NDerivTwo_R1(Func1, x, h),
		pow(10.0, -8.0));
  EXPECT_C_NEAR(-sin(x),
		NDerivTwo_C1(Func1, x, h),
		pow(10.0, -8.0));
  EXPECT_C_NEAR(-sin(x),
		NDerivTwo_R3(Func1, x, h),
		pow(10.0, -8.0));  
  
}

TEST(Erfc, real_erfc) {

  using namespace erfc_mori;

  // this function is forbidden
  //  erfc_add_Eh_q<int>(1, 2);
  
  double y, x, expect;
  ErfcCalcData calc_data;

  x = 1.0;
  expect =0.157299207050285130658779364917390740703933002034;
  //erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);
  
  x = 1.0;
  x = 1.0 / 100.0;
  expect = 0.98871658444415038308409047645193078905089904517;
  //  erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

  x = 3.0;
  expect = 0.0000220904969985854413727761295823203798477070873992;
  //  erfc_d(x, y, calc_data);
  // erfc(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

}
TEST(Erfc, complex_erfc) {

  using namespace erfc_mori;

  dcomplex y;
  ErfcCalcData calc_data;
  double eps = 10.0 * machine_eps();

  dcomplex x(1, -1);
  dcomplex y_expect( -0.31615128169794764488027108024367,
	       +0.190453469237834686284108861969162);

  Erfc(x, y, calc_data);

  EXPECT_DOUBLE_EQ(y.real(), y_expect.real());
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_DOUBLE_EQ(y.imag(), y_expect.imag());
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

  x = dcomplex(0.0157073173118206757532953533099,
	 0.9998766324816605986389071277312);
  y_expect = dcomplex(0.95184545524179913420177473658805,
		-1.64929108965086517748934245403332);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(     y.real(), y_expect.real(), eps);
  EXPECT_NEAR(     y.imag(), y_expect.imag(), eps);

  x = dcomplex(0.001564344650402308690101053194671668923139,
	 0.009876883405951377261900402476934372607584);
  y_expect = dcomplex(0.9982346553205423153337357292658472915601,
		-0.0111452046101524188315708507537751407281);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

}

TEST(ExpInt, STO_Int_Rplus) {

  EXPECT_C_EQ(4.09808073219042,  STOInt_Rplus(3, 1.1)+1.0);
  EXPECT_C_EQ(0.413223140495868, GTOInt_Rplus(3, 1.1));
  EXPECT_C_EQ(0.768167471819406, GTOInt_R(2, 1.1));
  EXPECT_C_EQ(0.0873211906359305 + 0.0197914200245872j,
	      STO_GTOInt_Rplus(3, 1.1, 1.3-0.2j));
  EXPECT_C_EQ(0.599365436693823j,
	      STO_GTOInt_R(3, dcomplex(0, 1.1), 1.2));  
  
  // -- see support/int_exp.py --
  //  dcomplex calc = 2.2*2.2*STOInt_Rplus(3+3+2, 1.1+1.1);
  //  dcomplex ref  = 161.644807242673;
  //  EXPECT_C_EQ(ref, calc);
  //  calc = 2.2*1.
  
}
TEST(ExpInt, fast_array) {

  dcomplex a(1.1, 0.3);
  dcomplex b(0.5, 0.2);

  vector<dcomplex> res;
  for(int maxn = 0; maxn < 4; maxn++) {
    STOInt_Rplus_array(maxn, a, &res);
    for(int n = 0; n <= maxn; n++)
      EXPECT_C_EQ(STOInt_Rplus(n, a), res[n]);
    GTOInt_Rplus_array(maxn, a, &res);
    for(int n = 0; n <= maxn; n++)
      EXPECT_C_EQ(GTOInt_Rplus(n, a), res[n]);
    STO_GTOInt_Rplus_array(maxn, a, b, &res);
    for(int n = 0; n <= maxn; n++)
      EXPECT_C_EQ(STO_GTOInt_Rplus(n, a, b), res[n]);        
  }
  
}
int main(int argc, char **args) {
  cout << "wa:" << endl;
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}

