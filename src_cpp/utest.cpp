/*
To be removed
*/

#include <iostream>
#include <complex>
#include <gtest/gtest.h>
#include "fact.hpp"
#include "erfc.hpp"
#include "prim.hpp"
#include "delta.hpp"
#include "cut_prim.hpp"
#include "op.hpp"
#include "lcomb.hpp"
#include "hatom.hpp"

using namespace erfc_mori;
using namespace cbasis;
using namespace fact;
using namespace std;

::testing::AssertionResult AssertComplexNear(const char* a_expr,
					     const char* b_expr,
					     const char* eps_expr,
					     complex<double> a,
					     complex<double> b,
					     double eps) {
  if(abs(a-b) < eps)
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
    << a_expr << " and " << b_expr << " are not near." << endl
    << a_expr << " : " << a << endl 
    << b_expr << " : " << b << endl
    << eps_expr << " : " << eps << endl
    << "|" << a_expr << "-" << b_expr<< "| : " << abs(a-b) << endl;

}

TEST(math, Factorial) {

  EXPECT_ANY_THROW(Factorial(-1));
  EXPECT_EQ(1, Factorial(0));
  EXPECT_EQ(1, Factorial(1));
  EXPECT_EQ(2, Factorial(2));
  EXPECT_EQ(6, Factorial(3));
  EXPECT_EQ(24, Factorial(4));
  EXPECT_EQ(120, Factorial(5));
  EXPECT_EQ(720, Factorial(6));

  EXPECT_ANY_THROW(DoubleFactorial(-1));
  EXPECT_EQ(1,   DoubleFactorial(0));
  EXPECT_EQ(1,   DoubleFactorial(1));
  EXPECT_EQ(2,   DoubleFactorial(2));
  EXPECT_EQ(3,   DoubleFactorial(3));
  EXPECT_EQ(8,   DoubleFactorial(4));
  EXPECT_EQ(15,  DoubleFactorial(5));
  EXPECT_EQ(48,  DoubleFactorial(6));
  EXPECT_EQ(105, DoubleFactorial(7));  
  
}
TEST(math, real) {


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
TEST(math, complex) {

  CD y;
  ErfcCalcData calc_data;
  double eps = 10.0 * machine_eps();

  CD x(1, -1);
  CD y_expect( -0.31615128169794764488027108024367,
	       +0.190453469237834686284108861969162);

  Erfc(x, y, calc_data);

  EXPECT_DOUBLE_EQ(y.real(), y_expect.real());
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_DOUBLE_EQ(y.imag(), y_expect.imag());
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);

  x = CD(0.0157073173118206757532953533099,
	 0.9998766324816605986389071277312);
  y_expect = CD(0.95184545524179913420177473658805,
		-1.64929108965086517748934245403332);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(     y.real(), y_expect.real(), eps);
  EXPECT_NEAR(     y.imag(), y_expect.imag(), eps);

  x = CD(0.001564344650402308690101053194671668923139,
	 0.009876883405951377261900402476934372607584);
  y_expect = CD(0.9982346553205423153337357292658472915601,
		-0.0111452046101524188315708507537751407281);

  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);
}
TEST(math, lower_gamma) {
  double eps = pow(10.0, -13.0);
  CD a(1.1, -0.3);

  EXPECT_ANY_THROW(LowerGamma(-1, a));

  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      1.0-exp(-a),
		      LowerGamma(1, a), 
		      eps);
  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      -a*exp(-a) + 1.0-exp(-a),
		      LowerGamma(2, a), 
		      eps);

}
TEST(Prim, Construct) {
  
  RSTO s0;
  RSTO s1(1.0, 2, 1.0);
  RGTO g1(1.0, 2, 1.0);
  RSTO s2(3, 1.1, Normalized);

  EXPECT_DOUBLE_EQ(0.0, s0.c());
  EXPECT_EQ(0, s0.n());
  EXPECT_DOUBLE_EQ(1.0, s0.z());
    
  EXPECT_EQ(2, s1.n());
  EXPECT_EQ(2, g1.n());
  //  EXPECT_EQ(3, s2.n());

  CSTO s3(1.0, 3, CD(1.0, -0.2));
  CGTO g3(1.0, 3, CD(1.0, -0.2));
  CSTO s4(s1);
  EXPECT_EQ(3, s3.n());
  EXPECT_EQ(3, g3.n());
  EXPECT_DOUBLE_EQ(1.0, s4.z().real());

  EXPECT_EQ(1, CSTO::exp_power);
}
TEST(Prim, accessor) {

  RSTO s1(1.1, 2, 0.2);
  double x0(3.0);
  EXPECT_DOUBLE_EQ(1.1 * x0 * x0 * exp(-0.2 * x0), 
		   s1.at(x0));

  s1.set_z(1.3);
  EXPECT_DOUBLE_EQ(1.3, s1.z());


}
TEST(Prim, CIP) {
  RSTO s1(2.5, 2, 1.1);
  CSTO s2(1.2, 3, CD(0.4, 0.2));
  RGTO g1(0.3, 1, 1.2);
  CGTO g2(0.4, 4, CD(0.1, -0.1));
  double eps = 10.0 * machine_eps();

  
  EXPECT_DOUBLE_EQ(2.9105687018397886, CIP(s1, s1));
  
  EXPECT_DOUBLE_EQ(20.98619989895233, CIP(CSTO(s1), s2).real());
  EXPECT_DOUBLE_EQ(20.98619989895233, CIP(s2, CSTO(s1)).real());
  EXPECT_DOUBLE_EQ(-21.40636768181864, CIP(CSTO(s1), s2).imag());
  EXPECT_DOUBLE_EQ(-21.40636768181864, CIP(s2, CSTO(s1)).imag());
  
  EXPECT_NEAR(0.0764996359892135, CIP(s1, g1), eps);
  EXPECT_NEAR(0.0764996359892135, CIP(g1, s1), eps);

  CSTO c_s1(s1);
  eps *= 100000;
  EXPECT_NEAR(5.562595882704702, CIP(c_s1, g2).real(),      eps);
  EXPECT_NEAR(5.562595882704702, CIP(g2, CSTO(s1)).real(),  eps);
  EXPECT_NEAR(+22.587241177071004, CIP(CSTO(s1), g2).imag(),eps);
  EXPECT_NEAR(+22.587241177071004, CIP(g2, CSTO(s1)).imag(),eps);

  EXPECT_NEAR(0.05270913901892936, CIP(CGTO(g1), g2).real(), eps);
  EXPECT_NEAR(0.05270913901892936, CIP(g2, CGTO(g1)).real(), eps);
  EXPECT_NEAR(0.012359047425198447, CIP(CGTO(g1), g2).imag(), eps);
  EXPECT_NEAR(0.012359047425198447, CIP(g2, CGTO(g1)).imag(), eps);
  
}
TEST(Prim, Normalized) {
  
  RSTO n_s1(2, 1.1, Normalized);
  RSTO s1(1.0, 2, 1.1);
  RSTO n_s2(3, 1.2, Normalized);
  RSTO s2(1.0, 3, 1.2);
  
  EXPECT_DOUBLE_EQ(1.0, CIP(n_s1, n_s1));

  Op<RSTO> op = OpDDr<RSTO>();

  EXPECT_NEAR( CIP(n_s1, op(n_s2)),
	       CIP(s1,   op(s2)) /
	       sqrt(CIP(s1, s1) * CIP(s2, s2)),
	       0.000000000001);

  RGTO n_g1(2, 1.1, Normalized);
  RGTO g1(1.0, 2, 1.1);
  RGTO n_g2(3, 1.2, Normalized);
  RGTO g2(1.0, 3, 1.2);
  
  EXPECT_DOUBLE_EQ(1.0, CIP(n_g1, n_g1));

  Op<RGTO> op_g = OpDDr<RGTO>();

  EXPECT_NEAR( CIP(n_g1, op_g(n_g2)),
	       CIP(g1,   op_g(g2)) /
	       sqrt(CIP(g1, g1) * CIP(g2, g2)),
	       0.000000000001);  
}
TEST(Prim, stream) {
  RSTO n_s1(2, 1.1, Normalized);
  cout << n_s1 << endl;

}
TEST(Prim, Operate) {

  CSTO s1(1.2, 3, complex<double>(1.5));
  CSTO s2 = OperateRm(2, s1);
  CSTO s3 = OperateCst(1.3, s1);

  EXPECT_EQ(5, s2.n());
  EXPECT_DOUBLE_EQ(1.2 * 1.3, s3.c().real());
}
TEST(Prim, ComplexConj) {

  CSTO s(CD(1.1, 1.2), 2, CD(2.1, -2.2));

  CSTO cs(s.ComplexConjugate());
  
  EXPECT_DOUBLE_EQ(1.1, cs.c().real());
  EXPECT_DOUBLE_EQ(-1.2, cs.c().imag());
  EXPECT_EQ(2, cs.n());
  EXPECT_DOUBLE_EQ(2.1, cs.z().real());
  EXPECT_DOUBLE_EQ(2.2, cs.z().imag());

  // if below line are tried to compile, it fails.
  RSTO rsto(1.1 ,2, 3.1);
  RSTO c_rsto(rsto.ComplexConjugate());
  EXPECT_DOUBLE_EQ(1.1, c_rsto.c());
  EXPECT_EQ(2, c_rsto.n());
  EXPECT_DOUBLE_EQ(3.1, c_rsto.z());
}
TEST(CutPrim, First) {

  double eps = pow(10.0, -10.0);

  CSTO csto(1.2, 2, 2.5);
  CutCSTO cut_csto(1.2, 2, 2.5, 10.0);
  
  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      csto.at(1.2), 
		      cut_csto.at(1.2), 
		      eps);

  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      0.0, 
		      cut_csto.at(10.1), 
		      eps);

  CutCSTO cut_csto2(1.1, 3, 1.5, 10.0);
  double sol = 0.0386718749998404;
  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      sol, 
		      CIP(cut_csto, cut_csto2), 
		      eps);

  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      sol, 
		      CIP(cut_csto2, cut_csto), 
		      eps);

  CSTO csto2(1.1, 3, 1.5);
  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      CIP(cut_csto2, cut_csto),
		      CIP(csto2, cut_csto),		      
		      eps);  

  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      CIP(cut_csto2, cut_csto),
		      CIP(cut_csto, csto2),
		      eps);  
}
TEST(Delta, STO) {
  
  double eps = pow(10.0, -10.0);
  DiracDelta<CD> d(2.3);
  CSTO sto(1.1, 2, 2.0);
  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      sto.at(2.3),
		      CIP(d, sto),
		      eps);
}
TEST(LinearComb, Construct) {

  LinearComb<RSTO> rstos;
  rstos += 1.3 * RSTO(1.0, 2, 1.5);
  EXPECT_DOUBLE_EQ(1.3, rstos.coef_i(0));

  rstos[0] = 1.1 * RSTO(1.0, 3, 1.2);
  EXPECT_DOUBLE_EQ(1.1, rstos.coef_i(0));

  LinearComb<RSTO> f(RSTO(2.1, 3, 0.1));
  f += rstos; 
  f += RSTO(1.0, 1, 1.44);
  EXPECT_EQ(3, f.size());
  EXPECT_DOUBLE_EQ(1.2, f.prim_i(1).z());
  EXPECT_DOUBLE_EQ(1.0, f.coef_i(2));
  EXPECT_DOUBLE_EQ(1.44,f.prim_i(2).z());

  LinearComb<RSTO> f_cloned = f.Clone();
  EXPECT_DOUBLE_EQ(1.44, f_cloned.prim_i(2).z());
  
}
TEST(LinearComb, CIP) {

  CSTO s1(1.1, 2, 3.1);
  CSTO s2(1.2, 3, 0.1);
  CSTO s3(1.2, 3, 0.1);

  LinearComb<CSTO> f1;
  LinearComb<CSTO> f2;
  f1 += 1.2 * s1;
  f1 += 1.3 * s2;
  f2 += 0.3 * s3;

  CD expe = 0.3 * 1.2 * CIP(s1, s3) + 0.3 * 1.3 * CIP(s2, s3);
  CD calc = CIP(f1, f2);
  EXPECT_DOUBLE_EQ(expe.real(), calc.real());
  EXPECT_DOUBLE_EQ(expe.imag(), calc.imag());

  EXPECT_DOUBLE_EQ(CIP(s1, f2).real(), 0.3 * CIP(s1, s3).real());
  EXPECT_DOUBLE_EQ(CIP(f2, s1).real(), 0.3 * CIP(s1, s3).real());

}
TEST(LinearComb, op) {

  LinearComb<RSTO> f1;
  f1 += 1.2 * RSTO(1.1, 2, 0.2);
  f1 += 1.3 * RSTO(0.9, 3, 0.3);

  Op<RSTO> op_s = OpRM<RSTO>(2);
  LinearComb<RSTO> f2 = op_s(f1);
  EXPECT_EQ(2, f2.size());
  EXPECT_EQ(4, f2.prim_i(0).n());
  EXPECT_EQ(5, f2.prim_i(1).n());
  EXPECT_DOUBLE_EQ(1.3, f2.coef_i(1));

  Op<RSTO> op_rm3 = OpRM<RSTO>(3);
  LinearComb<RSTO> f3 = op_rm3(RSTO(1.0, 2, 1.0));
  EXPECT_EQ(1, f3.size());
  EXPECT_EQ(5, f3.prim_i(0).n());
  
  Op<RSTO> op_ddr = OpDDr<RSTO>();
  LinearComb<RSTO> df1 = op_ddr(RSTO(1.1, 2, 1.2));

  EXPECT_EQ(2, df1.size());

  LinearComb<RSTO> df2 = OpDDr2<RSTO>()(RSTO(1.1, 2, 1.2));
   
  EXPECT_EQ(3, df2.size());  
  
}
TEST(LinearComb, AtX) {

  RSTO s1(1.1, 1, 0.2);
  RSTO s2(1.2, 1, 0.3);
  LinearComb<RSTO> sto;
  sto += 0.3 * s1;
  sto += 0.4 * s2;
  
  double x0(3.0);
  EXPECT_DOUBLE_EQ(x0 * (0.3 * 1.1 * exp(-0.2 * x0) +
			 0.4 * 1.2 * exp(-0.3*x0)),
		   sto.at(x0));
  
}
TEST(LinearComb, setter) {

  LinearComb<RSTO> stos;
  stos.resize(2);

  EXPECT_EQ(2, stos.size());

  stos.set_coef_i(0, 1.1);
  stos.set_coef_i(1, 2.2);

  stos.ref_prim_i(0) = RSTO(0.1, 2, 0.2);
  stos.ref_prim_i(1) = RSTO(0.2, 3, 0.25);

  EXPECT_DOUBLE_EQ(1.1, stos.coef_i(0));
  EXPECT_DOUBLE_EQ(0.1, stos.prim_i(0).c());
  EXPECT_DOUBLE_EQ(0.25, stos.prim_i(1).z());

}
TEST(LinearComb, dbasis_STO) {

  double dz = 0.001;
  RSTO s1(2, 0.3, Normalized);
  RSTO s1_p(2, 0.3 + dz, Normalized);
  RSTO s1_m(2, 0.3 - dz, Normalized);
  
  LinearComb<RSTO> d_s1;
  d_s1.SetD1Normalized(s1);

  LinearComb<RSTO> dd_s1;
  dd_s1.SetD2Normalized(s1);

  double x0 = 3.3;
  cout << s1.at(x0) << endl;
  cout << d_s1.at(x0) << endl;
	    
  EXPECT_NEAR( (s1_p.at(x0)-s1_m.at(x0)) / (2.0 * dz),
	       d_s1.at(x0),  0.00001);
  EXPECT_NEAR( (s1_p.at(x0) + s1_m.at(x0) - 2.0 * s1.at(x0))
	       /(dz * dz),  dd_s1.at(x0), 0.0001);
  
  typedef std::complex<double> CD;
  CSTO g1(2, CD(0.3), Normalized);
  CSTO g1_p(2, CD(0.3 + dz), Normalized);
  CSTO g1_m(2, CD(0.3 - dz), Normalized);
  
  LinearComb<CSTO> d_g1; 
  d_g1.SetD1Normalized(g1);
  LinearComb<CSTO> dd_g1;
  dd_g1.SetD2Normalized(g1);

  CD cx0(3.3, 0.0);
  EXPECT_NEAR( ((g1_p.at(cx0) - g1_m.at(cx0)) / (2.0 * dz)).real(),
	       (d_g1.at(cx0)).real(),  0.00001);
  EXPECT_NEAR( ((g1_p.at(cx0) + g1_m.at(cx0))
		-2.0 * g1.at(cx0)).real()
	       /(dz * dz),
	       dd_g1.at(cx0).real(), 0.0001);
}
TEST(LinearComb, dbasis_GTO) {

  double dz = 0.0001;
  typedef std::complex<double> CD;
  CGTO g1(2, CD(0.3), Normalized);
  CGTO g1_p(2, CD(0.3 + dz), Normalized);
  CGTO g1_m(2, CD(0.3 - dz), Normalized);
  
  LinearComb<CGTO> d_g1  = D1Normalized(g1);
  LinearComb<CGTO> dd_g1 = D2Normalized(g1);

  CD cx0(3.3, 0.0);
  EXPECT_NEAR( ((g1_p.at(cx0) - g1_m.at(cx0)) /
		(2.0 * dz)).real(),
	       d_g1.at(cx0).real(),  0.00001);
  EXPECT_NEAR( ((g1_p.at(cx0) + g1_m.at(cx0))
		-2.0 * g1.at(cx0)).real()
	       /(dz * dz),
	       dd_g1.at(cx0).real(), 0.0001);
}
TEST(LinearComb, ComplexConjugate) {

  CGTO g1(CD(1.1, 0.3), 1, CD(0.2, 0.5));
  CGTO g2(1.2, 1, 0.3);
  LinearComb<CGTO> gto;
  gto += 0.3 * g1;
  gto += 0.4 * g2;

  LinearComb<CGTO> c_gto = gto.ComplexConjugate();
  EXPECT_EQ(2, c_gto.size());
  EXPECT_DOUBLE_EQ(gto.at(0.3).real(), +c_gto.at(0.3).real());
  EXPECT_DOUBLE_EQ(gto.at(0.3).imag(), -c_gto.at(0.3).imag());

}
TEST(LinearComb, WithDelta) {

  double eps = pow(10.0, -10.0);
  double r0(0.5);

  LinearComb<DiracDelta<CD> > delta;
  delta += 1.2 * DiracDelta<CD>(r0);

  LinearComb<CSTO> stos;
  stos += 1.1 * CSTO(1.3, 2, 1.2);

  EXPECT_PRED_FORMAT3(AssertComplexNear, 
		      1.2 * stos.at(r0),
		      CIP(stos, delta),
		      eps);  

}
TEST(Op, Rm) {

  Op<RSTO> op = OpRM<RSTO>(3);
  RSTO s1(2, 3.0, Normalized);
  LinearComb<RSTO> s2 = op(s1);
  LinearComb<RSTO> s3 = op(s2);
  EXPECT_EQ(2 + 3 + 3, s3.prim_i(0).n());

}
TEST(Op, Cst) {

  Op<RSTO> op = OpCst<RSTO>(2.2);
  RSTO s1 (1.1, 2, 3.0);
  LinearComb<RSTO> s2 = op(s1);
  EXPECT_DOUBLE_EQ(2.2 * 1.1, s2.prim_i(0).c());

}
TEST(Op, DDr) {

  // d/dr 2r^3 exp(-4r)
  // = (6r^2 -8r^3)exp(-4r)
  // d2/dr2 2^3 exp(-4r)
  // = (12r -24r^2 -24r^2 + 32r^3) exp(-4r)
  // = (12r -48r^2 + 32r^3) exp(-4r)
   
  Op<RSTO> op = OpDDr<RSTO>();
  RSTO s1(2.0, 3, 4.0);

  LinearComb<RSTO> ds = OperateDDr(s1);

  LinearComb<RSTO> s2 = op(s1);
  LinearComb<RSTO> s3_1 = OpDDr<RSTO>()(s2);
  LinearComb<RSTO> s3_2 = OpDDr2<RSTO>()(s1);
  double r = 2.2;

  EXPECT_DOUBLE_EQ( (6*r*r - 8*r*r*r) * exp(-4*r), 
		    ds.at(r));  

  EXPECT_DOUBLE_EQ( (6*r*r - 8*r*r*r) * exp(-4*r), 
		    s2.at(r));
  EXPECT_DOUBLE_EQ( (12.0 - 48*r + 32*r*r)*r*exp(-4*r), 
		    s3_1.at(r));
  EXPECT_DOUBLE_EQ(s3_1.at(r), s3_2.at(r));
}
TEST(Op, algebra) {

  Op<RSTO> op;
  op.Add(1.3, OpCst<RSTO>(2.2));;
  op.Add(2.4, OpRM<RSTO>(3));
  op.ScalarProduct(1.1);

  RSTO s1(2, 0.5, Normalized);

  double x = 2.3;
  EXPECT_DOUBLE_EQ( op(s1).at(x),
		    1.1 * 1.3 * 2.2 * s1.at(x) +
		    1.1 * 2.4 * OpRM<RSTO>(3)(s1).at(x));
		    

}
TEST(HAtom, EigenState) {

  HLikeAtom<CD> hatom00(1, 1.0, 0);
  HLikeAtom<CD> hatom10(2, 1.0, 0);
  HLikeAtom<CD> hatom11(2, 1.0, 1);
  HLikeAtom<CD> hatom_3d(3, 1.0, 2);

  LinearComb<CSTO> f00 = hatom00.EigenState();
  LinearComb<CSTO> f10 = hatom10.EigenState();
  LinearComb<CSTO> f11 = hatom11.EigenState();
  LinearComb<CSTO> psi_3d = hatom_3d.EigenState();
  EXPECT_EQ(1, f00.size());

  EXPECT_EQ(0.0, abs(CIP(f00, f10)));
  EXPECT_EQ(0.0, abs(CIP(f00, f10)));

  EXPECT_DOUBLE_EQ(1.0, CIP(f00, f00).real());
  EXPECT_DOUBLE_EQ(1.0, CIP(f10, f10).real());
  EXPECT_DOUBLE_EQ(1.0, CIP(f11, f11).real());
  EXPECT_DOUBLE_EQ(1.0, CIP(psi_3d, psi_3d).real());

  double eps = machine_eps() * 10;

  Op<CSTO> hop = hatom_3d.Hamiltonian<CSTO>();
  Op<CSTO> hop_m05 = hatom_3d.HMinusEnergy<CSTO>(0.5);

  EXPECT_NEAR
    (hatom_3d.EigenEnergy(),
     CIP(psi_3d,
	 hop(psi_3d.prim_i(0))).real(),
     eps);

  /*
  EXPECT_NEAR(hatom00.EigenEnergy(),
	      CIP(f00, Op(hatom00.Hamiltonian<CSTO>(), f00)).real(), eps);
  EXPECT_NEAR(0.0,
	      CIP(f00, Op(hatom10.Hamiltonian<CSTO>(), f10)).real(), eps);
  EXPECT_NEAR(hatom10.EigenEnergy(),
	      CIP(f10, Op(hatom10.Hamiltonian<CSTO>(), f10)).real(), eps);
  EXPECT_NEAR(hatom11.EigenEnergy(),
	      CIP(f11, Op(hatom11.Hamiltonian<CSTO>(), f11)).real(), eps);
  */
}
TEST(HAtom, DipoleInitLength) {

  HLikeAtom<CD> hatom(2, 1.0, 1);
  LinearComb<CSTO> mu_phi_00 = hatom.DipoleInitLength(2);
  EXPECT_EQ(1, mu_phi_00.size());
  
}



