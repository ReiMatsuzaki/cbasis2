#include <iostream>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <boost/timer.hpp>
#include <gtest/gtest.h>

#include "gtest_plus.hpp"

#include "exp_func.hpp"
#include "normal_exp.hpp"
#include "cut_exp.hpp"
#include "delta.hpp"

#include "cip_exp.hpp"

#include "hatom_h.hpp"

// #include "hatom_h.hpp"

using namespace std;
using namespace cbasis;
typedef std::complex<double> CD;

TEST(Func, NTermFunc) {

  NTermFunc<3, RSTO>::type stos = func_add_func(RSTO(1.1,2,1.1), 
						func_add_func(RSTO(1.0, 3, 0.2), 
							      RSTO(2.0, 2, 0.1)));
}
TEST(Func, ExpFuncConstruct) {
  
  RSTO s0(1.0, 1, 1.0);
  RSTO s1(1.0, 2, 1.0);
  RGTO g1(1.0, 2, 1.0);

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


  RSTO s0_copy(s0);
  EXPECT_DOUBLE_EQ(s0.at(0.1), s0_copy.at(0.1));
  RSTO s0_ass = s0;
  EXPECT_DOUBLE_EQ(s0.at(0.1), s0_ass.at(0.1));
}
TEST(Func, accessor) {

  RSTO s1(1.1, 2, 0.2);
  double x0(3.0);
  EXPECT_DOUBLE_EQ(1.1 * x0 * x0 * exp(-0.2 * x0), 
		   s1.at(x0));

  s1.set_z(1.3);
  EXPECT_DOUBLE_EQ(1.3, s1.z());


}
TEST(Func, DerivExpFunc) {

  RSTO s(0.1, 2, 0.3);
  RSTO ds = s.DerivParamOne();
  RSTO dds = s.DerivParamTwo();

  EXPECT_DOUBLE_EQ(-0.1, ds.c());
  EXPECT_EQ(3, ds.n());

  EXPECT_DOUBLE_EQ(0.1, dds.c());
  EXPECT_EQ(4, dds.n());

}
TEST(Func, stream) {
  RSTO n_s1(2.1, 2, 1.1);
  cout << n_s1 << endl;
}
TEST(Func, ComplexConj) {

  CSTO s(CD(1.1, 1.2), 2, CD(2.1, -2.2));

  CSTO cs(s); cs.SetComplexConjugate();
  
  EXPECT_DOUBLE_EQ(1.1, cs.c().real());
  EXPECT_DOUBLE_EQ(-1.2, cs.c().imag());
  EXPECT_EQ(2, cs.n());
  EXPECT_DOUBLE_EQ(2.1, cs.z().real());
  EXPECT_DOUBLE_EQ(2.2, cs.z().imag());

  // if below line are tried to compile, it fails.
  RSTO rsto(1.0 ,2, 3.1);
  RSTO c_rsto(rsto); c_rsto.SetComplexConjugate();
  EXPECT_DOUBLE_EQ(1.0, c_rsto.c());
  EXPECT_EQ(2, c_rsto.n());
  EXPECT_DOUBLE_EQ(3.1, c_rsto.z());
}
TEST(Func, Product) {

  CSTO s1(1.2, 3, std::complex<double>(1.5));
  s1.SetScalarProd(2.0);
  s1.SetRmProd(2);

  EXPECT_EQ(5, s1.n());
  EXPECT_DOUBLE_EQ(2.4, s1.c().real());
}
TEST(Func, CutExpConstruct) {

  CSTO csto(1.2, 2, 2.5);
  CutCSTO cut_csto(1.2, 2, 2.5, 10.0);

  EXPECT_C_EQ(csto.at(1.2), cut_csto.at(1.2));
  EXPECT_C_EQ(0.0, cut_csto.at(10.1));


}
TEST(Func, LinFuncConstruct) {

  RSTO s1(1.0, 2, 1.2);
  RSTO s2(1.0, 3, 1.4);
  LinFunc<RSTO> stos; 
  stos.Add(2.5, s1);
  stos.Add(1.2, s2);

  double r0(0.3);

  EXPECT_DOUBLE_EQ(2.5*s1.at(r0) + 1.2*s2.at(r0),
		   stos.at(r0));

}
TEST(Func, LinFuncSetter) {
  LinFunc<RSTO> func;
  //RSTO s1(1.1, 2, 1.3); func.Add(1.2, s1);
  func.Add(1.2, RSTO(1.1, 2, 1.3));
  RSTO s2(1.3, 3, 1.4); func.Add(1.2, s2);
  func.SetScalarProd(1.4);
  EXPECT_DOUBLE_EQ(1.2*1.4, func.begin()->first);

  func.SetNormalize();
  EXPECT_DOUBLE_EQ(1.0, CIP(func, func));
}
TEST(Func, NormalExpFunc) {

  double dz = 0.0003;
  NormalRSTO s1(2, 0.3);
  NormalRSTO s1_p(2, 0.3 + dz);
  NormalRSTO s1_m(2, 0.3 - dz);

  EXPECT_DOUBLE_EQ(1.0, CIP(s1, s1));

  NormalRSTO::FuncDerivOne d_s1 = s1.DerivParamOne();
  NormalRSTO::FuncDerivTwo dd_s1= s1.DerivParamTwo();

  double x0 = 3.3;
	    
  EXPECT_NEAR( (s1_p.at(x0)-s1_m.at(x0)) / (2.0 * dz),
	       CIP(RDelta(x0), d_s1),  0.00001);

  EXPECT_NEAR( (s1_p.at(x0) + s1_m.at(x0) - 2.0 * s1.at(x0)) / (dz*dz),
	       CIP(RDelta(x0), dd_s1), 0.00001);


  s1.set_z(0.4);
  EXPECT_DOUBLE_EQ(1.0, CIP(s1, s1));

}
TEST(CIP, ExpFunc) {
  RSTO s1(2.5, 2, 1.1);
  CSTO s2(1.2, 3, CD(0.4, 0.2));
  RGTO g1(0.3, 1, 1.2);
  CGTO g2(0.4, 4, CD(0.1, -0.1));
  double eps = pow(10.0, -9.0);

  
  EXPECT_DOUBLE_EQ(2.9105687018397886, CIP(s1, s1));

  CD sol(20.98619989895233, -21.40636768181864);
  EXPECT_C_NEAR(sol, CIP(CSTO(s1), s2), eps);

  sol = CD(20.98619989895233, -21.40636768181864);
  EXPECT_C_EQ(sol, CIP(s2, CSTO(s1)));
  
  EXPECT_NEAR(0.0764996359892135, CIP(s1, g1), eps);
  EXPECT_NEAR(0.0764996359892135, CIP(g1, s1), eps);

  CSTO c_s1(s1);
  EXPECT_NEAR(5.562595882704702, CIP(c_s1, g2).real(),      eps);
  EXPECT_NEAR(5.562595882704702, CIP(g2, CSTO(s1)).real(),  eps);
  EXPECT_NEAR(+22.587241177071004, CIP(CSTO(s1), g2).imag(),eps);
  EXPECT_NEAR(+22.587241177071004, CIP(g2, CSTO(s1)).imag(),eps);

  EXPECT_NEAR(0.05270913901892936, CIP(CGTO(g1), g2).real(), eps);
  EXPECT_NEAR(0.05270913901892936, CIP(g2, CGTO(g1)).real(), eps);
  EXPECT_NEAR(0.012359047425198447, CIP(CGTO(g1), g2).imag(), eps);
  EXPECT_NEAR(0.012359047425198447, CIP(g2, CGTO(g1)).imag(), eps);


}
TEST(CIP, CutExpFunc) {

  CutCSTO cut_csto(1.2, 2, 2.5, 10.0);
  CutCSTO cut_csto2(1.1, 3, 1.5, 10.0);
  double sol = 0.0386718749998404;
  EXPECT_C_EQ(sol, CIP(cut_csto, cut_csto2));
  EXPECT_C_EQ(sol, CIP(cut_csto2, cut_csto));

  CSTO csto2(1.1, 3, 1.5);
  EXPECT_C_EQ(CIP(csto2, cut_csto), 
	      CIP(cut_csto, csto2));
  EXPECT_C_EQ(CIP(cut_csto2, cut_csto), 
	      CIP(cut_csto, csto2));

  double r0 = 4.0;
  CutCSTO cut_csto3(1.0, 2, 0.2, r0);
  CSTO csto3(1.0, 2, 0.2);
  EXPECT_C_EQ(csto3.at(r0), cut_csto3.at(r0));
  CDelta d3(3.0);
  EXPECT_C_EQ(cut_csto.at(3.0), CIP(d3, cut_csto));

}
TEST(CIP, Normalized) {
  
  RSTO s1(1.0, 2, 1.1);
  RGTO g1(1.0, 2, 1.1);
  CSTO s2(1.0, 3, CD(2.1, -0.3));
  
  EXPECT_DOUBLE_EQ(CIP(s1, s1), CNorm(s1)*CNorm(s1));
  EXPECT_DOUBLE_EQ(CIP(s1, s1), CNorm2(s1));
  EXPECT_DOUBLE_EQ(1.0, CNorm(CNormalize(s1)));
  EXPECT_DOUBLE_EQ(1.0, CNorm(CNormalize(g1)));


/*
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
*/
}
TEST(CIP, Delta) {

  double r0(2.4);
  RDelta d0(r0);
  RSTO  sto(1.1, 2, 3.2);

  EXPECT_DOUBLE_EQ(sto.at(2.4), CIP(d0, sto));
  EXPECT_DOUBLE_EQ(sto.at(2.4), CIP(sto, d0));
  
}
TEST(CIP, LinFunc) {

  LinFunc<CSTO> f;
  CD c0(1.0, 1.2);
  CSTO u0(1.2, 2, CD(0.4, -0.4));
  f.Add(c0, u0);
  CD c1(0.3, 0.2);
  CSTO u1(1.0, 1, CD(0.1, -0.2));
  f.Add(c1, u1);

  LinFunc<CSTO> g;
  CD d0(0.1, 1.2);
  CSTO v0(0.3, 1, CD(0.1, -0.5));
  g.Add(d0, v0);
  CD d1(0.1, 0.5);
  CSTO v1(0.8, 1, CD(0.8, -0.0));
  g.Add(d1, v1);  

  EXPECT_C_EQ(CIP(f, g),
	      c0*d0*CIP(u0, v0)+c1*d0*CIP(u1, v0)+c0*d1*CIP(u0, v1)+c1*d1*CIP(u1, v1));
  
}
/*

TEST(CIP, OpRm) {


  RSTO s1(2.0, 3, 4.0);
  RSTO s2(1.1, 2, 2.2);


  EXPECT_DOUBLE_EQ(CIP(s1, OP_lin(RRm(2), s2)), 
		   CIP(s1, RRm(2), s2));

  RGTO g1(2.0, 4, 1.0);
  RGTO g2(2.0, 4, 1.0);
  
  EXPECT_DOUBLE_EQ(CIP(g1, OP_lin(RRm(3), g2)), 
		   CIP(g1, RRm(3), g2));

}
TEST(CIP, OpD2) {

  RSTO s1(2.0, 3, 4.0);
  RSTO s2(1.1, 2, 2.2);


  EXPECT_DOUBLE_EQ(CIP(s1, OP_lin(RD2(), s2)), 
		   CIP(s1, RD2(), s2));

  RGTO g1(2.0, 4, 1.0);
  RGTO g2(2.0, 4, 1.0);
  
  EXPECT_DOUBLE_EQ(CIP(g1, OP_lin(RD2(), g2)), 
		   CIP(g1, RD2(), g2));
}
TEST(CIP, time) {
  RSTO s1(1.2, 5, 0.3);
  LinFunc<RSTO> ss; 
  for(int i = 0; i < 100; i++)
    ss.Add(1.2, s1);

  int num(100);

  boost::timer t0;
  for(int i = 0; i < num; i++) 
    CIP(s1, OpD2(), s1);  
  cout << "t[CIP(A,O,B)] = " << t0.elapsed() << endl;

  boost::timer t1;
  for(int i = 0; i < num; i++) 
    CIP(ss, OP_lin(OpD2(), ss));
  cout << "t[CIP(A,O[B])] = " << t1.elapsed() << endl;

}
*/ 
TEST(CIP, Op_algebra) {

  RSTO s1(1.2, 5, 3.0);
  RSTO s2(1.1, 2, 5.0);

  EXPECT_DOUBLE_EQ(
		   CIP(s1, op_add_op(RRm(2), RD1()), s2),
		   CIP(s1, RRm(2), s2) + CIP(s1, RD1(), s2)
		   );

  EXPECT_DOUBLE_EQ(
		   CIP(s1, scalar_mult_op(0.3, RD1()), s2),
		   0.3 * CIP(s1, RD1(), s2)
		   );

  EXPECT_DOUBLE_EQ(
		   CIP(s1, op_add_op(RRm(2), scalar_mult_op(0.3, RD1())), s2),
		   CIP(s1, RRm(2), s2) + 0.3*CIP(s1, RD1(), s2)
		   );

}
TEST(CIP, Func_algebra) {

  RSTO s1(0.1, 2, 0.4);
  RSTO s2(0.2, 3, 0.3);
  RSTO s3(0.3, 2, 0.1);

  EXPECT_DOUBLE_EQ(CIP(s1, func_add_func(s2, s3)),
		   CIP(s1, s2) + CIP(s1, s3));

  EXPECT_DOUBLE_EQ(CIP(func_add_func(s1, s2), s3),
		   CIP(s1, s3) + CIP(s2, s3));

}
TEST(CIP, symmetry) {
  using namespace cbasis::h_atom_rad;

  RSTO s1(1.1, 2, 0.3);
  RSTO s2(1.2, 2, 0.4);
  RSTO s3(1.3, 2, 0.5);

  EXPECT_DOUBLE_EQ(CIP(s1, RRm(2), s2),
		   CIP(s2, RRm(2), s1));

  EXPECT_DOUBLE_EQ(CIP(s1, RD2(), s2),
		   CIP(s2, RD2(), s1));

  EXPECT_DOUBLE_EQ(CIP(s1, HOp<double, 0>()(), func_add_func(s2, s3)),
		   -0.5*CIP(s1, RD2(), s2) -0.5*CIP(s1, RD2(), s3)
		   -1.0*CIP(s1, RRm(-1), s2) -1.0*CIP(s1, RRm(-1), s3));

  EXPECT_DOUBLE_EQ(CIP(s1, HOp<double, 0>()(), s2),
		   CIP(s2, HOp<double, 0>()(), s1));

}
TEST(HAtom, hop) {
  using namespace cbasis::h_atom_rad;
  EXPECT_DOUBLE_EQ(1.0, CNorm2(EigenFunc<double, 1, 0>()()));
  EXPECT_DOUBLE_EQ(0.0, CIP(EigenFunc<double, 1, 0>()(),
			    EigenFunc<double, 2, 0>()()));
  EXPECT_DOUBLE_EQ(0.0, CIP(EigenFunc<double, 2, 0>()(),
			    EigenFunc<double, 1, 0>()()));

  EXPECT_DOUBLE_EQ(-0.5, CIP(EigenFunc<double, 1, 0>()(),
			     HOp<double, 0>()(),
			     EigenFunc<double, 1, 0>()()));

  EXPECT_DOUBLE_EQ(-0.125, CIP(EigenFunc<double, 2, 0>()(),
			       HOp<double, 0>()(),
			       EigenFunc<double, 2, 0>()()));

  EXPECT_DOUBLE_EQ(0.0, CIP(EigenFunc<double, 2, 0>()(),
			    HOp<double, 0>()(),
			    EigenFunc<double, 1, 0>()()));

  EXPECT_DOUBLE_EQ(0.0, CIP(EigenFunc<double, 1, 0>()(),
			    HOp<double, 0>()(),
			    EigenFunc<double, 2, 0>()()));

  EXPECT_DOUBLE_EQ(-0.125, CIP(EigenFunc<double, 2, 1>()(),
			       HOp<double, 1>()(),
			       EigenFunc<double, 2, 1>()()));
			    
}
TEST(HAtom, length) {
  using namespace cbasis::h_atom_rad;
  typedef Length<double, 1, 0, 1> LengthForm;
  LengthForm::type s = LengthForm()();
  EXPECT_EQ(2, s.n());
}


/*

TEST(OP, RRm) {


  RSTO s1(2.3, 2, 1.1);
  RSTO r2_s1 = OP_lin(RRm(2), s1);

  EXPECT_EQ(4, r2_s1.n());

}
TEST(OP, RD1) {

  // d/dr 2r^3 exp(-4r)
  // = (6r^2 -8r^3)exp(-4r)
  
  RSTO s1(2.0, 3, 4.0);
  FuncAdd<RSTO, RSTO> s2 = OP_lin(RD1(), s1);
  double r = 2.2;
  EXPECT_DOUBLE_EQ((6.0*r*r-8*r*r*r)*exp(-4*r),
		   s2.at(r));

  // d/dr 2r^3 exp(-4r^2)
  // = (6r^2 -16r^4) *exp(-4rr)
  RGTO g1(2.0, 3, 4.0);
  FuncAdd<RGTO, RGTO> g2 = OP_lin(RD1(), g1);
  EXPECT_DOUBLE_EQ((6.0*r*r-16*r*r*r*r)*exp(-4*r*r),
		   g2.at(r));

}
TEST(OP, RD2) {

  RSTO s(2.0, 3, 4.0);
  
  //  cout << s1 << endl;
  //  cout << s2 << endl;

  double r0(2.1);
  EXPECT_DOUBLE_EQ(OP_lin(RD1(), OP_lin(RD1(), s)).at(r0),
		   OP_lin(RD2(), s).at(r0));
  
}
TEST(OP_lin, Add) {


  RSTO s(1.1, 2, 1.2);
  double r0(2.1);
  EXPECT_DOUBLE_EQ(OP_lin(RRm(2), s).at(r0) + OP_lin(RD2(), s).at(r0),
		   OP_lin(op_add_op(RRm(2), RD2()), s).at(r0));


}
TEST(OP_lin, ScalarProd) {

  RSTO s(1.1, 2, 1.2);

  double r0(0.4);
  EXPECT_DOUBLE_EQ(
		   0.3 * (OP_lin(RRm(2), s).at(r0) + OP_lin(RD2(), s).at(r0)),
		   OP_lin(scalar_mult_op(0.3, op_add_op(RRm(2), RD2())), s).at(r0));
}

TEST(HAtom, eigenstate) {
  HLikeAtom<double> hatom;
  
  HPsi<1,0,double> f10(hatom);
  HPsi<2,0,double> f20(hatom);
  //  HPsi<2,1,double> f21;
  EXPECT_DOUBLE_EQ(1.0,
		   CIP(f10.value, f10.value )
		   );
  EXPECT_DOUBLE_EQ(0.0,
		   CIP(f20.value, f10.value )
		   );

  HOp<0,double> hop(hatom);
  EXPECT_DOUBLE_EQ(EigenEnergy<1>(hatom),
		   CIP(f10.value, hop.value, f10.value)
		   );

  EXPECT_DOUBLE_EQ(0.0,
		   CIP(f20.value, hop.value, f10.value)
		   );

  HminusEOp<0, double> Lop(hatom, 0.7);
  EXPECT_DOUBLE_EQ(-1.2,
		   CIP(f10.value, Lop.value, f10.value)
		   );

  HLength<1, 0, 1, double> mu_phi(hatom);
  EXPECT_EQ(2, mu_phi.value.n());

}
TEST(HAtom, p_wave) {

  HLikeAtom<double> hatom;
  HPsi<2,1,double> f(hatom);
  HOp<1,double> hop(hatom);

  EXPECT_DOUBLE_EQ(EigenEnergy<2>(hatom),
		   CIP(f.value, hop.value, f.value)
		   );

}
*/
int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}

/*
int practice_fusion () {

  using boost::fusion::vector;
  using boost::fusion::at_c;

  vector<int, double> x(1, 0.2);
  std::cout << at_c<0>(x) << endl;
  return 0;
}
int main ()
{
  ExpFunc<double, 1> stoA(1.0, 2, 2.5);
  ExpFunc<double, 1> stoB(2.0, 3, 2.5);
  OpRm r2(2);

  double ab = CIP(stoA, stoB);
  double abc = CIP(stoA, r2, stoB);

  cout << ab << endl;
  cout << abc<< endl;

  std::vector<double> cs;
  cs.push_back(1.1); cs.push_back(1.2); cs.push_back(0.7);
  OpLin<double , OpRm, OpRm, OpRm> op_lin(cs, r2, r2, r2);
  op_lin.print();
  //  cout << op_lin << endl;

  cout << endl;
  cout << CIP(stoA, op_lin, stoB) << endl;

  std::stringstream ss;
  
  return 0;
}
*/
