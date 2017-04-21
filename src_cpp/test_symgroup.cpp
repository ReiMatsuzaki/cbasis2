#include <iostream>
#include <gtest/gtest.h>
#include "../utils/gtest_plus.hpp"
#include "../utils/macros.hpp"
#include "symgroup.hpp"

using namespace std;
using namespace Eigen;
using namespace cbasis;

dcomplex random_complex() {

  double r0 = (double)rand() / ((double)RAND_MAX+1);
  double i0 = (double)rand() / ((double)RAND_MAX+1);
  return dcomplex(r0, i0);

}
void RandomPrim(PrimGTO* g) {
  int maxn(4);
  g->nx = rand() % maxn;
  g->ny = rand() % maxn;
  g->nz = rand() % maxn;
  g->x = random_complex();
  g->y = random_complex();
  g->z = random_complex();
}
bool IsSameOp(SymOp a, SymOp b, int num = 5) {

  PrimGTO x0, ay, by;
  for(int i = 0; i < num; i++) {
    RandomPrim(&x0);

    bool prim_a, prim_b;
    int  sig_a,  sig_b;
    a->getOp(x0, &ay, &sig_a, &prim_a);
    b->getOp(x0, &by, &sig_b, &prim_b);
    if(!prim_a || !prim_b) {
      string msg; SUB_LOCATION(msg);
      msg += "result is not primitive";
      throw runtime_error(msg);
    }

    if(!IsNear(ay, by))
      return false;
    if(sig_a != sig_b)
      return false;
  }
  return true;
}
void ExpectOpEq(SymOp a, SymOp b, int num = 5) {
  EXPECT_TRUE(IsSameOp(a, b, num))
    << a->str() << endl
    << b->str() << endl;
}
void ExpectSymmetryGroup(SymmetryGroup g) {

  typedef _SymmetryGroup::ItSymOp It;
  It i0 = g->sym_op_.begin();
  It end = g->sym_op_.end();
  for(It it = i0; it != end; ++it) {
    for(It jt = i0; jt != end; ++jt) {
      SymOp ij = prod(*it, *jt);
      bool res(false);
      for(It kt = i0; kt != end; ++kt) {
	if(IsSameOp(ij, *kt)) {
	  res = true;
	  break;
	}
      }
      EXPECT_TRUE(res);
    }
  }
}

TEST(first, first) {
  EXPECT_EQ(2, 1+1);
}
TEST(SymOp, Id) {

  PrimGTO a(2, 1, 3, 0.1, 0.2, 0.3);
  PrimGTO b(a);
  PrimGTO c(2, 1, 3, 0.0, 0.1, 0.3);
  ISymOp *id = new Id();
  //  bool is_prim;
  //  int sig;
  
  EXPECT_EQ(1, id->Op(a, b));
  EXPECT_EQ(0, id->Op(a, c));

  delete id;
  
}
TEST(SymOp, C2) {

  ISymOp *cx2 = new Cyclic(CoordX, 2);
  ISymOp *cy2 = new Cyclic(CoordY, 2);
  ISymOp *cz2 = new Cyclic(CoordZ, 2);

  PrimGTO s(0, 0, 0, 0.1, 0.2, 0.3);
  PrimGTO px(1, 0, 0, 0.1, 0.2, 0.3);
  PrimGTO px_x(1, 0, 0, 0.1, -0.2, -0.3);
  PrimGTO px_y(1, 0, 0, -0.1, +0.2, -0.3);
  PrimGTO px_z(1, 0, 0, -0.1, -0.2, 0.3);

  EXPECT_EQ(+1, cx2->Op(px, px_x));
  EXPECT_EQ(-1, cy2->Op(px, px_y));
  EXPECT_EQ(-1, cz2->Op(px, px_z));

  PrimGTO py(0, 1, 0, 0.1, 0.2, 0.3);
  PrimGTO pz(0, 0, 1, 0.1, 0.2, 0.3);
  
  delete cx2;
  delete cy2;
  delete cz2;

}
TEST(SymOp, C2_C4) {

  Coord axis_list[3] = {CoordX, CoordY, CoordZ};

  for(int i = 0; i < 3; i++) {
    Coord axis = axis_list[i];
    SymOp C2 = cyclic(axis, 2);
    SymOp C4 = cyclic(axis, 4);

    ExpectOpEq(id(), mult(C2, 2));
    ExpectOpEq(id(), prod(C2, C2));
    ExpectOpEq(id(), mult(C4, 4));
    ExpectOpEq(C2, mult(C4, 2));
    ExpectOpEq(C2, prod(C4, C4));

  }
}
TEST(SymOp, Sig) {
  Coord axis_list[3] = {CoordX, CoordY, CoordZ};
  for(int i = 0; i < 3; i++) {
    Coord coord = axis_list[i];
    SymOp SIG  = reflect(coord);
    ExpectOpEq(id(), prod(SIG, SIG));
  }
}
TEST(SymOp, SigC2) {

  Coord axis_list[3] = {CoordX, CoordY, CoordZ};
  for(int i = 0; i < 3; i++) {
    Coord axis = axis_list[i];
    SymOp C2   = cyclic(axis, 2);
    SymOp SIG  = reflect(axis);
    ExpectOpEq(inv(), prod(C2, SIG));
    ExpectOpEq(inv(), prod(SIG, C2));
  }  
}
TEST(SymOp, D2h) {
  SymOp C2z = cyclic(CoordZ, 2);
  SymOp C2x = cyclic(CoordX, 2);
  SymOp C2y = cyclic(CoordY, 2);
  SymOp SIGxy = reflect(CoordZ);
  SymOp SIGzx = reflect(CoordY);
  SymOp SIGyz = reflect(CoordX);

  ExpectOpEq(C2y, prod(C2z, C2x));
  ExpectOpEq(C2x, prod(C2z, C2y));
  ExpectOpEq(inv(), prod(C2z, SIGxy));
  ExpectOpEq(C2z, prod(C2x, C2y));
  ExpectOpEq(SIGzx, prod(C2x, SIGxy));
  ExpectOpEq(SIGyz, prod(C2y, SIGxy));

}

TEST(SymGroup, C1) {

  SymmetryGroup C1 = SymmetryGroup_C1();
  ExpectSymmetryGroup(C1);  
    
  EXPECT_TRUE(C1->Non0_Scalar(0, 0));
  EXPECT_TRUE(C1->Non0_Z(0, 0));
  EXPECT_TRUE(C1->Non0_4(0, 0, 0, 0));
}
TEST(SymGroup, Cs) {

  SymmetryGroup Cs = SymmetryGroup_Cs();
  ExpectSymmetryGroup(Cs);  

  EXPECT_TRUE(Cs->Non0_Scalar(0, 0));
  EXPECT_TRUE(Cs->Non0_Scalar(1, 1));
  EXPECT_FALSE(Cs->Non0_Scalar(0, 1));
  EXPECT_FALSE(Cs->Non0_Scalar(1, 0));

  EXPECT_TRUE(!Cs->Non0_Z(0, 0));
  EXPECT_TRUE(!Cs->Non0_Z(1, 1));
  EXPECT_TRUE( Cs->Non0_Z(0, 1));
  EXPECT_TRUE( Cs->Non0_Z(1, 0));

  EXPECT_TRUE( Cs->Non0_4(1, 0, 1, 0));

  vector<PrimGTO> gs(2);
  gs[0] = PrimGTO(0, 0, 1, 0.0, 0.0, +0.2);
  gs[1] = PrimGTO(0, 0, 1, 0.0, 0.0, -0.2);
  
  MatrixXi sym_Ii;
  MatrixXi sig_Ii;
  Cs->CalcSymMatrix(gs, sym_Ii, sig_Ii);
  
  MatrixXi ref_sym_Ii(2, 2), ref_sig_Ii(2, 2);
  ref_sym_Ii << 0, 1, 1, 0;
  ref_sig_Ii << 1, 1,-1, -1;
  EXPECT_MATXI_EQ(ref_sym_Ii, sym_Ii);
  EXPECT_MATXI_EQ(ref_sig_Ii, sig_Ii);
  
  EXPECT_EQ(0, Cs->GetIrrep("A'"));
  EXPECT_EQ(1, Cs->GetIrrep("A''"));
  EXPECT_ANY_THROW(Cs->GetIrrep("A''''"));

}
TEST(SymGroup, C2h) {

  SymmetryGroup C2h = SymmetryGroup_C2h();
  ExpectSymmetryGroup(C2h);  

  int Ag = 0;   
  int Bg = 1; 
  int Au = 2; 
  int Bu = 3; 

  EXPECT_TRUE(C2h->prod_table_(Ag, Ag, Ag));
  EXPECT_TRUE(C2h->prod_table_(Ag, Bg, Bg));
  EXPECT_TRUE(C2h->prod_table_(Ag, Au, Au));
  EXPECT_TRUE(C2h->prod_table_(Ag, Bu, Bu));
			                  
  EXPECT_TRUE(C2h->prod_table_(Bg, Ag, Bg));
  EXPECT_TRUE(C2h->prod_table_(Bg, Bg, Ag));
  EXPECT_TRUE(C2h->prod_table_(Bg, Au, Bu));
  EXPECT_TRUE(C2h->prod_table_(Bg, Bu, Au));
			                  
  EXPECT_TRUE(C2h->prod_table_(Au, Ag, Au));
  EXPECT_TRUE(C2h->prod_table_(Au, Bg, Bu));
  EXPECT_TRUE(C2h->prod_table_(Au, Au, Ag));
  EXPECT_TRUE(C2h->prod_table_(Au, Bu, Bg));
			                  
  EXPECT_TRUE(C2h->prod_table_(Bu, Ag, Bu));
  EXPECT_TRUE(C2h->prod_table_(Bu, Bg, Au));
  EXPECT_TRUE(C2h->prod_table_(Bu, Au, Bg));
  EXPECT_TRUE(C2h->prod_table_(Bu, Bu, Ag));
}
TEST(SymGroup, C2v) {
  SymmetryGroup C2v = SymmetryGroup_C2v();
  ExpectSymmetryGroup(C2v);

  int A1 = C2v->GetIrrep("A1");
  int A2 = C2v->GetIrrep("A2");
  int B1 = C2v->GetIrrep("B1");
  int B2 = C2v->GetIrrep("B2");

  EXPECT_TRUE(C2v->prod_table_(A1, A1, A1));
  EXPECT_TRUE(C2v->prod_table_(A1, A2, A2));
  EXPECT_TRUE(C2v->prod_table_(A1, B1, B1));
  EXPECT_TRUE(C2v->prod_table_(A1, B2, B2));

  EXPECT_TRUE(C2v->prod_table_(A2, A1, A2));
  EXPECT_TRUE(C2v->prod_table_(A2, A2, A1));
  EXPECT_TRUE(C2v->prod_table_(A2, B1, B2));
  EXPECT_TRUE(C2v->prod_table_(A2, B2, B1));

  EXPECT_TRUE(C2v->prod_table_(B1, A1, B1));
  EXPECT_TRUE(C2v->prod_table_(B1, A2, B2));
  EXPECT_TRUE(C2v->prod_table_(B1, B1, A1));
  EXPECT_TRUE(C2v->prod_table_(B1, B2, A2));

  EXPECT_TRUE(C2v->prod_table_(B2, A1, B2));
  EXPECT_TRUE(C2v->prod_table_(B2, A2, B1));
  EXPECT_TRUE(C2v->prod_table_(B2, B1, A2));
  EXPECT_TRUE(C2v->prod_table_(B2, B2, A1));

}
TEST(SymGroup, D2h) {
  
  SymmetryGroup D2h = SymmetryGroup_D2h();
  ExpectSymmetryGroup(D2h);
  //  cout << D2h.str() << endl;
  int Ag = 0;    
  int B1g = 1;
  int B2g = 2;  
  int B3g = 3;
  int Au = 4;
  int B1u = 5;   
  int B2u = 6;   
  int B3u = 7;   
  
  EXPECT_TRUE(D2h->prod_table_(Ag, Au, Au));
  EXPECT_TRUE(D2h->prod_table_(B2g, B3g, B1g));
  EXPECT_TRUE(D2h->prod_table_(Au, B2u, B2g));
  EXPECT_TRUE(D2h->prod_table_(B1u, B1g, Au));
  EXPECT_TRUE(D2h->prod_table_(B3u, B2u, B1g));
  
  EXPECT_FALSE(D2h->prod_table_(Au, B2u, B3g));
}
TEST(SymGroup, C4) {
  SymmetryGroup C4 = SymmetryGroup_C4();
  //  cout << C4.str() << endl;
  ExpectSymmetryGroup(C4);

  Irrep A(0);
  Irrep B(1);
  Irrep E(2);

  EXPECT_TRUE(C4->prod_table_(A, A, A));
  EXPECT_TRUE(C4->prod_table_(A, B, B));
  EXPECT_TRUE(C4->prod_table_(A, E, E));

  EXPECT_TRUE(C4->prod_table_(B, A, B));
  EXPECT_TRUE(C4->prod_table_(B, B, A));
  EXPECT_TRUE(C4->prod_table_(B, E, E));

  EXPECT_TRUE(C4->prod_table_(E, A, E));
  EXPECT_TRUE(C4->prod_table_(E, B, E));
  EXPECT_TRUE(C4->prod_table_(E, E, A));
  EXPECT_TRUE(C4->prod_table_(E, E, B));

  EXPECT_TRUE(!C4->prod_table_(E, E, E));
  
}
TEST(SymGroup, Equality) {

  vector<SymmetryGroup> syms;
  syms.push_back(SymmetryGroup_C1());
  syms.push_back(SymmetryGroup_Cs());
  syms.push_back(SymmetryGroup_C2h());
  syms.push_back(SymmetryGroup_C2v());
  syms.push_back(SymmetryGroup_D2h());
  syms.push_back(SymmetryGroup_C4());
  
  typedef vector<SymmetryGroup>::iterator It;
  for(It it = syms.begin(); it != syms.end(); ++it)
    for(It jt = syms.begin(); jt != syms.end(); ++jt) {

      if(it == jt)
	EXPECT_TRUE((*it)->IsSame(*jt));
      else
	EXPECT_FALSE((*it)->IsSame(*jt));

    }
    

}
TEST(SymGroup, SymPosList) {

  vector<Vector3cd> xs(1);
  xs[0] = Vector3cd(0.1, 0.2, 0.3);
  
  SymmetryGroup cs = SymmetryGroup_Cs();

  vector<Vector3cd> ys;
  cs->CalcSymPosList(xs, &ys);

  EXPECT_EQ(2, ys.size());
  EXPECT_C_EQ(0.1, ys[1][0]);
  EXPECT_C_EQ(0.2, ys[1][1]);
  EXPECT_C_EQ(-0.3, ys[1][2]);
}
int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
