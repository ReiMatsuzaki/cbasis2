/*
  Unit test for cints.cpp, 
*/

#include <gtest/gtest.h>
#include <Eigen/Core>

#include "exp_func.hpp"
#include "cip_exp.hpp"

#include "utils.hpp"
#include "cints.hpp"
#include "angmoment.hpp"
#include "gto3d.hpp"
#include "gto3dset.hpp"
#include "molint.hpp"
#include "eigen_plus.hpp"
#include "spec_func.hpp"

using namespace std;
using namespace cbasis;

int lm(int l, int m) {
  return lm_index(l, m);
}

TEST(MOLINT, overlap) {

  dcomplex old = overlap(C(1.2,0.2), 2, 1, 0,
			 1.0,  0.0,   0.0,
			 C(0.7,-0.5), 1, 1, 2,
			 0.0, 0.3, 0.0);
  dcomplex calc = gto_overlap(2, 1, 0,
			   1.0, 0.0, 0.0,
			   dcomplex(1.2, 0.2),
			   1, 1, 2,
			   0.0, 0.3, 0.0,
			      dcomplex(0.7, -0.5));
  EXPECT_C_EQ(old, calc);

  old = overlap(C(1.2,0.2), 2, 1, 0,
		0.0,  0.0,   0.0,
		C(0.7,-0.5), 2, 1, 0,
		0.0, 0.0, 0.0);
  calc = gto_overlap(2, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(1.2, 0.2),
		     2, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

}
TEST(MOLINT, kinetic) {

  dcomplex old = kinetic(C(1.2,0.2), 0, 0, 0,
			 0.0,  0.0,   0.0,
			 C(0.7,-0.5), 0, 0, 0,
			 0.0, 0.0, 0.0);
  dcomplex calc = gto_kinetic(0, 0, 0,
			      0.0, 0.0, 0.0,
			      dcomplex(1.2, 0.2),
			      0, 0, 0,
			      0.0, 0.0, 0.0,
			      dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

  old = kinetic(C(1.2,0.2), 0, 1, 0,
		0.0,  0.0,   0.0,
		C(0.7,-0.5), 0, 1, 0,
		0.0, 0.0, 0.0);
  calc = gto_kinetic(0, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(1.2, 0.2),
		     0, 1, 0,
		     0.0, 0.0, 0.0,
		     dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

  old = kinetic(dcomplex(1.2,0.2), 2, 0, 3,
		0.2,  0.3,   0.4,
		dcomplex(0.7,-0.5), 1, 1, 2,
		0.0, 0.0, 0.0);
  calc = gto_kinetic(2, 0, 3,
		     0.2, 0.3, 0.4,
		     dcomplex(1.2, 0.2),
		     1, 1, 2,
		     0.0, 0.0, 0.0,
		     dcomplex(0.7, -0.5));
  EXPECT_TRUE(abs(calc)>0.0001);  
  EXPECT_C_EQ(old, calc);

}
TEST(MOLINT, IncGamma) {

  /*
    python code :
    from scipy.special import hyp1f1
    def inc_gamma(m, z):
        return 1.0/(2*m+1) * hyp1f1(m+0.5, m+1.5, -z)
    print inc_gamma(0, 50.0-5.50j)
    print inc_gamma(1, 50.0-5.50j)
    print inc_gamma(2, 50.0-5.50j)
   */

  dcomplex* us = new dcomplex[11];
  IncompleteGamma(10, dcomplex(0.1, -0.2), us);
  EXPECT_C_EQ(dcomplex(0.96392503253589557,
		       +0.06262986651323893),		       
	      us[0]);
  EXPECT_C_EQ(dcomplex(0.100986690286, 0.0166268474338),
	      us[4]);

  IncompleteGamma(10, dcomplex(21.0, 15.5), us);
  EXPECT_C_EQ(dcomplex(-3.59002512901*pow(10.0, -6),
		       -0.000190939738972),
	      us[2]);
  
  IncompleteGamma(10, dcomplex(50.0, -5.5), us);

  double step(0.00001);
  EXPECT_C_EQ(dcomplex(0.124767690392,
		       0.00684158939256),
	      us[0]);  
  EXPECT_C_EQ(dcomplex(0.0012253247264,
		       0.00020320161383),
	      us[1]);
  EXPECT_C_EQ(dcomplex(3.5657718077638774*step,
		       1.0018397403429266*step),
	      us[2]);

  IncompleteGamma(10, 0.0, us);
  EXPECT_C_EQ(1.0, us[0]);
  EXPECT_C_EQ(1.0/3.0, us[1]);
  EXPECT_C_EQ(1.0/5.0, us[2]);

  delete us;
  
}
TEST(MOLINT, nuclear_attraction) {

  dcomplex zetaA(1.2, -0.3);
  dcomplex zetaB(0.6, 0.0);  
  dcomplex old = nuclear_attraction(0.0,  0.0,  0.0,
				    1.0,
				    0, 0, 0, zetaA,
				    0.0, 0.0, 0.0,
				    1.0,
				    0, 0, 0, zetaB,
				    0.0, 0.0, 0.0);
  dcomplex calc = gto_nuclear_attraction(0, 0, 0, 
					 0.0, 0.0, 0.0,
					 zetaA,
					 0, 0, 0,
					 0.0, 0.0, 0.0,
					 zetaB,
					 0.0, 0.0, 0.0);  
  double eps(pow(10.0, -8));
  EXPECT_C_NEAR(old, calc, eps*10.0);

  old = nuclear_attraction(0.0,  0.0,  0.0,
			   1.0,
			   1, 0, 0, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 0, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(1, 0, 0, 
				0.0, 0.0, 0.0,
				zetaA,
				0, 0, 0,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_NEAR(old, calc, eps);  

  old = nuclear_attraction(0.3,  0.0,  0.0,
			   1.0,
			   0, 0, 0, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 0, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(0, 0, 0, 
				0.3, 0.0, 0.0,
				zetaA,
				0, 0, 0,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_NEAR(old, calc, eps);  

  old = nuclear_attraction(0.3,  0.0,  0.0,
			   1.0,
			   1, 0, 0, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 0, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(1, 0, 0, 
				0.3, 0.0, 0.0,
				zetaA,
				0, 0, 0,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_NEAR(old, calc, eps);  

  old = nuclear_attraction(0.3,  0.0,  0.0,
			   1.0,
			   1, 2, 0, zetaA,
			   0.0, dcomplex(0.1, 0.01), 0.0,
			   1.0,
			   0, 2, 0, zetaB,
			   0.3, 0.0, 0.0);

  calc = gto_nuclear_attraction(1, 2, 0, 
				0.3, 0.0, 0.0,
				zetaA,
				0, 2, 0,
				0.0, dcomplex(0.1, 0.01), 0.0,
				zetaB,
				0.3, 0.0, 0.0);
  EXPECT_C_NEAR(old, calc, eps);  

  old = nuclear_attraction(0.0,  0.0,  0.0,
			   1.0,
			   0, 0, 1, zetaA,
			   0.0, 0.0, 0.0,
			   1.0,
			   0, 0, 1, zetaB,
			   0.0, 0.0, 0.0);

  calc = gto_nuclear_attraction(0, 0, 1, 
				0.0, 0.0, 0.0,
				zetaA,
				0, 0, 1,
				0.0, 0.0, 0.0,
				zetaB,
				0.0, 0.0, 0.0);
  EXPECT_C_NEAR(old, calc, eps);  

}

TEST(GTOs, offset) {

  GTOs gtos;
  gtos.AddSphericalGTOs(1, 0.0, 0.0, 0.0, 1.2);
  gtos.AddSphericalGTOs(2, 0.0, 0.0, 0.0, 1.2);
  gtos.AddSphericalGTOs(0, 0.0, 0.0, 0.0, 1.2);
  EXPECT_EQ(0, gtos.offset_ish[0]);
  EXPECT_EQ(3, gtos.offset_ish[1]);
  EXPECT_EQ(8, gtos.offset_ish[2]);

}
TEST(GTOs, d_coef) {

  dcomplex zetaP(1.1, 0.3);
  dcomplex wPk(0.666);
  dcomplex xi(0.111);
  dcomplex xj(0.125);
  dcomplex* buf = new dcomplex[1000];
  MultArray<dcomplex, 3> ary = calc_d_coef(2, 3, 4, zetaP, wPk, xi, xj, buf);

  int idx(0);
  for(int ni = 0; ni <= 2; ni++)
    for(int nj = 0; nj <= 3; nj++)
      for(int n = 0; n <=4; n++) {
	dcomplex calc = ary.get(ni, nj, n);
	dcomplex ref0 = coef_d(zetaP, wPk, xi, xj, ni, nj, n);
	EXPECT_C_EQ(ref0, calc)<< ni << nj << n;	  
	idx++;
      }
  delete[] buf;
}
TEST(GTOs, size) {

  GTOs gtos;
  
  for(int L = 2; L < 3; L++)
    for(int n = 0; n < 1; n++ ) {

      dcomplex zeta(0.3);
      dcomplex x(0.0);
      dcomplex y(0.0);
      dcomplex z(0.0);

      gtos.AddSphericalGTOs(L, x, y, z, zeta);
  }

  EXPECT_EQ(6, gtos.size_prim());
  EXPECT_EQ(5, gtos.size_basis());

}
TEST(GTOs, CalcMat) {
  SphericalGTOSet gto_ref;
  GTOs gtos;

  for(int L = 0; L < 3; L++)
    for(int n = 0; n < 1; n++ ) {
      
      dcomplex zeta(0.1+n*0.1, 0.0);
      dcomplex x(0.4*(n+1));
      dcomplex y(-0.1+0.1*(n+1));
      dcomplex z(0.3+0.2*(n+1));

      gtos.AddSphericalGTOs(L, x, y, z, zeta);
      gto_ref.AddBasis(    L, x, y, z, zeta);
      
  }

  gtos.AddAtom(1.0, 0.0, 0.0, 0.0);

  // below code contain serious error.
  dcomplex* smat_ref  = gto_ref.SMat();
  dcomplex* tmat_ref  = gto_ref.TMat();
  dcomplex* zmat_ref  = gto_ref.XyzMat(0, 0, 1);
  dcomplex* vmat_ref = gto_ref.VMat(1.0, 0.0, 0.0, 0.0);
 
  int num = gtos.size_basis();  
  MatrixSet res = gtos.Calc();

  double eps(pow(10.0, -10.0));

  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      //      std::cout << i << j << ", " << std::endl;
      int idx = i * num + j;
      //      std::cout << i << j << std::endl;
      EXPECT_C_NEAR(smat_ref[idx],
		    res.get("s")[idx],
		    eps) 
		<< "(" << i << ", " << j << ")" << std::endl;

      EXPECT_C_NEAR(tmat_ref[idx], res.get("t")[idx], eps) 
	<< "(" << i << ", " << j << ")" << std::endl;
      EXPECT_C_NEAR(zmat_ref[idx],
		    res.get("z")[idx],
		    eps);
      EXPECT_C_NEAR(vmat_ref[idx],
		    res.get("v")[idx], eps)
	<< "(" << i << ", " << j << ")" << std::endl;

    }

  delete[] smat_ref;
  delete[] tmat_ref;
  delete[] vmat_ref;
  delete[] zmat_ref;

}
TEST(GTOs, AddGTO) {

  GTOs gtos1;
  GTOs gtos2;
  for(int L = 0; L < 3; L++)
    for(int n = 0; n < 1; n++ ) {
      
      dcomplex zeta(0.1+n*0.1, 0.0);
      dcomplex x(0.4*(n+1));
      dcomplex y(-0.1+0.1*(n+1));
      dcomplex z(0.3+0.2*(n+1));

      gtos1.AddSphericalGTOs(L, x, y, z, zeta);
      for(int M = -L; M <= L; M++) {
	gtos2.AddOneSphericalGTO(L, M, x, y, z, zeta);
      }
  }

  gtos1.AddAtom(1.0, dcomplex(1.1, 0.1), 0.0, dcomplex(0.2, -0.2));
  gtos2.AddAtom(1.0, dcomplex(1.1, 0.1), 0.0, dcomplex(0.2, -0.2));

  MatrixSet res1 = gtos1.Calc();
  MatrixSet res2 = gtos2.Calc();
 
  int num = gtos1.size_basis();  
  int num2 = gtos2.size_basis();  
  EXPECT_EQ(num, num2);
  
  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      int idx = i * num + j;
      EXPECT_C_EQ(res1.get("s")[idx], res2.get("s")[idx]) << i << j;
      EXPECT_C_EQ(res1.get("v")[idx], res2.get("v")[idx]) << i << j;
      EXPECT_C_EQ(res1.get("t")[idx], res2.get("t")[idx]) << i << j;
      EXPECT_C_EQ(res1.get("z")[idx], res2.get("z")[idx]) << i << j;
    }

}
TEST(GTOs, H_atom) {

  GTOs gtos;
  for(int n = -5; n < 5; n++ ) {
    dcomplex zeta(pow(2.0, n), 0.0);
    gtos.AddSphericalGTOs(0, 0.0, 0.0, 0.0, zeta);
  }
  gtos.AddAtom(1.0, 0.0, 0.0, 0.0);
  MatrixSet res = gtos.Calc();
  int num = gtos.size_basis();
  
  Eigen::MatrixXcd H(num, num);
  Eigen::MatrixXcd S(num, num);
  dcomplex* s_ptr = res.get("s");
  dcomplex* t_ptr = res.get("t");
  dcomplex* v_ptr = res.get("v");
  for(int i = 0; i < num; i++) {
    for(int j = 0; j < num; j++) {
      int idx(i*num+j);
      H(i, j) = t_ptr[idx] + v_ptr[idx];
      S(i, j) = s_ptr[idx];
    }
  }
  Eigen::MatrixXcd CS(5, 5);
  Eigen::VectorXcd eig(5, 1);
  generalizedComplexEigenSolve(H, S, &CS, &eig);
}
TEST(GTOs, mat) {

  Eigen::MatrixXcd m(5, 5);
  Eigen::MatrixXcd s(5, 5);
  for(int i = 0; i < 5; i++) {
    for(int j = 0; j < 5; j++) {
      m(i, j) = 1.0 * (i + j + 1);
      s(i, j) = i==j ? 1.0 : 0.0;
    }
  }
  Eigen::MatrixXcd CS(5, 5);
  Eigen::VectorXcd res(5, 1);
  generalizedComplexEigenSolve(m, s, &CS, &res);
  //  std::cout << CS << std::endl;
  std::cout << res<< std::endl;

}
TEST(GTOs, at_r_ylm) {

  // -- common --
  dcomplex zeta1(0.1, 0.01);
  dcomplex zeta2(0.1, 0.02);

  // -- GTO --
  CGTO g1(1.0, 1, zeta1); dcomplex cnorm1 = CIP(g1, g1);
  g1 = CGTO(1.0/sqrt(cnorm1), 1, zeta1);
  CGTO g2(1.0, 1, zeta2); dcomplex cnorm2 = CIP(g2, g2);
  g2 = CGTO(1.0/sqrt(cnorm2), 1, zeta2);

  // -- GTO set --
  GTOs gtos;  
  gtos.AddSphericalGTOs(0, 0.0, 0.0, 0.0, zeta1);
  gtos.AddSphericalGTOs(0, 0.0, 0.0, 0.0, zeta2);

  dcomplex* cs = new dcomplex[2];
  cs[0] = 1.1; cs[1] = 1.3;
  int r_num = 2;
  dcomplex* rs = new dcomplex[r_num];
  for(int i = 0; i < r_num; i++)
    rs[i] = 1.1 + i * 0.1;

  dcomplex* vs = new dcomplex[r_num];
  gtos.AtR_Ylm(0, 0, rs, r_num, cs, vs);

  EXPECT_C_EQ(g1.c(), gtos.coef_ylm_ish[0]);

  for(int i = 0; i < r_num; i++) {
    std::cout << std::setprecision(15) << rs[i] << ", " << vs[i] << std::endl;
    EXPECT_C_EQ(g1.at(rs[i]) * cs[0] + g2.at(rs[i]) * cs[1],
		vs[i]);
  }

  delete[] rs;
  delete[] vs;
  delete[] cs;
  
}
TEST(GTOs, location) {

  GTOs gtos1;  
  GTOs gtos2;  
  dcomplex x(1.1, 0.1);
  dcomplex y(0.2, -0.03);
  dcomplex z(-0.3);
  
  dcomplex zeta(1.1, 0.2);
  for(int L = 0; L < 3; L++) {
    gtos1.AddSphericalGTOs(L, 0.0, 0.0, 0.0, zeta);
    gtos2.AddSphericalGTOs(L, x, y, z, zeta);
    for(int M = -L; M <= L; M++) {
      gtos1.AddOneSphericalGTO(L, M, 0.0, 0.0, 0.0, zeta);
      gtos2.AddOneSphericalGTO(L, M, x, y, z, zeta);
    }
  }
  gtos1.AddAtom(1.1, 0.0, 0.0, 0.0);
  gtos2.AddAtom(1.1, x, y, z);

  MatrixSet s1 = gtos1.Calc();
  MatrixSet s2 = gtos2.Calc();

  int num = gtos1.size_basis();

  double eps(pow(10.0, -10.0));
  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      int idx = i*num+j;
      EXPECT_C_NEAR(s1.get("s")[idx], s2.get("s")[idx], eps);
      EXPECT_C_NEAR(s1.get("v")[idx], s2.get("v")[idx], eps);
      EXPECT_C_NEAR(s1.get("t")[idx], s2.get("t")[idx], eps);
    }
}
TEST(GTOs, CalcZMatOther) {

  GTOs gtos1;
  GTOs gtos2;

  dcomplex x(1.0);
  dcomplex y(0.3, 0.1);
  dcomplex z(2.0);
  dcomplex zeta(1.0, 0.1);
  gtos1.AddSphericalGTOs(0, x, y, z, zeta);
  gtos1.AddSphericalGTOs(2, x, y, z, zeta);
  
  gtos2.AddSphericalGTOs(2, x, y, z, zeta);

  MatrixSet m1 = gtos1.Calc();
  MatrixSet m2 = gtos1.CalcZMatOther(gtos2);

  double eps(pow(10.0, -10.0));

  for(int i = 0; i < 5; i++)
    for(int j = 0; j < 5; j++)
      EXPECT_C_NEAR(m1.get("z")[1+i + (j+1)*6],
		    m2.get("z")[1+i +     j*6],
		    eps) << i << j;
}

TEST(CINTS, First) {
  EXPECT_TRUE(true);
}
TEST(CINTS, lgamma) {
  C val(-0.008197780565406, -0.05732294041672);
  C calc(lgamma(C(1.0, 0.1)));
  double eps(0.00000001);
  EXPECT_C_NEAR(val, calc, eps);
}
TEST(CINTS,  overlap) {
  double eps(0.000000001);
  C s = overlap(C(1.1,0.2), 0, 0, 0,
		0.1,  0.2,   0.3,
		C(1.3,0.5), 0, 0, 0,
		0.1, -0.22, -0.3);
  EXPECT_C_NEAR(s, C(0.88910293518825, -0.5004029971544), eps);
  //  EXPECT_C_EQ(s, C(0.88910293518825, -0.5004029971544));

  C s1 = overlap(C(1.1,0.2), 2, 1, 0,
		 0.1,  C(0.2, 0.03),   0.3,
		 C(1.3,0.02), 0, 1, 3,
		 0.1, -0.22, C(-0.3,-0.7));
  EXPECT_C_NEAR(s1, C(0.012403718856620667, 0.003134556944411788), eps);
  
  C s2 = overlap(C(0.2,0.3), 2, 1, 1,
		 0.0, 0.0, 0.0,
		 C(0.2,-0.7), 2,1,2,
		 0.0, C(0.2,-0.1), 0.3);
  EXPECT_C_NEAR(s2, C(-14.3092, -3.67256), 0.00005);
  
}
TEST(CINTS,  kinetic) {
  double eps(0.000001);
  C s = kinetic(C(1.1,0.2), 0, 0, 0,
		0.1,  0.2,   0.3,
		C(1.3,0.5), 0, 0, 0,
		0.1, -0.22, -0.3);
  C se(1.422908003563741695, -0.476445156336799588);
  EXPECT_C_NEAR(s, se, eps);


/*
  C s1 = kinetic(C(1.1,0.2), 2, 1, 0,
		 0.1,  C(0.2, 0.03),   0.3,
		 C(1.3,0.02), 0, 1, 3,
		 0.1, -0.22, C(-0.3,-0.7));
  //  EXPECT_NEAR_COMPLEX(s1, C(0.0143520521644942813,0.0108902064709103336), eps);
  
  C s2 = kinetic( C(0.1+0.3, 0.3), 1, 0, 1,
		  C(0.3, 0.1), 0.7, 0.0,
		  C(0.3, -0.1), 1, 1, 0,
		  C(0.3, 0.1), C(0.0,-0.0), 0.3);
  //  EXPECT_NEAR_COMPLEX(s2, C(0.404777, -0.50371), 0.000005);
  
  C s3 = kinetic( C(0.1+0.3, 0.3), 2, 1, 1,
		  C(0.3, 0.1), 0.7, 0.0,
		  C(0.3, -0.1), 2, 2, 1,
		  C(0.3, 0.1), C(0.0,-0.0), 0.3);
  //  EXPECT_NEAR_COMPLEX(s3, C(0.404777, -0.50371), 0.000005);
  
  C s3inv = kinetic( C(0.3, -0.1), 2, 2, 1,
		     C(0.3, 0.1), C(0.0,-0.0), 0.3,
		     C(0.1+0.3, 0.3), 2, 1, 1,
		     C(0.3, 0.1), 0.7, 0.0);

  //  EXPECT_NEAR_COMPLEX(s3inv, C(0.404777, -0.50371), eps);
*/
}
TEST(CINTS, product_center_1D) {


  /*
  C a1(1.1,0.2), a2(1.3,0.5), a3(1.1,0.2), a4(1.3,0.02);
  
  C rys   = product_center_1D_rys(1.1, 1.2,1.3,1.4);
  C cints = product_center_1D(1.1, 1.2,1.3,1.4);
  EXPECT_NEAR_COMPLEX(rys, cints, eps);

  
  rys   = product_center_1D_Rys(a1, 0.1, a2, 0.3);
  cints = product_center_1D(a1, 0.1, a2, 0.3);
  
  EXPECT_NEAR_COMPLEX(rys, cints, eps);
  */
}
TEST(CINTS,  nuclear_attraction) {
  double eps(0.000000001);
  C sc = nuclear_attraction(0.1,  C(0.2, 0.03), 0.3, 1.0, 
			    2, 1, 0, C(1.1,0.2),
			    0.1, -0.22, C(-0.3,-0.7), 1.0,
			    0, 1, 3, C(1.3,0.02),
			    -0.1, C(-0.1,0.1), 0.3);
  C se(-0.0111712963403, -0.0039848461450);
  EXPECT_C_NEAR(sc, se, eps);
}
TEST(CartGTO, SetScalarProd) {

  CartGTO<C, d3> g(1.0, i3(1, 1, 2), d3(0.0, 0.0, 0.0), C(1.0, 0.2));
  g.SetScalarProd(1.2);
  EXPECT_C_EQ(1.2, g.c());

}
TEST(CartGTO, SetNormalize) {
  
  double eps(0.000000001);
  CartGTO<C, d3> g(1.0, i3(1, 1, 2), d3(0.3, 0.1, 0.2), C(1.0, 0.2));
  g.SetNormalize();
  EXPECT_C_NEAR(1.0, CIP(g, g), eps);

}
TEST(CartGTO, Overlap) {
  double eps(0.00005);
  double c1(1.1);
  double c2(1.2);
  CartGTO<C, c3> g1(c1, i3(2, 1, 1), c3(0.0, 0.0, 0.0),        C(0.2, 0.3));
  CartGTO<C, c3> g2(c2, i3(2, 1, 2), c3(0.0, C(0.2,-0.1), 0.3), C(0.2, -0.7));
  C s = CIP(g1, g2);
  
  EXPECT_C_NEAR(s, c1*c2*C(-14.3092, -3.67256), eps);  
}
TEST(CargGTO, Kinetic) {

  double eps(0.000000001);
  double c1(1.1);
  double c2(1.2);
  CartGTO<C, c3> g1(c1,
		   i3(0, 0, 0),
		   c3(0.1, 0.2, 0.3),
		   C(1.1, 0.2));
  CartGTO<C, c3> g2(c2, 
		   i3(0, 0, 0),
		   c3(0.1, -0.22, -0.3),
		   C(1.3, 0.5));

  C s = CIP(g1, OpKE<C, c3>(), g2);
  C se(1.422908003563741695, -0.476445156336799588);
  EXPECT_C_NEAR(s, c1*c2*se, eps);    
  

}
TEST(CartGTO, NuclearAttraction) {
  double eps(0.000000001);
  double c1(0.33);
  double c2(0.11);
  double q1(1.1);
  CartGTO<C, c3> g1(c1,
		   i3(2, 1, 0),
		   c3(0.1,  C(0.2, 0.03), 0.3), 
		   C(1.1,0.2));
  CartGTO<C, c3> g2(c2, 
		   i3(0, 1, 3),
		   c3(0.1, -0.22, C(-0.3,-0.7)),
		   C(1.3,0.02));
  OpNA<C, c3> op_na(q1, c3(-0.1, C(-0.1,0.1), 0.3));
    
  C s = CIP(g1, op_na, g2);
  C se(-0.0111712963403, -0.0039848461450);
  
  EXPECT_C_NEAR(s, q1*c1*c2*se, eps);    
}

TEST(Angmoment, cg_coef) {
  int j1 = 2;
  int m1 = 1;
  EXPECT_DOUBLE_EQ(-1.0/sqrt(2*j1+1),
		   cg_coef(j1,j1,m1,-m1, 0,0));
  EXPECT_NEAR(pow(fact(2*j1),2)/(6 * sqrt(fact(8)*1.0)),
	      cg_coef(j1,j1,m1,-m1,2*j1,0),
	      0.0000000001);
  std::cout << cg_coef(1,1,1,0,0,0) << std::endl;
}
TEST(Angmoment, lm_pair) {

  EXPECT_EQ(0, lm_index(0, 0));

  EXPECT_EQ(1, lm_index(1, -1));
  EXPECT_EQ(2, lm_index(1,  0));
  EXPECT_EQ(3, lm_index(1, +1));

  EXPECT_EQ(4, lm_index(2, -2));
  EXPECT_EQ(5, lm_index(2, -1));
  EXPECT_EQ(6, lm_index(2,  0));
  EXPECT_EQ(7, lm_index(2, +1));
  EXPECT_EQ(8, lm_index(2, +2));

  EXPECT_EQ(9,  lm_index(3, -3));
  EXPECT_EQ(10, lm_index(3, -2));
  EXPECT_EQ(11, lm_index(3, -1));
  EXPECT_EQ(12, lm_index(3,  0));
  EXPECT_EQ(13, lm_index(3, +1));
  EXPECT_EQ(14, lm_index(3, +2));
  EXPECT_EQ(15, lm_index(3, +3));
}
TEST(Angmoment, num_lm_pair) {

  EXPECT_EQ(1, num_lm_pair(0));
  EXPECT_EQ(1+3, num_lm_pair(1));
  EXPECT_EQ(1+3+5, num_lm_pair(2));
  EXPECT_EQ(1+3+5+7, num_lm_pair(3));

}
TEST(Angmoment, mod_spherical_bessel) {

  C xs[3];
  xs[0] = C(1.1, 0.3);
  xs[1] = C(0.1, 0.002);
  xs[2] = C(10.0, 2.0);
  
  for(int i = 0; i < 3; i++) {
    C x(xs[i]);
    C* vs = new C[11];
    ModSphericalBessel(x, 10, vs);
  
    EXPECT_C_EQ(sinh(x)/x,                 vs[0]) << x;
    EXPECT_C_EQ((x*cosh(x)-sinh(x))/(x*x), vs[1]) << x;
    EXPECT_C_EQ(((x*x+3.0)*sinh(x)-3.0*x*cosh(x))/(x*x*x), vs[2]) << x;
    EXPECT_C_EQ(((x*x*x+15.0*x)*cosh(x)-(6.0*x*x+15.0)*sinh(x))/pow(x, 4.0), vs[3])
      << x;
    
    delete[] vs;
  }

  xs[0] = C(0.01, 0.000002);
  xs[1] = C(0.00005, 0.002);
  xs[2] = C(0.0001, 0.002);
  for(int i = 0; i < 3; i++) {
    C* vs = new C[11];
    ModSphericalBessel(xs[i], 10, vs);
    for(int n = 2; n <= 8; n++) {
      EXPECT_C_EQ(1.0, (vs[n-1]-vs[n+1])*xs[i]/ ((2*n+1.0)*vs[n])) << xs[i];    
    }
    delete[] vs;
  }
    
  C x(0.0, 0.0);
  C* vs = new C[11];
  ModSphericalBessel(x, 10, vs);
  EXPECT_C_EQ(1.0, vs[0]);
  for(int i = 1; i < 10; i++)
    EXPECT_C_EQ(0.0, vs[i]);
  delete[] vs;

}
TEST(Angmoment, associated_legendre) {

  C x( 0.6, 0.0);
  int lmax(10);
  C* vs = new C[num_lm_pair(lmax)];
  AssociatedLegendre(x, lmax, vs);  
  EXPECT_C_EQ(1.0, vs[lm(0, +0)]);

  EXPECT_C_EQ(-sqrt(1.0-x*x), vs[lm(1, +1)]);
  EXPECT_C_EQ(x,              vs[lm(1, +0)]);  

  EXPECT_C_EQ(0.5*(3.0*x*x-1.0),  vs[lm(2, 0)]);
  EXPECT_C_EQ(-3.0*x*sqrt(1.0-x*x), vs[lm(2, 1)]);
  EXPECT_C_EQ(3.0*(1.0-x*x),      vs[lm(2,+2)]);

  for(int L = 0; L < 4; L++) 
    for(int M = 1; M <= L; M++) 
      EXPECT_C_EQ(vs[lm(L, -M)], pow(-1.0, M)*fact(L-M)*1.0/fact(L+M)*vs[lm(L, M)]);
  delete[] vs;
}
TEST(Angmoment, real_spherical_harm) {

  C t = C(1.1, 0.2);
  C p = C(1.1, 0.1);

  int max_l(10);
  C* vs = new C[num_lm_pair(max_l)];
  RealSphericalHarmonics(t, p, max_l, vs);
  
  C a = 1.0/sqrt(4.0*M_PI);
  EXPECT_C_EQ(a,                         vs[lm(0, 0)]);
  EXPECT_C_EQ(a*sqrt(3.0)*cos(t),        vs[lm(1, 0)]);
  EXPECT_C_EQ(a*sqrt(3.0)*sin(t)*cos(p), vs[lm(1,+1)]);
  EXPECT_C_EQ(a*sqrt(3.0)*sin(t)*sin(p), vs[lm(1,-1)]);
  
  EXPECT_C_EQ(a*sqrt(5.0/4.0)*(3.0*cos(t)*cos(t)-1.0),   vs[lm(2, 0)]);
  EXPECT_C_EQ(a*sqrt(15.0)*sin(t)*cos(t)*cos(p),         vs[lm(2,+1)]);
  EXPECT_C_EQ(a*sqrt(15.0)*sin(t)*cos(t)*sin(p),         vs[lm(2,-1)]);
  EXPECT_C_EQ(a*sqrt(15.0/4.0)*sin(t)*sin(t)*cos(2.0*p), vs[lm(2,+2)]);
  EXPECT_C_EQ(a*sqrt(15.0/4.0)*sin(t)*sin(t)*sin(2.0*p), vs[lm(2,-2)]);

  delete[] vs;
}
TEST(Angmoment, gto_00_r) {

  double pi(M_PI);  
  int num_r(5);
  C* rs = new C[num_r];
  C* vs = new C[num_r];
  
  for(int i = 0; i < num_r; i++) {
    rs[i] = 0.1 + 0.2*i;
  }
  C zeta(1.0, 0.0);
  int L(0);
  int M(0);
  C* work = new C[num_lm_pair(L) + L + 1];
  gto_00_r(0.0, 0.0, 0.0, L, M, zeta, rs, num_r, work, vs);
  for(int i = 0; i < num_r; i++) {
    C r(rs[i]);
    C calc(sqrt(4.0*pi)*exp(-zeta*r*r));    
    EXPECT_C_EQ(vs[i], calc) << i;    
  }
  delete[] vs;
  delete[] rs;
  delete[] work;
}

TEST(SphericalGTO, OrthNormality) {
  double eps(0.000000001);
  c3 xyz(0.1, 0.2, 0.3);
  C zeta(1.1, -0.3);
  
  int num = 9;
  SphericalGTO<C, C>** gs = new SphericalGTO<C, C>*[num];
  
  // LinFunc<CartGTO<C, c3> >* gs;
  //gs = new LinFunc<CartGTO<C, c3> >[num];
  gs[0] = new SphericalGTO<C, C>(0, 0, xyz, zeta);
  gs[1] = new SphericalGTO<C, C>(1, 0, xyz, zeta);
  gs[2] = new SphericalGTO<C, C>(1, 1, xyz, zeta);
  gs[3] = new SphericalGTO<C, C>(1,-1, xyz, zeta);
  gs[4] = new SphericalGTO<C, C>(2, 0, xyz, zeta);
  gs[5] = new SphericalGTO<C, C>(2, 1, xyz, zeta);
  gs[6] = new SphericalGTO<C, C>(2,-1, xyz, zeta);
  gs[7] = new SphericalGTO<C, C>(2, 2, xyz, zeta);
  gs[8] = new SphericalGTO<C, C>(2,-2, xyz, zeta);

  
  for(int i = 0; i < num; i++) {
    EXPECT_C_NEAR(C(1), CIP(*gs[i], *gs[i]), eps) << i;
    for(int j = 0; j < i; j++) {
      EXPECT_C_NEAR(C(0), CIP(*gs[i], *gs[j]), eps) << i << j;
      EXPECT_C_NEAR(C(0), CIP(*gs[j], *gs[i]), eps) << i << j;
    }
  }
  for(int i = 0; i < num; i++)
    delete gs[i];
  delete [] gs;
}
TEST(SphericalGTO, at) {

  c3 xyz(0.1, 0.2, 0.3); C zeta(1.1, -0.3);

  SphericalGTO<C, C> g(1, 0, xyz, zeta);
  EXPECT_C_EQ(0.0, g.at(c3(0.33, 0.44, 0.3)));

  SphericalGTO<C, C> g11(1, 1, xyz, zeta);
  EXPECT_C_EQ(0.0, g11.at(c3(0.1, 0.44, 0.44)));

  SphericalGTO<C, C> g1m1(1, -1, xyz, zeta);
  EXPECT_C_EQ(0.0, g1m1.at(c3(0.11, 0.2, 0.44)));
  
}
TEST(SphericalGTO, Exception) {
  
  LinFunc<CartGTO<C, c3> > func;
  try {
    SetSphericalGTO<C, C>(1, 2, c3(1.1, 1.2, 1.3), 1.2, &func);
  } catch(const std::exception& e) {
    std::cerr << e.what() << std::endl;    
  }

}

int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
