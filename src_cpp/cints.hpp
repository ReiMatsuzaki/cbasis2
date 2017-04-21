/**********************************************************************
cints.hpp
Simple math function using complex number.
This source was converted from cints.c from PyQuntum.
The comment of original source is followed :
 **********************************************************************/  

/*************************************************************************
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **************************************************************************/
/* My routines */

#include <complex>

using namespace std;
typedef complex<double> C;

// z -> Log_e Gamma[z]
 C lgamma(C);

 C fB(int i, int l1, int l2, C px, C ax, C bx, 
	  int r, C g);
 C Bfunc(int i, int r, C g);
 C contr_coulomb(int ia, C *aexps, C *acoefs, C *anorms,
			    C xa, C ya, C za, int la, int ma, int na, 
			    int ib, C *bexps, C *bcoefs, C *bnorms,
			    C xb, C yb, C zb, int lb, int mb, int nb, 
			    int ic, C *cexps, C *ccoefs, C *cnorms,
			    C xc, C yc, C zc, int lc, int mc, int nc, 
			    int id, C *dexps, C *dcoefs, C *dnorms,
			    C xd, C yd, C zd, int ld, int md, int nd);

// (ab|cd) = Int a(1)b(1) v(1,2) c(2)d(2)
// This result can be proof numerically by check_ee_property();
 C coulomb_repulsion(C xa, C ya, C za, C norma,
				int la, int ma, int na, C alphaa,
				C xb, C yb, C zb, C normb,
				int lb, int mb, int nb, C alphab,
				C xc, C yc, C zc, C normc,
				int lc, int mc, int nc, C alphac,
				C xd, C yd, C zd, C normd,
				int ld, int md, int nd, C alphad);

 C *B_array(int l1, int l2, int l3, int l4, C p, C a,
		C b, C q, C c, C d,
		C g1, C g2, C delta);

 C B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
		     int l3, int l4, C Px, C Ax, C Bx,
		     C Qx, C Cx, C Dx, C gamma1,
		     C gamma2, C delta);
 C kinetic(C alpha1, int l1, int m1, int n1,
		      C xa, C ya, C za,
		      C alpha2, int l2, int m2, int n2,
		      C xb, C yb, C zb);
 C overlap(C alpha1, int l1, int m1, int n1,
		      C xa, C ya, C za,
		      C alpha2, int l2, int m2, int n2,
		      C xb, C yb, C zb);
 C overlap_1D(int l1, int l2, C PAx,
			 C PBx, C gamma);
 C nuclear_attraction(C x1, C y1, C z1, C norm1,
		      int l1, int m1, int n1, C alpha1,
		      C x2, C y2, C z2, C norm2,
		      int l2, int m2, int n2, C alpha2,
		      C x3, C y3, C z3);
 C nuclear_attraction0(C x1, C y1, C z1, 
		      int l1, int m1, int n1, C alpha1,
		      C x2, C y2, C z2,
		      int l2, int m2, int n2, C alpha2,
		      C x3, C y3, C z3);
 C A_term(int i, int r, int u, int l1, int l2,
		     C PAx, C PBx, C CPx, C gamma);
 C *A_array(int l1, int l2, C PA, C PB,
		       C CP, C g);

 int fact(int n);
 int fact2(int n);
 C dist2(C x1, C y1, C z1, 
		    C x2, C y2, C z2);
 C dist(C x1, C y1, C z1, 
		   C x2, C y2, C z2);
 C binomial_prefactor(int s, int ia, int ib, C xpa, C xpb);
 int binomial(int a, int b);

 C Fgamma(C m, C x);
 C gamm_inc(C a, C x);

 int ijkl2intindex(int i, int j, int k, int l);

 int fact_ratio2(int a, int b);

 C product_center_1D(C alphaa, C xa, 
			 C alphab, C xb);

 C three_center_1D(C xi, int ai, C alphai,
			      C xj, int aj, C alphaj,
			      C xk, int ak, C alphak);

/* Routines from Numerical Recipes */
 void gser(C *gamser, C a, C x, C *gln);
 void gcf(C *gammcf, C a, C x, C *gln);


