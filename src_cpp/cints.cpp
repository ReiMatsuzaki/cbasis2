/**********************************************************************
cints.cpp
Simple math function using complex number.
This source was converted from cints.c from PyQuntum.
The comment of original source is followed :
 **********************************************************************/  

/**********************************************************************
 * cints.c  C implementation of simple math functions in pyutil.
 *
 * The equations herein are based upon
 * 'Gaussian Expansion Methods for Molecular Orbitals.' H. Taketa,
 * S. Huzinaga, and K. O-ohata. H. Phys. Soc. Japan, 21, 2313, 1966.
 * [THO paper].
 *
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
 **********************************************************************/

#include "cints.hpp"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define SMALL 0.00000001

// ---- gto_at ----
 

// safaty pow
C spow(C z, int n) {
  if(n == 0) 
    return 1.0;
  else
    return std::pow(z, n);
}

C lgamma(C z) {
    C c[7];
    C x,y ,tmp, ser, v;
    int i;

    if (real(z)<=0) return 0;

    c[0]=2.5066282746310005;
    c[1]=76.18009172947146;
    c[2]=-86.50532032941677;
    c[3]=24.01409824083091;
    c[4]=-1.231739572450155;
    c[5]=0.1208650973866179e-2;
    c[6]=-0.5395239384953e-5;

    x   = z;
    y   = x;
    tmp = x+5.5;
    tmp = (x+0.5)*log(tmp)-tmp;
    ser = 1.000000000190015;
    for (i=1; i<7; i++) {
        y   += 1.0;
        ser += c[i]/y;
        }
    v = tmp+log(c[0]*ser/x);
    return v;
    }

C fB(int i, int l1, int l2, C px, C ax, C bx, 
		 int r, C g){
  return binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,r,g);
}

C Bfunc(int i, int r, C g){
  return (fact_ratio2(i,r)*1.0)*pow(4.0*g,r-i);
}

C contr_coulomb(int lena, C *aexps, C *acoefs,
			    C *anorms, C xa, C ya, C za,
			    int la, int ma, int na, 
			    int lenb, C *bexps, C *bcoefs,
			    C *bnorms, C xb, C yb, C zb,
			    int lb, int mb, int nb, 
			    int lenc, C *cexps, C *ccoefs,
			    C *cnorms, C xc, C yc, C zc,
			    int lc, int mc, int nc, 
			    int lend, C *dexps, C *dcoefs,
			    C *dnorms, C xd, C yd, C zd,
			    int ld, int md, int nd){

  int i,j,k,l;
  C Jij = 0.,incr=0.;

  for (i=0; i<lena; i++)
    for (j=0; j<lenb; j++)
      for (k=0; k<lenc; k++)
	for (l=0; l<lend; l++){
	  incr = coulomb_repulsion(xa,ya,za,anorms[i],la,ma,na,aexps[i],
			      xb,yb,zb,bnorms[j],lb,mb,nb,bexps[j],
			      xc,yc,zc,cnorms[k],lc,mc,nc,cexps[k],
			      xd,yd,zd,dnorms[l],ld,md,nd,dexps[l]);
	  
	  Jij += acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*incr;
	}
  return Jij;
}

C coulomb_repulsion(C xa, C ya, C za, C norma,
				int la, int ma, int na, C alphaa,
				C xb, C yb, C zb, C normb,
				int lb, int mb, int nb, C alphab,
				C xc, C yc, C zc, C normc,
				int lc, int mc, int nc, C alphac,
				C xd, C yd, C zd, C normd,
				int ld, int md, int nd, C alphad){

  C rab2, rcd2,rpq2,xp,yp,zp,xq,yq,zq,gamma1,gamma2,delta,sum;
  C *Bx, *By, *Bz;
  int II,J,K;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  rcd2 = dist2(xc,yc,zc,xd,yd,zd);
  xp = product_center_1D(alphaa,xa,alphab,xb);
  yp = product_center_1D(alphaa,ya,alphab,yb);
  zp = product_center_1D(alphaa,za,alphab,zb);
  xq = product_center_1D(alphac,xc,alphad,xd);
  yq = product_center_1D(alphac,yc,alphad,yd);
  zq = product_center_1D(alphac,zc,alphad,zd);
  rpq2 = dist2(xp,yp,zp,xq,yq,zq);
  gamma1 = alphaa+alphab;
  gamma2 = alphac+alphad;
  delta = (1./gamma1+1./gamma2)/4.;

  Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta);
  By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta);
  Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta);

  sum = 0.;
  for (II=0; II<la+lb+lc+ld+1;II++)
    for (J=0; J<ma+mb+mc+md+1;J++)
      for (K=0; K<na+nb+nc+nd+1;K++)
	sum += Bx[II]*By[J]*Bz[K]*Fgamma(II+J+K,0.25*rpq2/delta);

  free(Bx);
  free(By);
  free(Bz);  
  
  return 2.*pow(M_PI,2.5)/(gamma1*gamma2*sqrt(gamma1+gamma2))
    *exp(-alphaa*alphab*rab2/gamma1) 
    *exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd;
}

C *B_array(int l1, int l2, int l3, int l4, C p, C a,
		C b, C q, C c, C d,
		C g1, C g2, C delta){
  int Imax,i1,i2,r1,r2,u,II,i;
  C *B;
  Imax = l1+l2+l3+l4+1;
  B = (C *)malloc(Imax*sizeof(C));
  for (i=0; i<Imax; i++) B[i] = 0.;

  for (i1=0; i1<l1+l2+1; i1++)
    for (i2=0; i2<l3+l4+1; i2++)
      for (r1=0; r1<i1/2+1; r1++)
	for (r2=0; r2<i2/2+1; r2++)
	  for (u=0; u<(i1+i2)/2-r1-r2+1; u++){
	    II = i1+i2-2*(r1+r2)-u;
	    B[II] = B[II] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
				 p,a,b,q,c,d,g1,g2,delta);
	  }

  return B;
}

C B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
	      int l3, int l4, C Px, C Ax, C Bx,
	      C Qx, C Cx, C Dx, C gamma1,
	      C gamma2, C delta){
  /* THO eq. 2.22 */
  return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)
    *pow(-1,i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)
    *pow(-1,u)*(1.0*fact_ratio2(i1+i2-2*(r1+r2),u))
    *spow(Qx-Px,i1+i2-2*(r1+r2)-2*u)
    /spow(delta,i1+i2-2*(r1+r2)-u);
}


C kinetic(C alpha1, int l1, int m1, int n1,
	       C xa, C ya, C za,
	       C alpha2, int l2, int m2, int n2,
	       C xb, C yb, C zb){

  C term0,term1,term2;
  term0 = alpha2*(2.0*(l2+m2+n2)+3)*
    overlap(alpha1,l1,m1,n1,xa,ya,za,
		   alpha2,l2,m2,n2,xb,yb,zb);
  term1 = -2.0*pow(alpha2,2)*
    (overlap(alpha1,l1,m1,n1,xa,ya,za,
		    alpha2,l2+2,m2,n2,xb,yb,zb)
     + overlap(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2+2,n2,xb,yb,zb)
     + overlap(alpha1,l1,m1,n1,xa,ya,za,
		      alpha2,l2,m2,n2+2,xb,yb,zb));
  term2 = -0.5*(l2*(l2-1)*1.0*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2-2,m2,n2,xb,yb,zb) +
		m2*(m2-1)*(1.0)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2-2,n2,xb,yb,zb) +
		n2*(n2-1)*(1.0)*overlap(alpha1,l1,m1,n1,xa,ya,za,
					 alpha2,l2,m2,n2-2,xb,yb,zb));
  return term0+term1+term2;
}

C overlap(C alpha1, int l1, int m1, int n1,
	  C xa, C ya, C za,
	  C alpha2, int l2, int m2, int n2,
	  C xb, C yb, C zb){
  /*Taken from THO eq. 2.12*/
  C rab2,gamma,xp,yp,zp,pre,wx,wy,wz;

  rab2 = dist2(xa,ya,za,xb,yb,zb);
  gamma = alpha1+alpha2;
  xp = product_center_1D(alpha1,xa,alpha2,xb);
  yp = product_center_1D(alpha1,ya,alpha2,yb);
  zp = product_center_1D(alpha1,za,alpha2,zb);

  pre = pow(M_PI/gamma,1.5)*exp(-alpha1*alpha2*rab2/gamma);

  wx = overlap_1D(l1,l2,xp-xa,xp-xb,gamma);
  wy = overlap_1D(m1,m2,yp-ya,yp-yb,gamma);
  wz = overlap_1D(n1,n2,zp-za,zp-zb,gamma);
  return pre*wx*wy*wz;
}

C overlap_1D(int l1, int l2, C PAx,
			 C PBx, C gamma){
  /*Taken from THO eq. 2.12*/
  int i;
  C sum;
  sum = 0.;
  for (i=0; i<(1+floor(0.5*(l1+l2))); i++)
    sum += 1.0 *binomial_prefactor(2*i,l1,l2,PAx,PBx)* 
      (1.0*fact2(2*i-1)/pow(2.0*gamma,i));
  return sum;
}
    
C nuclear_attraction(C x1, C y1, C z1, C norm1,
		     int l1, int m1, int n1, C alpha1,
		     C x2, C y2, C z2, C norm2,
		     int l2, int m2, int n2, C alpha2,
		     C x3, C y3, C z3){
  int II,J,K;
  C gamma,xp,yp,zp,sum,rab2,rcp2;
  C *Ax,*Ay,*Az;

  gamma = alpha1+alpha2;

  xp = product_center_1D(alpha1,x1,alpha2,x2);
  yp = product_center_1D(alpha1,y1,alpha2,y2);
  zp = product_center_1D(alpha1,z1,alpha2,z2);

  rab2 = dist2(x1,y1,z1,x2,y2,z2);
  rcp2 = dist2(x3,y3,z3,xp,yp,zp);

  Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma);
  Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma);
  Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma);

  sum = 0.;
  for (II=0; II<l1+l2+1; II++)
    for (J=0; J<m1+m2+1; J++)
      for (K=0; K<n1+n2+1; K++)
	sum += Ax[II]*Ay[J]*Az[K]*Fgamma(II+J+K,rcp2*gamma);

  free(Ax);
  free(Ay);
  free(Az);
  return -norm1*norm2*
    2.0*M_PI/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum;
}
C nuclear_attraction0(C x1, C y1, C z1,
		     int l1, int m1, int n1, C alpha1,
		     C x2, C y2, C z2,
		     int l2, int m2, int n2, C alpha2,
		     C x3, C y3, C z3){
  int II,J,K;
  C gamma,xp,yp,zp,sum,rab2,rcp2;
  C *Ax,*Ay,*Az;

  gamma = alpha1+alpha2;

  xp = product_center_1D(alpha1,x1,alpha2,x2);
  yp = product_center_1D(alpha1,y1,alpha2,y2);
  zp = product_center_1D(alpha1,z1,alpha2,z2);

  rab2 = dist2(x1,y1,z1,x2,y2,z2);
  rcp2 = dist2(x3,y3,z3,xp,yp,zp);

  Ax = A_array(l1,l2,xp-x1,xp-x2,xp-x3,gamma);
  Ay = A_array(m1,m2,yp-y1,yp-y2,yp-y3,gamma);
  Az = A_array(n1,n2,zp-z1,zp-z2,zp-z3,gamma);

  sum = 0.;
  for (II=0; II<l1+l2+1; II++)
    for (J=0; J<m1+m2+1; J++)
      for (K=0; K<n1+n2+1; K++)
	sum += Ax[II]*Ay[J]*Az[K]*Fgamma(II+J+K,rcp2*gamma);

  free(Ax);
  free(Ay);
  free(Az);
  return -2.0*M_PI/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum;
}
        
C A_term(int i, int r, int u, int l1, int l2,
		     C PAx, C PBx, C CPx, C gamma){
  /* THO eq. 2.18 */
  return pow(-1,i)*binomial_prefactor(i,l1,l2,PAx,PBx)*
    pow(-1,u)*(1.0*fact(i))*spow(CPx,i-2*r-2*u)*
    pow(0.25/gamma,r+u)/(1.0*fact(r))/(1.0*fact(u))/(1.0*fact(i-2*r-2*u));
}

C *A_array(int l1, int l2, C PA, C PB,
		C CP, C g){
  /* THO eq. 2.18 and 3.1 */
  int Imax,i,r,u,II;
  C *A;

  Imax = l1+l2+1;
  A = (C *)malloc(Imax*sizeof(C));
  for (i=0; i<Imax; i++) A[i] = 0.;
  for (i=0; i<Imax; i++)
    for (r=0; r<floor(i/2)+1;r++)
      for (u=0; u<floor((i-2*r)/2.)+1; u++){
	II = i-2*r-u;
	A[II] += A_term(i,r,u,l1,l2,PA,PB,CP,g);
      }
  return A;
}

int fact(int n){
  if (n <= 1) return 1;
  return n*fact(n-1);
}

int fact2(int n){ /* C factorial function = 1*3*5*...*n */
  if (n <= 1) return 1;
  return n*fact2(n-2);
}

C dist2(C x1, C y1, C z1,
		    C x2, C y2, C z2){
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}
C dist(C x1, C y1, C z1,
		   C x2, C y2, C z2){
  return sqrt(dist2(x1,y1,z1,x2,y2,z2));
}

// sum(t=0..t|s-ia<=t, t<=ib) (ia)C(s-t) (ib)C(t) xpa^(ia-s+t) xpb^()
C binomial_prefactor(int s, int ia, int ib, C xpa, C xpb){
  int t;
  C sum=0.;
  for (t=0; t<s+1; t++)
    if ((s-ia <= t) && (t <= ib)) 
      sum += (1.0*binomial(ia,s-t)*binomial(ib,t))*spow(xpa,ia-s+t)*spow(xpb,ib-t);
  return sum;
} 

// aCb = a!/(b! (a-b)!)
int binomial(int a, int b){return fact(a)/(fact(b)*fact(a-b));}

// check!
C Fgamma(C m, C x){
  C val;
  if (abs(x) < SMALL) x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5*pow(x,-m-0.5)*val; 
}

// check!
C gamm_inc(C a, C x){ /* Taken from NR routine gammap */
  C gamser,gammcf,gln;
  
  //  assert (abs(x) >= 0.);
  //  assert (real(a) > 0.);
  if (abs(x) < abs(a+1.0)) {
    gser(&gamser,a,x,&gln);
    return exp(gln)*gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return exp(gln)*(1.0-gammcf);
  }
}
 
void gser(C *gamser, C a, C x, C *gln){
  int n;
  C sum,del,ap;

  *gln=lgamma(a);
  /*
  if (real(x) <= 0.0) {
    assert(real(x)>=0.);
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ap = 1.0 + ap;
      del *= x/ap;
      sum += del;
      if (abs(del) < abs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser");
    return;
  }
  */
  ap=a;
  del=sum=1.0/a;
  for (n=1;n<=ITMAX;n++) {
    ap = 1.0 + ap;
    del *= x/ap;
    sum += del;
    if (abs(del) < abs(sum)*EPS) {
      *gamser=sum*exp(-x+a*log(x)-(*gln));
      return;
    }
  }
  printf("a too large, ITMAX too small in routine gser");
  return;
}
 
void gcf(C *gammcf, C a, C x, C *gln){
  int i;
  C an,b,c,d,del,h;
  
  *gln=lgamma(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -(1.0*i)*(1.0*i-a);
    b += 2.0;
    d=an*d+b;
    if (abs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (abs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (abs(del-1.0) < EPS) break;
  }
  assert(i<=ITMAX);
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

int ijkl2intindex(int i, int j, int k, int l){
  int tmp,ij,kl;
  if (i<j) return ijkl2intindex(j,i,k,l);
  if (k<l) return ijkl2intindex(i,j,l,k);
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij<kl){
    tmp = ij; 
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

int ijkl2intindex_old(int i, int j, int k, int l){
  int tmp,ij,kl;
  if (i<j){
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k<l){
    tmp = k;
    k = l;
    l = tmp;
  }
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij<kl){
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

int fact_ratio2(int a, int b){ return fact(a)/fact(b)/fact(a-2*b); }

C product_center_1D(C alphaa, C xa, 
			 C alphab, C xb){
  return (alphaa*xa+alphab*xb)/(alphaa+alphab);
}

C three_center_1D(C xi, int ai, C alphai,
			      C xj, int aj, C alphaj,
			      C xk, int ak, C alphak){

  C gamma, dx, px, xpi,xpj,xpk,intgl;
  int q,r,s,n;
  
  gamma = alphai+alphaj+alphak;
  dx = exp(-alphai*alphaj*spow(xi-xj,2)/gamma) *
    exp(-alphai*alphak*spow(xi-xk,2)/gamma) *
    exp(-alphaj*alphak*spow(xj-xk,2)/gamma);
  px = (alphai*xi+alphaj*xj+alphak*xk)/gamma;
    
  xpi = px-xi;
  xpj = px-xj;
  xpk = px-xk;
  intgl = 0;
  for (q=0; q<ai+1; q++){
    for (r=0; r<aj+1; r++){
      for (s=0; s<ak+1; s++){
	n = (q+r+s)/2;
	if ((q+r+s)%2 == 0) {
	  intgl += 1.0*binomial(ai,q)*binomial(aj,r)*binomial(ak,s)*
	    spow(xpi,ai-q)*spow(xpj,aj-r)*spow(xpk,ak-s)*
	    (1.0*fact2(2*n-1)/spow(2.0*gamma,n)*sqrt(M_PI/gamma));
	}
      }
    }
  }
  return dx*intgl;
}

#undef ITMAX
#undef EPS
#undef FPMIN
