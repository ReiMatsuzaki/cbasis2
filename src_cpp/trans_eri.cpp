#include "trans_eri.hpp"

using namespace Eigen;

namespace cbasis {

  void TransformERI_Slow(IB2EInt* ao, MO mo, IB2EInt* res) {
    
    /**
       Transform AO basis to MO basis for ERI
     */

    res->Init(ao->size());

    BMat& C = mo->C;

    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;

    // ==== set zero =====
    ao->Reset();
    while(ao->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) { 
      res->Set(ib,jb,kb,lb,i,j,k,l,0.0);
    }
    
    // ==== calculate ====
    ao->Reset();
    while(ao->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) { 
      MatrixXcd& Ci = C[make_pair(ib, ib)]; 
      MatrixXcd& Cj = C[make_pair(jb, jb)];
      MatrixXcd& Ck = C[make_pair(kb, kb)];
      MatrixXcd& Cl = C[make_pair(lb, lb)];
      int ni(Ci.rows());
      int nj(Cj.rows());
      int nk(Ck.rows());
      int nl(Cl.rows());
      // ii,jj,kk,ll are index for MO.
      for(int ii = 0; ii < ni; ii++)
	for(int jj = 0; jj < nj; jj++)
	  for(int kk = 0; kk < nk; kk++)
	    for(int ll = 0; ll < nl; ll++) {
	      dcomplex c = (Ci(i, ii) *
			    Cj(j, jj) *
			    Ck(k, kk) *
			    Cl(l, ll));
	      res->At(ib, jb, kb, lb, ii, jj, kk, ll) += c * v;
	    }
    }

  }
  void TransformERI(IB2EInt* ao, MO mo, IB2EInt* res) {

    BMat& C = mo->C;

    int n_irrep(C.size());
    vector<int> num_irrep(n_irrep);
    for(int irrep = 0; irrep < n_irrep; irrep++ ) {
      MatrixXcd& Ci = C[make_pair(irrep, irrep)];
      num_irrep[irrep] = Ci.rows();
    }

    /*
      int ib,jb,kb,lb,i,j,k,l,t;
      dcomplex v;
      ao->Reset();
    */

  }
  void CalcJK_MO(IB2EInt* eri_ao, BMat& C, int A0, int a0, int B0, int b0,
		 dcomplex cJ, dcomplex cK, BMat* bmat) {
    
    /**
       ijkl : ao index
       abcd : mo index
       AB   : A,B symmetry
       Compute 
       {{ <ab|a0,b0>  | a<-A, b<-B } | A,B}
       {{ <a,b0|a0,b> | a<-A, b<-B}} | A,B}
     */

    for(int A = 0; A < (int)C.size(); A++) {
      pair<Irrep, Irrep> II(A, A);
      MatrixXcd& CA = C[II];
      int num = CA.rows();
      (*bmat)[II] = MatrixXcd::Zero(num, num);
    }

    VectorXcd cA0 = C[make_pair(A0, A0)].col(a0);
    VectorXcd cB0 = C[make_pair(B0, B0)].col(b0);

    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;
    eri_ao->Reset();
    while(eri_ao->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {

      if(ib == jb && kb == A0 && lb == B0) {
	pair<Irrep, Irrep> AB(ib, jb);
	(*bmat)[AB](i, j) += cJ * cA0(k) * cB0(l) * v;
      }

      if(ib == lb && jb == B0 && kb == A0) {
	pair<Irrep, Irrep> AB(ib, lb);
	(*bmat)[AB](i, l) += cK * cA0(k) * cB0(j) * v;
      }
    }

    for(BMat::iterator it = bmat->begin(); it != bmat->end(); ++it) {
      MatrixXcd HA = it->second; // AO basis
      MatrixXcd CC = C[it->first];
      it->second = CC.transpose() * HA * CC;
    }

  }
}
