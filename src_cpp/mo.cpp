#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "mo.hpp"
#include "one_int.hpp"
#include "two_int.hpp"
#include "../utils/eigen_plus.hpp"

using namespace Eigen;

namespace cbasis {
  _MO::_MO() {}

  typedef pair<Irrep, dcomplex> IrrepComplex;
  struct Compare_IrrepEig {
    bool operator()(const IrrepComplex& a, const IrrepComplex& b) {
      return real(a.second) < real(b.second);
    }
  };
  vector<int> CalcOccNum(const BVec& eigs, int num_sym, int total_occ_orb) {
  /**
     Give occupation number from orbital energies for each symmetry.
     
     Inputs 
     ------
     const BVec& eigs : orbital energies for each symmetries.
     int num_sym      : number of symmetry
     int num_orb      : total number of occupied orbitals
     
     Outputs
     -------
     vector<int> res : number of occupied orbital for each symmetry
   */

    // -- build vector<IrrepIndexEig> --
    vector<IrrepComplex> irrep_eig_list;
    for(BVec::const_iterator it = eigs.begin(); it != eigs.end(); ++it) {

      Irrep irrep = it->first;
      const VectorXcd& eig = it->second;

      for(int index = 0; index < eig.size(); index++) {
	IrrepComplex ic(irrep, eig[index]);
	irrep_eig_list.push_back(ic);
      }
    }

    // -- sort --
    sort(irrep_eig_list.begin(), irrep_eig_list.end(), Compare_IrrepEig());

    // -- extract results --
    vector<int> res(num_sym, 0);
    for(int i = 0; i < total_occ_orb; i++) {
      res[irrep_eig_list[i].first]++;
    }
    return res;
  }
  vector<Irrep> CalcIrrepList(const BMat& bmat) {
    vector<Irrep> irrep_list;
    int num_block(bmat.size());
    for(Irrep irrep = 0; irrep < num_block; irrep++) {
      if(bmat.has_block(irrep, irrep))
	irrep_list.push_back(irrep);
    }
    return irrep_list;
  }
  vector<Irrep> CalcIrrepList(BMatSet mat_set) {
    //    cout << "T: " << endl;
    //    cout << mat_set->GetBlockMatrix("t");
    vector<Irrep> res= CalcIrrepList(mat_set->GetBlockMatrix("t"));
    //    cout << "size of res = " << res.size() << endl;
    return res;
    /*
    vector<Irrep> irrep_list;
    int num_block(mat_set->block_num());
    for(Irrep irrep = 0; irrep < num_block; irrep++) {
      if(mat_set->Exist("t", irrep, irrep)) {
	irrep_list.push_back(irrep);
      }
    } 

    return irrep_list;
    */   

  }
  MO NewMO(SymmetryGroup sym, BMat& _H, BMat& _C, BVec& _eigs, int num_ele) {

    MO mo;
    int nocc = num_ele/2;
    mo->sym = sym;
    mo->H = _H;
    mo->C = _C;
    mo->eigs = _eigs;
    mo->num_occ_irrep = CalcOccNum(mo->eigs, sym->num_class(), nocc);      
    mo->irrep_list = CalcIrrepList(_C);
    return mo;
    
  }
  MO CalcOneEle(SymmetryGroup sym, BMatSet mat_set, int) {

    MO mo(new _MO);

    // ---- get non0 symmetry ----
    mo->irrep_list = CalcIrrepList(mat_set);
    typedef vector<Irrep>::iterator It;
    for(It it = mo->irrep_list.begin(), end = mo->irrep_list.end();
	it != end; ++it) {
      Irrep irrep(*it);
      pair<Irrep, Irrep> ii(irrep, irrep);
      mo->H[ii] = (mat_set->GetMatrix("t", irrep, irrep) +
		   mat_set->GetMatrix("v", irrep, irrep));
      mo->S[ii] = mat_set->GetMatrix("s", irrep, irrep);
      mo->C[ii] = MatrixXcd::Zero(1, 1);
      mo->eigs[irrep] = VectorXcd::Zero(1);
      generalizedComplexEigenSolve(mo->H[ii], mo->S[ii], &mo->C[ii],
				   &mo->eigs[irrep]);
    }

    // ---- other value ----
    mo->sym = sym;
    mo->num_occ_irrep = CalcOccNum(mo->eigs, sym->num_class(), 1);
    for(It it = mo->irrep_list.begin(), end = mo->irrep_list.end();
	it != end; ++it) {
      if(mo->num_occ_irrep[*it] == 1)
	mo->energy = mo->eigs[*it][0];
    }
    return mo;
  }
  void AddJK(B2EInt eri, BMat& C, int I0, int i0,
	     dcomplex coef_J, dcomplex coef_K, BMat& H, bool use_real) {

    /*
      Add coef_J J + coef_K K to matrix H.
      J and K is build from I0 symmetry i0 th MO.
      eri is assumed from CalcERIComplex or CalcERIHermite.
    */
    
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;
    eri->Reset();
    const MatrixXcd& CC = C[make_pair(I0, I0)];

    if(use_real) {
      while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      
	if(ib == jb && kb == I0 && lb == I0) {
	  // Add J
	  pair<Irrep, Irrep> ij(ib, jb);
	  H[ij](i, j) += (coef_J * CC(k, i0) * CC(l, i0) * v).real();
	}
	
	if(ib == lb && jb == I0 && kb == I0) {
	  // Add K
	  pair<Irrep, Irrep> il(ib, lb);
	  H[il](i, l) += (coef_K * CC(j, i0) * CC(k, i0) * v).real();
	}
      }
    } else {

      while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      
	if(ib == jb && kb == I0 && lb == I0) {
	  // Add J
	  pair<Irrep, Irrep> ij(ib, jb);
	  H[ij](i, j) += coef_J * CC(k, i0) * CC(l, i0) * v;
	}
	
	if(ib == lb && jb == I0 && kb == I0) {
	  // Add K
	  pair<Irrep, Irrep> il(ib, lb);
	  H[il](i, l) += coef_K * CC(j, i0) * CC(k, i0) * v;
	}
      }
      
    }

  }
  void AddJK_Slow(B2EInt eri, BMat& C, int I0, int i0,
		  dcomplex coef_J, dcomplex coef_K, BMat& H) {

    VectorXcd c0 = C[make_pair(I0, I0)].col(i0);
    for(BMat::iterator it = H.begin(); it != H.end(); ++it) {

      pair<Irrep, Irrep> II(it->first);
      Irrep irrep(II.first);
      Irrep jrrep(II.second);
      MatrixXcd& HII = it->second;
      int n(HII.rows());
      int m(HII.cols());
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  dcomplex cumsum(0);
	  for(int a = 0; a < c0.size(); a++) {
	    for(int b = 0; b < c0.size(); b++) {
	      if(eri->Exist(irrep, jrrep, I0, I0, i, j, a, b)) 
		cumsum += coef_J*c0[a] * c0[b] *
		  eri->At(irrep, irrep, I0, I0, i, j, a, b);
	      if(eri->Exist(irrep, I0, I0, jrrep, i, a, b, j))
		cumsum += coef_K*c0[a] * c0[b] *
		  eri->At(irrep, I0, I0, jrrep, i, a, b, j);
	    }
	  }
	  HII(i, j) += cumsum;
	}
      }
    }
  }
  void AddJ(B2EInt eri, const VectorXcd& Ca, Irrep ir_a, dcomplex coef, BMat& J, bool use_real) {

    /**
       Add 
       (u_i| J_a |v_j) = (u_i(1)^* v_j(1) phi_a(2)^* phi_a(2))
       .               = (u_i(1)^* v_j(1) w_k(2)^* w_l(2)) C_ka C_la
       to mat J.

       In this function complex conjugate(*) is ignored.
       
       eri is calculated from
       .     eri = CalcERI(g_u, g_a, g_v, g_a, method)
     */
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;

    eri->Reset();

    while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      if(ib == jb && kb == ir_a && lb == ir_a) {
	pair<Irrep, Irrep> ij(ib, jb);
	if(use_real)
	  J[ij](i, j) += (v * Ca(k) * Ca(l) * coef).real();
	else
	  J[ij](i, j) += v * Ca(k) * Ca(l) * coef;
      }
    }

  }
  void AddK(B2EInt eri, const VectorXcd& Ca, Irrep ir_a, dcomplex coef, BMat& K, bool use_real) {
    
    /**
       Add
       (u_i | K_a | v_l) = (u_i(1)* phi_a(1) phi_a(2)* v_l(2) )
       .                 = (u_i(1)* w_j(1) w_k(2)* v_l(2)) C_ja C_ka
     */
    
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;

    eri->Reset();
    while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {

      if(ib == lb && jb == ir_a && kb == ir_a) {
	pair<Irrep, Irrep> il(ib, lb);
	if(use_real)
	  K[il](i, l) += (v * Ca(j) * Ca(k) * coef).real();
	else
	  K[il](i, l) += v * Ca(j) * Ca(k) * coef;
      }
    }    

  }
  /*
  void CalcJ(B2EInt eri, const Eigen::VectorXcd& Ca, Irrep ir_a, BMat *J) {
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;

    eri->Reset();

    while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      if(ib == jb && kb == ir_a && lb == ir_a) {	
	if(not J->has_block(ib, jb)) {
	  (*J) = 
	}
	pair<Irrep, Irrep> ij(ib, jb);
	J[ij](i, j) += v * Ca(k) * Ca(l);
      }
    }
  }
  */
  void CalcK(B2EInt eri, const Eigen::VectorXcd& Ca, Irrep ir_a, BMat *K) {
  }
  MO CalcRHF(SymGTOs gtos, int nele, ERIMethod method, const SCFOptions& opts,
	     bool *is_conv) {

    BMatSet mat_set = CalcMat_Complex(gtos, true);
    B2EInt eri      = CalcERI_Complex(gtos, method);
    SymmetryGroup sym = gtos->sym_group();

    MO mo0 = CalcOneEle(sym, mat_set, 0);
    BVec E0 = mo0->eigs;
    BMat C0 = mo0->C;    

    return CalcRHF(gtos->sym_group(), mat_set, eri, nele, E0, C0,
		   opts, is_conv);

  }
  MO CalcRHF(SymmetryGroup sym, BMatSet mat_set, B2EInt eri, int nele, 
	     const BVec& E0, const BMat& C0, const SCFOptions& opts,
	     bool *is_conv) {
    
    if(nele == 1) {
      *is_conv = true;
      return CalcOneEle(sym, mat_set, opts.debug_lvl);
    }
    if(nele % 2 == 1) {
      string msg; SUB_LOCATION(msg); msg += "nele must be even integer.";
      throw runtime_error(msg);
    }

    bool use_real = opts.use_real;
    int max_iter = opts.max_iter;
    double eps = opts.tol;
    int debug_lvl = opts.debug_lvl;
    
    int nocc = nele/2;
    *is_conv = false;
    MO mo(new _MO);
    mo->sym = sym;

    // ---- get non0 symmetry ----
    mo->irrep_list = CalcIrrepList(mat_set);

    // ---- initilize ----
    BMat FOld;    
    vector<int> num_irrep(sym->num_class(), 0);
    typedef vector<Irrep>::iterator It;
    for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it ) {
      Irrep irrep = *it;      
      pair<Irrep, Irrep> ii(irrep, irrep);
      mo->H[ii] = (mat_set->GetMatrix("t", irrep, irrep) +
		   mat_set->GetMatrix("v", irrep, irrep));
      mo->S[ii] = mat_set->GetMatrix("s", irrep, irrep);
      mo->F[ii] = mo->H[ii];
      int n(mo->H[ii].rows());
      num_irrep[irrep] = n;
      //      mo->C[ii] = MatrixXcd::Zero(n, n);
      mo->C[ii] = C0[ii];
      mo->P[ii] = MatrixXcd::Zero(n, n);
      FOld[ii] = MatrixXcd::Zero( n, n);
      //      mo->eigs[*it] = VectorXcd::Zero(n);
      mo->eigs[*it] = E0[*it];
    }

    // ---- SCF calculation ----
    for(int iter = 0; iter < max_iter; iter++) {
      
      // -- number of occupied orbitals --
      mo->num_occ_irrep = CalcOccNum(mo->eigs, sym->num_class(), nocc);      

      // -- density matrix --
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	int n_irrep(num_irrep[*it]);
	MatrixXcd& P_ii = mo->P[ii];
	MatrixXcd& C_ii = mo->C[ii];
	if(use_real) {
	  for(int i = 0; i < mo->num_occ_irrep[*it]; i++) 
	    for(int k = 0; k < n_irrep; k++)
	      for(int l = 0; l < n_irrep; l++)
		P_ii(k, l) = 2.0 * C_ii(k, i).real() * C_ii(l, i).real();
	} else {
	  for(int i = 0; i < mo->num_occ_irrep[*it]; i++) 
	    for(int k = 0; k < n_irrep; k++)
	      for(int l = 0; l < n_irrep; l++)
		P_ii(k, l) = 2.0 * C_ii(k, i) * C_ii(l, i);
	}
	
      }
      
      // -- update Fock matrix --
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	Irrep irrep(*it);
	pair<Irrep, Irrep> ii(make_pair(irrep, irrep));
	FOld[ii].swap(mo->F[ii]);
	mo->F[ii] = mo->H[ii].real().cast<dcomplex>();
      }
      VectorXcd ca = mo->C[make_pair(0,0)].col(0);
      AddJ(eri, ca, 0, 2.0,  mo->F, use_real);
      AddK(eri, ca, 0, -1.0, mo->F, use_real);

      // -- diag F matrix --
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	if(debug_lvl > 1) {
	  cout << "(irrep, im[F], im[S], im[H]) = " << *it
	       << " " << (mo->F[ii]-mo->F[ii].conjugate()).norm()
	       << " " << (mo->S[ii]-mo->S[ii].conjugate()).norm()
	       << " " << (mo->H[ii]-mo->H[ii].conjugate()).norm()
	       << endl;
	}
	generalizedComplexEigenSolve(mo->F[ii], mo->S[ii],
				     &mo->C[ii], &mo->eigs[*it]);
      }

      // -- print current status --
      if(debug_lvl > 0) {
	for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	  Irrep irrep = *it;
	  cout << irrep << ": " << mo->num_occ_irrep[irrep];
	  for(int i = 0; i < mo->num_occ_irrep[irrep]; i++) {
	    cout << mo->eigs[irrep][i] << ", ";
	  }
	  cout << endl;	
	}
	cout << endl;       
      }

      // -- continue if first
      if(iter == 0)
	continue;
      
      // -- Convergence check --
      bool conv(true);
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	if(abs((FOld[ii] - mo->F[ii]).mean()) > eps) {
	  conv = false;
	}
      }
      FOld = mo->F; // copy

      // -- calculate total energy --
      mo->energy = 0.0;
      for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
	pair<Irrep, Irrep> ii(make_pair(*it, *it));
	MatrixXcd& P(mo->P[ii]);
	MatrixXcd& H(mo->H[ii]);
	MatrixXcd& F(mo->F[ii]);
	int n(num_irrep[*it]);
	for(int i = 0; i < n; i++)
	  for(int j = 0; j < n; j++) 
	    mo->energy += 0.5 * P(i, j) * (H(i, j) + F(i, j));
      }
      
      // -- if convergence, break loop. --
      if(conv) {
	*is_conv = true;
	break;
      }
    }
    return mo;
  }
  void CalcSEHamiltonian(MO mo, B2EInt eri, Irrep I0, int i0, BMat* h_stex, int method) {

    typedef vector<Irrep>::iterator It;
    BMat res;
    
    // set res as Fock matrix
    for(It it = mo->irrep_list.begin(); it != mo->irrep_list.end(); ++it) {
      Irrep irrep(*it);
      pair<Irrep, Irrep> ii(make_pair(irrep, irrep));
      res[ii] = mo->H[ii];
    }

    // loop ERI and add to J and K
    if(method == 0) {
      int ib,jb,kb,lb,i,j,k,l,t;
      dcomplex v;
      eri->Reset();
      while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
	if(ib == jb && kb == I0 && lb == I0) {
	  // Add J
	  pair<Irrep, Irrep> ij(ib, jb), kl(I0, I0);
	  MatrixXcd& C = mo->C[kl];
	  res[ij](i, j) += C(k, i0) * C(l, i0) * v;
	}
	
	if(ib == lb && jb == I0 && kb == I0) {
	  // Add K
	  pair<Irrep, Irrep> il(ib, lb), jk(I0, I0);
	  MatrixXcd& C = mo->C[jk];
	  res[il](i, l) += C(j, i0) * C(k, i0) * v;
	}
      }  
    } else if(method == 2) {
      int ib,jb,kb,lb,i,j,k,l,t;
      dcomplex v;
      eri->Reset();
      while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
	if(ib == jb && kb == I0 && lb == I0) {
	  // Add J
	  pair<Irrep, Irrep> ij(ib, jb), kl(I0, I0);
	  MatrixXcd& C = mo->C[kl];
	  res[ij](i, j) += C(k, i0) * C(l, i0) * v;
	}
	/*
	if(ib == lb && jb == I0 && kb == I0) {
	  // Add K
	  pair<Irrep, Irrep> il(ib, lb), jk(I0, I0);
	  MatrixXcd& C = mo->C[jk];
	  res[il](i, l) += C(j, i0) * C(k, i0) * v;
	}
	*/
      }  
    } else if(method == 3) {
      int ib,jb,kb,lb,i,j,k,l,t;
      dcomplex v;
      eri->Reset();
      while(eri->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
	if(ib == jb && kb == I0 && lb == I0) {
	  // Add J
	  pair<Irrep, Irrep> ij(ib, jb), kl(I0, I0);
	  MatrixXcd& C = mo->C[kl];
	  res[ij](i, j) += C(k, i0) * C(l, i0) * v;
	}
	if(ib == lb && jb == I0 && kb == I0) {
	  // Add K
	  pair<Irrep, Irrep> il(ib, lb), jk(I0, I0);
	  MatrixXcd& C = mo->C[jk];
	  res[il](i, l) -= C(j, i0) * C(k, i0) * v;
	}
      }  
    } else if(method == 1) {
      vector<int> num_isym;
      SymmetryGroup sym = mo->sym;
      BMat J0, K0;
      for(int ib = 0; ib < sym->num_class(); ib++) {
	int num = mo->H[make_pair(ib, ib)].rows();
	J0[make_pair(ib ,ib)] = MatrixXcd::Zero(num, num);
	K0[make_pair(ib ,ib)] = MatrixXcd::Zero(num, num);
	num_isym.push_back(num);
      }

      for(int ijb = 0; ijb < sym->num_class(); ijb++)
      for(int klb = 0; klb < sym->num_class(); klb++) {
	if(klb == I0) {
	  MatrixXcd& Jij = J0[make_pair(ijb, ijb)];
	  MatrixXcd& Kij = K0[make_pair(ijb, ijb)];
	  MatrixXcd& Ckl = mo->C[make_pair(klb, klb)];
	  for(int i  = 0; i  < num_isym[ijb]; i++)
	  for(int j  = 0; j  < num_isym[ijb]; j++) {
	    dcomplex tmp_j(0), tmp_k(0);
	    for(int k  = 0; k  < num_isym[klb]; k++)
	    for(int l  = 0; l  < num_isym[klb]; l++) {
	      tmp_j += Ckl(k, i0) * Ckl(l, i0) * eri->At(ijb, ijb, klb, klb,
							 i,   j,   k,   l);
	      tmp_k += Ckl(k, i0) * Ckl(l, i0) * eri->At(ijb, klb, klb, ijb,
							 i,   l,   k,   j);
	    }
	    Jij(i, j) = tmp_j;
	    Kij(i, j) = tmp_k;
	  }
	  res[make_pair(ijb, ijb)] += Jij + Kij;
	}
      }
    }

    // swap
    h_stex->swap(res);
  }
  dcomplex CalcAlpha(MO mo, BMatSet mat_set, Irrep I0, int i0, BMat& h_stex,
		     double w, Coord coord, int method) {

    typedef vector<Irrep>::const_iterator It;
    dcomplex ene;
    if(method == 0) {
      // -- Koopsman's theorem --
      ene = w + mo->eigs[I0](i0);
    } else if(method == 1){
      // -- calulate ionization potential using HF energy --
      double E_ion = -2.0;
      ene = w + mo->energy - E_ion;
    } else if(method == 2) {
      // -- experimental value from "McQuarrie and Simon" p.302 --
      ene = w - 0.903724375;
    }

    // -- set matrix name for direction --
    string mat_name;
    if(coord == CoordX)
      mat_name = "x";
    else if(coord == CoordY)
      mat_name = "y";
    else if(coord == CoordZ)
      mat_name = "z";

    dcomplex a(0);
    It it; 
    It end= mo->irrep_list.end();
    for(it = mo->irrep_list.begin(); it != end; ++it) {      
      Irrep irrep = *it;
      if(mat_set->Exist(mat_name, irrep, I0)) {
	pair<Irrep, Irrep> ii(irrep, irrep);
	const MatrixXcd& S = mat_set->GetMatrix("s", irrep, irrep);
	const MatrixXcd& H = h_stex[ii];
	const MatrixXcd& Z = mat_set->GetMatrix(mat_name, irrep, I0);
	
	MatrixXcd L = S * ene - H;
	VectorXcd m = Z * mo->C[make_pair(I0, I0)].col(i0);
	VectorXcd c = L.colPivHouseholderQr().solve(m);
	int num(m.size());
	for(int i = 0; i < num; i++) {
	  a += m[i] * c[i];
	}
      }
    }
    return a;
  }
  double PITotalCrossSection(dcomplex alpha, double w, int num_occ_ele) {
    double au2mb(pow(5.291772, 2));
    double c(137.035999258);
    double c0 = -4.0 * M_PI * w / c * imag(alpha) * au2mb;
    if(num_occ_ele == 1)
      return c0;
    else
      return 2.0 * c0;
  }
}
