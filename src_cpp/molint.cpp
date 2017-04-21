#include <iostream>
#include <iomanip>
#include "molint.hpp"
#include "macros.hpp"
#include "angmoment.hpp"
#include "spec_func.hpp"

namespace cbasis {

  using namespace std;
  typedef MultArray<dcomplex, 3> A3dc;

  // ==== Utility ====
  double div(int a, int b) {
    return a*1.0/b;
  }

  dcomplex gto_int(dcomplex z, int n) {
    
    dcomplex res;
    if(n % 2 == 0) {
      
      int nn = n/2;
      res = DoubleFactorial(2*nn-1) * sqrt(M_PI) /
	(dcomplex(pow(2, nn+1)) * pow(sqrt(z), 2*nn+1));
    
    } else {
      
      int nn = (n-1)/2;
      res = dcomplex(Factorial(nn)) / (dcomplex(2) * pow(z, nn+1));
    
    }
    return res;
  } 

  dcomplex gto_overlap(int nAx, int nAy, int nAz, 
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB) {

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx    = (zetaA * wAx + zetaB * wBx) / zetaP;
    dcomplex wPy    = (zetaA * wAy + zetaB * wBy) / zetaP;
    dcomplex wPz    = (zetaA * wAz + zetaB * wBz) / zetaP;
		       
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*dist2(wAx-wBx, wAy-wBy, wAz-wBz));

    return (eAB * pow(M_PI/zetaP, 1.5) *
	    coef_d(zetaP, wPx, wAx, wBx, nAx, nBx, 0) *
	    coef_d(zetaP, wPy, wAy, wBy, nAy, nBy, 0) *
	    coef_d(zetaP, wPz, wAz, wBz, nAz, nBz, 0));

  }

  dcomplex gto_kinetic(int nAx, int nAy, int nAz,
		       dcomplex wAx, dcomplex wAy, dcomplex wAz,
		       dcomplex zetaA,
		       int nBx, int nBy, int nBz, 
		       dcomplex wBx, dcomplex wBy, dcomplex wBz,
		       dcomplex zetaB) {
    dcomplex res(0.0);
    dcomplex s000 = gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx, nBy, nBz, wBx, wBy, wBz, zetaB);
    dcomplex sp00 = gto_overlap(nAx,   nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx+2, nBy, nBz, wBx, wBy, wBz, zetaB);
    dcomplex s0p0 = gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx, nBy+2, nBz, wBx, wBy, wBz, zetaB);
    dcomplex s00p = gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
				nBx, nBy, nBz+2, wBx, wBy, wBz, zetaB);
    res += -2.0*zetaB*(2.0*nBx+2*nBy+2*nBz+3.0) * s000;
    res += 4.0*zetaB*zetaB*(sp00 + s0p0 + s00p);
    if(nBx > 1)
      res += 1.0*nBx*(nBx-1)* gto_overlap(nAx,   nAy, nAz, wAx, wAy, wAz, zetaA,
					  nBx-2, nBy, nBz, wBx, wBy, wBz, zetaB);
    
    if(nBy > 1)
      res += 1.0*nBy*(nBy-1)* gto_overlap(nAx, nAy,   nAz, wAx, wAy, wAz, zetaA,
					  nBx, nBy-2, nBz, wBx, wBy, wBz, zetaB);     

    if(nBz > 1)
      res += 1.0*nBz*(nBz-1)* gto_overlap(nAx, nAy, nAz,   wAx, wAy, wAz, zetaA,
					  nBx, nBy, nBz-2, wBx, wBy, wBz, zetaB);
    

    return -0.5*res;
  }

  dcomplex gto_moment_z(int nAx, int nAy, int nAz,
			dcomplex wAx, dcomplex wAy, dcomplex wAz,
			dcomplex zetaA,
			int nBx, int nBy, int nBz, 
			dcomplex wBx, dcomplex wBy, dcomplex wBz,
			dcomplex zetaB) {
    
    return (gto_overlap(nAx, nAy, nAz+1, wAx, wAy, wAz, zetaA,
			nBx, nBy, nBz,   wBx, wBy, wBz, zetaB)
	    + wAz *
	    gto_overlap(nAx, nAy, nAz, wAx, wAy, wAz, zetaA,
			nBx, nBy, nBz, wBx, wBy, wBz, zetaB));
    
  }
  
  dcomplex gto_nuclear_attraction(int nAx, int nAy, int nAz, 
				  dcomplex wAx, dcomplex wAy, dcomplex wAz,
				  dcomplex zetaA,
				  int nBx, int nBy, int nBz, 
				  dcomplex wBx, dcomplex wBy, dcomplex wBz,
				  dcomplex zetaB,
				  dcomplex wCx, dcomplex wCy, dcomplex wCz) {

    dcomplex zetaP = zetaA + zetaB;
    dcomplex wPx    = (zetaA * wAx + zetaB * wBx) / zetaP;
    dcomplex wPy    = (zetaA * wAy + zetaB * wBy) / zetaP;
    dcomplex wPz    = (zetaA * wAz + zetaB * wBz) / zetaP;
		       
    dcomplex eAB = exp(-zetaA*zetaB/zetaP*dist2(wAx-wBx, wAy-wBy, wAz-wBz));

    dcomplex* dxs = new dcomplex[nAx+nBx+1];
    for(int nx = 0; nx <= nAx+nBx; nx++)
      dxs[nx] = coef_d(zetaP, wPx, wAx, wBx, nAx, nBx, nx);

    dcomplex* dys = new dcomplex[nAy+nBy+1];
    for(int ny = 0; ny <= nAy+nBy; ny++)
      dys[ny] = coef_d(zetaP, wPy, wAy, wBy, nAy, nBy, ny);

    dcomplex* dzs = new dcomplex[nAz+nBz+1];
    for(int nz = 0; nz <= nAz+nBz; nz++) 
      dzs[nz] = coef_d(zetaP, wPz, wAz, wBz, nAz, nBz, nz);
    
    dcomplex* Fjs = new dcomplex[nAx+nAy+nAz+nBx+nBy+nBz+1];
    IncompleteGamma(nAx+nAy+nAz+nBx+nBy+nBz,
		    zetaP * dist2(wPx-wCx, wPy-wCy, wPz-wCz), Fjs);

    dcomplex cumsum(0);    
    for(int nx = 0; nx <= nAx+nBx; nx++)
      for(int ny = 0; ny <= nAy+nBy; ny++)
	for(int nz = 0; nz <= nAz+nBz; nz++) {
	  cumsum += (-2.0*M_PI/zetaP *
		     dxs[nx] * dys[ny] * dzs[nz] * 
		     coef_R(zetaP, wPx, wPy, wPz, wCx, wCy, wCz,
			    nx, ny, nz, 0, Fjs));
	}
    delete dxs;
    delete dys;
    delete dzs;

    /*
    dcomplex cumsum(0);    
    for(int nx = 0; nx <= nAx+nBx; nx++)
      for(int ny = 0; ny <= nAy+nBy; ny++)
	for(int nz = 0; nz <= nAz+nBz; nz++) {
	  cumsum += (-2.0*M_PI/zetaP *
		     coef_d(zetaP, wPx, wAx, wBx, nAx, nBx, nx) * 
		     coef_d(zetaP, wPy, wAy, wBy, nAy, nBy, ny) * 
		     coef_d(zetaP, wPz, wAz, wBz, nAz, nBz, nz) * 
		     coef_R(zetaP, wPx, wPy, wPz, wCx, wCy, wCz,
			    nx, ny, nz, 0));
	}
	*/
    return eAB * cumsum;
  //  dcomplex coef_d(int max_mA, int max_mB, dcomplex* ds, dcomplex zetap,
  }

  MatrixSet::MatrixSet(int nbi, int nbj) {
    nbasis_i_ = nbi;
    nbasis_j_ = nbj;
  }
  MatrixSet::~MatrixSet() {

    /*
    for(std::map<std::string, dcomplex*>::iterator it = mat_map_.begin();
	it != mat_map_.end(); ++it) {
      if(it->second != NULL) {
	std::cout <<it->first << std::endl;
	delete[] (it->second);
      }
    }    
    */
  }
  int MatrixSet::size_basis_i() const { return nbasis_i_;}
  int MatrixSet::size_basis_j() const { return nbasis_j_;}
  void MatrixSet::set(std::string label, dcomplex* ptr) {
    mat_map_[label] = ptr;;
  }
  dcomplex* MatrixSet::get(std::string label) {
    return mat_map_[label];
  }

  GTOs::GTOs() {
    offset_ish.push_back(0);
  }
  void GTOs::Add(dcomplex _zeta,
		 dcomplex x, dcomplex y, dcomplex z,
		 vector<int> _nx, vector<int> _ny, vector<int> _nz, 
		 vector<vector<dcomplex> > _coef, int L, vector<int> ms) {
    // new shell index
    int ish = this->size_sh();

    // offset for the new shell
    int offset0 = offset_ish[ish];

    int basis_size0 = _coef.size();
    offset_ish.push_back(offset0 + basis_size0);

    zeta_ish.push_back(_zeta);
    x_ish.push_back(x);
    y_ish.push_back(y);
    z_ish.push_back(z);
    nx_ish_iprim.push_back(_nx);
    ny_ish_iprim.push_back(_ny);
    nz_ish_iprim.push_back(_nz);
    coef_ish_ibasis_iprim.push_back(_coef);

    l_ish.push_back(L);
    m_ish_ibasis.push_back(ms);

    dcomplex c2 = gto_int(_zeta*2.0, 2*L+2);
    coef_ylm_ish.push_back(1.0/sqrt(c2));
    
  }
  void GTOs::AddSphericalGTOs(int L, dcomplex x, dcomplex y, dcomplex z,
			      dcomplex _zeta) {
    
    if(L == 0) {
      vector<int> n0; n0.push_back(0);
      vector<dcomplex> c1; c1.push_back(1.0);      
      vector<vector<dcomplex> > c; c.push_back(c1);
      vector<int> ms(1); ms[0] = 0;
      this->Add(_zeta, x, y, z, n0, n0, n0, c, 0, ms);     
    }
    else if(L == 1) {
      vector<int> nx(3), ny(3), nz(3);
      nx[0] = 1; ny[0] = 0; nz[0] = 0;
      nx[1] = 0; ny[1] = 1; nz[1] = 0;
      nx[2] = 0; ny[2] = 0; nz[2] = 1;

      vector<int> ms(3);
      vector<vector<dcomplex> > cs(3, vector<dcomplex>(3));
      cs[0][0] = 0.0; cs[0][1] = 1.0; cs[0][2] = 0.0; ms[0] = -1; // M =-1
      cs[1][0] = 0.0; cs[1][1] = 0.0; cs[1][2] = 1.0; ms[1] = +0; // M = 0
      cs[2][0] = 1.0; cs[2][1] = 0.0; cs[2][2] = 0.0; ms[2] = +1; // M = 1

      this->Add(_zeta, x, y, z, nx, ny, nz, cs, L, ms);     
    }
    else if(L == 2) {
      vector<int> nx(6), ny(6), nz(6);
      nx[0] = 2; ny[0] = 0; nz[0] = 0; 
      nx[1] = 0; ny[1] = 2; nz[1] = 0; 
      nx[2] = 0; ny[2] = 0; nz[2] = 2; 
      nx[3] = 1; ny[3] = 1; nz[3] = 0; 
      nx[4] = 1; ny[4] = 0; nz[4] = 1;
      nx[5] = 0; ny[5] = 1; nz[5] = 1;
      
      vector<int> ms(5);
      vector<vector<dcomplex> > cs(5, vector<dcomplex>(6));
      // M = -2
      ms[0] = -2;
      cs[0][0] = 0.0; cs[0][1] = 0.0; cs[0][2] = 0.0;
      cs[0][3] = 1.0; cs[0][4] = 0.0; cs[0][5] = 0.0; 

      // M = -1
      ms[1] = -1;
      cs[1][0] = 0.0; cs[1][1] = 0.0; cs[1][2] = 0.0;
      cs[1][3] = 0.0; cs[1][4] = 0.0; cs[1][5] = 1.0; 

      // M = 0
      ms[2] = 0;
      cs[2][0] = -1.0; cs[2][1] = -1.0; cs[2][2] = 2.0;
      cs[2][3] = +0.0; cs[2][4] = +0.0; cs[2][5] = 0.0; 

      // M = 1
      ms[3] = +1;
      cs[3][0] = 0.0; cs[3][1] = 0.0; cs[3][2] = 0.0;
      cs[3][3] = 0.0; cs[3][4] = 1.0; cs[3][5] = 0.0; 

      // M = 2
      ms[4] = 2;
      cs[4][0] = 1.0; cs[4][1] =-1.0; cs[4][2] = 0.0;
      cs[4][3] = 0.0; cs[4][4] = 0.0; cs[4][5] = 0.0; 

      this->Add(_zeta, x, y, z, nx, ny, nz, cs, L, ms);     

    } else {
      std::string msg;
      SUB_LOCATION(msg);      
      msg += "Not implemented yet for this L";
      throw std::runtime_error(msg);      
    }
    
    this->Normalize();

  }
  void GTOs::AddOneSphericalGTO(int L, int M, dcomplex x, dcomplex y, dcomplex z,
				dcomplex _zeta) {
    
    
    vector<int> nx(1), ny(1), nz(1), ms(1);
    vector<vector<dcomplex> > cs(1, vector<dcomplex>(1));
    ms[0] = M;
    if(L == 0 && M == 0) {
      nx[0] = 0; ny[0] = 0; nz[0] = 0; cs[0][0] = 1.0; 
    }
    else if(L == 1 && M == -1) {
      nx[0] = 0; ny[0] = 1; nz[0] = 0; cs[0][0] = 1.0;
    }
    else if(L == 1 && M == 0) {
      nx[0] = 0; ny[0] = 0; nz[0] = 1; cs[0][0] = 1.0;
    }
    else if(L == 1 && M == 1) {
      nx[0] = 1; ny[0] = 0; nz[0] = 0; cs[0][0] = 1.0;
    }
    else if(L == 2) {
      if(M == -2) {
	nx[0] = 1; ny[0] = 1; nz[0] = 0; cs[0][0]= 1.0; 
      }
      else if(M == -1) {
	nx[0] = 0; ny[0] = 1; nz[0] = 1; cs[0][0] = 1.0;
      }
      else if(M == 0) {
	nx.push_back(2); ny.push_back(0); nz.push_back(0);
	nx.push_back(0); ny.push_back(2); nz.push_back(0);
	nx.push_back(0); ny.push_back(0); nz.push_back(2);
	cs[0].push_back(-1.0);
	cs[0].push_back(-1.0);
	cs[0].push_back(2.0);
      }
      else if(M == +1) {
	nx[0] = 1; ny[0] = 0; nz[0] = 1; cs[0][0] = 1.0;
      }
      else if(M == +2) {
	nx.push_back(2); ny.push_back(0); nz.push_back(0);
	nx.push_back(0); ny.push_back(2); nz.push_back(0);
	cs[0].push_back(+1.0);
	cs[0].push_back(-1.0);
      }

    } else {
      std::string msg;
      SUB_LOCATION(msg);      
      msg += "Not implemented yet for this L";
      throw std::runtime_error(msg);      
    }

    this->Add(_zeta, x, y, z, nx, ny, nz, cs, L, ms);     
    this->Normalize();

  }
  void GTOs::AddAtom(dcomplex q, dcomplex x, dcomplex y, dcomplex z) {
    x_iat.push_back(x);
    y_iat.push_back(y);
    z_iat.push_back(z);
    q_iat.push_back(q);
  }
  void GTOs::Normalize() {

    for(int ish = 0; ish < this->size_sh(); ish++) {
      for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	dcomplex cumsum(0.0);
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	  for(int jprim = 0; jprim < this->size_prim_ish(ish); jprim++) {
	    dcomplex s = gto_overlap(nx_ish_iprim[ish][iprim],
				     ny_ish_iprim[ish][iprim],
				     nz_ish_iprim[ish][iprim],
				     x_ish[ish], y_ish[ish], z_ish[ish],
				     zeta_ish[ish],
				     nx_ish_iprim[ish][jprim],
				     ny_ish_iprim[ish][jprim],
				     nz_ish_iprim[ish][jprim],
				     x_ish[ish], y_ish[ish], z_ish[ish],
				     zeta_ish[ish]);
	    s *= coef_ish_ibasis_iprim[ish][ibasis][iprim];
	    s *= coef_ish_ibasis_iprim[ish][ibasis][jprim];
	    cumsum += s;
	  }
	}
	dcomplex scale = 1.0/sqrt(cumsum);	
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) 
	  coef_ish_ibasis_iprim[ish][ibasis][iprim] *= scale;	
      }
    }    
  }
  MatrixSet GTOs::Calc() {
    dcomplex *s, *t, *dz, *v;
    this->CalcMat(&s, &t, &dz, &v);
    MatrixSet matrix_set(this->size_basis(), this->size_basis());
    matrix_set.set("s", s);
    matrix_set.set("t", t);
    matrix_set.set("v", v);
    matrix_set.set("z", dz);
    return matrix_set;
  }
  void GTOs::CalcMat(dcomplex** s, dcomplex** t, dcomplex** dz, dcomplex** v) {
    
    int nb = this->size_basis();    
    int np = this->size_prim();
    dcomplex* S = new dcomplex[nb*nb];
    dcomplex* T = new dcomplex[nb*nb];
    dcomplex* Dz = new dcomplex[nb*nb];
    dcomplex* V = new dcomplex[nb*nb];    
    dcomplex* smat_prim = new dcomplex[np*np];
    dcomplex* zmat_prim = new dcomplex[np*np];    
    dcomplex* vmat_prim = new dcomplex[np*np];
    dcomplex* tmat_prim = new dcomplex[np*np];
    dcomplex* dsx_buff = new dcomplex[5*7*10];
    dcomplex* dsy_buff = new dcomplex[5*7*10];
    dcomplex* dsz_buff = new dcomplex[5*7*10];
    dcomplex** Fjs_iat = new dcomplex*[this->size_atom()];
    
    for(int idx=0; idx<nb*nb; idx++) {
      S[idx] = 7.7; T[idx] = 7.7; Dz[idx] = 7.7; V[idx] = 7.7;
    }

    // compute maximum of nk (k=x,y,z) for each shell.
    int maxn(0);
    int* maxn_ish = new int[this->size_sh()];
    for(int ish = 0; ish < this->size_sh(); ish++) {
      int sum_ni(0);
      for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++) {
	  int ni = (nx_ish_iprim[ish][iprim] +
		    ny_ish_iprim[ish][iprim] +
		    nz_ish_iprim[ish][iprim]); 
	  sum_ni = sum_ni < ni ? ni : sum_ni;
      }
      maxn_ish[ish] = sum_ni;
      if(maxn < sum_ni)
	maxn = sum_ni;
    }
    for(int iat = 0; iat < this->size_atom(); iat++) {
      Fjs_iat[iat] = new dcomplex[2*maxn+1];
    }
    
    for(int ish = 0; ish < this->size_sh(); ish++) {
      int jsh0 = ish;
      for(int jsh = jsh0; jsh < this->size_sh(); jsh++) {
	dcomplex zetai = zeta_ish[ish]; dcomplex zetaj = zeta_ish[jsh];	
	dcomplex zetaP = zetai + zetaj;
	dcomplex xi(x_ish[ish]); dcomplex yi(y_ish[ish]); dcomplex zi(z_ish[ish]);
	dcomplex xj(x_ish[jsh]); dcomplex yj(y_ish[jsh]); dcomplex zj(z_ish[jsh]);
	dcomplex wPx    = (zetai*xi + zetaj*xj)/zetaP;
	dcomplex wPy    = (zetai*yi + zetaj*yj)/zetaP;
	dcomplex wPz    = (zetai*zi + zetaj*zj)/zetaP;	
	dcomplex d2 = dist2(xi-xj,yi-yj,zi-zj);	
	dcomplex eAB = exp(-zetai*zetaj/zetaP*d2);
	dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);

	int mi = maxn_ish[ish]; int mj = maxn_ish[jsh];
	A3dc dxmap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPx,xi,xj,dsx_buff);
	A3dc dymap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPy,yi,yj,dsy_buff);
	A3dc dzmap = calc_d_coef(mi,mj+2,mi+mj,zetaP,wPz,zi,zj,dsz_buff);
	
	for(int iat = 0; iat < this->size_atom(); iat++) {
	  dcomplex cx(x_iat[iat]); dcomplex cy(y_iat[iat]); dcomplex cz(z_iat[iat]);
	  IncompleteGamma(maxn_ish[ish] + maxn_ish[jsh],
			  zetaP * dist2(wPx-cx,wPy-cy,wPz-cz),
			  Fjs_iat[iat]);
	}

	int npi = this->size_prim_ish(ish);
	int npj = this->size_prim_ish(jsh);
	for(int iprim = 0; iprim < npi; iprim++) {
	  int jprim0 = ish == jsh ? iprim : 0;
	  for(int jprim = jprim0; jprim < npj; jprim++) {
	    int nxi = nx_ish_iprim[ish][iprim]; int nxj = nx_ish_iprim[jsh][jprim];
	    int nyi = ny_ish_iprim[ish][iprim]; int nyj = ny_ish_iprim[jsh][jprim];
	    int nzi = nz_ish_iprim[ish][iprim]; int nzj = nz_ish_iprim[jsh][jprim];

	    // ---- S mat ----
	    dcomplex dx00 = dxmap.get(nxi, nxj, 0);
	    dcomplex dy00 = dymap.get(nyi, nyj, 0);
	    dcomplex dz00 = dzmap.get(nzi, nzj, 0);

	    // ---- z mat ----
	    dcomplex dz01 = dzmap.get(nzi, nzj+1, 0);
	    
	    // ---- t mat ----
	    dcomplex dx02 = dxmap.get(nxi, nxj+2, 0);
	    dcomplex dy02 = dymap.get(nyi, nyj+2, 0);
	    dcomplex dz02 = dzmap.get(nzi, nzj+2, 0);
	    dcomplex t_ele(0.0);
	    t_ele += -2.0*zetaj * (2*nxj+2*nyj+2*nzj+3.0) * dx00*dy00*dz00;
	    t_ele += 4.0*zetaj*zetaj*(dx02*dy00*dz00+dx00*dy02*dz00+dx00*dy00*dz02);
	    if(nxj > 1) {
	      dcomplex dx = dxmap.get(nxi, nxj-2, 0);
	      t_ele += 1.0*nxj*(nxj-1) * dx * dy00 * dz00;
	    }
	    if(nyj > 1) {
	      dcomplex dy = dymap.get(nyi, nyj-2, 0);
	      t_ele += 1.0*nyj*(nyj-1) * dx00 * dy * dz00;
	    }
	    if(nzj > 1) {
	      dcomplex dz = dzmap.get(nzi, nzj-2, 0);
	      t_ele += 1.0*nzj*(nzj-1) * dx00 * dy00 * dz;
	    }
	    
	    dcomplex v_ele(0);
	    for(int nx = 0; nx <= nxi + nxj; nx++)
	      for(int ny = 0; ny <= nyi + nyj; ny++)
		for(int nz = 0; nz <= nzi + nzj; nz++)
		  for(int iat = 0; iat < this->size_atom(); iat++) {
		    v_ele += (q_iat[iat] *
			      dxmap.get(nxi, nxj, nx) *
			      dymap.get(nyi, nyj, ny) *
			      dzmap.get(nzi, nzj, nz) *
			      coef_R(zetaP, wPx, wPy, wPz,
				     x_iat[iat], y_iat[iat], z_iat[iat],
				     nx, ny, nz, 0, Fjs_iat[iat]));
		  
		}
	    	    
	    int idx = iprim * npj + jprim;
	    smat_prim[idx] = ce * dx00 * dy00 * dz00;
	    zmat_prim[idx] = ce*(dx00*dy00*dz01 + zj*dx00*dy00*dz00);
	    tmat_prim[idx] = -0.5 * ce * t_ele;
	    vmat_prim[idx] = -2.0*M_PI/zetaP * eAB * v_ele;
	    if(ish == jsh && iprim != jprim) {
	      int idx1 = jprim * npj + iprim;
	      smat_prim[idx1] = smat_prim[idx];
	      zmat_prim[idx1] =zmat_prim[idx];
	      tmat_prim[idx1] = tmat_prim[idx];
	      vmat_prim[idx1] = vmat_prim[idx];
	    }
	  }
	}

	for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	  int jbasis0 = ish == jsh ? ibasis : 0;
	  for(int jbasis = jbasis0; jbasis < this->size_basis_ish(jsh); jbasis++) {
	    dcomplex s(0.0); dcomplex dz(0.0); dcomplex t(0.0);	dcomplex v(0.0);
	    
	    for(int iprim = 0; iprim < npi; iprim++) {
	      
	      for(int jprim = 0; jprim < npj; jprim++) {
		int idx_prim0(iprim * npj + jprim);
		dcomplex cc(coef_ish_ibasis_iprim[ish][ibasis][iprim]*
			    coef_ish_ibasis_iprim[jsh][jbasis][jprim]);
		s += cc*smat_prim[idx_prim0];
		t += cc*tmat_prim[idx_prim0];
		dz += cc*zmat_prim[idx_prim0];
		v += cc*vmat_prim[idx_prim0];
		idx_prim0++;
	      }
	    }
	    int i = offset_ish[ish] + ibasis;
	    int j = offset_ish[jsh] + jbasis;
	    int idx = i+j*nb;	
	    S[idx] = s; Dz[idx] = dz; T[idx] = t; V[idx] = v;
	    if(ish != jsh || ibasis != jbasis) {
	      int idx1 = j + i*nb;
	      S[idx1] = s; Dz[idx1] = dz; T[idx1] = t; V[idx1] = v;
	    }
	  }
	}

      }
    }
    *s = S; *dz = Dz; *t = T; *v = V;    
    delete smat_prim;
    delete tmat_prim;
    delete zmat_prim;
    delete vmat_prim;
    delete dsx_buff;
    delete dsy_buff;
    delete dsz_buff;
    for(int iat = 0; iat < this->size_atom(); iat++)
      delete Fjs_iat[iat];
    delete[] Fjs_iat;
  }
  MatrixSet GTOs::CalcZMatOther(const GTOs& o) {

    dcomplex* Z = new dcomplex[this->size_basis() * o.size_basis()];
    dcomplex* z_prim = new dcomplex[this->size_prim() * o.size_prim()];
    dcomplex* dsx_buff = new dcomplex[5*7*10];
    dcomplex* dsy_buff = new dcomplex[5*7*10];
    dcomplex* dsz_buff = new dcomplex[5*7*10];

    for(int ish = 0; ish < this->size_sh(); ish++)
      for(int jsh = 0; jsh < o.size_sh(); jsh++) {

	// -- compute d coeff --
	dcomplex zetai = this->zeta_ish[ish];
	dcomplex zetaj =  o.zeta_ish[jsh];
	dcomplex zetaP = zetai + zetaj;
	dcomplex xi(this->x_ish[ish]); dcomplex xj(o.x_ish[jsh]);
	dcomplex yi(this->y_ish[ish]); dcomplex yj(o.y_ish[jsh]);
	dcomplex zi(this->z_ish[ish]); dcomplex zj(o.z_ish[jsh]);
	dcomplex wPx = (zetai*xi + zetaj*xj)/zetaP;
	dcomplex wPy = (zetai*yi + zetaj*yj)/zetaP;
	dcomplex wPz = (zetai*zi + zetaj*zj)/zetaP;
	dcomplex eAB = exp(-zetai*zetaj/zetaP*dist2(xi-xj,yi-yj,zi-zj));
	dcomplex ce = eAB * pow(M_PI/zetaP, 1.5);
	A3dc dxmap = calc_d_coef(this->l_ish[ish], o.l_ish[jsh], 0,
				 zetaP, wPx, xi, xj, dsx_buff);
	A3dc dymap = calc_d_coef(this->l_ish[ish], o.l_ish[jsh], 0,
				 zetaP, wPy, yi, yj, dsy_buff);
	A3dc dzmap = calc_d_coef(this->l_ish[ish], o.l_ish[jsh]+1, 0,
				 zetaP, wPz, zi, zj, dsz_buff);

	// -- compute primitive matrix --
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++)
	  for(int jprim = 0; jprim < o.size_prim_ish(jsh); jprim++) {
	    dcomplex dx00 = dxmap.get(this->nx_ish_iprim[ish][iprim],
				      o.nx_ish_iprim[jsh][jprim], 0);
	    dcomplex dy00 = dymap.get(this->ny_ish_iprim[ish][iprim],
				      o.ny_ish_iprim[jsh][jprim], 0);
	    dcomplex dz00 = dzmap.get(this->nz_ish_iprim[ish][iprim],
				      o.nz_ish_iprim[jsh][jprim], 0);	    
	    dcomplex dz01 = dzmap.get(this->nz_ish_iprim[ish][iprim],
				      o.nz_ish_iprim[jsh][jprim] + 1, 0);
	    z_prim[iprim * o.size_prim_ish(jsh) + jprim] =
	      ce * (dx00*dy00*dz01 + zj*dx00*dy00*dz00);
	  }
	
	// -- transform --
	for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++)
	  for(int jbasis = 0; jbasis < o.size_basis_ish(jsh); jbasis++) {
	    dcomplex cumsum(0.0);
	    for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++)
	      for(int jprim = 0; jprim < o.size_prim_ish(jsh); jprim++) {
		cumsum += (this->coef_ish_ibasis_iprim[ish][ibasis][iprim] * 
			   o.coef_ish_ibasis_iprim[jsh][jbasis][jprim] * 
			   z_prim[iprim * o.size_prim_ish(jsh) + jprim]);
	      }
	    int i = this->offset_ish[ish] + ibasis;
	    int j = o.offset_ish[jsh] + jbasis;
	    int idx = i+j*this->size_basis();	
	    Z[idx] = cumsum;
	  }
      }
    MatrixSet matrix_set(this->size_basis(), o.size_basis());
    matrix_set.set("z", Z);
    return matrix_set;
  }

  
  void GTOs::AtR_Ylm(int L, int M, dcomplex* rs, int num_r, dcomplex* cs,
		     dcomplex* res) {

    int lmax(0);
    for(int ish = 0; ish < this->size_sh(); ish++)
      lmax = (lmax < this->l_ish[ish]) ? this->l_ish[ish] : lmax;
    dcomplex* ylm = new dcomplex[num_lm_pair(L)];
    dcomplex* il  = new dcomplex[L+1];
    double eps(0.0000001);

    for(int i = 0; i < num_r; i++)
      res[i] = 0.0;

    for(int ish = 0; ish < this->size_sh(); ish++) {
      dcomplex x = x_ish[ish];
      dcomplex y = y_ish[ish];
      dcomplex z = z_ish[ish];
      dcomplex a2 = dist2(x, y, z);
      dcomplex xxyy = dist2(x, y, 0.0);
      dcomplex a = sqrt(a2);
      dcomplex theta = acos(z / a);
      dcomplex phi   = acos(x / sqrt(xxyy));
      dcomplex zeta = zeta_ish[ish];

      if(abs(a2) < eps) {
	// ---- on center ----
	for(int ibasis = 0; ibasis < size_prim_ish(ish); ibasis++) {
	  int idx = offset_ish[ish] + ibasis;
	  dcomplex coef = coef_ylm_ish[ish] * cs[idx];
	  if(l_ish[ish] == L && m_ish_ibasis[ish][ibasis] == M) {
	    for(int i = 0; i < num_r; i++) {
	      dcomplex r(rs[i]);
	      res[i] += coef * pow(r, L+1) * exp(-zeta * r * r);
	    }
	  }
	}
      } else {
	// ---- no center ----
	if(l_ish[ish] != 0) {
	  std::string msg; SUB_LOCATION(msg);
	  msg += "non center non (l,m)=(0,0) is not implemented.";
	  throw std::runtime_error(msg);
	} else {
	  // --- (l,m) == (0,0)
	  RealSphericalHarmonics(theta, phi, L, ylm);
	  for(int ibasis = 0; ibasis < size_basis_ish(ish); ibasis++) {
	    int idx = offset_ish[ish] + ibasis;	    
	    dcomplex coef = coef_ylm_ish[ish] * cs[idx];
	    for(int i = 0; i < num_r; i++) {
	      dcomplex r(rs[i]);
	      ModSphericalBessel(2.0*zeta*a*r, L, il);
	      res[i] += coef * (4.0*M_PI * exp(-zeta*(r*r+a2)) * il[L] *
			 pow(-1.0, M) * ylm[lm_index(L, -M)]);
	    }
	  }
	}
      }
    }
  }

  void show_complex(dcomplex x, int n) {
    double eps(pow(10.0, -10.0));
    if(abs(x.imag()) < eps) 
      cout << setw(2*n) << right << x.real();
    else if(abs(x.real()) < eps) 
      cout << setw(2*n) << right << x.imag() << "j";
    else if(x.imag()>0.0) 
      cout << setw(n) << right << x.real() << "+" << x.imag() << "j";
    else if(x.imag()<0.0)
      cout << setw(n) << right << x.real() << "-" << -x.imag() << "j";
  }
  void GTOs::Show() const {
    cout << " ==== GTOs ====" << endl;
    cout << "ish zeta x   y   z   L M" << endl;
    for(int ish = 0; ish < this->size_sh(); ish++) {
      
      for(int ibasis = 0; ibasis < this->size_basis_ish(ish); ibasis++) {
	cout << ish << " " << zeta_ish[ish] << " " << x_ish[ish] << " ";
	cout << y_ish[ish] << " " << z_ish[ish] << " " << l_ish[ish] << " ";
	cout << m_ish_ibasis[ish][ibasis];
	for(int iprim = 0; iprim < this->size_prim_ish(ish); iprim++)
	  cout << coef_ish_ibasis_iprim[ish][ibasis][iprim];
	cout << endl;
      }
    }
  }

}
