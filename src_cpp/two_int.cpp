#include <iostream>
#include <stdexcept>
#include "../utils/typedef.hpp"
#include "two_int.hpp"

using namespace std;

namespace cbasis {

  // ==== Definition ====
  typedef vector<SubSymGTOs>::iterator SubIt;
  typedef vector<Reduction>::iterator RdsIt;
  typedef vector<SubSymGTOs>::const_iterator cSubIt;
  typedef vector<Reduction>::const_iterator cRdsIt;
  typedef MultArray<dcomplex, 1> A1dc;
  typedef MultArray<dcomplex, 2> A2dc;
  typedef MultArray<dcomplex, 3> A3dc;
  typedef MultArray<dcomplex, 4> A4dc;

  // ==== Slow routines ====
  dcomplex ERIEle(CartGTO& i, CartGTO& j, CartGTO& k, CartGTO& l) {
    dcomplex zetaP = i.zeta + j.zeta;
    dcomplex zetaPp= k.zeta + l.zeta;
    dcomplex wPx = (i.zeta*i.x + j.zeta*j.x)/zetaP;
    dcomplex wPy = (i.zeta*i.y + j.zeta*j.y)/zetaP;
    dcomplex wPz = (i.zeta*i.z + j.zeta*j.z)/zetaP;
    dcomplex wPpx= (k.zeta*k.x + l.zeta*l.x)/zetaPp;
    dcomplex wPpy= (k.zeta*k.y + l.zeta*l.y)/zetaPp;
    dcomplex wPpz= (k.zeta*k.z + l.zeta*l.z)/zetaPp;
    dcomplex d2 = dist2(i.x-j.x, i.y-j.y, i.z-j.z);
    dcomplex d2p= dist2(k.x-l.x, k.y-l.y, k.z-l.z);

    static const int num(1000);
    MultArray<dcomplex, 3> dx(num),  dy(num),  dz(num);
    MultArray<dcomplex, 3> dxp(num), dyp(num), dzp(num);
    calc_d_coef(i.nx, j.nx, i.nx+j.nx, zetaP,  wPx,  i.x, j.x, dx);
    calc_d_coef(i.ny, j.ny, i.ny+j.ny, zetaP,  wPy,  i.y, j.y, dy);
    calc_d_coef(i.nz, j.nz, i.nz+j.nz, zetaP,  wPz,  i.z, j.z, dz);
    calc_d_coef(k.nx, l.nx, k.nx+l.nx, zetaPp, wPpx, k.x, l.x, dxp);
    calc_d_coef(k.ny, l.ny, k.ny+l.ny, zetaPp, wPpy, k.y, l.y, dyp);
    calc_d_coef(k.nz, l.nz, k.nz+l.nz, zetaPp, wPpz, k.z, l.z, dzp);
    
    double delta(0.0000000000001);
    dcomplex arg = zetaP * zetaPp / (zetaP + zetaPp) * dist2(wPx-wPpx, wPy-wPpy, wPz-wPpz);
    int mm = i.nx+j.nx+k.nx+l.nx+i.ny+j.ny+k.ny+l.ny+i.nz+j.nz+k.nz+l.nz;
    dcomplex Fj_or_Gj[num];
    dcomplex c(0);
    if(real(arg)+delta > 0.0) {
      IncompleteGamma(mm, arg, Fj_or_Gj);
      c = exp(-i.zeta*j.zeta/zetaP*d2 -k.zeta*l.zeta/zetaPp*d2p);
    } else {
      ExpIncompleteGamma(mm, -arg, Fj_or_Gj);
      c = exp(-i.zeta*j.zeta/zetaP*d2 -k.zeta*l.zeta/zetaP*d2p -arg);
    }

    MultArray<dcomplex, 3> Rrs(num);
    Rrs.SetRange(0, mm, 0, mm, 0, mm);
    for(int nx = 0; nx <= mm; nx++)
      for(int ny = 0; ny <= mm; ny++)
	for(int nz = 0; nz <= mm; nz++)
	  if(nx + ny + nz <= mm) {
	    Rrs(nx, ny, nz) = c * coef_R(zetaP*zetaPp/(zetaP+zetaPp),
					 wPx, wPy, wPz, wPpx, wPpy, wPpz,
					 nx, ny, nz, 0, Fj_or_Gj);
	  }
    dcomplex cumsum(0);
    for(int Nx  = 0; Nx  <= i.nx + j.nx; Nx++)
    for(int Nxp = 0; Nxp <= k.nx + l.nx; Nxp++)
    for(int Ny  = 0; Ny  <= i.ny + j.ny; Ny++)
    for(int Nyp = 0; Nyp <= k.ny + l.ny; Nyp++)
    for(int Nz  = 0; Nz  <= i.nz + j.nz; Nz++)
    for(int Nzp = 0; Nzp <= k.nz + l.nz; Nzp++) {
      dcomplex r0;
      r0 = Rrs(Nx+Nxp, Ny+Nyp, Nz+Nzp);
      cumsum += (dx(i.nx, j.nx, Nx) * dxp(k.nx, l.nx, Nxp) *
		 dy(i.ny, j.ny, Ny) * dyp(k.ny, l.ny, Nyp) *
		 dz(i.nz, j.nz, Nz) * dzp(k.nz, l.nz, Nzp) *
		 r0 * pow(-1.0, Nxp+Nyp+Nzp));
    }
    
    dcomplex lambda = 2.0*pow(M_PI, 2.5)/(zetaP * zetaPp * sqrt(zetaP + zetaPp));
    return lambda * cumsum;
  }
  
  // ==== R_coef ====
  void calc_R_coef_eri0(dcomplex zarg,
		       dcomplex wPx,  dcomplex wPy,  dcomplex wPz,
		       dcomplex wPpx, dcomplex wPpy, dcomplex wPpz,
		       int max_n, dcomplex *Fjs, dcomplex mult_coef,
		       A3dc& res) {
    res.SetRange(0, max_n, 0, max_n, 0, max_n);
    for(int nx = 0; nx <= max_n; nx++)
      for(int ny = 0; ny <= max_n; ny++)
	for(int nz = 0; nz <= max_n; nz++) {
	  if(nx + ny + nz <= max_n) {
	    dcomplex v = coef_R(zarg,
			      wPx, wPy, wPz,
			      wPpx, wPpy, wPpz,
			      nx, ny, nz, 0, Fjs);
	    res(nx, ny, nz) = mult_coef * v;
	  }
	}

  }
  MultArray<dcomplex, 4> coef_R_map_val(1000);
  MultArray<bool,     4> coef_R_map_has(1000);
  dcomplex coef_R_with_memo(dcomplex zetaP,
			    dcomplex wPx, dcomplex wPy, dcomplex wPz,
			    dcomplex cx,  dcomplex cy,  dcomplex cz,
			    int mx, int my, int mz, int j, dcomplex* Fjs) {
    /* Compute function coef_R with memorize;     */

    dcomplex& ref = coef_R_map_val(mx, my, mz, j);
    bool& has = coef_R_map_has(mx, my, mz, j);
    if(has) 
      return ref;

    if(mx == 0 && my == 0 && mz == 0) {
      ref = pow(-2.0*zetaP, j) * Fjs[j];
      has = true;
      return ref;
    }
    
    if(mx > 0) {
      ref = 0.0;
      if(mx > 1) 
	ref += (mx-1.0) * coef_R_with_memo(zetaP, wPx, wPy, wPz, cx, cy, cz,
					    mx-2, my, mz, j+1, Fjs);
      ref += (wPx-cx) * coef_R_with_memo(zetaP, wPx, wPy, wPz, cx, cy, cz,
					  mx-1, my, mz, j+1, Fjs);
      has = true;
      return ref;
    }
    if(my > 0) {
      ref = 0.0;
      if(my > 1) 
	ref += (my-1.0) * coef_R_with_memo(zetaP, wPx, wPy, wPz, cx, cy, cz,
					    mx, my-2, mz, j+1, Fjs);
      ref += (wPy-cy) *   coef_R_with_memo(zetaP, wPx, wPy, wPz, cx, cy, cz,
					    mx, my-1, mz, j+1, Fjs);
      has = true;
      return ref;
    }
    if(mz > 0) {
      ref = 0.0;
      if(mz > 1) 
	ref += (mz-1.0) * coef_R_with_memo(zetaP, wPx, wPy, wPz, cx, cy, cz,
					    mx, my, mz-2, j+1, Fjs);
      ref += (wPz-cz) *   coef_R_with_memo(zetaP, wPx, wPy, wPz, cx, cy, cz,
					    mx, my, mz-1, j+1, Fjs);
      has = true;
      return ref;
    }

    std::string msg;
    SUB_LOCATION(msg);
    msg += " one of mx, my, mz is negative integer.";
    throw std::runtime_error(msg);    
    return 0.0;    

  }
  void calc_R_coef_eri1(dcomplex zarg,
			dcomplex wPx, dcomplex wPy, dcomplex wPz,
			dcomplex wPpx,dcomplex wPpy,dcomplex wPpz,
			int max_n, dcomplex *Fjs, dcomplex mult_coef, A3dc& res) {

    coef_R_map_val.SetRange(0, max_n, 0, max_n, 0, max_n, 0, max_n);
    coef_R_map_has.SetRange(0, max_n, 0, max_n, 0, max_n, 0, max_n);

    coef_R_map_has.SetValue(false);
      
    for(int nx = 0; nx <= max_n; nx++)
      for(int ny = 0; ny <= max_n; ny++)
	for(int nz = 0; nz <= max_n; nz++) {
	  if(nx + ny + nz <= max_n) {
	    dcomplex v = coef_R_with_memo(zarg,
					  wPx, wPy, wPz,
					  wPpx, wPpy, wPpz,
					  nx, ny, nz, 0, Fjs);
	    res(nx, ny, nz) = mult_coef * v;
	  }
	}
  }
  void coef_R_eri_switch(dcomplex zarg,
			 dcomplex wPx, dcomplex wPy, dcomplex wPz,
			 dcomplex wPpx,dcomplex wPpy,dcomplex wPpz,
			 int max_n, dcomplex *Fjs, dcomplex mult_coef, A3dc& res,
			 ERIMethod method) {

    res.SetRange(0, max_n, 0, max_n, 0, max_n);


    if(method.coef_R_memo == 0) {
      calc_R_coef_eri0(zarg, wPx, wPy, wPz,
		       wPpx, wPpy, wPpz, max_n, Fjs, mult_coef, res);
    } else { 
      calc_R_coef_eri1(zarg, wPx, wPy, wPz,
		       wPpx, wPpy, wPpz, max_n, Fjs, mult_coef, res);
    }

  }

  // ==== SymGTOs ====
  bool ExistNon0(SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub) {
    SymmetryGroup sym = isub->sym_group();
    bool non0(false);
    for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds)
    for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end(); ++jrds)
    for(RdsIt krds = ksub->rds.begin(); krds != ksub->rds.end(); ++krds)
    for(RdsIt lrds = lsub->rds.begin(); lrds != lsub->rds.end(); ++lrds) {
      if(sym->Non0_4(irds->irrep, jrds->irrep, krds->irrep, lrds->irrep))
	non0 = true;
    }
    return non0;
  }

  struct ERI_buf {
    A3dc dx, dy, dz, dxp, dyp, dzp;
    A1dc Fjs;
    A3dc Rrs;
    //    dcomplex eij, ekl, lambda;
    dcomplex lambda;
    ERI_buf(int n) :
      dx(n), dy(n), dz(n), dxp(n), dyp(n), dzp(n), Fjs(n), Rrs(n) {}
  };
  void CalcCoef(dcomplex xi, dcomplex yi, dcomplex zi,
	       int mi, dcomplex zetai,
	       dcomplex xj, dcomplex yj, dcomplex zj,
	       int mj, dcomplex& zetaj,
	       dcomplex xk, dcomplex yk, dcomplex zk,
	       int mk, dcomplex& zetak,
	       dcomplex xl, dcomplex yl, dcomplex zl,
		int ml, dcomplex zetal, ERI_buf& buf, ERIMethod method) {
    dcomplex zetaP = zetai + zetaj;
    dcomplex zetaPp= zetak + zetal;
    dcomplex wPx((zetai*xi + zetaj*xj)/zetaP);
    dcomplex wPy((zetai*yi + zetaj*yj)/zetaP);
    dcomplex wPz((zetai*zi + zetaj*zj)/zetaP);
    
    calc_d_coef(mi, mj, mi+mj, zetaP,  wPx,  xi, xj, buf.dx);
    calc_d_coef(mi, mj, mi+mj, zetaP,  wPy,  yi, yj, buf.dy);
    calc_d_coef(mi, mj, mi+mj, zetaP,  wPz,  zi, zj, buf.dz);

    dcomplex wPpx((zetak*xk + zetal*xl)/zetaPp);
    dcomplex wPpy((zetak*yk + zetal*yl)/zetaPp);
    dcomplex wPpz((zetak*zk + zetal*zl)/zetaPp);
    
    calc_d_coef(mk, ml, mk+ml, zetaPp, wPpx, xk, xl, buf.dxp);
    calc_d_coef(mk, ml, mk+ml, zetaPp, wPpy, yk, yl, buf.dyp);
    calc_d_coef(mk, ml, mk+ml, zetaPp, wPpz, zk, zl, buf.dzp);

    dcomplex zarg(zetaP * zetaPp / (zetaP + zetaPp));
    dcomplex argIncGamma(zarg * dist2(wPx-wPpx, wPy-wPpy, wPz-wPpz));
    double delta(0.0000000000001);
    if(real(argIncGamma)+delta > 0.0) {
      IncompleteGamma(mi+mj+mk+ml, argIncGamma, &buf.Fjs(0));      
      int mm = mi + mj + mk + ml;
      dcomplex eij = exp(-zetai * zetaj / zetaP *  dist2(xi-xj, yi-yj, zi-zj));
      dcomplex ekl = (exp(-zetak * zetal / zetaPp * dist2(xk-xl, yk-yl, zk-zl)));
      coef_R_eri_switch(zarg, wPx, wPy, wPz, wPpx, wPpy, wPpz, mm,
			&buf.Fjs(0), eij * ekl, buf.Rrs, method);
    } else {      
      int mm = mi + mj + mk + ml;
      ExpIncompleteGamma(mm, -argIncGamma, &buf.Fjs(0)); 
      dcomplex arg_other = 
	-zetai * zetaj / zetaP *  dist2(xi-xj, yi-yj, zi-zj)
	-zetak * zetal / zetaPp * dist2(xk-xl, yk-yl, zk-zl)
	-argIncGamma;
      coef_R_eri_switch(zarg, wPx, wPy, wPz, wPpx, wPpy, wPpz, mm,
			&buf.Fjs(0), exp(arg_other), buf.Rrs, method);
    }

  }
  dcomplex CalcPrimOne(int nxi, int nyi, int nzi,
		   int nxj, int nyj, int nzj,
		   int nxk, int nyk, int nzk,
		   int nxl, int nyl, int nzl, ERI_buf& buf) {
    dcomplex cumsum(0);
    for(int Nx  = 0; Nx  <= nxi + nxj; Nx++)
      for(int Nxp = 0; Nxp <= nxk + nxl; Nxp++)
	for(int Ny  = 0; Ny  <= nyi + nyj; Ny++)
	  for(int Nyp = 0; Nyp <= nyk + nyl; Nyp++)
	    for(int Nz  = 0; Nz  <= nzi + nzj; Nz++)
	      for(int Nzp = 0; Nzp <= nzk + nzl; Nzp++) {
		dcomplex r0;
		r0 = buf.Rrs(Nx+Nxp, Ny+Nyp, Nz+Nzp);
		cumsum += (buf.dx(nxi, nxj, Nx) * buf.dxp(nxk, nxl, Nxp) *
			   buf.dy(nyi, nyj, Ny) * buf.dyp(nyk, nyl, Nyp) *
			   buf.dz(nzi, nzj, Nz) * buf.dzp(nzk, nzl, Nzp) *
			   r0 * pow(-1.0, Nxp+Nyp+Nzp));
	      }
    return buf.lambda * cumsum;;
	
  }
  dcomplex OneERI(dcomplex xi, dcomplex yi, dcomplex zi,
		  int nxi, int nyi, int nzi, dcomplex zetai,
		  dcomplex xj, dcomplex yj, dcomplex zj,
		  int nxj, int nyj, int nzj, dcomplex& zetaj,
		  dcomplex xk, dcomplex yk, dcomplex zk,
		  int nxk, int nyk, int nzk, dcomplex& zetak,
		  dcomplex xl, dcomplex yl, dcomplex zl,
		  int nxl, int nyl, int nzl, dcomplex zetal,
		  ERIMethod method) {
    
    static dcomplex Fjs[100];
    static A3dc Rrs(100);
    static int n(1000);
    static A3dc dxmap(n), dymap(n), dzmap(n), dxmap_p(n), dymap_p(n), dzmap_p(n);    

    dcomplex zetaP = zetai + zetaj; dcomplex zetaPp= zetak + zetal;
    dcomplex lambda = 2.0*pow(M_PI, 2.5)/(zetaP * zetaPp * sqrt(zetaP + zetaPp));

    dcomplex wPx((zetai*xi + zetaj*xj)/zetaP);
    dcomplex wPy((zetai*yi + zetaj*yj)/zetaP);
    dcomplex wPz((zetai*zi + zetaj*zj)/zetaP);
    dcomplex eij(exp(-zetai * zetaj / zetaP *  dist2(xi-xj, yi-yj, zi-zj)));
    int mi = nxi + nyi + nzi;
    int mj = nxj + nyj + nzj;
    int mk = nxk + nyk + nzk;
    int ml = nxl + nyl + nzl;

    calc_d_coef(mi, mj, mi+mj, zetaP,  wPx,  xi, xj, dxmap);
    calc_d_coef(mi, mj, mi+mj, zetaP,  wPy,  yi, yj, dymap);
    calc_d_coef(mi, mj, mi+mj, zetaP,  wPz,  zi, zj, dzmap);

    dcomplex wPpx((zetak*xk + zetal*xl)/zetaPp);
    dcomplex wPpy((zetak*yk + zetal*yl)/zetaPp);
    dcomplex wPpz((zetak*zk + zetal*zl)/zetaPp);
    dcomplex ekl(exp(-zetak * zetal / zetaPp * dist2(xk-xl, yk-yl, zk-zl)));
    calc_d_coef(mk, ml, mk+ml, zetaPp, wPpx, xk, xl, dxmap_p);
    calc_d_coef(mk, ml, mk+ml, zetaPp, wPpy, yk, yl, dymap_p);
    calc_d_coef(mk, ml, mk+ml, zetaPp, wPpz, zk, zl, dzmap_p);
      
    dcomplex zarg(zetaP * zetaPp / (zetaP + zetaPp));
    dcomplex argIncGamma(zarg * dist2(wPx-wPpx, wPy-wPpy, wPz-wPpz));
    IncompleteGamma(mi+mj+mk+ml, argIncGamma, Fjs);      

    int mm = mi + mj + mk + ml;
    coef_R_eri_switch(zarg, wPx, wPy, wPz, wPpx, wPpy, wPpz, mm, Fjs, 1.0, Rrs,
		      method);
      

    // -- determine 4 primitive GTOs for integration (ij|kl) --
	  
    dcomplex cumsum(0);
    for(int Nx  = 0; Nx  <= nxi + nxj; Nx++)
    for(int Nxp = 0; Nxp <= nxk + nxl; Nxp++)
    for(int Ny  = 0; Ny  <= nyi + nyj; Ny++)
    for(int Nyp = 0; Nyp <= nyk + nyl; Nyp++)
    for(int Nz  = 0; Nz  <= nzi + nzj; Nz++)
    for(int Nzp = 0; Nzp <= nzk + nzl; Nzp++) {
      dcomplex r0;
      r0 = Rrs(Nx+Nxp, Ny+Nyp, Nz+Nzp);
      cumsum += (dxmap(nxi, nxj, Nx) * dxmap_p(nxk, nxl, Nxp) *
		 dymap(nyi, nyj, Ny) * dymap_p(nyk, nyl, Nyp) *
		 dzmap(nzi, nzj, Nz) * dzmap_p(nzk, nzl, Nzp) *
		 r0 * pow(-1.0, Nxp+Nyp+Nzp));
    }
    return eij * ekl * lambda * cumsum;
  }
  int one_dim(int ip, int ni, int jp, int nj, int kp, int nk, int lp) {
    return ip + jp * ni + kp * ni * nj + lp * ni * nj * nk;
  }
  void CheckEqERI(SymmetryGroup sym_group, 
		  SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		  int ip, int jp, int kp, int lp,
		  int ni, int nj, int nk,
		  int *mark_I, bool *is_zero, bool *is_youngest) {
    *is_zero = false;
    *is_youngest = true;;
    mark_I[0] = 1;
    for(int I = 1; I < sym_group->num_class(); I++) {
      int ipt = isub->ip_jg_kp(I, ip);
      int jpt = jsub->ip_jg_kp(I, jp);
      int kpt = ksub->ip_jg_kp(I, kp);
      int lpt = lsub->ip_jg_kp(I, lp);
      int sig_ipt = isub->sign_ip_jg_kp(I, ip);
      int sig_jpt = jsub->sign_ip_jg_kp(I, jp);
      int sig_kpt = ksub->sign_ip_jg_kp(I, kp);
      int sig_lpt = lsub->sign_ip_jg_kp(I, lp);
      int sig = sig_ipt * sig_jpt * sig_kpt * sig_lpt;

      if(sig == 0) {
	mark_I[I] = 0;
	
      } else if(ipt == ip && jpt == jp && kp == kpt && lp == lpt) {

	if(sig == 1) {
	  mark_I[I] = 0;
	} else {
	  *is_zero = true;
	  return;
	}

      } else {
	if(one_dim(ip, ni, jp, nj, kp, nk, lp) >
	   one_dim(ipt,ni, jpt,nj, kpt,nk, lpt)) {
	  *is_youngest = false;
	}
	mark_I[I] = sig;
      }
    }

  }

  // ==== Primitive ====
  // -- very simple --
  void CalcPrimERI0(SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		    dcomplex zetai, dcomplex zetaj, dcomplex zetak, dcomplex zetal,
		    A4dc& prim, ERIMethod method) {

    static ERI_buf buf(1000);
    dcomplex zetaP = zetai + zetaj; dcomplex zetaPp = zetak + zetal;

    int nati, natj, natk, natl;
    nati = isub->size_at(); natj = jsub->size_at();
    natk = ksub->size_at(); natl = lsub->size_at();

    int npni, npnj, npnk, npnl;
    npni = isub->size_pn(); npnj = jsub->size_pn();
    npnk = ksub->size_pn(); npnl = lsub->size_pn();

    prim.SetRange(0, nati*npni, 0, natj*npnj, 0, natk*npnk, 0, natl*npnl);

    int mi, mj, mk, ml;
    mi = isub->maxn; mj = jsub->maxn; mk = ksub->maxn; ml = lsub->maxn;
    buf.lambda = 2.0*pow(M_PI, 2.5)/(zetaP * zetaPp * sqrt(zetaP + zetaPp));    

    prim.SetValue(0.0);

    for(int iat = 0; iat < nati; iat++) 
    for(int jat = 0; jat < natj; jat++)       
    for(int kat = 0; kat < natk; kat++) 
    for(int lat = 0; lat < natl; lat++) {
      CalcCoef(isub->x(iat), isub->y(iat), isub->z(iat), mi, zetai,
	       jsub->x(jat), jsub->y(jat), jsub->z(jat), mj, zetaj, 
	       ksub->x(kat), ksub->y(kat), ksub->z(kat), mk, zetak, 
	       lsub->x(lat), lsub->y(lat), lsub->z(lat), ml, zetal, buf, method);
      
      for(int ipn = 0; ipn < npni; ipn++) 
      for(int jpn = 0; jpn < npnj; jpn++) 
      for(int kpn = 0; kpn < npnk; kpn++) 
      for(int lpn = 0; lpn < npnl; lpn++) {
	int ip = isub->ip_iat_ipn(iat, ipn); int jp = jsub->ip_iat_ipn(jat, jpn);
	int kp = ksub->ip_iat_ipn(kat, kpn); int lp = lsub->ip_iat_ipn(lat, lpn);
	dcomplex v;
	v = CalcPrimOne(isub->nx(ipn), isub->ny(ipn), isub->nz(ipn),
			jsub->nx(jpn), jsub->ny(jpn), jsub->nz(jpn),
			ksub->nx(kpn), ksub->ny(kpn), ksub->nz(kpn),
			lsub->nx(lpn), lsub->ny(lpn), lsub->nz(lpn), buf);
	prim(ip, jp, kp, lp) = v;
      }
    }
  }
  // -- Symmetry considerration --
  void CalcPrimERI1(SymmetryGroup sym,
		    SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		    dcomplex zetai, dcomplex zetaj, dcomplex zetak, dcomplex zetal,
		    A4dc& prim, ERIMethod method) {

    static ERI_buf buf(1000);
    dcomplex zetaP = zetai + zetaj; dcomplex zetaPp = zetak + zetal;

    int nati, natj, natk, natl;
    nati = isub->size_at(); natj = jsub->size_at();
    natk = ksub->size_at(); natl = lsub->size_at();

    int npni, npnj, npnk, npnl;
    npni = isub->size_pn(); npnj = jsub->size_pn();
    npnk = ksub->size_pn(); npnl = lsub->size_pn();

    prim.SetRange(0, nati*npni, 0, natj*npnj, 0, natk*npnk, 0, natl*npnl);

    int mi, mj, mk, ml;
    mi = isub->maxn; mj = jsub->maxn; mk = ksub->maxn; ml = lsub->maxn;
    buf.lambda = 2.0*pow(M_PI, 2.5)/(zetaP * zetaPp * sqrt(zetaP + zetaPp));    

    int numI(sym->order());

    prim.SetValue(0.0);

    for(int iat = 0; iat < nati; iat++) 
    for(int jat = 0; jat < natj; jat++)       
    for(int kat = 0; kat < natk; kat++) 
    for(int lat = 0; lat < natl; lat++) {

      // -- check 0 or non 0 --
      bool find_non0(false);
      for(int ipn = 0; ipn < npni; ipn++) 
      for(int jpn = 0; jpn < npnj; jpn++) 
      for(int kpn = 0; kpn < npnk; kpn++) 
      for(int lpn = 0; lpn < npnl; lpn++) {
	int ip = isub->ip_iat_ipn(iat, ipn); int jp = jsub->ip_iat_ipn(jat, jpn);
	int kp = ksub->ip_iat_ipn(kat, kpn); int lp = lsub->ip_iat_ipn(lat, lpn);
	int mark_I[10]; bool is_zero, is_youngest;
	CheckEqERI(isub->sym_group(), isub, jsub, ksub, lsub, ip, jp, kp, lp,
		   nati*npni, natj*npnj, natk*npnk, mark_I, &is_zero, &is_youngest);
	if(!is_zero && is_youngest)
	  find_non0 = true;
      }

      // -- compute if found non0 --
      if(find_non0) {
      
	CalcCoef(isub->x(iat), isub->y(iat), isub->z(iat), mi, zetai,
		 jsub->x(jat), jsub->y(jat), jsub->z(jat), mj, zetaj, 
		 ksub->x(kat), ksub->y(kat), ksub->z(kat), mk, zetak, 
		 lsub->x(lat), lsub->y(lat), lsub->z(lat), ml, zetal, buf, method);

	for(int ipn = 0; ipn < npni; ipn++) 
	for(int jpn = 0; jpn < npnj; jpn++) 
	for(int kpn = 0; kpn < npnk; kpn++) 
	for(int lpn = 0; lpn < npnl; lpn++) {
	  int ip = isub->ip_iat_ipn(iat, ipn); int jp = jsub->ip_iat_ipn(jat, jpn);
	  int kp = ksub->ip_iat_ipn(kat, kpn); int lp = lsub->ip_iat_ipn(lat, lpn);

	  int mark_I[10]; bool is_zero, is_youngest;
	  CheckEqERI(isub->sym_group(), isub, jsub, ksub, lsub, ip, jp, kp, lp,
		     nati*npni, natj*npnj, natk*npnk, mark_I, &is_zero, &is_youngest);
	  if(!is_zero && is_youngest) {
	    dcomplex v;
	    v = CalcPrimOne(isub->nx(ipn), isub->ny(ipn), isub->nz(ipn),
			    jsub->nx(jpn), jsub->ny(jpn), jsub->nz(jpn),
			    ksub->nx(kpn), ksub->ny(kpn), ksub->nz(kpn),
			    lsub->nx(lpn), lsub->ny(lpn), lsub->nz(lpn), buf);
	    for(int I = 0; I < numI; I++) {
	      if(mark_I[I] != 0) {
		int ipt = isub->ip_jg_kp(I, ip);
		int jpt = jsub->ip_jg_kp(I, jp);
		int kpt = ksub->ip_jg_kp(I, kp);
		int lpt = lsub->ip_jg_kp(I, lp);
		prim(ipt, jpt, kpt, lpt) = dcomplex(mark_I[I]) * v;
	      }
	    }
	  }
	}
      }
    }

  }
  // -- Simple --
  void CalcPrimERI2(SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		    int iz, int jz, int kz, int lz,
		    A4dc& prim, ERIMethod method) {
    static A3dc dxmap(1000);
    static A3dc dymap(1000);
    static A3dc dzmap(1000);
    static A3dc dxmap_p(1000);
    static A3dc dymap_p(1000);
    static A3dc dzmap_p(1000);

    dcomplex zetai, zetaj, zetak, zetal, zetaP, zetaPp;
    zetai = isub->zeta_iz[iz]; zetaj = jsub->zeta_iz[jz];
    zetak = ksub->zeta_iz[kz]; zetal = lsub->zeta_iz[lz];
    zetaP = zetai + zetaj;     zetaPp= zetak + zetal;

    int nati, natj, natk, natl;
    nati = isub->size_at(); natj = jsub->size_at();
    natk = ksub->size_at(); natl = lsub->size_at();

    int npni, npnj, npnk, npnl;
    npni = isub->size_pn(); npnj = jsub->size_pn();
    npnk = ksub->size_pn(); npnl = lsub->size_pn();

    int mi, mj, mk, ml;
    mi = isub->maxn; mj = jsub->maxn; mk = ksub->maxn; ml = lsub->maxn;
    dcomplex lambda = 2.0*pow(M_PI, 2.5)/(zetaP * zetaPp * sqrt(zetaP + zetaPp));    

    prim.SetRange(0, nati*npni, 0, natj*npnj, 0, natk*npnk, 0, natl*npnl);

    static dcomplex Fjs[100];
    static A3dc Rrs(100);
    int idx(0);
    for(int iat = 0; iat < nati; iat++) {
    for(int jat = 0; jat < natj; jat++) {      
    for(int kat = 0; kat < natk; kat++) {
    for(int lat = 0; lat < natl; lat++) {      
      dcomplex xi(isub->x(iat)), yi(isub->y(iat)), zi(isub->z(iat));
      dcomplex xj(jsub->x(jat)), yj(jsub->y(jat)), zj(jsub->z(jat));
      dcomplex wPx((zetai*xi + zetaj*xj)/zetaP);
      dcomplex wPy((zetai*yi + zetaj*yj)/zetaP);
      dcomplex wPz((zetai*zi + zetaj*zj)/zetaP);
      dcomplex eij(exp(-zetai * zetaj / zetaP *  dist2(xi-xj, yi-yj, zi-zj)));
      calc_d_coef(mi, mj, mi+mj, zetaP,  wPx,  xi, xj, dxmap);
      calc_d_coef(mi, mj, mi+mj, zetaP,  wPy,  yi, yj, dymap);
      calc_d_coef(mi, mj, mi+mj, zetaP,  wPz,  zi, zj, dzmap);

      dcomplex xk(ksub->x(kat)), yk(ksub->y(kat)), zk(ksub->z(kat));
      dcomplex xl(lsub->x(lat)), yl(lsub->y(lat)), zl(lsub->z(lat));
      dcomplex wPpx((zetak*xk + zetal*xl)/zetaPp);
      dcomplex wPpy((zetak*yk + zetal*yl)/zetaPp);
      dcomplex wPpz((zetak*zk + zetal*zl)/zetaPp);
      dcomplex ekl(exp(-zetak * zetal / zetaPp * dist2(xk-xl, yk-yl, zk-zl)));
	calc_d_coef(mk, ml, mk+ml, zetaPp, wPpx, xk, xl, dxmap_p);
	calc_d_coef(mk, ml, mk+ml, zetaPp, wPpy, yk, yl, dymap_p);
	calc_d_coef(mk, ml, mk+ml, zetaPp, wPpz, zk, zl, dzmap_p);

	dcomplex zarg(zetaP * zetaPp / (zetaP + zetaPp));
	dcomplex argIncGamma(zarg * dist2(wPx-wPpx, wPy-wPpy, wPz-wPpz));
	IncompleteGamma(mi+mj+mk+ml, argIncGamma, Fjs);      

	int mm = mi + mj + mk + ml;
	coef_R_eri_switch(zarg, wPx, wPy, wPz, wPpx, wPpy, wPpz, mm, Fjs, 1.0, Rrs, method);


	for(int ipn = 0; ipn < npni; ipn++) 
	for(int jpn = 0; jpn < npnj; jpn++) 
	for(int kpn = 0; kpn < npnk; kpn++) 
	for(int lpn = 0; lpn < npnl; lpn++) {      
	  int ip = isub->ip_iat_ipn(iat, ipn);
	  int jp = jsub->ip_iat_ipn(jat, jpn);
	  int kp = ksub->ip_iat_ipn(kat, kpn);
	  int lp = lsub->ip_iat_ipn(lat, lpn);

	  // -- determine 4 primitive GTOs for integration (ij|kl) --
	  int nxi, nxj, nxk, nxl, nyi, nyj, nyk, nyl, nzi, nzj, nzk, nzl;
	  nxi = isub->nx(ipn); nyi = isub->ny(ipn); nzi = isub->nz(ipn);
	  nxj = jsub->nx(jpn); nyj = jsub->ny(jpn); nzj = jsub->nz(jpn);
	  nxk = ksub->nx(kpn); nyk = ksub->ny(kpn); nzk = ksub->nz(kpn);
	  nxl = lsub->nx(lpn); nyl = lsub->ny(lpn); nzl = lsub->nz(lpn);
	  
	  dcomplex cumsum(0);
	  for(int Nx  = 0; Nx  <= nxi + nxj; Nx++)
	  for(int Nxp = 0; Nxp <= nxk + nxl; Nxp++)
	  for(int Ny  = 0; Ny  <= nyi + nyj; Ny++)
	  for(int Nyp = 0; Nyp <= nyk + nyl; Nyp++)
	  for(int Nz  = 0; Nz  <= nzi + nzj; Nz++)
	  for(int Nzp = 0; Nzp <= nzk + nzl; Nzp++) {
	    dcomplex r0;
	    r0 = Rrs(Nx+Nxp, Ny+Nyp, Nz+Nzp);
	    cumsum += (dxmap(nxi, nxj, Nx) * dxmap_p(nxk, nxl, Nxp) *
		       dymap(nyi, nyj, Ny) * dymap_p(nyk, nyl, Nyp) *
		       dzmap(nzi, nzj, Nz) * dzmap_p(nzk, nzl, Nzp) *
		       r0 * pow(-1.0, Nxp+Nyp+Nzp));
		   
	  }
	  prim(ip, jp, kp, lp) = eij * ekl * lambda * cumsum;
	  ++idx;
	}
      }
    }}}
  }  
  
  // ==== Transform ====
  void CalcTransERI(SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		   int iz, int jz, int kz, int lz,
		    A4dc& prim, B2EInt eri,
		    bool use_perm=false, bool perm_0=false, bool perm_1=false) {

    int nati, natj, natk, natl;
    nati = isub->size_at(); natj = jsub->size_at();
    natk = ksub->size_at(); natl = lsub->size_at();

    int npni, npnj, npnk, npnl;
    npni = isub->size_pn(); npnj = jsub->size_pn();
    npnk = ksub->size_pn(); npnl = lsub->size_pn();

    for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds)
    for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end(); ++jrds)
    for(RdsIt krds = ksub->rds.begin(); krds != ksub->rds.end(); ++krds)
    for(RdsIt lrds = lsub->rds.begin(); lrds != lsub->rds.end(); ++lrds) {
      
      

      if(isub->sym_group()->Non0_4(irds->irrep, jrds->irrep, krds->irrep, lrds->irrep)) {
	dcomplex cz(irds->coef_iz(iz) * 
		    jrds->coef_iz(jz) * 
		    krds->coef_iz(kz) * 
		    lrds->coef_iz(lz));
	dcomplex cumsum(0);
	for(int iat = 0; iat < nati; iat++) {
	for(int ipn = 0; ipn < npni; ipn++) {
        for(int jat = 0; jat < natj; jat++) {
	for(int jpn = 0; jpn < npnj; jpn++) {
        for(int kat = 0; kat < natk; kat++) {
	for(int kpn = 0; kpn < npnk; kpn++) {
        for(int lat = 0; lat < natl; lat++) {
        for(int lpn = 0; lpn < npnl; lpn++) {
	  int ip = isub->ip_iat_ipn(iat, ipn);
	  dcomplex ci(irds->coef_iat_ipn(iat, ipn));
	  int jp = jsub->ip_iat_ipn(jat, jpn);
	  dcomplex cj(jrds->coef_iat_ipn(jat, jpn));
	  int kp = ksub->ip_iat_ipn(kat, kpn);
	  dcomplex ck(krds->coef_iat_ipn(kat, kpn));
	  int lp = lsub->ip_iat_ipn(lat, lpn);
	  dcomplex cc = ci * cj * ck * lrds->coef_iat_ipn(lat, lpn) * cz;
	  cumsum += cc * prim(ip, jp, kp, lp);
	}}}}}}}}
	
	Irrep ir(irds->irrep), jr(jrds->irrep), kr(krds->irrep), lr(lrds->irrep);
	int idx(irds->offset+iz), jdx(jrds->offset+jz), kdx(krds->offset+kz);
	int ldx(lrds->offset+lz);
	
	eri->Set(  ir, jr, kr, lr, idx, jdx, kdx, ldx, cumsum);
	if(use_perm)
	  eri->Set(kr, lr, ir, jr, kdx, ldx, idx, jdx, cumsum);
	if(perm_0)
	  eri->Set(jr, ir, kr, lr, jdx, idx, kdx, ldx, cumsum);
	if(perm_1)
	  eri->Set(ir, jr, lr, kr, idx, jdx, ldx, kdx, cumsum);
      }
    }

  }
  void CalcTransERI0(SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		     int iz, int jz, int kz, int lz,
		     A4dc& prim, B2EInt eri) {

    int nati, natj, natk, natl;
    nati = isub->size_at(); natj = jsub->size_at();
    natk = ksub->size_at(); natl = lsub->size_at();

    int npni, npnj, npnk, npnl;
    npni = isub->size_pn(); npnj = jsub->size_pn();
    npnk = ksub->size_pn(); npnl = lsub->size_pn();
    

    for(RdsIt irds = isub->rds.begin(); irds != isub->rds.end(); ++irds)
    for(RdsIt jrds = jsub->rds.begin(); jrds != jsub->rds.end(); ++jrds)
    for(RdsIt krds = ksub->rds.begin(); krds != ksub->rds.end(); ++krds)
    for(RdsIt lrds = lsub->rds.begin(); lrds != lsub->rds.end(); ++lrds) {
      
      dcomplex cz(irds->coef_iz(iz) * 
		  jrds->coef_iz(jz) * 
		  krds->coef_iz(kz) * 
		  lrds->coef_iz(lz));
      dcomplex cumsum(0);
      int idx(0);
      for(int iat = 0; iat < nati; iat++) {
      for(int ipn = 0; ipn < npni; ipn++) {
      for(int jat = 0; jat < natj; jat++) {
      for(int jpn = 0; jpn < npnj; jpn++) {
      for(int kat = 0; kat < natk; kat++) {
      for(int kpn = 0; kpn < npnk; kpn++) {
      for(int lat = 0; lat < natl; lat++) {
      for(int lpn = 0; lpn < npnl; lpn++) {
	int ip = isub->ip_iat_ipn(iat, ipn);
	dcomplex ci(irds->coef_iat_ipn(iat, ipn));
	int jp = jsub->ip_iat_ipn(jat, jpn);
	dcomplex cj(jrds->coef_iat_ipn(jat, jpn));
	int kp = ksub->ip_iat_ipn(kat, kpn);
	dcomplex ck(krds->coef_iat_ipn(kat, kpn));
	int lp = lsub->ip_iat_ipn(lat, lpn);
	dcomplex cc = ci * cj * ck * lrds->coef_iat_ipn(lat, lpn) * cz;
	cumsum += cc * prim(ip, jp, kp, lp); ++idx;
      }}}}}}}}
      eri->Set(irds->irrep, jrds->irrep, krds->irrep, lrds->irrep,
	       irds->offset+iz, jrds->offset+jz, krds->offset+kz, lrds->offset+lz,
	       cumsum);
    }
  }

  // ==== calc for Sub  ====
  void CalcERI0(SymmetryGroup sym, SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		A4dc& prim, ERIMethod method, B2EInt eri) {
    for(int iz = 0; iz < isub->size_zeta(); ++iz)
    for(int jz = 0; jz < jsub->size_zeta(); ++jz)
    for(int kz = 0; kz < ksub->size_zeta(); ++kz)
    for(int lz = 0; lz < lsub->size_zeta(); ++lz) {
      if(method.symmetry == 0) {
	CalcPrimERI0(isub, jsub, ksub, lsub,
		     isub->zeta_iz[iz], jsub->zeta_iz[jz],
		     ksub->zeta_iz[kz], lsub->zeta_iz[lz], prim, method);
	CalcTransERI0(isub, jsub, ksub, lsub, iz, jz, kz, lz, prim, eri);
      } else {
	CalcPrimERI1(sym, isub, jsub, ksub, lsub,
		     isub->zeta_iz[iz], jsub->zeta_iz[jz],
		     ksub->zeta_iz[kz], lsub->zeta_iz[lz],
		     prim, method);
	CalcTransERI0(isub, jsub, ksub, lsub, iz, jz, kz, lz, prim, eri);
      }
    }
    
  }
  // -- permutation symmetry --
  void CalcERI1(SymGTOs& gi, SymGTOs& gj,SymGTOs& gk,SymGTOs& gl,
		SubIt isub, SubIt jsub, SubIt ksub, SubIt lsub,
		A4dc& prim, ERIMethod method, B2EInt eri) {
    if(!ExistNon0(isub, jsub, ksub, lsub))
      return;

    int n_ij = (distance(gi->subs().begin(), isub) +
		distance(gj->subs().begin(), jsub) * gi->subs().size());
    int n_kl = (distance(gk->subs().begin(), ksub) +
		distance(gl->subs().begin(), lsub) * gk->subs().size());    

    if(isub == ksub && jsub == lsub) {
      for(int iz = 0; iz < isub->size_zeta(); ++iz)
      for(int jz = 0; jz < jsub->size_zeta(); ++jz)
      for(int kz = 0; kz < ksub->size_zeta(); ++kz)
      for(int lz = 0; lz < lsub->size_zeta(); ++lz) {
	int nnnij = iz + isub->size_zeta() * jz;
	int nnnkl = kz + ksub->size_zeta() * lz;
	if(nnnij >= nnnkl) {
	  CalcPrimERI1(gi->sym_group(), isub, jsub, ksub, lsub,
		       isub->zeta_iz[iz], jsub->zeta_iz[jz],
		       ksub->zeta_iz[kz], lsub->zeta_iz[lz],
		       prim, method);
	  if(nnnij == nnnkl)
	    CalcTransERI(isub, jsub, ksub, lsub, kz, lz, iz, jz, prim, eri, false);
	  else
	    CalcTransERI(isub, jsub, ksub, lsub, kz, lz, iz, jz, prim, eri, true);
	}
      }
    } else if(n_ij > n_kl){
      for(int iz = 0; iz < isub->size_zeta(); ++iz)
      for(int jz = 0; jz < jsub->size_zeta(); ++jz)
      for(int kz = 0; kz < ksub->size_zeta(); ++kz)
      for(int lz = 0; lz < lsub->size_zeta(); ++lz) {
	CalcPrimERI1(gi->sym_group(), isub, jsub, ksub, lsub,
		     isub->zeta_iz[iz], jsub->zeta_iz[jz],
		     ksub->zeta_iz[kz], lsub->zeta_iz[lz],
		     prim, method);
	CalcTransERI(isub, jsub, ksub, lsub, iz, jz, kz, lz, prim, eri, true);
      }      
    }

  }
  void CalcERI_ijkk(SymGTOs& g0, SymGTOs& g1,SymGTOs& g2,
		    SubIt isub, SubIt jsub, SubIt ksub,
		    A4dc& prim, ERIMethod method, B2EInt eri) {
    if(!ExistNon0(isub, jsub, ksub, ksub))
      return;
    SubIt lsub = ksub;
    for(int iz = 0; iz < isub->size_zeta(); ++iz)
    for(int jz = 0; jz < jsub->size_zeta(); ++jz)
    for(int kz = 0; kz < ksub->size_zeta(); ++kz)
    for(int lz = 0; lz < lsub->size_zeta(); ++lz) {
      if(kz > lz) {
	CalcPrimERI1(g0->sym_group(), isub, jsub, ksub, lsub,
		     isub->zeta_iz[iz], jsub->zeta_iz[jz],
		     ksub->zeta_iz[kz], lsub->zeta_iz[lz],
		     prim, method);
	CalcTransERI(isub, jsub, ksub, lsub, kz, lz, iz, jz, prim, eri,
		     false, false, true);

      } else if(kz == lz) {
	CalcPrimERI1(g0->sym_group(), isub, jsub, ksub, lsub,
		     isub->zeta_iz[iz], jsub->zeta_iz[jz],
		     ksub->zeta_iz[kz], lsub->zeta_iz[lz],
		     prim, method);
	CalcTransERI(isub, jsub, ksub, lsub, kz, lz, iz, jz, prim, eri,
		     false, false, false);
      }
    }
  }


  // ==== Interface ====
  B2EInt CalcERI_Complex(SymGTOs i, ERIMethod method) {
    return CalcERI(i, i, i, i, method);
  }
  B2EInt CalcERI_Hermite(SymGTOs i, ERIMethod method) {
    SymGTOs ci = i->Conj();
    return CalcERI(ci, i, ci, i, method);
  }
  B2EInt CalcERI(SymGTOs gi, SymGTOs gj, SymGTOs gk, SymGTOs gl,
		 ERIMethod method) {

    if(not gi->setupq)
      gi->SetUp();
    
    if(not gj->setupq)
      gj->SetUp();

    if(not gk->setupq)
      gk->SetUp();

    if(not gl->setupq)
      gl->SetUp();
    
    B2EInt eri(new B2EIntMem);
    eri->Init(gi->size_basis() * gj->size_basis() * gk->size_basis() * gl->size_basis());
    A4dc prim(gi->max_num_prim() * gj->max_num_prim() *
	      gk->max_num_prim() * gl->max_num_prim());
    
    for(SubIt isub = gi->subs().begin(); isub != gi->subs().end(); ++isub) 
      for(SubIt jsub = gj->subs().begin(); jsub != gj->subs().end(); ++jsub)
	for(SubIt ksub = gk->subs().begin(); ksub != gk->subs().end(); ++ksub)
	  for(SubIt lsub = gl->subs().begin(); lsub != gl->subs().end(); ++lsub)
	      if(method.perm == 0) 
		CalcERI0(gi->sym_group(),isub, jsub, ksub, lsub, prim, method, eri);
	      else
		CalcERI1(gi, gj, gk, gl, isub, jsub, ksub, lsub, prim, method, eri);

    return eri;

  }
  /*
  B2EInt CalcERI_0122(SymGTOs g0, SymGTOs g1, SymGTOs g2) {

    if(g0->setupq)
      g0->SetUp();
    if(g1->setupq)
      g1->SetUp();
    if(g2->setupq)
      g2->SetUp();    

    SymGTOs g3 = g2;
    B2EInt eri(new B2EInt);
    eri->Init(g0->size_basis()*g1->size_basis()*g2->size_basis()*g3->size_basis());
    A4dc prim(g0->max_num_prim() * g2->max_num_prim() *
	      g1->max_num_prim() * g3->max_num_prim());    
    for(SubIt isub = g0->subs.begin(); isub != g0->subs.end(); ++isub) 
	for(SubIt jsub = g1->subs.begin(); jsub != g1->subs.end(); ++jsub)
	  for(SubIt ksub = g2->subs.begin(); ksub != g2->subs.end(); ++ksub) {

	    // (ijkk) type

	    // (ijkl) (k!=l) type
	    for(SubIt lsub = g2->subs.begin(); lsub != ksub; ++lsub) {
	    }
	  }      
  }
  */

  // =====  Transform ERI ====
  
  
}
