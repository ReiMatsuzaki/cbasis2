#include <fstream>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Core>
#include <Eigen/QR>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../utils/read_json.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/angmoment.hpp"
#include "../src_cpp/mo.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/two_int.hpp"
#include "../src_cpp/symmol_read_json.hpp"
#include "../src_cpp/ivo.hpp"

using namespace std;
using namespace Eigen;
using namespace cbasis;
using namespace picojson;
using namespace boost::assign;
using boost::format;

// -- const --
double au2ev = 27.2114;
double c_light = 137.035999139;
double au2mb = 5.291772 * 5.291772;

// -- Input --
string comment;
string in_json, out_json, cs_csv;
object obj;

// -- Output --
ofstream f_pad;

// -- target --
SymmetryGroup sym;
Molecule mole;
dcomplex Z;
Molecule mole0;
vector<double> w_list;
int ne;
string calc_type;
//bool use_stex;
enum ECalcTerm {ECalcTerm_One, ECalcTerm_Full };
ECalcTerm calc_term;
ERIMethod eri_method;

// -- Initial state --
SymGTOs basis0;
Irrep irrep0;
int i0;
dcomplex E0;
VectorXcd c0;

// -- psi1 --
SymGTOs basis1;

// -- psi0/chi0 --
vector<int> Ls;
map<int, SymGTOs> basis_psi0_L;
map<int, SymGTOs> basis_c_psi0_L;
map<int, SymGTOs> basis_chi0_L;

// -- solver --
string solve_driv_type;
LinearSolver linear_solver;

// -- write_psi --
bool write_psi_q;
string write_psi_filename;
VectorXd write_psi_rs;
int   write_psi_lmax;

// -- intermediate --
// ---- for psi1 ----
vector<Irrep> irreps;
vector<int>   ms_xyz;
BMat S1, T1, V1, L1;
BVec s1, sD1, c1, cD1;

// ---- for psi0_L --
map<int, BMat> S0L, T0L, V0L, L0L;
map<int, BMat> V0L1, HV0L1;
map<int, BVec> s0L, sD0L, c0L, Hc0L;
map<int, BVec> s0L_chi; // driven term for psi_p

dcomplex muphi_psi1[10], muphi_psi1_v[10];
MultArray<dcomplex, 2> impsi0_muphi(100),  impsi0_muphi_v(100);
MultArray<dcomplex, 2> impsi0_v_psi1(100), impsi0_v_psi1_v(100);
MultArray<dcomplex, 2> impsi0_chi(100);

// -- results --
MatrixXd result;
int num_header = 14;
int idx_cs(0), idx_cs_sigu(1), idx_cs_piu(2), idx_beta(3);
int idx_cs_alpha(4), idx_cs_sigu_alpha(5), idx_cs_piu_alpha(6);
int idx_cs_v(7), idx_cs_sigu_v(8), idx_cs_piu_v(9), idx_beta_v(10);
int idx_cs_alpha_v(11), idx_cs_sigu_alpha_v(12), idx_cs_piu_alpha_v(13);

void CalcSTEXHamiltonian(SymGTOs g1, SymGTOs g0, ERIMethod m, VectorXcd& c0,
			 BMat *H) {

  BMat S, V;
  
  B2EInt eri_J = CalcERI(g1, g1, g0, g0, m);
  B2EInt eri_K = CalcERI(g1, g0, g0, g1, m);

  CalcSTVMat(g1, g1, &S, H, &V);
  H->Add(1.0, V);
  AddJ(eri_J, c0, irrep0, 1.0, *H);
  AddK(eri_K, c0, irrep0, 1.0, *H);
  
}

void s2y(MultArray<dcomplex, 1>& sM,
	 MultArray<dcomplex, 1> *yM, int L) {
  /**
     compute spherical harmonics values from solid spherical harmonics values.
     
     case M > 0
     <0 | Y_LM> = (-)^M / sqrt(2) * (<0|S_LM + j S_L-M>)
     .          = (-)^M / sqrt(2) * (<0|S_LM >+ j <0|S_L-M>)
     <0 | Y_L-M> =   1.0 / sqrt(2) * (<0|S_LM > - j <0|S_L-M>)
   */

  dcomplex ii(0, 1);
  (*yM)(0) = sM(0);
  for(int M = 1; M <= L; M++) {
    (*yM)(M) = pow(-1, M) * (sM(M) + ii*sM(-M)) / sqrt(2.0);
    (*yM)(-M)=              (sM(M) - ii*sM(-M)) / sqrt(2.0);
  }
}
void s2y_bra(MultArray<dcomplex, 1>& sM,
	 MultArray<dcomplex, 1> *yM, int L) {
  /**
     <Y_LM | 0> = <0^* |Y_LM^*> = 
   */

  dcomplex ii(0, 1);
  (*yM)(0) = sM(0);
  for(int M = 1; M <= L; M++) {
    (*yM)(M) = pow(-1, M) * (sM(M) - ii*sM(-M)) / sqrt(2.0);
    (*yM)(M) =              (sM(M) + ii*sM(-M)) / sqrt(2.0);
  }

}
void ss0_2_yy0(MultArray<dcomplex, 1>& sM,
	       MultArray<dcomplex, 1> *yM, int L) {
  /**
     <S_L0 | S_L'0 | 0> = <Y_L0 | Y_L'0 | 0>
     <S_L+M | S_L'0 | 0> = 0
     <S_L+M | S_L'M | 0> 
     .= 2^(-0.5) { (-)^M <Y_L+M|S_L'M|0> -j(-)^M<Y_L-M|S_L'M|0> }
     .= 1/2 {    <Y_L+M|Y_L'+M|0> + j <Y_L+M|Y_L'-M|0>
     .        -j (-)^M <Y_L-M|Y_L'M|0> -1(-)^M}
     <S_L-M | S_L'0 | 0> 
     .= 2^(-0.5)      { <Y_L+M|S_L'0|0> +j <Y_L-M|S_L'0|0> }
     <S_L+M | S_L' | 0> 
     .= 2^(-0.5)(-)^M { <Y_L+M|S_L'0|0> -j <Y_L-M|S_L'0|0> }
     <S_L+M | S_L'0 | 0> = 2^(-0.5)(-)^M {  <Y_L+M | S_L'0 | 0> 
     .                                    -j<Y_L-M | S_L'0 | 0> }


     <Y_+M | Y_M | 0> 
     .  = 2^(-1) { <S_+M | S_M | 0> + <S_-M | S_-M | 0> }
     <Y_-M | Y_-M | 0> 
     .  = 2^(-1) { <S_+M | S_M | 0> + <S_-M | S_-M | 0> }

   */
}
double CoulombShift(dcomplex eta, int L) {
  // -- eta = Z1*Z2/k
  dcomplex ii(0,1);
  dcomplex argment = 1.0 + L + ii * eta;
  double re = argment.real();
  double im = argment.imag();
  gsl_sf_result lnr, arg;
  gsl_sf_lngamma_complex_e(re, im, &lnr, &arg);
  return arg.val;
}
dcomplex CoefAzeta0(dcomplex w, MultArray<dcomplex, 2>& Alm,
		    int zeta,
		    vector<int>& L_list, vector<int>& Ms) {
  /**
     see
     J. C. Tully, R. S. Berry, and B. J. Dalton, 
     Physical Review 176, 95 (1968)
  */
  
  dcomplex c0 = 4.0 * M_PI * M_PI / (c_light*w);
  dcomplex cumsum(0);
  BOOST_FOREACH(int L1, L_list) {
    BOOST_FOREACH(int L2, L_list) {
      dcomplex c1 = sqrt(1.0 * (2*L1+1)*(2*L2+1)) / (4*M_PI*(2*zeta+1));
	
      BOOST_FOREACH(int M1, Ms) {
	BOOST_FOREACH(int M2, Ms) {
	  int mzeta = -M1+M2;
	  dcomplex c3 = (cg_coef(1,  1,-M1,M2,  zeta,mzeta) *
			 cg_coef(1,  1,0,0,     zeta,0) *
			 cg_coef(L2,L1,+M2,-M1, zeta,mzeta) *
			 cg_coef(L2,L1,0,0,     zeta,0));
	  dcomplex c4 = Alm(L1,M1) * conj(Alm(L2,M2));
	  cumsum += c0*c1*c3*c4;
	}
      }
    }
  }
  
  return cumsum;
}
double Coef_FixedMole(dcomplex w, MultArray<dcomplex, 2>& Alm,
		      vector<int>& L_list, int L, int M ) {
  /**
     Gives coefficient of differential cross section for 
     fixed molecule axis with paralled electronic field.
     see note @ 2017/3/15

     Inputs
     ------
     Alm(L, M) = w i sqrt(2/3) (i^-L) e^{i eta_l} <phi^-, mu phi_i>

     d sigma/d Omega = sum_L C_L P_L

   */

  dcomplex cumsum(0);
  BOOST_FOREACH(int L1, L_list) {    
    BOOST_FOREACH(int L2, L_list) {
      dcomplex c1 = sqrt(1.0*(2*L1+1)*(2*L2+1));
      dcomplex c2 = cg_coef(L1,L2,0,0, L,0) * cg_coef(L1,L2,M,-M, L,0);
      dcomplex c3 = Alm(L1, M) * conj(Alm(L2, M));
      cumsum += c1 * c2 * c3;
    }
  }
  return (cumsum).real();
}
void Parse() {
  PrintTimeStamp("Parse", NULL);  
  try {
    ifstream f(in_json.c_str());
    if(f.fail()) {
      throw runtime_error("failed to open input json file");
    }
    value json; f >> json;
    if(not json.is<object>()) 
      throw runtime_error("invalid json file");
    obj = json.get<object>();
    comment = ReadJson<string>(obj, "comment");
    calc_type = ReadJson<string>(obj, "calc_type");
    
    // -- write wave function --    
    if(obj.find("write_psi") != obj.end())  {
      write_psi_q = true;
      CheckObject<object>(obj, "write_psi");      
      object& write_psi_obj = obj["write_psi"].get<object>();
      write_psi_filename = ReadJson<string>(write_psi_obj, "file");
      write_psi_rs = ReadJson<VectorXd>(write_psi_obj, "rs");
      write_psi_lmax = ReadJson<int>(write_psi_obj, "lmax");
    } else {
      write_psi_q = false;
    }

    // -- number of calculation term --
    string str_calc_term = ReadJsonWithDefault<string>(obj, "calc_term", "full");
    if(str_calc_term == "one") {
      calc_term = ECalcTerm_One;
    } else if (str_calc_term == "full") {
      calc_term = ECalcTerm_Full;
    } else {
      throw runtime_error("calc_term must be \"one\" or \"full\"");
    }

    // solver for driven equation
    solve_driv_type = ReadJsonWithDefault
      <string>(obj, "solve_driv_type", "linear_solve");

    if(solve_driv_type == "linear_solve") {
      linear_solver = ReadJsonWithDefault
	<LinearSolver>(obj, "linear_solver", LinearSolver());
    } else if(solve_driv_type == "eigen_value") {

    } else {
      THROW_ERROR("unsupported solve_driv_type");
    }
    ne = ReadJson<int>(obj, "num_ele");

    if(calc_type == "STEX" or calc_type == "RPA" or calc_type == "propagator") {
      eri_method = ReadJson<ERIMethod>(obj, "eri_method");
    }
    sym     = ReadJson<SymmetryGroup>(obj, "sym");
    out_json = ReadJson<string>(obj, "out_json");
    cs_csv   = ReadJson<string>(obj, "cs_csv");
    mole = NewMolecule(sym); ReadJson_Molecule(obj, "molecule", mole);
    Z = ReadJson<dcomplex>(obj, "Z");
    int lmax = ReadJsonWithDefault<int>(obj, "lmax", 3);
    for(int L = 1; L <= lmax; L += 2) {
      Ls.push_back(L);
    }
    cout << "Ls: ";
    BOOST_FOREACH(int L, Ls) {
      cout << L << " " ;
    }
    cout << endl;
    cout << "Ls0 = " << Ls[0] << endl;
    cout << "Ls1 = " << Ls[Ls.size()-1] << endl;
    impsi0_muphi.SetRange( Ls[0], Ls[Ls.size()-1], -1,1);
    impsi0_muphi_v.SetRange( Ls[0], Ls[Ls.size()-1], -1,1);

    impsi0_v_psi1.SetRange(Ls[0], Ls[Ls.size()-1], -1,1);
    impsi0_v_psi1_v.SetRange(Ls[0], Ls[Ls.size()-1], -1,1);
    impsi0_chi.SetRange(   Ls[0], Ls[Ls.size()-1], -1,1);

    // -- write coeff for molecule fixed electronic field --
    string pad_csv = ReadJsonWithDefault<string>(obj, "pad_csv", "");
    cout << "pad_csv: " << pad_csv << endl;
    if(pad_csv != "") {
      cout << 1 << endl;
      f_pad.open(pad_csv.c_str(), ios::out);
      cout << 2 << endl;
      if(f_pad.fail()) {
	throw runtime_error("failed to open file.");
      }
      int lmax = 2*Ls[Ls.size()-1];
      f_pad << "w";
      for(int M = 0; M <= 1; M++)
	for(int L = 0; L <= lmax; L++)
	  f_pad << format(",C_l_%d_%d,C_v_%d_%d") % L % M % L % M;
      f_pad << endl;
    }


    cout << "mole" << endl;
    mole0 = NewMolecule(sym);
    Atom cen0 = NewAtom("CEN0", Z); cen0->Add(0, 0, 0);
    mole0->Add(cen0);
    
    basis0 = NewSymGTOs(mole); ReadJson_SymGTOs_Subs(obj, "basis0", basis0);
    basis1 = NewSymGTOs(mole); ReadJson_SymGTOs_Subs(obj, "basis1", basis1);
    ReadJson_Orbital(obj, "orbital0", sym, &E0, &c0, &irrep0, &i0);

    VectorXi Ms(3); Ms << -1, 0, 1;
    VectorXcd zeta_chi(1); zeta_chi << ReadJson<dcomplex>(obj, "zeta_chi");
    BOOST_FOREACH(int L, Ls) {
      string name;
      if(L == 1) 
	name = "p";
      else if(L == 3)
	name = "f";
      else if(L == 5)
	name = "h";
      else
	throw runtime_error("unsupported L");
      cout << "start raeding psi0 for " << name << endl;
      VectorXcd zeta0 = ReadJson<VectorXcd>(obj, "zeta0_" + name);
      basis_psi0_L[L] = NewSymGTOs(mole0);
      basis_psi0_L[L]->NewSub("CEN0").SolidSH_Ms(L, Ms, zeta0);
      basis_chi0_L[L] = NewSymGTOs(mole0);
      basis_chi0_L[L]->NewSub("CEN0").SolidSH_Ms(L, Ms, zeta_chi);
      cout << "end raeding psi0 for " << name << endl;
    }

    VectorXd _ws = ReadJson<VectorXd>(obj, "ws");
    for(int i = 0; i < _ws.size(); i++) {
      w_list.push_back(_ws(i));
    }
    result = MatrixXd::Zero(w_list.size(), num_header);

    
  } catch(exception& e) {
    cerr << "error on parsing json" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  cout << "basis0 and basis1 setup" << endl;
  try {
    basis0->SetUp();    
    basis1->SetUp();
  } catch(exception& e) {
    cerr << "error on SetUp basis: " << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  cout << "basis_psi0_L setup" << endl;
  try {
    BOOST_FOREACH(int L, Ls) {
      basis_psi0_L[L]->SetUp();
    }
  } catch(exception& e) {
    cerr << "error on SetUp basis: " << endl;
    cerr << e.what() << endl;
    cerr << basis_psi0_L[1]->str() << endl;
    cerr << basis_psi0_L[3]->str() << endl;
    exit(1);
  }
  cout << "basis_chi setup" << endl;
  try {
    BOOST_FOREACH(int L, Ls) {
      basis_chi0_L[L]->SetUp();
    }
  } catch(exception& e) {
    cerr << "error on SetUp basis_chi0_L" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  cout << "basis_c_psi0 conj" << endl;
  BOOST_FOREACH(int L, Ls) {
    basis_c_psi0_L[L] = basis_psi0_L[L]->Conj();
  }
  irreps += sym->irrep_x(),sym->irrep_y(),sym->irrep_z();
  ms_xyz += +1, -1, 0;
}
void PrintIn() {
  PrintTimeStamp("PrintIn", NULL);
  cout << "comment: " << comment << endl;
  cout << "in_json: " << in_json << endl;
  if(write_psi_q) {
    cout << "write_psi_filename: " << write_psi_filename << endl;
    VectorXd& rs = write_psi_rs;
    cout << format("write_psi_rs: [%5.3f, %5.3f, ..., %5.3f]\n") % rs[0] % rs[1] % rs[rs.size()-1];
    cout << "write_psi_lmax: " << write_psi_lmax << endl;
  } else {
    cout << "write_psi: no\n";
  }
  cout << "cs_csv: "  << cs_csv;
  cout << "out_json: " << out_json << endl;
  cout << "calc_type: " << calc_type << endl;
  cout << "calc_term: " << (calc_term == ECalcTerm_One ? "one" : "full") << endl;
  cout << "solve_driv_type: " << solve_driv_type << endl;  
  cout << "linear_solver: " << linear_solver.str() << endl;
  cout << "ERIMethod_use_symmetry: " << eri_method.symmetry << endl;
  cout << "ERIMethod_use_memo: " << eri_method.coef_R_memo << endl;
  cout << "ERIMethod_use_perm: " << eri_method.perm << endl;  
  cout << "Ne: " << ne << endl;
  cout << "E0: " << E0 << endl;
  cout << "Z: " << Z << endl;
  cout << "symmetry: " << sym->name() << endl;
  cout << "molecule:" << endl << mole->show() << endl;
  cout << "basis0:" << endl << basis0->show() << endl;  
  cout << "basis1:" << endl << basis1->show() << endl;
  BOOST_FOREACH(int L, Ls) {
    cout << "basis0_:" << L << endl << basis_psi0_L[L]->show() << endl;
  }
}
void CalcMat() {

  S1.set_name("S1"); T1.set_name("T1"); V1.set_name("V1");
  
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z();

  PrintTimeStamp("psi1", NULL);
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);
  
  PrintTimeStamp("psi1/init", NULL);
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  s1(x) = X1i(x, 0) * c0; sD1(x) = DX1i(x, 0) * c0; 
  s1(y) = Y1i(y, 0) * c0; sD1(y) = DY1i(y, 0) * c0; 
  s1(z) = Z1i(z, 0) * c0; sD1(z) = DZ1i(z, 0) * c0;  

  PrintTimeStamp("psi0,chi0", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];
    SymGTOs chi0 = basis_chi0_L[L];
    CalcSTVMat(psi0, psi0, &S0L[L], &T0L[L], &V0L[L]);

    BMat S0L_chi; CalcSMat(psi0, chi0, &S0L_chi);
    s0L_chi[L](x) = S0L_chi(x, x).col(0);
    s0L_chi[L](y) = S0L_chi(y, y).col(0);
    s0L_chi[L](z) = S0L_chi(z, z).col(0);

    BMat X0i, DX0i, Y0i, DY0i, Z0i, DZ0i;
    CalcDipMat(psi0, basis0, &X0i, &Y0i, &Z0i, &DX0i, &DY0i, &DZ0i);
    s0L[L](x) = X0i(x,0) * c0; sD0L[L](x) = DX0i(x,0) * c0;
    s0L[L](y) = Y0i(y,0) * c0; sD0L[L](y) = DY0i(y,0) * c0;
    s0L[L](z) = Z0i(z,0) * c0; sD0L[L](z) = DZ0i(z,0) * c0;

    BMat V01_full, HV01_full, V01_0th, HV01_0th;
    CalcVMat(psi0,   mole,  basis1, &V01_full);
    CalcVMat(c_psi0, mole,  basis1, &HV01_full);
    CalcVMat(psi0,   mole0, basis1, &V01_0th);
    CalcVMat(c_psi0, mole0, basis1, &HV01_0th);

    V0L1[L](x,x) = V01_full(x,x) - V01_0th(x,x);
    V0L1[L](y,y) = V01_full(y,y) - V01_0th(y,y);
    V0L1[L](z,z) = V01_full(z,z) - V01_0th(z,z);
    
    HV0L1[L](x,x) = HV01_full(x,x) - HV01_0th(x,x);
    HV0L1[L](y,y) = HV01_full(y,y) - HV01_0th(y,y);
    HV0L1[L](z,z) = HV01_full(z,z) - HV01_0th(z,z);
  }
  
}
void CalcMatSTEX() {
  PrintTimeStamp("MatSTEX_1", NULL);

  B2EInt eri_J_11 = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K_11 = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  AddJ(eri_J_11, c0, irrep0, 1.0, V1);
  AddK(eri_K_11, c0, irrep0, 1.0, V1);

  /*
  BMat H_STEX;
  CalcSTEXHamiltonian(basis1, basis0, eri_method, c0, &H_STEX);
  Copy(H_STEX, V1);
  V1.Add(-1.0, T1);
  */

  PrintTimeStamp("MatSTEX_01", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];
    B2EInt eri_JC = CalcERI(psi0,   basis1, basis0, basis0, eri_method);
    B2EInt eri_JH = CalcERI(c_psi0, basis1, basis0, basis0, eri_method);
    B2EInt eri_KC = CalcERI(psi0,   basis0, basis0, basis1, eri_method);
    B2EInt eri_KH = CalcERI(c_psi0, basis0, basis0, basis1, eri_method);
    AddJ(eri_JC, c0, irrep0, 1.0, V0L1[L]);
    AddK(eri_KC, c0, irrep0, 1.0, V0L1[L]);
    AddJ(eri_JH, c0, irrep0, 1.0, HV0L1[L]);
    AddK(eri_KH, c0, irrep0, 1.0, HV0L1[L]);
  }
}
void CalcMatRPA() {
  PrintTimeStamp("MatRPA_1", NULL);

  // ==== psi1 ====
  B2EInt eri_J_11 = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K_11 = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  
  // -- compute A matrix 
  // -- A = T + V_ne + J + K - epsilon0  
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);  
  BMat A("A"); Copy(T1, A);    
  try {
    A.Add(1.0, V1);
  } catch(exception& e) {
    cout << e.what() << endl;
    THROW_ERROR("ERROR on A+=V1");
  } 
  AddJ(eri_J_11, c0, irrep0, 1.0, A);
  AddK(eri_K_11, c0, irrep0, 1.0, A);
  try {
    A.Add(-E0, S1);
  } catch(exception& e) {
    cout << e.what() << endl;
    THROW_ERROR("error on A+=-E0 S1");
  }

  // -- compute B matrix --
  // -- B = K
  BMat B("B"); Copy(T1, B); B.SetZero();  
  AddK(eri_K_11, c0, irrep0, 1.0, B);

  // -- compute V0 --
  BMat V1_0th("V1_0th");
  CalcVMat(basis1, mole0, basis1, &V1_0th);
  
  // -- compute interaction for psi1 --
  // -- V = (H-T) = sqrt(sqrt(A-B)(A+B)sqrt(A-B)) - T + epsilon0
  BMat ApB("ApB"); Copy(A, ApB); ApB.Add(+1.0, B);  
  BMat AmB("AmB"); Copy(A, AmB); AmB.Add(-1.0, B);
  BMat sqrt_AmB("sqrt_AmB"); BMatSqrt(AmB, sqrt_AmB);
  BMat H2("H2");  Multi3(sqrt_AmB, ApB, sqrt_AmB, H2);
  BMatSqrt(H2, V1);
  V1.Add(+E0, S1);
  V1.Add(-1.0, T1);
  
  // ==== psi0/psi1 ====  
  PrintTimeStamp("MatRPA_2", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];

    BMat T01, HT01, S01, HS01, V01, HV01;
    BMat V01_0th, HV01_0th;
    CalcSTMat(psi0,   basis1, &S01,  &T01);
    CalcSTMat(c_psi0, basis1, &HS01, &HT01);
    CalcVMat(psi0,    mole,  basis1, &V01);
    CalcVMat(c_psi0,  mole,  basis1, &HV01);
    CalcVMat(psi0,    mole0, basis1, &V01_0th);
    CalcVMat(c_psi0,  mole0, basis1, &HV01_0th);

    // -- ERI --
    B2EInt eri_JC = CalcERI(psi0,   basis1, basis0, basis0, eri_method);
    B2EInt eri_JH = CalcERI(c_psi0, basis1, basis0, basis0, eri_method);
    B2EInt eri_KC = CalcERI(psi0,   basis0, basis0, basis1, eri_method);
    B2EInt eri_KH = CalcERI(c_psi0, basis0, basis0, basis1, eri_method);

    // -- A matrix --
    // -- A = T + V_ne + J + K - epsilon0  
    BMat CA; Copy(T01,  CA); CA.Add(1.0, V01);
    BMat HA; Copy(HT01, HA); HA.Add(1.0, HV01);    
    AddJ(eri_JC, c0, irrep0, 1.0, CA); AddK(eri_KC, c0, irrep0, 1.0, CA);
    AddJ(eri_JH, c0, irrep0, 1.0, HA); AddK(eri_KH, c0, irrep0, 1.0, HA);
    CA.Add(-E0, S01);
    HA.Add(-E0, HS01);

    // -- B matrix --
    BMat CB; Copy(T1, CB); CB.SetZero();
    BMat HB; Copy(T1, HB); HB.SetZero();
    AddK(eri_KC, c0, irrep0, 1.0, CB);
    AddK(eri_KH, c0, irrep0, 1.0, HB);

    // -- interaction --
    // -- V = H-H_0 = sqrt(sqrt(A-B)(A+B)sqrt(A-B)) + epsilon - T - V0
    BMat C_ApB; Copy(CA, C_ApB); C_ApB.Add(+1.0, CB);
    BMat C_AmB; Copy(CA, C_AmB); C_AmB.Add(-1.0, CB);
    BMat C_sqrt; BMatSqrt(C_AmB, C_sqrt);
    BMat C_H2; Multi3(C_sqrt, C_ApB, C_sqrt, C_H2);
    BMatSqrt(C_H2, V0L1[L]);
    V0L1[L].Add(E0, S01);
    V0L1[L].Add(-1.0, T01);
    V0L1[L].Add(-1.0, V01_0th);    
    BMat H_ApB; Copy(HA, H_ApB); H_ApB.Add(+1.0, HB);
    BMat H_AmB; Copy(HA, H_AmB); H_AmB.Add(-1.0, HB);
    BMat H_sqrt; BMatSqrt(H_AmB, H_sqrt);
    BMat H_H2; Multi3(H_sqrt, H_ApB, H_sqrt, H_H2);
    BMatSqrt(H_H2, HV0L1[L]);
    HV0L1[L].Add(E0,   HS01);
    HV0L1[L].Add(-1.0, HT01);
    HV0L1[L].Add(-1.0, HV01_0th);

  }
  
}
void CalcDriv_linear_solve(int iw) {
  PrintTimeStamp("CalcDriv_linear_solve", NULL);
  double w = w_list[iw];

  dcomplex ene = E0 + w;
  
  // -- Compute psi1 --
  BOOST_FOREACH(Irrep i, irreps) {
    L1(i,i) = S1(i,i) * ene - T1(i,i) - V1(i,i);
    ColPivHouseholderQR<MatrixXcd> pivx = L1(i,i).colPivHouseholderQr();
    c1(i)  = pivx.solve(s1(i));
    cD1(i) = pivx.solve(sD1(i));
  }

  // -- Compute psi0_p --
  BOOST_FOREACH(int L, Ls) {
    BOOST_FOREACH(Irrep i, irreps) {
      L0L[L](i,i) = S0L[L](i,i) * ene - T0L[L](i,i) - V0L[L](i,i);
      c0L[L](i)   = L0L[L](i,i).colPivHouseholderQr().solve(s0L_chi[L](i));
      Hc0L[L](i)  = c0L[L](i).conjugate();
    }
  }
  
}
void CalcDriv_eigen_value(int iw) {
  PrintTimeStamp("CalcDriv_eigen_value", NULL);
  double w = w_list[iw];
  dcomplex ene = E0 + w;
  
  BOOST_FOREACH(Irrep i, irreps) {

    // -- psi 1 --
    cout << "calculate psi1" << endl;  
    MatrixXcd H = T1(i,i)+V1(i,i);
    SymGenComplexEigenSolver solver(H, S1(i,i));
    const VectorXcd& eig = solver.eigenvalues();    
    const MatrixXcd& U1 = solver.eigenvectors();
    CtAC(U1, T1(i,i)); CtAC(U1, V1(i,i)); CtAC(U1, S1(i,i));
    Ctx(U1, s1(i));    Ctx(U1, sD1(i));
    c1(i) =  VectorXcd::Zero(s1(i).size());      
    cD1(i) = VectorXcd::Zero(s1(i).size());          
    for(int j = 0; j < U1.rows(); j++) {
      c1(i)(j) = s1(i)(j) / (ene - eig(j));
      cD1(i)(j) = sD1(i)(j) / (ene - eig(j));
    }      
        
    BOOST_FOREACH(int L, Ls) {
      // -- psi0 --
      cout << format("calculate psi0. L = %d\n") % L;
      MatrixXcd H0 = T0L[L](i,i) + V0L[L](i,i);
      MatrixXcd S0 = S0L[L](i,i);
      SymGenComplexEigenSolver solver(H0, S0);
      const VectorXcd& eig0 = solver.eigenvalues();
      cout << eig0 << endl;
      const MatrixXcd& U0   = solver.eigenvectors();
      CtAC(U0, T0L[L](i,i)); CtAC(U0, V0L[L](i,i)); CtAC(U0, S0L[L](i,i));
      Ctx(U0, s0L_chi[L](i));
      int num_Li = s0L_chi[L](i).size();
      c0L[L](i) = VectorXcd::Zero(num_Li);
      for(int j = 0; j < num_Li; j++) {
	c0L[L](i)(j) = s0L_chi[L](i)(j) / (ene - eig0(j));
      }
      Hc0L[L](i) = c0L[L](i).conjugate();

      // -- psi0/psi1 --
      cout << format("convert V, HV. L = %d\n") % L;
      CtAD(U0, U1, V0L1[L](i,i));
      CaAD(U0, U1, HV0L1[L](i,i));

      // -- psi0/phi_I --
      Ctx(U0, s0L[L](i)); Ctx(U0, sD0L[L](i));
      
    }
    cout << "r(L=1)" << endl;
    cout << s0L[1] << endl;
  }

}
void CalcBraket() {

  /**
     compute some braket.
   */
  // <ImPsi|0> = <(y-Hy)/2j|0>
  //           = -1/2j*(<y|0>-<Hy|0>)
  //           = -1/2j*((Hy|0)-(y|0))
  //           = -1/2j*((Hy|0)-(y|0))  
  
  PrintTimeStamp("calc_braket", NULL);
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  dcomplex m2(0, 2);

  BOOST_FOREACH(Irrep i, irreps) {
    muphi_psi1[i] = TDot(c1(i), s1(i)) / 3.0 * (1.0*ne);
    muphi_psi1_v[i] = TDot(cD1(i), sD1(i)) / 3.0 * (1.0*ne);
  }
  
  
  BOOST_FOREACH(int L, Ls) {
    BVec& c0 = c0L[L];
    BVec& Hc0 = Hc0L[L];

    // -- <ImPsi0 | mu | PhiInit> --
    impsi0_muphi(L, +1) = TDot(c0(x), s0L[L](x)).imag();
    impsi0_muphi(L, -1) = TDot(c0(y), s0L[L](y)).imag();
    impsi0_muphi(L,  0) = TDot(c0(z), s0L[L](z)).imag();
    
    impsi0_muphi_v(L, +1) = TDot(c0(x), sD0L[L](x)).imag();
    impsi0_muphi_v(L, -1) = TDot(c0(y), sD0L[L](y)).imag();
    impsi0_muphi_v(L,  0) = TDot(c0(z), sD0L[L](z)).imag();

    BMat& V  = V0L1[L];
    BMat& HV = HV0L1[L];
    
    impsi0_v_psi1(L, +1) = (TDot(c0(x),   V( x,x)*c1(x))
			    -TDot(Hc0(x), HV(x,x)*c1(x)))/m2;
    impsi0_v_psi1(L, -1) = (TDot(c0(y),   V( y,y)*c1(y))
			    -TDot(Hc0(y), HV(y,y)*c1(y)))/m2;
    impsi0_v_psi1(L,  0) = (TDot(c0(z),   V( z,z)*c1(z))
			    -TDot(Hc0(z), HV(z,z)*c1(z)))/m2;

    impsi0_v_psi1_v(L, +1) = (TDot(c0(x),   V( x,x)*cD1(x))
			      -TDot(Hc0(x), HV(x,x)*cD1(x)))/m2;
    impsi0_v_psi1_v(L, -1) = (TDot(c0(y),   V( y,y)*cD1(y))
			      -TDot(Hc0(y), HV(y,y)*cD1(y)))/m2;
    impsi0_v_psi1_v(L,  0) = (TDot(c0(z),   V( z,z)*cD1(z))
			      -TDot(Hc0(z), HV(z,z)*cD1(z)))/m2;
    
    impsi0_chi(L, +1) = TDot(c0(x), s0L_chi[L](x)).imag();
    impsi0_chi(L, -1) = TDot(c0(y), s0L_chi[L](y)).imag();
    impsi0_chi(L,  0) = TDot(c0(z), s0L_chi[L](z)).imag();    
    
  }

}
void CalcMain_alpha(int iw) {

  /**
     Inputs
     -------
     iw       : index for ws
     c1, s1   : for length form
     cD1, sD1 : for velocity form
   */
  
  PrintTimeStamp("calc_alpha", NULL);

  double w = w_list[iw];
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();

  dcomplex alpha = muphi_psi1[x] + muphi_psi1[y] + muphi_psi1[z];

  double cs = 4.0*M_PI*w/ c_light * alpha.imag()* au2mb;
  double cs_sigu = 4.0*M_PI*w/ c_light * muphi_psi1[z].imag()* au2mb;
  double cs_piu  = 4.0*M_PI*w/ c_light * (muphi_psi1[x]+muphi_psi1[y]).imag()* au2mb;  
  result(iw, idx_cs_alpha) = cs;
  result(iw, idx_cs_sigu_alpha) = cs_sigu;
  result(iw, idx_cs_piu_alpha) = cs_piu;  

  dcomplex alpha_d = muphi_psi1_v[x] + muphi_psi1_v[y] + muphi_psi1_v[z];
  double cs_v      = 4.0*M_PI/(w*c_light) * alpha_d.imag() * au2mb;
  double cs_sigu_v = 4.0*M_PI/(w*c_light) * muphi_psi1_v[z].imag() * au2mb;
  double cs_piu_v  = 4.0*M_PI/(w*c_light) * (muphi_psi1_v[x]+muphi_psi1_v[y]).imag() * au2mb;
  result(iw, idx_cs_alpha_v)      = cs_v;
  result(iw, idx_cs_sigu_alpha_v) = cs_sigu_v;
  result(iw, idx_cs_piu_alpha_v)  = cs_piu_v;

  cout << format("Cs(total,alpha): %20.10f, %20.10f\n") % cs % cs_v;
  cout << format("Cs(sigu,alpha):  %20.10f, %20.10f\n") % cs_sigu % cs_sigu_v;
  cout << format("Cs(piu,alpha):   %20.10f, %20.10f\n") % cs_piu % cs_piu_v;
  
}
void CalcMain(int iw) {
  // -- element of dl_ss is solid spherical, but
  // -- in the special case, it is equivalent to spherical haromonics
  
  // -- total --
  vector<int> Ms_sigu; Ms_sigu += 0;
  vector<int> Ms_piu; Ms_piu += 1,-1;
  vector<int> Ms; Ms += -1,0,+1;
  
  double w = w_list[iw];
  dcomplex k = sqrt(2.0*(w + E0));
  MultArray<dcomplex, 2> Alm(100);   Alm.SetRange(1,3,  -1,1);
  MultArray<dcomplex, 2> Alm_v(100); Alm_v.SetRange(1,3,  -1,1);
  dcomplex c = (sqrt(k/2.0) * dcomplex(0,2) / sqrt(k) * sqrt(3.0/(4.0*M_PI))
		* sqrt(1.0*ne));
  BOOST_FOREACH(int L, Ls) {
    dcomplex ii(0, 1);    
    BOOST_FOREACH(int M, Ms) {
      cout << L << M << impsi0_chi(L,M) 
	   << impsi0_muphi(L,M) << impsi0_v_psi1(L,M)
	   << impsi0_muphi_v(L,M) << impsi0_v_psi1_v(L,M)
	   <<endl;
	
      dcomplex sign = impsi0_chi(L,M)/abs(impsi0_chi(L,M));
      dcomplex etal = exp(ii*CoulombShift(-Z/k, L));
      dcomplex coef_L = sign * w * ii * sqrt(2.0/3.0) * pow(ii, -L) * etal * c; 
      dcomplex psi0_other, psi0_other_v;
      if(calc_term == ECalcTerm_One) {
	psi0_other   = impsi0_muphi(L,M);
	psi0_other_v = impsi0_muphi_v(L,M);
      } else if(calc_term == ECalcTerm_Full) {
	psi0_other   = impsi0_muphi(L,M)   + impsi0_v_psi1(L,M);
	psi0_other_v = impsi0_muphi_v(L,M) + impsi0_v_psi1_v(L,M);
      }
      Alm(L, M)   = coef_L * psi0_other   / sqrt(sign * impsi0_chi(L,M));
      Alm_v(L, M) = coef_L * psi0_other_v / sqrt(sign * impsi0_chi(L,M));
    }
  }

  // Length form, molecular axis averaged
  dcomplex A00 = CoefAzeta0(w, Alm, 0, Ls, Ms);
  dcomplex A20 = CoefAzeta0(w, Alm, 2, Ls, Ms);
  dcomplex A00_sigu = CoefAzeta0(w, Alm, 0, Ls, Ms_sigu);
  dcomplex A00_piu = CoefAzeta0(w, Alm, 0, Ls, Ms_piu);
  double cs   =    (4.0 * M_PI * A00      * au2mb).real(); 
  double cs_sigu = (4.0 * M_PI * A00_sigu * au2mb).real();
  double cs_piu =  (4.0 * M_PI * A00_piu  * au2mb).real();
  double beta = (A20/A00).real();    
  result(iw, idx_cs) = cs;
  result(iw, idx_beta) = beta;
  result(iw, idx_cs_sigu) = cs_sigu;
  result(iw, idx_cs_piu) = cs_piu;

  cout << "coefficient A_lm\n";
  BOOST_FOREACH(int L, Ls) {
    BOOST_FOREACH(int M, Ms) {
      if(L==0 && M!=0) 
	break;
      cout << format("%d, %d, %20.10f, %20.10f\n") % L % M
	% Alm(L, M).real() % Alm(L, M).imag();
    }
  }

  // Velocity form, molecular axis averaged
  dcomplex A00_v = CoefAzeta0(w, Alm_v, 0, Ls, Ms);
  dcomplex A20_v = CoefAzeta0(w, Alm_v, 2, Ls, Ms);
  dcomplex A00_piu_v  = CoefAzeta0(w, Alm_v, 0, Ls, Ms_piu);
  dcomplex A00_sigu_v = CoefAzeta0(w, Alm_v, 0, Ls, Ms_sigu);
  double cs_v      = (4 * M_PI * A00_v     /(w*w) * au2mb).real();
  double cs_piu_v  = (4 * M_PI * A00_piu_v /(w*w) * au2mb).real();
  double cs_sigu_v = (4 * M_PI * A00_sigu_v/(w*w) * au2mb).real();
  double beta_v    = (A20_v/A00_v).real();
  result(iw, idx_cs_v) = cs_v;
  result(iw, idx_beta_v) = beta_v; 
  result(iw, idx_cs_sigu_v) = cs_sigu_v;
  result(iw, idx_cs_piu_v)  = cs_piu_v;
  //  cout << "cross sections (length form, velocity form)" << endl;
  cout << format("Cs(total): %20.10f, %20.10f\n") % cs % cs_v;
  cout << format("Cs(sig_u): %20.10f, %20.10f\n") % cs_sigu % cs_sigu_v;
  cout << format("Cs(pi_u) : %20.10f, %20.10f\n") % cs_piu % cs_piu_v;
  cout << format("beta     : %20.10f, %20.10f\n") % beta % beta_v;

  // fixed molecule
  double coeff =  4.0*M_PI*M_PI / (c_light*w) * 1.0/(4.0*M_PI) * au2mb;
  cout << format("coef: %20.10f\n") % coeff;
  int lmax = 2*Ls[Ls.size()-1];

  if(f_pad.is_open())
    f_pad << format("%10.5f") % w;
  for(int M = 0; M <= 1; M++) {
    cout << format("Coef(M=%d) = \n") % M;
    for(int L = 0; L <= lmax; L++) {      
      double c_l = Coef_FixedMole(w, Alm,   Ls, L, M)      *coeff;
      double c_v = Coef_FixedMole(w, Alm_v, Ls, L, M)/(w*w)*coeff;
      cout << format("%d: %20.10f, %20.10f\n") % L % c_l % c_v; 
      if(f_pad.is_open())
	f_pad << format(",%20.10f,%20.10f") % c_l % c_v;
    }    
    double c0_l = Coef_FixedMole(w, Alm,   Ls, 0, M)       * coeff * 4.0*M_PI/3.0;
    double c0_v = Coef_FixedMole(w, Alm_v, Ls, 0, M)/(w*w) * coeff * 4.0*M_PI/3.0;
    if(M==1) {
      c0_l *= 2.0; c0_v *= 2.0;
    }
    cout << format("Cs(M=%d,fixed): %20.10f, %20.10f\n") % M % c0_l % c0_v;
  }
  if(f_pad.is_open())
    f_pad << endl;

}
void PrintOut() {

  PrintTimeStamp("PrintOut", NULL);
  int num(w_list.size());
  
  ofstream f(cs_csv.c_str(), ios::out);
  f << "w,ene,cs,cs_sigu,cs_piu,beta,cs_alpha,cs_sigu_alpha,cs_piu_alpha,cs_v,cs_sigu_v,cs_piu_v,beta_v,cs_alpha_v,cs_sigu_alpha_v,cs_piu_alpha_v" << endl;
  for(int i = 0; i < num; i++) {
    f << w_list[i] << "," << w_list[i] + E0.real();
    for(int idx = 0; idx < num_header; idx++) {
      f << "," << result(i, idx);
    }
    f << endl;
  }
  f.close();
  
  if(write_psi_q) {
    VectorXd& rs = write_psi_rs;
    int irr = sym->irrep_z();
    ofstream f(write_psi_filename.c_str(), ios::out);
    f << "r";
    for(int L = 0; L <= write_psi_lmax; L++) {
      f << format(",re_%d,im_%d") % L % L;
    }
    f << endl;
    
    vector<VectorXcd> ys(write_psi_lmax+1), dys(write_psi_lmax+1);
    for(int L = 0; L <= write_psi_lmax; L++) {
      basis1->AtR_Ylm(L, 0, irr, c1(irr), rs.cast<dcomplex>(), &ys[L], &dys[L]);
    }
    
    for(int i = 0; i < rs.size(); i++) {
      f << format("%f") % rs[i];
      for(int L = 0; L <= write_psi_lmax; L++) {
	dcomplex y = ys[L][i];
	f << format(",%f,%f") % y.real() % y.imag();
      }
      f << endl;
    }
  }
}


void Calc_one_lin() {
  PrintTimeStamp("Calc_one_lin", NULL);
  PrintTimeStamp("CalcMat", NULL);
  PrintTimeStamp("psi1", NULL);
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);
  
  PrintTimeStamp("psi1/init", NULL); 
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z();
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  s1(x) = X1i(x, 0) * c0; sD1(x) = DX1i(x, 0) * c0; 
  s1(y) = Y1i(y, 0) * c0; sD1(y) = DY1i(y, 0) * c0; 
  s1(z) = Z1i(z, 0) * c0; sD1(z) = DZ1i(z, 0) * c0;  

  PrintTimeStamp("psi0,chi0", NULL);
  BOOST_FOREACH(int L, Ls) {
    cout << "L = " << L << endl;
    SymGTOs psi0   = basis_psi0_L[L];
    SymGTOs c_psi0 = basis_c_psi0_L[L];
    SymGTOs chi0 = basis_chi0_L[L];
    CalcSTVMat(psi0, psi0, &S0L[L], &T0L[L], &V0L[L]);

    BMat S0L_chi; CalcSMat(psi0, chi0, &S0L_chi);
    
    BMat X0i, DX0i, Y0i, DY0i, Z0i, DZ0i;
    CalcDipMat(psi0, basis0, &X0i, &Y0i, &Z0i, &DX0i, &DY0i, &DZ0i);
   
    BMat V01_full, HV01_full, V01_0th, HV01_0th;
    CalcVMat(psi0,   mole,  basis1, &V01_full);
    CalcVMat(c_psi0, mole,  basis1, &HV01_full);
    CalcVMat(psi0,   mole0, basis1, &V01_0th);
    CalcVMat(c_psi0, mole0, basis1, &HV01_0th);

    BOOST_FOREACH(Irrep i, irreps) {
      s0L_chi[L](i) = S0L_chi(i, i).col(0);      
      V0L1[L](i,i) = V01_full(i,i) - V01_0th(i,i);
      HV0L1[L](i,i) = HV01_full(i,i) - HV01_0th(i,i);
    }
    s0L[L](x) = X0i(x,0) * c0; sD0L[L](x) = DX0i(x,0) * c0;
    s0L[L](y) = Y0i(y,0) * c0; sD0L[L](y) = DX0i(y,0) * c0;
    s0L[L](z) = Z0i(z,0) * c0; sD0L[L](z) = DX0i(z,0) * c0;
    
  }
  
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    PrintTimeStamp("CalcDriv(linear_solve)", NULL);
    //    double w = w_list[iw];    
    //    dcomplex ene = E0 + w;

    CalcDriv_linear_solve(iw);
    CalcMain_alpha(iw);    
    CalcBraket();    
    CalcMain(iw);
    
    /*  
    // -- Compute psi1 --
    BOOST_FOREACH(Irrep i, irreps) {
      L1(i,i) = S1(i,i) * ene - T1(i,i) - V1(i,i);
      ColPivHouseholderQR<MatrixXcd> pivx = L1(i,i).colPivHouseholderQr();
      c1(i)  = pivx.solve(s1(i));
      cD1(i) = pivx.solve(sD1(i));
    }

    // -- Compute psi0_p --
    BOOST_FOREACH(int L, Ls) {
      BOOST_FOREACH(Irrep i, irreps) {
	L0L[L](i,i) = S0L[L](i,i) * ene - T0L[L](i,i) - V0L[L](i,i);
	c0L[L](i)   = L0L[L](i,i).colPivHouseholderQr().solve(s0L_chi[L](i));
	Hc0L[L](i)  = c0L[L](i).conjugate();
      }
    }
    */

    // -- alpha --
    /*
    dcomplex alpha[3];
    for(int idx = 0; idx < 3; idx++) {
      int i = irreps[idx];
      alpha[idx] = TDot(c1(i), s1(i))/3.0*(1.0*ne);
    }
    dcomplex alpha_all = alpha[0] + alpha[1] + alpha[2];
    double cs = 4.0*M_PI*w/ c_light * alpha_all.imag()* au2mb;
    double cs_sigu = 4.0*M_PI*w/ c_light * alpha[2].imag()* au2mb;
    double cs_piu  = 4.0*M_PI*w/ c_light * (alpha[0]+alpha[1]).imag()* au2mb;  
    result(iw, idx_cs_alpha) = cs;
    result(iw, idx_cs_sigu_alpha) = cs_sigu;
    result(iw, idx_cs_piu_alpha) = cs_piu;  
    
    for(int idx = 0; idx < 3; idx++) {
      int i = irreps[idx];
      alpha[idx] = TDot(cD1(i), sD1(i))/3.0*(1.0*ne);
    }
    alpha_all = alpha[2] + alpha[0] + alpha[1];
    double cs_v = 4.0*M_PI/(w*c_light) * alpha_all.imag()* au2mb;
    double cs_sigu_v = 4.0*M_PI/(w*c_light) * alpha[2].imag()* au2mb;
    double cs_piu_v  = 4.0*M_PI/(w*c_light) * (alpha[1]+alpha[0]).imag()* au2mb;  
    result(iw, idx_cs_alpha_v)      = cs_v;
    result(iw, idx_cs_sigu_alpha_v) = cs_sigu_v;
    result(iw, idx_cs_piu_alpha_v)  = cs_piu_v;

    cout << format("Cs(total,alpha): %20.10f, %20.10f\n") % cs % cs_v;
    cout << format("Cs(sigu,alpha):  %20.10f, %20.10f\n") % cs_sigu % cs_sigu_v;
    cout << format("Cs(piu,alpha):   %20.10f, %20.10f\n") % cs_piu % cs_piu_v;
    */
    
    /*
    // -- braket --
    BOOST_FOREACH(int L, Ls) {
      for(int idx = 0; idx < 3; idx++) {
	int i = irreps[idx];
	int m = ms_xyz[idx];
	VectorXcd& c1i = c1(i);
	VectorXcd& cD1i = cD1(i);
	VectorXcd& c0 = c0L[L](i);
	VectorXcd& Hc0 = Hc0L[L](i);
	VectorXcd& s0 = s0L[L](i);
	VectorXcd& sD0 = sD0L[L](i);
	MatrixXcd& V  = V0L1[L](i, i);
	MatrixXcd& HV = HV0L1[L](i, i);
	impsi0_muphi(L, m) = TDot(c0, s0).imag();
	impsi0_muphi_v(L, m) = TDot(c0, sD0).imag();
	dcomplex m2(0, 2);
	impsi0_v_psi1(L, m)   = (TDot(c0, V*c1i)  -TDot(Hc0, HV*c1i)) /m2;
	impsi0_v_psi1_v(L, m) = (TDot(c0, V*cD1i) -TDot(Hc0, HV*cD1i))/m2;
	impsi0_chi(L, m) = TDot(c0, s0L_chi[L](i)).imag();
      }
    }
    */

  }
  
}
void Calc_RPA_Eigen() {
  
  // ==== build RPA Hamiltonian ====
  BVec eig1;
  BMat U1;

  // -- A = T + V_ne + J + K - epsilon0  
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);  
  BMat A("A"); Copy(T1, A);    
  A.Add(1.0, V1);
  B2EInt eri_J_11 = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K_11 = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  AddJ(eri_J_11, c0, irrep0, 1.0, A);
  AddK(eri_K_11, c0, irrep0, 1.0, A);
  A.Add(-E0, S1);

  // -- B = K --
  BMat B("B"); Copy(T1, B); B.SetZero();  
  AddK(eri_K_11, c0, irrep0, 1.0, B);

  // -- H^2= sqrt(sqrt(A-B)(A+B)sqrt(A-B))
  BMat ApB("ApB"); Copy(A, ApB); ApB.Add(+1.0, B);  
  BMat AmB("AmB"); Copy(A, AmB); AmB.Add(-1.0, B);
  BMat sqrt_AmB("sqrt_AmB"); BMatSqrt(AmB, sqrt_AmB);
  BMat H2("H2");  Multi3(sqrt_AmB, ApB, sqrt_AmB, H2);

  // -- solve eigen value problem --
  BMat hpg, hmg;
  BOOST_FOREACH(Irrep i, irreps) {
    SymGenComplexEigenSolver solver(H2(i,i), S1(i,i));
    eig1(i) = solver.eigenvalues();
    int n = eig1(i).size();
    for(int I = 0; I < n; I++) {
      eig1(i)(I) = sqrt(eig1(i)(I));
    }
    U1(i, i) = solver.eigenvectors();
    
    hpg(i, i) = MatrixXcd::Zero(n, n);
    hmg(i, i) = MatrixXcd::Zero(n, n);
    
    for(int I = 0; I < n; I++) {
      for(int j = 0; j < n; j++) {
	for(int k = 0; k < n; k++) {
	  hpg(i,i)(j, I) += sqrt_AmB(i,i)(j, k) * U1(i,i)(k, I) / sqrt(eig1(i)(I));
	}
      }
    }
    
    for(int I = 0; I < n; I++) {
      VectorXcd UI = U1(i, i).col(I);
      VectorXcd AmB_inv_sqrt_UI = sqrt_AmB(i,i).colPivHouseholderQr().solve(UI);
      for(int j = 0; j < n; j++) {
	hmg(i,i)(j, I) = sqrt(eig1(i)(I)) * AmB_inv_sqrt_UI(j);
      }
    }
    
    for(int I = 0; I < n; I++) {
      cout << "norm: "
	   << TDot(hmg(i,i).col(I), hpg(i,i).col(I))
	   << endl;
    }
  }
  
  // -- compute driven term --
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z(); 
  s1(x) = X1i(x, 0) * c0; sD1(x) = DX1i(x, 0) * c0; 
  s1(y) = Y1i(y, 0) * c0; sD1(y) = DY1i(y, 0) * c0; 
  s1(z) = Z1i(z, 0) * c0; sD1(z) = DZ1i(z, 0) * c0;

  // -- transform to excited states basis --
  BVec s1_rpa, sD1_rpa;
  BOOST_FOREACH(Irrep i, irreps) {
    Ctx(hpg(i,i), s1(i),   &s1_rpa(i));
    Ctx(hpg(i,i), sD1(i), &sD1_rpa(i));
  }

  // -- start loop --
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw];
    cout << "w_eV: " << w * au2ev << endl;
    cout << "w_au: " << w << endl;
    cout << "E_eV: " << (w + E0) * au2ev << endl;
    cout << "E_au: " << w + E0 << endl;
    cout << "k_au: " << sqrt(2.0*(w + E0)) << endl;

    // -- calculate alpha --
    // <mu phi0, psi1> = sum_I <phi_0, mu phi_I> <phi_I, mu phi_0>/(w-dEi)
    // <phi_0, mu phi_I> = <phi_0, mu u_i> U_{iI}
    BOOST_FOREACH(Irrep i, irreps) {
      muphi_psi1[i] = 0.0;
      muphi_psi1_v[i] = 0.0;
      for(int I = 0; I < U1(i,i).rows(); I++) {
	muphi_psi1[i]   += s1_rpa(i)(I)*s1_rpa(i)(I)   / (w-eig1(i)(I));
	muphi_psi1_v[i] += sD1_rpa(i)(I)*sD1_rpa(i)(I) / (w-eig1(i)(I));
      } 
    }
    
    CalcMain_alpha(iw);    
    //    CalcMain(iw);
  } 
  
}
void Calc_RPA_Eigen_IVO() {

  PrintTimeStamp("start RPA/IVO", NULL);
  
  // -- STEX Hamiltonian  --
  PrintTimeStamp("STEX", NULL);
  CalcSTVMat(basis1, basis1, &S1, &T1, &V1);
  B2EInt eri_J_11 = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K_11 = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  BMat H_STEX;
  Copy(T1, H_STEX);
  H_STEX.Add(1.0, V1);
  AddJ(eri_J_11, c0, irrep0, 1.0, H_STEX); AddK(eri_K_11, c0, irrep0, 1.0, H_STEX);

  // -- calculate IVO --
  // I : index for IVO
  // i : index for AO
  PrintTimeStamp("IVO", NULL);
  BVec eig_ivo_I;
  BMat U_ivo_iI;
  BOOST_FOREACH(Irrep irr, irreps) {
    SymGenComplexEigenSolver solver(H_STEX(irr,irr), S1(irr,irr));
    eig_ivo_I(irr) = solver.eigenvalues();
    U_ivo_iI(irr, irr) = solver.eigenvectors();
  }

  // -- Build A and B in IVO basis --
  PrintTimeStamp("A and B", NULL);
  BMat A; 
  BMatDiag(eig_ivo_I, &A);
  A.Shift(-E0); 
  BMat B; Copy(A, B); B.SetZero();
  AddK(eri_K_11, c0, irrep0, 1.0, B); 
  BOOST_FOREACH(Irrep irr, irreps) {
    CtAC(U_ivo_iI(irr, irr), B(irr, irr)); 
  }

  // -- build Hamiltonian --
  PrintTimeStamp("H2", NULL);
  BMat ApB("ApB"); Copy(A, ApB); ApB.Add(+1.0, B);  
  BMat AmB("AmB"); Copy(A, AmB); AmB.Add(-1.0, B);
  BMat sqrt_AmB("sqrt_AmB"); BMatSqrt(AmB, sqrt_AmB);
  BMat inv_sqrt_AmB("inv_sqrt_AmB"); Copy(sqrt_AmB, inv_sqrt_AmB);  
  BOOST_FOREACH(Irrep irr, irreps) {
    LinearSolver lin;
    lin.Inv(sqrt_AmB(irr, irr), &inv_sqrt_AmB(irr, irr));
  }  
  BMat H2("H2");  Multi3(sqrt_AmB, ApB, sqrt_AmB, H2);

  // -- diagonalize Hamiltonian --
  PrintTimeStamp("diag H2", NULL);
  BVec w2s_rpa;
  BMat YpZ, YmZ;
  BOOST_FOREACH(Irrep irr, irreps) {
    int n = H2(irr, irr).rows();
    MatrixXcd S = MatrixXcd::Identity(n, n);
    SymGenComplexEigenSolver solver(H2(irr,irr), S);
    w2s_rpa(irr) = solver.eigenvalues();
    const MatrixXcd& U = solver.eigenvectors();
    YpZ(irr, irr) = MatrixXcd::Zero(n, n);
    YmZ(irr, irr) = MatrixXcd::Zero(n, n);
    for(int i = 0; i < n; i++) {
      for(int I = 0; I < n; I++) {
	dcomplex w = sqrt(w2s_rpa(irr)(I));
	dcomplex acc1(0);
	dcomplex acc2(0);
	for(int k = 0; k < n; k++) {
	  acc1 += 1.0/w * sqrt_AmB(irr,irr)(i, k) * U(k, I);
	  acc2 += w     * inv_sqrt_AmB(irr,irr)(i, k) * U(k, I);
	}
	YpZ(irr, irr)(i, I) = acc1;
	YmZ(irr, irr)(i, I) = acc1;
      }
    }
  }

  // -- check normalization --
  BOOST_FOREACH(Irrep irr, irreps) {
    VectorXcd Y = 0.5 * (YpZ(irr, irr) + YmZ(irr, irr));
    VectorXcd Z = 0.5 * (YpZ(irr, irr) - YmZ(irr, irr));
    cout << TDot(Y, Y) - TDot(Z, Z) << endl;
  }
  
}

void Solve_0th(map<int, BMat> *U0, map<int, BVec> *eig0,	       
	       map<int, BVec> *r, map<int, BVec> *d,
	       map<int, BVec> *s) {

  PrintTimeStamp("Solve_0th", NULL);
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  
  // -- HF for ground state --
  int ng = basis0->size_basis_isym(0);
  BMat U_g;
  U_g(0, 0) = MatrixXcd::Zero(ng, ng);
  U_g(0, 0).col(0) = c0;  

  BOOST_FOREACH(int L, Ls) {
    cout << format("L = %d\n") % L;    
    SymGTOs b0L = basis_psi0_L[L];

    // -- solve eigen system --
    BMat S, H, V;
    CalcSTVMat(b0L, b0L, &S, &H, &V);
    H.Add(1.0, V);
    BMatEigenSolve(H, S, &(*U0)[L], &(*eig0)[L]);

    // -- dipole moment --
    BMat X, DX, Y, DY, Z, DZ, r_ao, d_ao, r_mat, d_mat;
    CalcDipMat(b0L, basis0, &X, &Y, &Z, &DX, &DY, &DZ);
    r_ao(x,0) = X(x,0); d_ao(x,0) = DX(x,0);
    r_ao(y,0) = Y(y,0); d_ao(y,0) = DY(y,0);
    r_ao(z,0) = Z(z,0); d_ao(z,0) = DZ(z,0);
    BMatCtAD((*U0)[L], U_g, r_ao, &r_mat);
    BMatCtAD((*U0)[L], U_g, d_ao, &d_mat);
    BMatCol(r_mat, 0, 0, &(*r)[L]);
    BMatCol(d_mat, 0, 0, &(*d)[L]);

    // -- artificial driven term --
    BMat s_ao; 
    CalcSMat(b0L, basis_chi0_L[L], &s_ao);
    BMat Us;
    Us(x,x)=MatrixXcd::Ones(1,1);
    Us(y,y)=MatrixXcd::Ones(1,1);
    Us(z,z)=MatrixXcd::Ones(1,1);
    BMat s_mat;
    BMatCtAD((*U0)[L], Us, s_ao, &s_mat);
    BOOST_FOREACH(Irrep irr, irreps) {
      (*s)[L](irr) = s_mat(irr, irr).col(0);
    }
    
  }
  
  
}
void Solve_1st_One(BMat *U1, BVec *ws, BVec *r, BVec *d) {
  
  PrintTimeStamp("Solve_Psi1_One", NULL);
  // --  ground state --
  int ng = basis0->size_basis_isym(0);
  BMat U_g;
  U_g(0, 0) = MatrixXcd::Zero(ng, ng);
  U_g(0, 0).col(0) = c0;

  // -- AO basis --
  BMat S, H, V;
  CalcSTVMat(basis1, basis1, &S, &H, &V);
  H.Add(1.0, V);

  BMatEigenSolve(H, S, U1, ws);
  ws->Shift(-E0);

  // -- dipole --
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i, r_ao, d_ao;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  r_ao(x,0) = X1i(x,0); d_ao(x,0) = DX1i(x,0);
  r_ao(y,0) = Y1i(y,0); d_ao(y,0) = DY1i(y,0);
  r_ao(z,0) = Z1i(z,0); d_ao(z,0) = DZ1i(z,0);
  BMat rmat, dmat;
  BMatCtAD(*U1, U_g, r_ao, &rmat);
  BMatCtAD(*U1, U_g, d_ao, &dmat);

  BOOST_FOREACH(Irrep irr, irreps) {
    (*r)[irr] = rmat(irr, 0).col(0);
    (*d)[irr] = dmat(irr, 0).col(0);
  }
  
}
void Solve_1st_STEX(BMat *U1, BVec *w1, BVec *r1, BVec *d1) {
  
  PrintTimeStamp("Solve_1st_STEX", NULL);
  
  // --  ground state --
  int ng = basis0->size_basis_isym(0);
  BMat U_g;
  U_g(0, 0) = MatrixXcd::Zero(ng, ng);
  U_g(0, 0).col(0) = c0;

  // -- AO basis --
  BMat S, T, V;
  CalcSTVMat(basis1, basis1, &S, &T, &V);  
  B2EInt eri_J = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  BMat K; Copy(S, K); K.SetZero();
  AddK(eri_K, c0, 0, 1.0, K);
  BMat H_IVO;
  CalcSTEXHamiltonian(T, V, 0, c0, eri_J, eri_K, &H_IVO);

  BMatEigenSolve(H_IVO, S, U1, w1);
  w1->Shift(-E0);

  // -- dipole -
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i, r_ao, d_ao;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  r_ao(x,0) = X1i(x,0); d_ao(x,0) = DX1i(x,0);
  r_ao(y,0) = Y1i(y,0); d_ao(y,0) = DY1i(y,0);
  r_ao(z,0) = Z1i(z,0); d_ao(z,0) = DZ1i(z,0);
  BMat rmat, dmat;
  BMatCtAD(*U1, U_g, r_ao, &rmat);
  BMatCtAD(*U1, U_g, d_ao, &dmat);

  BOOST_FOREACH(Irrep irr, irreps) {
    (*r1)[irr] = rmat(irr, 0).col(0);
    (*d1)[irr] = dmat(irr, 0).col(0);
  }
    
}
void Solve_1st_RPA(CalcRPA *rpa, BMat* U_HF, BVec *r_RPA, BVec *d_RPA) {
  
  PrintTimeStamp("Solve_RPA", NULL);
  // -- HF for ground state --
  int ng = basis0->size_basis_isym(0);
  BMat U_g;
  U_g(0, 0) = MatrixXcd::Zero(ng, ng);
  U_g(0, 0).col(0) = c0;
  
  // -- AO basis --
  BMat S, T, V;
  CalcSTVMat(basis1, basis1, &S, &T, &V);  
  B2EInt eri_J = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  BMat K; Copy(S, K); K.SetZero();
  AddK(eri_K, c0, 0, 1.0, K);
  BMat H_IVO;
  CalcSTEXHamiltonian(T, V, 0, c0, eri_J, eri_K, &H_IVO);
  
  // -- HF --
  BMat H_HF;
  CalcFock(T, V, 0, c0, eri_J, eri_K, &H_HF);
  BVec eig_HF;
  BMatEigenSolve(H_HF, S, U_HF, &eig_HF);

  // -- RPA --
  rpa->CalcH_HF(H_IVO, E0, K, eig_HF, *U_HF);
  rpa->SolveEigen();
  
  // -- dipole --
  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i, r_ao, d_ao;
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  r_ao(x,0) = X1i(x,0); d_ao(x,0) = DX1i(x,0);
  r_ao(y,0) = Y1i(y,0); d_ao(y,0) = DY1i(y,0);
  r_ao(z,0) = Z1i(z,0); d_ao(z,0) = DZ1i(z,0);
  BMat r_HF, d_HF;
  BMatCtAD(*U_HF, U_g, r_ao, &r_HF);
  BMatCtAD(*U_HF, U_g, d_ao, &d_HF);
  rpa->CalcOneInt(r_HF, r_RPA, true);
  rpa->CalcOneInt(d_HF, d_RPA, false);
  
}
void CalcOStrength(const BVec& w, const BVec& r, const BVec& d, dcomplex *os_r, dcomplex *os_d) {

  *os_r = 0.0;
  *os_d = 0.0;
  BOOST_FOREACH(Irrep irr, irreps) {
    int n = w[irr].size();
    for(int I = 0; I < n; I++) {
      dcomplex wI = w(irr)(I);
      dcomplex rI = r(irr)(I);
      dcomplex dI = d(irr)(I);
      *os_r += 2.0/3.0*wI*rI*rI;
      *os_d += 2.0/3.0*wI*dI*dI;
    }
  }
}
void Calc_V_one(const map<int, BMat>& U0, const BMat& U1,
		map<int, BMat> *v, map<int, BMat> *u) {
  PrintTimeStamp("Calc_V_one", NULL);
  BOOST_FOREACH(int L, Ls) {

    SymGTOs b0L =  basis_psi0_L[L];    
    BMat V1, V0;
    CalcVMat(b0L,  mole,  basis1, &V1);
    CalcVMat(b0L,  mole0, basis1, &V0);
    V1.Add(-1.0, V0);
    BMatCtAD(U0.at(L), U1, V1, &((*v)[L]));

    SymGTOs cb0L =  basis_c_psi0_L[L];
    BMat CV1, CV0;
    CalcVMat(cb0L, mole,  basis1, &CV1);    
    CalcVMat(cb0L, mole0, basis1, &CV0);
    CV1.Add(-1.0, CV0);
    BMatCaAD(U0.at(L), U1, CV1, &((*u)[L]));
  
  }  
  
}
void Calc_V_STEX(const map<int, BMat>& U0, const BMat& U1,
		 map<int, BMat> *v, map<int, BMat> *u) {
  
  PrintTimeStamp("Calc_V_STEX", NULL);
  BOOST_FOREACH(int L, Ls) {

    SymGTOs b0L =  basis_psi0_L[L];    
    BMat V1, V0;
    CalcVMat(b0L,  mole,  basis1, &V1);
    CalcVMat(b0L,  mole0, basis1, &V0);    
    B2EInt eri_J = CalcERI(b0L, basis1, basis0, basis0, eri_method);
    B2EInt eri_K = CalcERI(b0L, basis0, basis0, basis1, eri_method);
    V1.Add(-1.0, V0);
    AddJ(eri_J, c0, irrep0, 1.0, V1);
    AddK(eri_K, c0, irrep0, 1.0, V1);
    BMatCtAD(U0.at(L), U1, V1, &((*v)[L]));

    SymGTOs cb0L =  basis_c_psi0_L[L];
    BMat CV1, CV0;
    CalcVMat(cb0L, mole,  basis1, &CV1);    
    CalcVMat(cb0L, mole0, basis1, &CV0);
    B2EInt c_eri_J = CalcERI(cb0L, basis1, basis0, basis0, eri_method);
    B2EInt c_eri_K = CalcERI(cb0L, basis0, basis0, basis1, eri_method);
    CV1.Add(-1.0, CV0);
    AddJ(c_eri_J, c0, irrep0, 1.0, CV1);
    AddK(c_eri_K, c0, irrep0, 1.0, CV1);
    BMatCaAD(U0.at(L), U1, CV1, &((*u)[L]));
  
  }
  
}
void Calc_V_RPA_STEX( const map<int, BMat>& U0, 
		const CalcRPA& rpa, const BMat U1_HF, 
		map<int, BMat> *v, map<int, BMat> *u) {

  PrintTimeStamp("Calc_V_RPA_STEX", NULL);
  BOOST_FOREACH(int L, Ls) {

    SymGTOs b0L =  basis_psi0_L[L];    
    BMat V1, V0;
    CalcVMat(b0L,  mole,  basis1, &V1);
    CalcVMat(b0L,  mole0, basis1, &V0);    
    B2EInt eri_J = CalcERI(b0L, basis1, basis0, basis0, eri_method);
    B2EInt eri_K = CalcERI(b0L, basis0, basis0, basis1, eri_method);
    V1.Add(-1.0, V0);
    AddJ(eri_J, c0, irrep0, 1.0, V1);
    AddK(eri_K, c0, irrep0, 1.0, V1);
    BMat VV;
    BMatCtAD(U0.at(L), U1_HF, V1, &VV);
    rpa.CalcOneIntRight(VV, &(*v)[L], true, 1.0, 1.0);

    SymGTOs cb0L =  basis_c_psi0_L[L];
    BMat CV1, CV0;
    CalcVMat(cb0L, mole,  basis1, &CV1);    
    CalcVMat(cb0L, mole0, basis1, &CV0);
    B2EInt c_eri_J = CalcERI(cb0L, basis1, basis0, basis0, eri_method);
    B2EInt c_eri_K = CalcERI(cb0L, basis0, basis0, basis1, eri_method);
    CV1.Add(-1.0, CV0);
    AddJ(c_eri_J, c0, irrep0, 1.0, CV1);
    AddK(c_eri_K, c0, irrep0, 1.0, CV1);
    BMat UU;
    BMatCaAD(U0.at(L), U1_HF, CV1, &UU);
    rpa.CalcOneIntRight(UU, &(*u)[L], true, 1.0, 1.0);
  
  }
  
}
void Calc_V_RPA_eigen(const map<int, BMat>& U0, const map<int, BVec>& eig0,		      
		      const CalcRPA& rpa, const BMat U1_HF, 
		      map<int, BMat> *v, map<int, BMat> *u) {

  PrintTimeStamp("Calc_V_RPA_eigen", NULL);
  BOOST_FOREACH(int L, Ls) {
    
    SymGTOs b0L = basis_psi0_L[L];
    BMat S_ao;
    CalcSMat(b0L, basis1, &S_ao);
    BMat S_eig;
    BMatCtAD(U0.at(L), U1_HF, S_ao, &S_eig);
    BMat vv;
    rpa.CalcOneIntRight(S_eig, &vv, true);

    SymGTOs cb0L = basis_c_psi0_L[L];
    BMat C_S_ao;
    CalcSMat(cb0L, basis1, &C_S_ao);
    BMat C_S_eig;
    BMatCaAD(U0.at(L), U1_HF, C_S_ao, &C_S_eig);
    BMat uu;
    rpa.CalcOneIntRight(C_S_eig, &uu, true);

    BOOST_FOREACH(Irrep irr, irreps) {
      int n0 = b0L->size_basis_isym(irr);
      int nI = rpa.w[irr].size();
      for(int i = 0; i < n0; i++)
	for(int I = 0; I < nI; I++) {
	  dcomplex wi0 = eig0.at(L)(irr)(i) - E0;
	  dcomplex wI1 = rpa.w(irr)(I);
	  vv(irr,irr)(i, I) *= wI1-wi0;
	  uu(irr,irr)(i, I) *= wI1-conj(wi0);
	}
    }
    (*v)[L] = vv;
    (*u)[L] = uu;
  }

}
void Calc_V_RPA_RPA(const map<int, BMat>& U0, 
		    const CalcRPA& rpa, const BMat U1_HF, 
		    map<int, BMat> *v, map<int, BMat> *u) {
  
  PrintTimeStamp("Calc_V_RPA_RPA", NULL);
  BOOST_FOREACH(int L, Ls) {

    SymGTOs b0L =  basis_psi0_L[L];    
    BMat V1, V0;
    CalcVMat(b0L,  mole,  basis1, &V1);
    CalcVMat(b0L,  mole0, basis1, &V0);    
    B2EInt eri_J = CalcERI(b0L, basis1, basis0, basis0, eri_method);
    B2EInt eri_K = CalcERI(b0L, basis0, basis0, basis1, eri_method);
    V1.Add(-1.0, V0);
    AddJ(eri_J, c0, irrep0, 1.0, V1);
    AddK(eri_K, c0, irrep0, 1.0, V1);
    BMat V_Y; BMatCtAD(U0.at(L), U1_HF, V1, &V_Y);
    BMat V_Z; Copy(V1, V_Z); V_Z.SetZero();
    AddK(eri_K, c0, irrep0, 1.0, V_Z);
    BMat VV_Y, VV_Z;
    rpa.CalcOneIntRight(V_Y, &VV_Y, true, 1.0, 0.0);
    rpa.CalcOneIntRight(V_Z, &VV_Z, true, 0.0, 1.0);
    Copy(VV_Y, (*v)[L]);
    v->at(L).Add(1.0, VV_Z);

    SymGTOs cb0L =  basis_c_psi0_L[L];
    BMat CV1, CV0;
    CalcVMat(cb0L, mole,  basis1, &CV1);    
    CalcVMat(cb0L, mole0, basis1, &CV0);
    B2EInt c_eri_J = CalcERI(cb0L, basis1, basis0, basis0, eri_method);
    B2EInt c_eri_K = CalcERI(cb0L, basis0, basis0, basis1, eri_method);
    CV1.Add(-1.0, CV0);
    AddJ(c_eri_J, c0, irrep0, 1.0, CV1);
    AddK(c_eri_K, c0, irrep0, 1.0, CV1);
    BMat CV_Y; BMatCaAD(U0.at(L), U1_HF, CV1, &CV_Y);
    BMat CV_Z; Copy(CV1, CV_Z); CV_Z.SetZero();
    AddK(c_eri_K, c0, irrep0, 1.0, CV_Z);
    BMat CVV_Y, CVV_Z;
    rpa.CalcOneIntRight(CV_Y, &CVV_Y, true, 1.0, 0.0);
    rpa.CalcOneIntRight(CV_Z, &CVV_Z, true, 0.0, 1.0);
    Copy(CVV_Y, (*u)[L]);
    u->at(L).Add(1.0, CVV_Z);    
  }
  
}
void CalcBraket_eigen(const map<int, BVec>& eig0,
		      const map<int, BVec>& r0, const map<int, BVec>& d0,
		      const map<int, BVec>& s0,
		      const BVec& w1s, const BVec& r1, const BVec& d1,
		      const map<int, BMat>& V, const map<int, BMat>& U,
		      dcomplex w) {

  PrintTimeStamp("CalcBraket_eigen", NULL);
  dcomplex os_r, os_d;
  CalcOStrength(w1s, r1, d1, &os_r, &os_d);
  cout << endl;
  cout << "oscillator strength" <<endl;
  cout << os_r << ", " << os_d << endl;
  BOOST_FOREACH(Irrep irr, irreps) {
    muphi_psi1[irr] = 0.0;
    muphi_psi1_v[irr] = 0.0;
    int n = w1s[irr].size();
    for(int I = 0; I < n; I++) {
      dcomplex dl = r1(irr)(I);
      dcomplex dv = d1(irr)(I);
      dcomplex wi = w1s(irr)[I];
      muphi_psi1[irr]   += dl*dl/(w-wi) / 3.0 * (1.0*ne);
      muphi_psi1_v[irr] += dv*dv/(w-wi) / 3.0 * (1.0*ne);
    }
  }

  dcomplex i2(0, 2);
  int ms[3]; ms[0]=+1; ms[1]=-1; ms[2]=0;
  BOOST_FOREACH(int L, Ls) {    
    for(int im = 0; im < 3; im++) {
      Irrep irr = irreps[im];
      int m = ms[im];
      impsi0_muphi(L, m) = 0.0; impsi0_muphi_v(L, m) = 0.0;
      impsi0_chi(L,m) = 0.0;
      impsi0_v_psi1(L,m) = 0.0; impsi0_v_psi1_v(L,m) = 0.0;
      int n0 = eig0.at(L)[irr].size();
      int n1 = w1s[irr].size();
      cout << "(L,m)" << L << m << endl;
      for(int i = 0; i < n0; i++) {
	dcomplex s0i(s0.at(L)(irr)(i)); 
	dcomplex c_s0i(conj(s0i));  
	dcomplex r0i(r0.at(L)(irr)(i));
	dcomplex d0i(d0.at(L)(irr)(i));
	dcomplex de0i = E0+w-eig0.at(L)(irr)(i);
	dcomplex c_de0i = conj(de0i); 
	impsi0_muphi(L,m)   += (s0i*r0i/de0i).imag(); 
	impsi0_muphi_v(L,m) += (s0i*d0i/de0i).imag(); 
	impsi0_chi(L,m)     += (s0i*s0i/de0i).imag(); 
	for(int I = 0; I < n1; I++) {
	  dcomplex ViI = V.at(L)(irr,irr)(i,I); 
	  dcomplex UiI = U.at(L)(irr,irr)(i,I); 
	  dcomplex rI = r1(irr)(I); 
	  dcomplex dI = d1(irr)(I); 
	  dcomplex deI = w-w1s(irr)(I);
	  impsi0_v_psi1(L,m) +=
	    (+s0i*ViI*rI/(de0i*deI) - c_s0i*UiI*rI/(c_de0i*deI))/i2;
	  impsi0_v_psi1_v(L,m)  +=
	    (+s0i*ViI*dI/(de0i*deI) - c_s0i*UiI*dI/(c_de0i*deI))/i2;
	}
      }
    }
  }
}
void PrintWloop(dcomplex w) {
  cout << "calculation at...." << endl;
  cout << "w_eV: " << w * au2ev << endl;
  cout << "w_au: " << w << endl;
  cout << "E_eV: " << (w + E0) * au2ev << endl;
  cout << "E_au: " << w + E0 << endl;
  cout << "k_au: " << sqrt(2.0*(w + E0)) << endl;
}
void Calc_One_Eigen() {

  map<int, BMat> U0;
  map<int, BVec> eig0, r0, d0, s0;
  Solve_0th(&U0, &eig0, &r0, &d0, &s0);
  
  BMat U1;
  BVec w1, r1, d1;
  Solve_1st_One(&U1, &w1, &r1, &d1);

  map<int, BMat> V, U;
  Calc_V_one(U0, U1, &V, &U);    

  // -- w loop --
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw];
    PrintWloop(w);
    CalcBraket_eigen(eig0, r0, d0, s0,		   
		     w1, r1, d1, V, U, w);
    CalcMain_alpha(iw);
    CalcMain(iw);
  }
  
}
void Calc_STEX_Eigen() {

  map<int, BMat> U0;
  map<int, BVec> eig0, r0, d0, s0;
  Solve_0th(&U0, &eig0, &r0, &d0, &s0);

  BMat U1;
  BVec w1, r1, d1;
  Solve_1st_STEX(&U1, &w1, &r1, &d1);

  map<int, BMat> V, U;
  Calc_V_STEX(U0, U1, &V, &U);

  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw];
    PrintWloop(w);
    CalcBraket_eigen(eig0, r0, d0, s0,
		     w1, r1, d1, V, U, w);
    CalcMain_alpha(iw);
    CalcMain(iw);
  }  
}
void Calc_RPA_HF() {

  string eval_v = "STEX";
  if(obj.find("RPA") != obj.end()) {
    object& obj_rpa = obj["RPA"].get<object>();
    eval_v = ReadJsonWithDefault<string>(obj_rpa, "eval_v", "STEX");    
  }
  cout << format("eval_v: %s\n") % eval_v;
  
  map<int, BMat> U0;
  map<int, BVec> eig0, r0, d0, s0;
  Solve_0th(&U0, &eig0, &r0, &d0, &s0);
  
  CalcRPA rpa;
  BMat U1_HF;
  BVec r1_RPA, d1_RPA;
  Solve_1st_RPA(&rpa, &U1_HF, &r1_RPA, &d1_RPA);

  cout << "RPA results" << endl;
  cout << "Y:" << endl;
  cout << rpa.Y << endl;
  cout << "Z:" << endl;
  cout << rpa.Z << endl;
  
  map<int, BMat> V, U;
  if(eval_v == "STEX") {
    Calc_V_RPA_STEX(U0, rpa, U1_HF, &V, &U);
  } else if(eval_v == "eigen") {
    Calc_V_RPA_eigen(U0, eig0, rpa, U1_HF, &V, &U);
  } else if(eval_v == "RPA") {
    Calc_V_RPA_RPA(U0, rpa, U1_HF, &V, &U);
  } else {
    THROW_ERROR("unsupported RPA.eval_v");
  }

  // -- w loop --
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw];
    PrintWloop(w);
    CalcBraket_eigen(eig0, r0, d0, s0,		   
		     rpa.w, r1_RPA, d1_RPA, V, U, w);
    CalcMain_alpha(iw);
    CalcMain(iw);
  }

  PrintOut();
}

// ==== propagator ====
void Calc_AB_HF(SymGTOs g, BMat *U, BMat *A, BMat *B) {

  BMat S, T, V;
  CalcSTVMat(g, g, &S, &T, &V);  
  B2EInt eri_J = CalcERI(g, g, basis0, basis0, eri_method);
  B2EInt eri_K = CalcERI(g, basis0, basis0, g, eri_method);
  BMat K; Copy(S, K); K.SetZero(); AddK(eri_K, c0, 0, 1.0, K);  
  BMat H_IVO; CalcSTEXHamiltonian(T, V, 0, c0, eri_J, eri_K, &H_IVO);
  
  // -- HF --
  BMat H_HF; CalcFock(T, V, 0, c0, eri_J, eri_K, &H_HF);  
  BVec eig_HF; BMatEigenSolve(H_HF, S, U, &eig_HF);

  // -- A, B (HF basis) --
  BMatCtAC(*U, H_IVO, A); A->Shift(-E0);
  BMatCtAC(*U, K, B);

}
void CalcDipMat(SymGTOs g1, SymGTOs g0, const BMat& U1, const BMat U0, BMat *r, BMat *d) {
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z();  
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i, r_ao, d_ao;
  CalcDipMat(g1, g0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  r_ao(x,0) = X1i(x,0); d_ao(x,0) = DX1i(x,0);
  r_ao(y,0) = Y1i(y,0); d_ao(y,0) = DY1i(y,0);
  r_ao(z,0) = Z1i(z,0); d_ao(z,0) = DZ1i(z,0);
  
  BMatCtAD(U1, U0, r_ao, r);
  BMatCtAD(U1, U0, d_ao, d);
}
void CalcDipMat(SymGTOs g1, SymGTOs g0, const BMat U0, BMat *r, BMat *d) {
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z();  
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i, r_ao, d_ao;
  CalcDipMat(g1, g0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  r_ao(x,0) = X1i(x,0); d_ao(x,0) = DX1i(x,0);
  r_ao(y,0) = Y1i(y,0); d_ao(y,0) = DY1i(y,0);
  r_ao(z,0) = Z1i(z,0); d_ao(z,0) = DZ1i(z,0);

  Multi(r_ao, U0, *r);
  Multi(d_ao, U0, *d);

}
void CalcPMat(SymGTOs g0, SymGTOs g1, const BMat& U_1, BMat *pe, BMat *pd) {
  
  /**
     Pe_{ba} = (E_{0b}V|e)_a 
     .       = <Phi_0 | E_{0b}V E_{a0} | Phi_0> 
     .       = A_{ba} - H^0_{ba}
     Pd_{ba} = (E_{0b}V|d)_a 
     .       = <Phi_0 | E_{0a}E_{0b}V  | Phi_0> 
     .       = K_{ba}
   */

  BMat V1, V0;
  CalcVMat(g0, mole,  g1, &V1);
  CalcVMat(g0, mole0, g1, &V0);
  B2EInt eri_J = CalcERI(g0, g1, basis0, basis0, eri_method);
  B2EInt eri_K = CalcERI(g0, basis0, basis0, g1, eri_method);
  BMat K; Copy(V1, K); K.SetZero(); AddK(eri_K, c0, 0, 1.0, K);  
  V1.Add(-1.0, V0);
  AddJ(eri_J, c0, irrep0, 1.0, V1);
  AddK(eri_K, c0, irrep0, 1.0, V1);
  
  Multi(V1, U_1, *pe); 
  Multi(K,  U_1, *pd); 
  
}
void CalcRPAMat(dcomplex w, const BMat& A, const BMat& B, BMat *M) {

  BOOST_FOREACH(Irrep irr, irreps) {
    int n1 = basis1->size_basis_isym(irr);
    MatrixXcd& MM = (*M)(irr, irr);
    const MatrixXcd& AA = A(irr, irr);
    const MatrixXcd& BB = B(irr, irr);
    if(MM.rows() != 2*n1 || MM.cols() != 2*n1) {
      MM = MatrixXcd::Zero(2*n1, 2*n1);
    }
    
    for(int a = 0; a < n1; a++) {
      for(int b = 0; b < n1; b++) {
	if(a==b) {
	  MM(a, b)       = +w -AA(a, b);
	  MM(a+n1, b+n1) = -w -AA(a, b);
	} else {
	  MM(a, b)       = -AA(a, b);
	  MM(a+n1, b+n1) = -AA(a, b);
	}
	MM(a, b+n1) = -BB(a, b);
	MM(a+n1, b) = -BB(a, b);
      }
    }
  }
  
}
void CalcQVec(const BMat& r, const BMat &d, BVec *qvec_l, BVec *qvec_v) {

  /**
     calculate [(e_jb|mu), (d_jb|mu)]
     (e_jb|mu) = +(b|mu|j)
     (d_jb|mu) = -(j|mu|b)

     Inputs 
     -------
     r(a, i) : <a|mu_length|i>
     d(a, i) : <a|mu_velocity|i>
   */

  BOOST_FOREACH(Irrep irr, irreps) {
    int n1 = basis1->size_basis_isym(irr);
    const MatrixXcd& rr = r(irr, irrep0);
    const MatrixXcd& dd = d(irr, irrep0);
    VectorXcd& ll = (*qvec_l)(irr);
    VectorXcd& vv = (*qvec_v)(irr);
    if(ll.size() != 2*n1)
      ll = VectorXcd::Zero(2*n1);
    if(vv.size() != 2*n1)
      vv = VectorXcd::Zero(2*n1);

    for(int a = 0; a < n1; a++) {
      ll(a)    = rr(a);
      ll(a+n1) = rr(a);
      vv(a)    = dd(a);
      vv(a+n1) = -dd(a);
      /*
      ll(a)    = rr(a);
      ll(a+n1) = -rr(a);
      vv(a)    = dd(a);
      vv(a+n1) = -(-dd(a));
      */
    }
  }
  
}
void CalcPVec(const BMat& r, const BMat& d, BVec *qvec_r, BVec *qvec_d) {
  /**
     calculate [(mu|e_jb), (mu|d_jb)]
     (mu|e_jb) = (j|mu|b)
     (mu|d_jb) = -(b|mu|j)

     Inputs 
     -------
     r(a, i) : <a|mu_length|i>
     d(a, i) : <a|mu_velocity|i>
   */

  BOOST_FOREACH(Irrep irr, irreps) {
    int n1 = basis1->size_basis_isym(irr);
    const MatrixXcd& rr = r(irr, irrep0);
    const MatrixXcd& dd = d(irr, irrep0);
    VectorXcd& ll = (*qvec_r)(irr);
    VectorXcd& vv = (*qvec_d)(irr);
    if(ll.size() != 2*n1)
      ll = VectorXcd::Zero(2*n1);
    if(vv.size() != 2*n1)
      vv = VectorXcd::Zero(2*n1);

    for(int a = 0; a < n1; a++) {
      /*
      ll(a)    = rr(a);
      ll(a+n1) = -rr(a);
      vv(a)    = dd(a);
      vv(a+n1) = dd(a);
      */
      ll(a)    = rr(a);
      ll(a+n1) = rr(a);
      vv(a)    = dd(a);
      vv(a+n1) = -dd(a);
    }

  }
 
}
void CalcPVec(const map<int, BMat>& Pe, const map<int, BMat>& Pd, const map<int, BVec>& C0,
	      map<int, BVec>* pvec) {

  BOOST_FOREACH(Irrep irr, irreps) {
    int n1 = basis1->size_basis_isym(irr);
    BOOST_FOREACH(int L, Ls) {
      int n0 = basis_psi0_L[L]->size_basis_isym(irr);
      (*pvec)[L](irr) = VectorXcd::Zero(2*n1);
      VectorXcd& P = (*pvec)[L](irr);
      for(int a = 0; a < n1; a++) {
	for(int p = 0; p < n0; p++) {
	  P(a)    += C0.at(L)(irr)(p) * Pe.at(L)(irr,irr)(p, a);
	  P(a+n1) += C0.at(L)(irr)(p) * Pd.at(L)(irr,irr)(p, a);
	}
      }
    }
  }
  
}
void Calc_Propagator_RPA() {
  PrintTimeStamp("propagator_RPA", NULL);
  
  // ---- Ground State ----
  int ng = basis0->size_basis_isym(0);
  BMat U_g;
  U_g(0, 0) = MatrixXcd::Zero(ng, ng);
  U_g(0, 0).col(0) = c0;

  // ---- 0 th psi ----  
  map<int, BMat> S0, T0, V0, r_0g, d_0g, s_0;
  BOOST_FOREACH(int L, Ls) {
    SymGTOs g0L = basis_psi0_L[L];
    CalcSTVMat(basis_psi0_L[L], basis_psi0_L[L], &S0[L], &T0[L], &V0[L]);
    CalcDipMat(basis_psi0_L[L], basis0,          U_g, &r_0g[L], &d_0g[L]);
    CalcSMat(basis_psi0_L[L],   basis_chi0_L[L], &s_0[L]);    
  }

  // ---- 1st psi ----
  BMat U_HF_1, A_1, B_1, r_1g, d_1g;
  Calc_AB_HF(basis1, &U_HF_1, &A_1, &B_1);
  CalcDipMat(basis1, basis0, U_HF_1, U_g, &r_1g, &d_1g);
  
  // ---- 0th psi <-> 1st psi (PVec) ----
  map<int, BMat> Pe, Pd, C_Pe, C_Pd;
  BOOST_FOREACH(int L, Ls) {
    CalcPMat(basis_psi0_L[L],   basis1, U_HF_1, &Pe[L],   &Pd[L]);
    CalcPMat(basis_c_psi0_L[L], basis1, U_HF_1, &C_Pe[L], &C_Pd[L]);
  }

  // ---- w loop ----
  PrintTimeStamp("start loop", NULL);
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    dcomplex w = w_list[iw];
    dcomplex ene = E0+w;
    PrintWloop(w);

    // -- 1st psi --
    BMat M; CalcRPAMat(w, A_1, B_1, &M);
    BVec QVec_l, QVec_v; CalcQVec(r_1g, d_1g, &QVec_l, &QVec_v);
    BVec PVec_l, PVec_v; CalcPVec(r_1g, d_1g, &PVec_l, &PVec_v);
    BVec X_l, X_v;
    BOOST_FOREACH(Irrep irr, irreps) {
      ColPivHouseholderQR<MatrixXcd> piv = M(irr, irr).colPivHouseholderQr();
      X_l(irr) = piv.solve(QVec_l(irr)); X_v(irr) = piv.solve(QVec_v(irr));      
    }
    
    // -- 0th psi --
    map<int, BVec> C0, C_C0;
    BOOST_FOREACH(Irrep irr, irreps) {
      BOOST_FOREACH(int L, Ls) {	
	MatrixXcd L0 =  ene*S0[L](irr,irr) - T0[L](irr,irr) - V0[L](irr,irr);
	C0[L](irr)   = L0.colPivHouseholderQr().solve(s_0[L](irr,irr).col(0));
	C_C0[L](irr) = C0[L](irr).conjugate();
      }
    }

    // -- 0th <-> 1st --
    map<int, BVec> PVec, C_PVec;
    CalcPVec(Pe, Pd, C0, &PVec); CalcPVec(C_Pe, C_Pd, C_C0, &C_PVec);

    for(int im = 0; im < 3; im++) {
      Irrep irr(irreps[im]); int m(ms_xyz[im]);
      // -- ground <-> 1st --
      muphi_psi1[irr]   = TDot(PVec_l(irr), X_l(irr)) / 3.0 * (1.0*ne);
      muphi_psi1_v[irr] = TDot(PVec_v(irr), X_v(irr)) / 3.0 * (1.0*ne);

      BOOST_FOREACH(int L, Ls) {
	// -- 0th <-> ground --
	dcomplex i2(0, 2);
	impsi0_muphi(L, m)   = TDot(C0[L](irr), r_0g[L](irr, irrep0).col(0)).imag();
	impsi0_muphi_v(L, m) = TDot(C0[L](irr), d_0g[L](irr, irrep0).col(0)).imag();
	impsi0_v_psi1(L, m) = (TDot(PVec[L](irr), X_l(irr)) - TDot(C_PVec[L](irr), X_l(irr)))/i2;
	impsi0_v_psi1_v(L, m) = (TDot(PVec[L](irr), X_v(irr)) - TDot(C_PVec[L](irr), X_v(irr)))/i2;
	impsi0_chi(L, m) = TDot(C0[L](irr), s_0[L](irr,irr).col(0)).imag();
      }
    }    
    CalcMain_alpha(iw);
    CalcMain(iw);    
  }
  PrintOut();
}
void Calc_Propagator_RPA_eigen() {
  PrintTimeStamp("propagator_RPA_eigen", NULL);

  int ng = basis0->size_basis_isym(0);
  BMat U_g; U_g(0, 0)=MatrixXcd::Zero(ng,ng); U_g(0,0).col(0) = c0;

  BMat U_HF, A, B, r, d;
  Calc_AB_HF(basis1, &U_HF, &A, &B);
  CalcDipMat(basis1, basis0, U_HF, U_g, &r, &d);

  CalcRPA rpa; rpa.A = A; rpa.B = B; rpa.CalcH(); rpa.SolveEigen();
  BMat Y = rpa.Y;
  BMat Z = rpa.Z;

  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    dcomplex w = w_list[iw];
    PrintWloop(w);

    for(int im = 0; im < 3; im++) {
      Irrep irr(irreps[im]); 
      dcomplex acc_l(0), acc_v(0);
      int nI = rpa.w(irr).size();
      for(int I = 0; I < nI; I++) {
	dcomplex wI = rpa.w(irr)(I);
	dcomplex rYI = TDot(Y(irr,irr).col(I), r(irr,0).col(0));
	dcomplex rZI = TDot(Z(irr,irr).col(I), r(irr,0).col(0));
	dcomplex dYI = TDot(Y(irr,irr).col(I), d(irr,0).col(0));
	dcomplex dZI = TDot(Z(irr,irr).col(I), d(irr,0).col(0));
	acc_l += (+rYI+rZI)*(+rYI-rZI)/(w-wI) + (+rZI+rYI)*(-rZI+rYI)/(-w-wI);
	acc_v += (+dYI-dZI)*(+dYI+dZI)/(w-wI) + (+dZI-dYI)*(-dZI-dYI)/(-w-wI);
      }
      muphi_psi1[irr]   = acc_l / 3.0 * (1.0*ne);
      muphi_psi1_v[irr] = acc_v / 3.0 * (1.0*ne);
    }
    CalcMain_alpha(iw);
  }
  PrintOut();
}
void Calc_Propagator_AO() {
  PrintTimeStamp("Prop_AO", NULL);
  // ---- Ground State ----
  int ng = basis0->size_basis_isym(0);
  BMat U_g;
  U_g(0, 0) = MatrixXcd::Zero(ng, ng);
  U_g(0, 0).col(0) = c0;

  // ---- first psi ----
  BMat S, T, V;
  CalcSTVMat(basis1, basis1, &S, &T, &V);  
  B2EInt eri_J = CalcERI(basis1, basis1, basis0, basis0, eri_method);
  B2EInt eri_K = CalcERI(basis1, basis0, basis0, basis1, eri_method);
  BMat B; Copy(S, B); B.SetZero(); AddK(eri_K, c0, 0, 1.0, B);  
  BMat A; CalcSTEXHamiltonian(T, V, 0, c0, eri_J, eri_K, &A);
  A.Add(-E0, S);

  // ---- dipole -----
  Irrep x = sym->irrep_x(); Irrep y = sym->irrep_y(); Irrep z = sym->irrep_z();  
  BVec r, d;
  BMat X1i, DX1i, Y1i, DY1i, Z1i, DZ1i;  
  CalcDipMat(basis1, basis0, &X1i, &Y1i, &Z1i, &DX1i, &DY1i, &DZ1i);
  r(x) = X1i(x,0)*c0; d(x) = DX1i(x,0)*c0;
  r(y) = Y1i(y,0)*c0; d(y) = DY1i(y,0)*c0;
  r(z) = Z1i(z,0)*c0; d(z) = DZ1i(z,0)*c0;

  // ---- build matrix ----
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    dcomplex w = w_list[iw];
  
    BOOST_FOREACH(Irrep irr, irreps) {
      int n1 = basis1->size_basis_isym(irr);
      MatrixXcd M(2*n1, 2*n1);
      VectorXcd qvec_l(2*n1), qvec_v(2*n1);
      for(int a = 0; a < n1; a++) {
	/*
	qvec_l(a)    = r(irr)(a);
	qvec_l(a+n1) = -r(irr)(a);
	qvec_v(a)    = d(irr)(a);
	qvec_v(a+n1) = d(irr)(a);
	*/
	qvec_l(a)    = r(irr)(a);
	qvec_l(a+n1) = r(irr)(a);
	qvec_v(a)    = d(irr)(a);
	qvec_v(a+n1) = -d(irr)(a);
	for(int b = 0; b < n1; b++) {
	  M(a, b)       = +w*S(irr,irr)(a,b) - A(irr,irr)(a,b);
	  M(a+n1, b+n1) = -w*S(irr,irr)(a,b) - A(irr,irr)(a,b);
	  M(a+n1, b) = -B(irr,irr)(a,b);
	  M(a, b+n1) = -B(irr,irr)(a,b);
	}
      }
      ColPivHouseholderQR<MatrixXcd> pivx = M.colPivHouseholderQr();
      VectorXcd xvec_l = pivx.solve(qvec_l);
      muphi_psi1[irr] = TDot(xvec_l, qvec_l) / 3.0 * (1.0*ne);
      MatrixXcd xvec_v = pivx.solve(qvec_v)  / 3.0 * (1.0*ne);
      muphi_psi1_v[irr] = TDot(xvec_v, qvec_v);
    }
    CalcMain_alpha(iw);
  }
  PrintOut();
}
int main(int argc, char *argv[]) {
  cout << ">>>> two_pot >>>>" << endl;
  if(argc == 1) {
    cerr << "need one argument" << endl;
    exit(1);
  }
  in_json = argv[1];
  Parse();
  PrintIn();
  if(calc_type == "propagator_ao") {
    Calc_Propagator_AO();
    return 0;
  }
  if(calc_type == "propagator" and solve_driv_type == "linear_solve") {
    Calc_Propagator_RPA();
    return 0;
  }
  if(calc_type == "propagator" and solve_driv_type == "eigen_value") {
    Calc_Propagator_RPA_eigen();
    return 0;
  }
  if(calc_type == "RPA" and solve_driv_type == "eigen_value") {
    Calc_RPA_HF();
    //Calc_RPA_Eigen_HF_old();
    return 0;
  }
  if(calc_type == "STEX" and solve_driv_type == "eigen_value") {
    Calc_STEX_Eigen();
    return 0;
  }
  if(calc_type == "one" and solve_driv_type == "eigen_value") {
    Calc_One_Eigen();
    return 0;
  }

  if(calc_type == "one") {
    CalcMat();
  } else if(calc_type == "STEX") {
    CalcMat();
    CalcMatSTEX();
  } else if(calc_type == "RPA") {
    CalcMat();
    CalcMatRPA();
  }
  PrintTimeStamp("Calc", NULL);
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw];
    cout << "w_eV: " << w * au2ev << endl;
    cout << "w_au: " << w << endl;
    cout << "E_eV: " << (w + E0) * au2ev << endl;
    cout << "E_au: " << w + E0 << endl;
    cout << "k_au: " << sqrt(2.0*(w + E0)) << endl;
    if(solve_driv_type == "linear_solve") {
      CalcDriv_linear_solve(iw);
    } else if(solve_driv_type == "eigen_value") {
      CalcDriv_eigen_value(iw);
    }
    CalcBraket();    
    CalcMain_alpha(iw);    
    CalcMain(iw);
  }
  
  PrintOut();
  cout << "<<<< two_pot <<<<" << endl;
  return 0;
}


