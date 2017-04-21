#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/read_json.hpp"

using namespace std;
using namespace Eigen;
using namespace cbasis;
using namespace picojson;

// -- Input --
string comment;
string in_json, out_json, cs_csv;
SymmetryGroup sym;
Molecule mole;
SymGTOs basis0, basis1;
Irrep irrep0;
int i0;
dcomplex E0;
VectorXcd c0;
vector<dcomplex> w_list;

// -- intermediate --
BMat S, T, V;
BMat X, Y, Z, DX, DY, DZ, L;
BVec sX, sY, sZ, sDX, sDY, sDZ;
BVec cX, cY, cZ, cDX, cDY, cDZ;

// -- result --
vector<dcomplex> alpha_x, alpha_y, alpha_z, alpha_dx, alpha_dy, alpha_dz;

void PrintHelp() {
  cout << "one_driv" << endl;
}
void Parse() {
  PrintTimeStamp("Parse", NULL);
  try {
    // -- open file --
    ifstream f(in_json.c_str());  
    value json; f >> json;
    if(not json.is<object>()) {
      throw runtime_error("input json file is not object");
    }
    object& obj = json.get<object>();

    comment = ReadJson<string>(obj, "comment");
    sym     = ReadJson<SymmetryGroup>(obj, "sym");
    out_json = ReadJson<string>(obj, "out_json");
    cs_csv   = ReadJson<string>(obj, "cs_csv");
    mole = NewMolecule(sym); ReadJson_Molecule(obj, "molecule", mole);
    basis0 = NewSymGTOs(mole); ReadJson_SymGTOs_Subs(obj, "basis0", basis0);
    basis1 = NewSymGTOs(mole); ReadJson_SymGTOs_Subs(obj, "basis1", basis1);
    ReadJson_Orbital(obj, "orbital0", sym, &E0, &c0, &irrep0, &i0);

    VectorXcd _ws = ReadJson<VectorXcd>(obj, "ws");
    for(int i = 0; i < _ws.size(); i++)
      w_list.push_back(_ws(i));

  } catch(exception& e) {
    cerr << "error on parsing json" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  
}
void SetUp() {
  PrintTimeStamp("Set Up", NULL);
  try {
    basis0->SetUp();
  } catch(exception& e) {
    cerr << "error on SetUp of basis0" << endl;
    cerr << e.what() << endl;
    cerr << "basis0:" << endl;
    cerr << basis0->str() << endl;
    exit(1);
  }
  try {
    basis1->SetUp();
  } catch(exception& e) {
    cerr << "error on SetUp of basis1" << endl;
    cerr << e.what() << endl;
    cerr << "basis1:" << endl;
    cerr << basis1->str() << endl;
    exit(1);
  }    
}
void PrintIn() {
  PrintTimeStamp("PrintIn", NULL);
  cout << "comment: " << comment << endl;
  cout << "in_json: " << in_json << endl;
  cout << "cs_csv: "  << cs_csv;
  cout << "out_json: " << out_json << endl;
  cout << "E0: " << E0 << endl;
  cout << "mole:" << endl << mole->show() << endl;
  cout << "basis0:" << endl << basis0->show() << endl;  
  cout << "basis1:" << endl << basis1->show() << endl; 
}
void Init() {
  PrintTimeStamp("Init", NULL);

  Irrep x = sym->irrep_x();
  Irrep y = sym->irrep_y();
  Irrep z = sym->irrep_z();

  InitBVec(basis1, &cX); InitBVec(basis1, &cY); InitBVec(basis1, &cZ);
  InitBVec(basis1, &cDX);InitBVec(basis1, &cDY);InitBVec(basis1, &cDZ);
  
  vector<Irrep> irrep_x_list; sym->Non0IrrepList(irrep0, x, &irrep_x_list);
  vector<Irrep> irrep_y_list; sym->Non0IrrepList(irrep0, y, &irrep_y_list);
  vector<Irrep> irrep_z_list; sym->Non0IrrepList(irrep0, z, &irrep_z_list);
  InitBVec(basis1, irrep_x_list, &sX); InitBVec(basis1, irrep_x_list, &sDX);
  InitBVec(basis1, irrep_y_list, &sY); InitBVec(basis1, irrep_y_list, &sDY);
  InitBVec(basis1, irrep_z_list, &sZ); InitBVec(basis1, irrep_z_list, &sDZ);
  InitBVec(basis1, irrep_x_list, &cX); InitBVec(basis1, irrep_x_list, &cDX);
  InitBVec(basis1, irrep_y_list, &cY); InitBVec(basis1, irrep_y_list, &cDY);
  InitBVec(basis1, irrep_z_list, &cZ); InitBVec(basis1, irrep_z_list, &cDZ);

  alpha_x.resize(w_list.size()); alpha_dx.resize(w_list.size()); 
  alpha_y.resize(w_list.size()); alpha_dy.resize(w_list.size()); 
  alpha_z.resize(w_list.size()); alpha_dz.resize(w_list.size()); 
}
void Calc() {
  PrintTimeStamp("Calc", NULL);
  CalcSTVMat(basis1, basis1, &S, &T, &V);
  CalcDipMat(basis1, basis0, &X, &Y, &Z, &DX, &DY, &DZ);
  for(int irrep = 0; irrep < sym->order(); irrep++) {
    if(X.has_block(irrep, irrep0)) {
      sX(irrep)  = X(irrep, irrep0)  * c0;
      sDX(irrep) = DX(irrep, irrep0) * c0;
    }
    if(Y.has_block(irrep, irrep0)) {
      sY(irrep)  = Y(irrep, irrep0)  * c0;
      sDY(irrep) = DY(irrep, irrep0) * c0;
    }
    if(Z.has_block(irrep, irrep0)) {
      sZ(irrep)  = Z(irrep, irrep0)  * c0;
      sDZ(irrep) = DZ(irrep, irrep0) * c0;      
    }
  }

  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    dcomplex w = w_list[iw];
    cout << "w = " << w << endl;
    for(int irrep = 0; irrep < sym->order(); irrep ++) {
      if(X.has_block(irrep, irrep0)) {
	MatrixXcd& Li = L(irrep, irrep);
	Li = S(irrep, irrep)*(E0+w) - T(irrep, irrep) - V(irrep, irrep);
	ColPivHouseholderQR<MatrixXcd> piv = Li.colPivHouseholderQr();
	cX(irrep)  = piv.solve(sX(irrep));
	cDX(irrep) = piv.solve(sDX(irrep));
	alpha_x[iw]  = TDot(cX(irrep), sX(irrep));
	alpha_dx[iw] = TDot(cDX(irrep), sDX(irrep));
	cout << "x : " << alpha_x[iw] << "  " << alpha_dx[iw] << endl;
      }
      if(Y.has_block(irrep, irrep0)) {
	MatrixXcd& Li = L(irrep, irrep);
	Li = S(irrep, irrep)*(E0+w) - T(irrep, irrep) - V(irrep, irrep);
	ColPivHouseholderQR<MatrixXcd> piv = Li.colPivHouseholderQr();
	cY(irrep)  = piv.solve(sY(irrep));
	cDY(irrep) = piv.solve(sDY(irrep));
	alpha_y[iw]  = TDot(cY(irrep),  sY(irrep));
	alpha_dy[iw] = TDot(cDY(irrep), sDY(irrep));
	cout << "y : " << alpha_y[iw] << "  " << alpha_dy[iw] << endl;
      }
      if(Z.has_block(irrep, irrep0)) {
	MatrixXcd& Li = L(irrep, irrep);
	Li = S(irrep, irrep)*(E0+w) - T(irrep, irrep) - V(irrep, irrep);
	ColPivHouseholderQR<MatrixXcd> piv = Li.colPivHouseholderQr();
	cZ(irrep)  = piv.solve(sZ(irrep));
	cDZ(irrep) = piv.solve(sDZ(irrep));
	alpha_z[iw]  = TDot(cZ(irrep),  sZ(irrep));
	alpha_dz[iw] = TDot(cDZ(irrep), sDZ(irrep));
	cout << "z : " << alpha_z[iw] << "  " << alpha_dz[iw] << endl;
      }
    }
  }
}
void PrintOut() {  
  PrintTimeStamp("PrintOut", NULL);
  
  double c = 137.035999139;
  double au2mb = 5.291772 * 5.291772;

  ofstream f(cs_csv.c_str(), ios::out);
  f << "w,E,x,y,z,totl,dx,dy,dz,totv" << endl;
  for(int iw = 0; iw < (int)w_list.size(); iw++) {
    double w = w_list[iw].real();
    double cs_x = 4.0*M_PI*w/c   * alpha_x[iw].imag() * au2mb;
    double cs_y = 4.0*M_PI*w/c   * alpha_y[iw].imag() * au2mb;
    double cs_z = 4.0*M_PI*w/c   * alpha_z[iw].imag() * au2mb;
    double cs_l = (cs_x+cs_y+cs_z);
    double cs_dx= 4.0*M_PI/(w*c) * alpha_dx[iw].imag() * au2mb;
    double cs_dy= 4.0*M_PI/(w*c) * alpha_dy[iw].imag() * au2mb;
    double cs_dz= 4.0*M_PI/(w*c) * alpha_dz[iw].imag() * au2mb;
    double cs_v = (cs_dx+cs_dy+cs_dz);
    f << w
      << "," << w+E0
      << "," << cs_x
      << "," << cs_y
      << "," << cs_z
      << "," << cs_l
      << "," << cs_dx
      << "," << cs_dy
      << "," << cs_dz
      << "," << cs_v
      << endl;

  }
}
int main(int argc, char *argv[]) {
  cout << ">>>> one_driv >>>>" << endl;
  if(argc == 1) {
    PrintHelp(); exit(1);
  }
  in_json = argv[1];
  Parse();
  SetUp();
  PrintIn();
  Init();
  Calc();
  PrintOut();
  cout << "<<<< one_driv <<<<" << endl;
  return 0;
}
