#include <fstream>
#include <iostream>
#include <boost/foreach.hpp>
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

/**
   Partialy Orthonormalized Plane Wave (POWP) method.
   |POPW> = |PW> - |phi0><phi0|PW>
   TDM = <phi0|D|POPW> 
   .   = <phi0|D|PW> -<phi0|D|phi0><phi0|PW>
   .   = sum_i c_i <g_i|D|PW>
 */

// -- Input --
string comment;
string in_json, out_json, cs_csv;
SymmetryGroup sym;
Molecule mole;
SymGTOs basis0;
Irrep irrep0;
int i0;
dcomplex E0;
VectorXcd c0;
vector<Vector3d> k_list;

// -- Intermediate --
BVec S, X, Y, Z;

// -- Results --
vector<dcomplex> tdm_x_list;
vector<dcomplex> tdm_y_list;
vector<dcomplex> tdm_z_list;
//vector<dcomplex> cs_list;

void PrintHelp() {
  
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
    ReadJson_Orbital(obj, "orbital0", sym, &E0, &c0, &irrep0, &i0);

    VectorXd k = ReadJson<VectorXd>(obj, "k", 3);
    k_list.push_back(Vector3d(k));
    
  } catch(exception& e) {
    cerr << "error on parsing json" << endl;
    cerr << e.what() << endl;
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
  cout << "k_list: "  << endl;

  BOOST_FOREACH(Vector3d& k, k_list) {
    cout << "[" << k[0] << ", " << k[1] << ", " << k[2] << "]" << endl;
  }
  cout << "mole:" << endl << mole->show() << endl;
  cout << "basis0:" << endl << basis0->show() << endl;

}
void Init() {
  PrintTimeStamp("Init", NULL);

  try {
    basis0->SetUp();
  } catch(exception& e) {
    cerr << "error on SetUp of basis0" << endl;
    cerr << e.what() << endl;
    cerr << "basis0:" << endl;
    cerr << basis0->str() << endl;
    exit(1);
  }  

  InitBVec(basis0, &S);
  InitBVec(basis0, &X);
  InitBVec(basis0, &Y);
  InitBVec(basis0, &Z);

  int n = k_list.size();
  tdm_x_list.resize(n);
  tdm_y_list.resize(n);
  tdm_z_list.resize(n);

}
void Calc() {
  PrintTimeStamp("Calc", NULL);

  int num = k_list.size();
  cout << "kx, ky, kz, (PW|X|phi0), (PW|Y|phi0), (PW|Z|phi0)" << endl;
  for(int i = 0; i < num; i++) {
    Vector3d k = k_list[i];
    Vector3cd kc = k.cast<dcomplex>();
    CalcPWVec(basis0, kc, &S, &X, &Y, &Z);
    tdm_x_list[i] = TDot(X[irrep0], c0);
    tdm_y_list[i] = TDot(Y[irrep0], c0);
    tdm_z_list[i] = TDot(Z[irrep0], c0);
    cout << k[0] << ", "<< k[1] << ", " << k[2] << ", "
	 << tdm_x_list[i] << ", " << tdm_y_list[i] << ", " << tdm_z_list[i]
	 << endl;
  }
  
}
void PrintOut() {
  PrintTimeStamp("PrintOut", NULL);
  
  ofstream f(cs_csv.c_str(), ios::out);
  f << "w,x,y,z" << endl;

}
int main(int argc, char *argv[]) {
  cout << ">>>> popw >>>>" << endl;
  if(argc == 1) {
    PrintHelp(); exit(1);
  }
  in_json = argv[1];
  Parse();
  PrintIn();
  Init();
  Calc();
  PrintOut();
  cout << "<<<< popw <<<<" << endl;

  return 0;
}
