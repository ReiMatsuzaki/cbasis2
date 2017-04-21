#include <fstream>
#include <iostream>
#include <stdexcept>
#include <boost/format.hpp>
#include "../external/picojson/picojson.h"
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../utils/read_json.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/two_int.hpp"
#include "../src_cpp/mo.hpp"
#include "../src_cpp/symmol_read_json.hpp"

using namespace std;
using namespace cbasis;

// -- Input --
string comment, in_json, in_eigvecs, in_eigvals, out_eigvecs, out_eigvals;
SymmetryGroup sym;
Molecule mole;
SymGTOs gtos;
int num_ele;
SCFOptions scf_opts;
ERIMethod eri_method;
BVec E0;
BMat C0;

// -- Results --
MO mo;
bool conv;
SCFOptions ReadJson_SCFOptions(picojson::object& obj0, string key) {
  SCFOptions opts;
  if(obj0.find(key) == obj0.end()) {
    throw runtime_error("not found key in ReadJson_SCFOptions. key = "+key);
  }
  if(not obj0[key].is<picojson::object>()) {
    throw runtime_error("in ReadJson_SCFOptions, json must be object");
  }

  picojson::object& obj = obj0[key].get<picojson::object>();
  
  if(obj.find("debug_lvl") != obj.end()) {
    opts.debug_lvl = ReadJson<int>(obj, "debug_lvl");
  }
  if(obj.find("max_iter") != obj.end()) {
    cout << "maxit" << endl;
    opts.max_iter = ReadJson<int>(obj, "max_iter");
  }
  if(obj.find("tol") != obj.end()) {
    opts.tol = ReadJson<double>(obj, "tol");
  }
  if(obj.find("use_real") != obj.end()) {
    opts.use_real = ReadJson<bool>(obj, "use_real");
  }
  return opts;
}
void Parse() {
  PrintTimeStamp("Parse", NULL);
  
  ifstream f;
  picojson::value json;

  f.open(in_json.c_str(), ios::in);
  if(f.fail()) {
    cerr << "opening input json file failed." << endl;
    cerr << "filename : " << in_json << endl;
    exit(1);
  }

  try {
    f >> json;
    if(not json.is<picojson::object>()) {
      throw runtime_error("json is not object");
    }
    picojson::object& obj = json.get<picojson::object>();
    comment = ReadJson<string>(obj, "comment");
    sym = ReadJson<SymmetryGroup>(obj, "sym");
    mole = NewMolecule(sym); ReadJson_Molecule(obj, "molecule", mole);
    gtos = NewSymGTOs(mole);
    ReadJson_SymGTOs_Subs(obj, "basis", gtos);
    num_ele   = ReadJson<int>(obj, "num_ele");
    scf_opts = ReadJson_SCFOptions(obj, "scf_opts");
    eri_method = ReadJson<ERIMethod>(obj, "eri_method");
    in_eigvecs = ReadJson<string>(obj, "in_eigvecs");
    in_eigvals = ReadJson<string>(obj, "in_eigvals");
    out_eigvecs = ReadJson<string>(obj, "out_eigvecs");
    out_eigvals = ReadJson<string>(obj, "out_eigvals");
    
  } catch(exception& e) {
    cerr << "error on parse json" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  f.close();

  try {
    gtos->SetUp();
  } catch(exception& e) {
    cerr << "error on setup of SymGTOs" << endl;
    cerr << e.what() << endl;
    cerr << "gtos:" << endl;
    cerr << gtos->str() << endl;
    exit(1);
  }

  try {
    E0.Read(in_eigvals);
  } catch(exception& e) {
    cerr << "error on reading in_eigvals" << endl;
    exit(1);
  }

  try {
    C0.Read(in_eigvecs);
  } catch(exception& e) {
    cerr << "error on reading in_eigvecs" << endl;
    exit(1);
  }
  
}
void PrintIn() {
  PrintTimeStamp("PrintIn", NULL);
  cout << "comment: " << comment << endl;
  cout << "in_json: " << in_json << endl;
  cout << "in_eigvecs: " << in_eigvecs << endl;
  cout << "in_eigvals: " << in_eigvals << endl;  
  cout << "out_eigvecs: " << out_eigvecs << endl;
  cout << "out_eigvals: " << out_eigvals << endl;
  cout << "ERIMethod_use_symmetry: " << eri_method.symmetry << endl;
  cout << "ERIMethod_use_memo: " << eri_method.coef_R_memo << endl;
  cout << "ERIMethod_use_perm: " << eri_method.perm << endl;
  cout << "SCFOptions_use_real: " << (scf_opts.use_real?"true":"false") << endl;
  cout << "SCFOptions_max_iter: " << scf_opts.max_iter << endl;
  cout << "SCFOptions_tol: " << scf_opts.tol << endl;
  cout << "symmetry: " << sym->name() << endl;
  cout << "molecule: " << endl << mole->show() << endl;
  cout << "num_ele: " << num_ele << endl;
  cout << "gtos: " << endl << gtos->show() << endl;
}
void CalcMat() {
  
  //PrintTimeStamp("Mat", NULL);
    /*
    try {
      E.Read(in_eigvals);
      C.Read(in_eigvecs);
    } catch(exception& e) {
      cerr << "error on reading initial eigvals or eigvecs" << endl;
      cerr << e.what() << endl;
      exit(1);
    }
    vector<int> occ_num = CalcOccNum(E, sym->order(), num_ele/2);    
    CalcSTVMat(gtos, gtos, &S, &T, &V);

    ERIMethod eri_method;
    eri = CalcERI_Complex(gtos, eri_method);
    */
    
}
void CalcMain() {
  
  PrintTimeStamp("Calc", NULL);
  BMatSet mat_set = CalcMat_Complex(gtos, true);
  B2EInt  eri     = CalcERI_Complex(gtos, eri_method);
  
  mo = CalcRHF(sym, mat_set, eri, num_ele, E0, C0, scf_opts, &conv);

}
void PrintOut() {

  PrintTimeStamp("PrintOut", NULL);

  try {
    mo->C.Write(out_eigvecs);
  } catch(exception& e) {
    cerr << "error on writing coef matrix" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  
  try {
    mo->eigs.Write(out_eigvals);
  } catch(exception& e) {
    cerr << "error on writing eig vector" << endl;
    cerr << e.what() << endl;
    exit(1);
  }
  
  cout << "convergence: " << (conv ? "Yes" : "No") << endl;
  cout << "orbital_energy: " << mo->eigs(0)(0) << endl;
  cout << "ele_energy: " << mo->energy << endl;
  cout << "mole_energy: " << mo->energy + mole->NucEnergy() << endl;

  /*
  picojson::object out;
  out["comment"] = picojson::value(comment);
  out["in_json"] = picojson::value(in_json);
  picojson::array ary;
  ary.push_back(picojson::value(mo->energy.real()));
  ary.push_back(picojson::value(mo->energy.imag()));
  out["E0"] = picojson::value(ary);
  ofstream ofs(out_json.c_str())
  */
}
int main (int argc, char *argv[]) {
  cout << ">>>> rhf >>>>" << endl;
  if(argc != 2) {
    cerr << "need input json file" << endl;
    exit(1);
  }
  in_json = argv[1];
  Parse();
  PrintIn();
  CalcMat();
  CalcMain();
  PrintOut();
  cout << "<<<< rhf <<<<" << endl;
}
