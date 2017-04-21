#include <fstream>
#include <boost/format.hpp>
#include "../external/picojson/picojson.h"
#include "../utils/eigen_plus.hpp"
#include "../utils/timestamp.hpp"
#include "../utils/read_json.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/symmol_read_json.hpp"


using namespace std;
using boost::format;
using namespace Eigen;
using namespace picojson;
using namespace cbasis;

void PrintHelp() {
  cout << "one_eig" << endl;
}
int main (int argc, char *argv[]) {
  cout << ">>>> one_eig >>>>" << endl;
  
  if(argc == 1) {
    PrintHelp();
    exit(1);
  }
  
  ifstream f;
  f.open(argv[1], ios::in);
  if(f.fail()) {
    cerr << "failed to open input file" << endl;
    cerr << "in_file: " << argv[1] << endl;
    exit(1);
  }
  
  picojson::value json;    
  // ==== parse json ====
  PrintTimeStamp("Parse", NULL);
  SymmetryGroup sym;
  SymGTOs gtos;
  Molecule mole;
  string comment, out_json, out_eigvecs, out_eigvals;
  int reduce_canonical_num;

  try {
    f >> json;
    picojson::object& obj = json.get<picojson::object>();

    // -- Comment --
    //comment = ReadJson<string>(obj, "comment");
    string key = "comment";
    comment = ReadJson<string>(obj, key);    

    // -- Basis --
    sym = ReadJson<SymmetryGroup>(obj, "sym");
    
    mole = NewMolecule(sym);
    ReadJson_Molecule(obj, "molecule", mole);

    gtos = NewSymGTOs(mole);
    ReadJson_SymGTOs_Subs(obj, "basis", gtos);
    if(obj.find("reduce_canonical_num") == obj.end())
      reduce_canonical_num =0;
    else
      reduce_canonical_num = ReadJson<int>(obj, "reduce_canonical_num");

    // -- out --
    out_json = ReadJson<string>(obj, "out_json");
    out_eigvecs =ReadJson<string>(obj, "out_eigvecs");
    out_eigvals =ReadJson<string>(obj, "out_eigvals");    
    
  } catch(std::exception& e) {
    string msg = "error on parsing json\n";
    msg += e.what();
    cerr << msg << endl;
    exit(1);
  }

  try {
    gtos->SetUp();
  } catch(exception& e) {
    cerr << "error on setup" << endl;
    cerr << "error message:" << endl;
    cerr << e.what() << endl;
    cerr << "gtos::" << endl;
    cerr << gtos->str() << endl;
    exit(1);
  }

  cout << "comment: " << comment << endl;
  cout << "in_json: " << argv[1] << endl;
  cout << "out_json: " << out_json << endl;
  cout << "out_eigvecs: " << out_eigvecs << endl;
  cout << "out_eigvals: " << out_eigvals << endl;
  cout << "reduce_canonical_num: " << reduce_canonical_num << endl;
  cout << "symmetry: " << sym->name() << endl;
  cout << "molecule: " << endl << mole->show() << endl;
  cout << "gtos: " << endl << gtos->show() << endl;

  // ==== calculation ====
  PrintTimeStamp("Calc", NULL);
  BMat S, T, V, C;
  BVec E;
  
  CalcSTVMat(gtos, gtos, &S, &T, &V);

  cout << "Eigensystem of overlap" << endl;
  for(int irrep = 0; irrep < sym->order(); irrep++) {    
    if(gtos->size_basis_isym(irrep) != 0) {
      MatrixXcd& s = S(irrep, irrep);
      ComplexEigenSolver<MatrixXcd> s_solver(s);
      const MatrixXcd& SVecs = s_solver.eigenvectors();
      const VectorXcd& SVals = s_solver.eigenvalues();

      int num = s.rows(); int numcol = 5;

      string irrep_name = sym->GetIrrepName(irrep);
      cout << "Irrep = " << irrep_name << endl;
      for(int nblock = 0; nblock < num/numcol; nblock++) {
	cout << "Eigvals | ";
	for(int n = nblock*numcol; n < (nblock+1)*numcol; n++) {
	  cout << format("%10.5e ") % abs(SVals[n]);
	}
	cout << endl;
	for(int m = 0; m < num; m++) {
	  cout << format("%7d | ") % m;
	  for(int n = nblock*numcol; n < (nblock+1)*numcol; n++) {
	    cout << format("%10.5f ") % abs(SVecs(m, n));
	  }
	  cout << endl;
	}
	cout << endl;
      }
    }
  }
  for(int irrep = 0; irrep < sym->order(); irrep++) {
    if(gtos->size_basis_isym(irrep) != 0) {
      MatrixXcd h = T(irrep, irrep) + V(irrep, irrep);
      MatrixXcd s = S(irrep, irrep);
      SymGenComplexEigenSolver solver(h, s, h.rows()-reduce_canonical_num);
      E(irrep) = solver.eigenvalues();
      C(irrep, irrep) = solver.eigenvectors();

      //      cout << "diff:" 
      //	   << (h * C(irrep, irrep).col(0) - E(irrep)(0)*s*C(irrep, irrep).col(0)).array().abs().maxCoeff()
      //	   << endl;

      }
      //      cout << "eig_of_s: " << s_solver.eigenvalues()[0] << endl;  
  }
  
  // ==== output ====
  PrintTimeStamp("Out", NULL);
  cout << "E0 = " << E(0)(0) << endl;

  cout << "writing E" << endl;
  E.Write(out_eigvals);
  cout << "writing C" << endl;
  C.Write(out_eigvecs);
  cout << "out.json" << endl;
  
  picojson::object out;
  out["comment"] = picojson::value(comment);
  out["in_json"] = picojson::value(argv[1]);
  out["out_eigvals"] = picojson::value(out_eigvals);
  out["out_eigvecs"] = picojson::value(out_eigvecs);
  out["E0"] = ToJson(E(0)(0));
  int irrep0(0);
  out["irrep0"] = ToJson(irrep0);
  picojson::value out_val(out);
  ofstream of(out_json.c_str(), ios::out);
  of << out_val.serialize();
  cout << "<<<< one_eig <<<<" << endl;
}
