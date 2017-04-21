#include <stdexcept>
#include <fstream>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/QR>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include "../utils/timestamp.hpp"
#include "../utils/eigen_plus.hpp"
#include "../utils/read_json.hpp"
#include "../r1basis/r1basis.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace picojson;
using namespace cbasis;

void DrivenTerm(string channel, string dipole, LC_STOs stos) {

  bool is_length = dipole == "length";
  bool is_velocity = dipole == "velocity";
  bool is_1skp = channel == "1s->kp";
  if(is_1skp && is_length) {
    stos->Add(2.0, 2, 1.0);
  } else if(is_1skp && is_velocity) {
    stos->Add(-2.0, 1, 1.0);
  } else {
    throw runtime_error("not supported combination of channel and dipole.");
  }
}

template<int MB, int MD>
class R1Driv {
public:
  typedef typename _EXPs<MB>::EXPs BasisEXPs;
  typedef typename _EXPs<MD>::LC_EXPs DrivLC_EXPs;
  object& obj;
  string in_file;
  string out_file;
  string comment;

  LinearSolver linear_solver;
  
  dcomplex Z;
  int L;
  dcomplex E;
  
  BasisEXPs basis;
  DrivLC_EXPs driv;

  bool write_matrix;

  bool write_psi;
  string write_psi_filename;
  VectorXd write_psi_rs;
  
  bool write_basis;
  string write_basis_filename;
  VectorXd write_basis_rs;

  R1Driv(object& _obj): obj(_obj), basis(new _EXPs<MB>()), driv(new _LC_EXPs<MD>()) {
    write_matrix = false;
    write_psi = false;
    write_basis = false;
  }
  void Parse() {
    PrintTimeStamp("Parse", NULL);

    // -- parse simple --            
    try {
      
      this->comment = ReadJson<string>(obj, "comment");
      this->linear_solver = ReadJsonWithDefault
	<LinearSolver>(obj, "linear_solver", LinearSolver());
      cout << format("comment: %s\n") % this->comment;
    } catch(std::exception& e) {
      cerr << "error on parsing json\n";
      cerr << e.what() << endl;
      exit(1);
    }

    // -- parse options --
    bool find_opt = (obj.find("opts") != obj.end());
    if(find_opt) {
      CheckObject<object>(obj, "opts");
      object& opt_obj = obj["opts"].get<object>();
      write_matrix = ReadJsonWithDefault<bool>(opt_obj, "write_matrix", false);
      cout << format("write_matrix: %s\n") % (write_matrix ? "yes" : "no");
    }
    
    // -- write wave func --
    this->write_psi = (obj.find("write_psi") != obj.end());
    if(this->write_psi) {
      
      try {
	CheckObject<object>(obj, "write_psi");
	object& wavefunc_obj = obj["write_psi"].get<object>();
	this->write_psi_filename = ReadJson<string>(wavefunc_obj, "file");
	cout << format("write_psi_file: %s\n") % this->write_psi_filename;
	this->write_psi_rs = ReadJson<VectorXd>(wavefunc_obj, "rs");
	cout << format("write_psi_rs\n");
      } catch(std::exception& e) {
	cerr << "error on parsing write_psi\n";
	cerr << e.what() << endl;
	exit(1);
      }
    }
    
    // -- write basis --
    this->write_basis = (obj.find("write_basis") != obj.end());
    if(this->write_basis) {
      
      try {
	CheckObject<object>(obj, "write_basis");
	object& write_basis_obj = obj["write_basis"].get<object>();
	this->write_basis_filename = ReadJson<string>(write_basis_obj, "file");
	cout << format("write_basis_file: %s\n") % this->write_basis_filename;
	this->write_basis_rs = ReadJson<VectorXd>(write_basis_obj, "rs");
	cout << format("write_basis_rs: \n");
      } catch(std::exception& e) {
	cerr << "error on pasing write_basis\n";
	cerr << e.what() << endl;
	exit(1);
      }
    }
    
    // -- parse target --
    try {      
      CheckObject<object>(obj, "target");
      object& target_obj = obj["target"].get<object>();
      string target_type = ReadJson<string>(target_obj, "type");
      if(target_type == "h_pi") {
	cout << "target_type: h_pi\n";
	string channel = ReadJson<string>(target_obj, "channel");
	cout << format("channel: %s\n") % channel;
	string dipole  = ReadJson<string>(target_obj, "dipole");
	cout << format("dipole: %s\n") % dipole;
	DrivenTerm(channel, dipole, this->driv);	
      } else if(target_type == "custom") {
	cout << format("target_type: %s\n") % target_type;

	// - driv -
	CheckObject<object>(target_obj, "driv");
	object& driv_obj = target_obj["driv"].get<object>();
	//string driv_type = ReadJso<string>(driv_obj, "type");
	CheckObject<picojson::array>(driv_obj, "value");
	picojson::array& vals = driv_obj["value"].get<picojson::array>();	  
	BOOST_FOREACH(value& val, vals) {
	  CheckValue<object>(val);
	  object& obj = val.get<object>();
	  dcomplex c = ReadJson<dcomplex>(obj, "c");
	  int      n = ReadJson<int>(obj, "n");
	  dcomplex z = ReadJson<dcomplex>(obj, "z");
	  this->driv->Add(c, n, z);
	}

	// - other -
	Z = ReadJson<dcomplex>(target_obj, "Z");
	cout << "Z: " << Z << endl;
	L = ReadJson<int>(target_obj, "L");
	cout << "L: " << L << endl;
	E = ReadJson<dcomplex>(target_obj, "E");
	cout << "E: " << E << endl;
	
      } else {
	throw runtime_error("not supported type in target");
      }
    } catch(std::exception& e) {
      cerr << "error on parsing target\n";
      cerr << e.what() << endl;
      exit(1);
    }

    // -- parse basis --
    try {
      string normalization_const = "none";
      CheckObject<object>(obj, "basis");
      object& basis_obj = obj["basis"].get<object>();
      if(basis_obj.find("opts") != basis_obj.end()) {
	CheckObject<object>(basis_obj, "opts");
	object& opt_obj = basis_obj["opts"].get<object>();
	normalization_const = ReadJson<string>(opt_obj, "normalization_const");
	if(normalization_const != "none" &&
	   normalization_const != "continue" &&
	   normalization_const != "sqrt" ) {
	  throw runtime_error("normalization_const <- {none, continue, sqrt}" );
	}
      }
      CheckObject<picojson::array>(basis_obj, "value");
      picojson::array& basis_vals = basis_obj["value"].get<picojson::array>();
      BOOST_FOREACH(value& val, basis_vals) {
	try {
	  CheckValue<object>(val);
	} catch(std::exception& e) {
	  string msg; msg = "error on element of basis/value\n";
	  msg += e.what();
	  throw runtime_error(msg);
	}
	object& obj0 = val.get<object>();
	int n = ReadJson<int>(obj0, "n");
	cout << "n: " << n << endl;
	if(obj0.find("z") != obj0.end()) {
	  dcomplex z = ReadJson<dcomplex>(obj0, "z");
	  this->basis->AddPrim(n, z);
	  cout << "z: " << z << endl;
	} else if(obj0.find("czs") != obj0.end()) {
	  MatrixXcd czs = ReadJson<MatrixXcd>(obj0, "czs");
	  cout << format("normalization_const: %s\n") % normalization_const;
	  cout << "czs: " << endl << czs << endl;
	  
	  typename _EXPs<MB>::LC_EXPs lc(new _LC_EXPs<MB>());
	  for(int i = 0; i < czs.rows(); i++) {
	    dcomplex c(czs(i, 0));
	    dcomplex z(czs(i, 1));
	    dcomplex nterm(1);
	    if(normalization_const == "continue") {
	      nterm = NormalizationTermContinue<MB>(n, z);
	    } else if(normalization_const == "sqrt") {
	      nterm = 1.0/sqrt(EXPInt<MB, MB>(2*n, z, z));
	    } 
	    lc->Add(c*nterm, n, z);
	  }
	  this->basis->AddLC(lc);
	} else {
	  throw runtime_error("key \"z\" or \"czs\" are necessary in value in basis");
	}
      }
      
    } catch(std::exception& e) {
      cerr << "error on parsing basis\n";
      cerr << e.what() << endl;
      exit(1);
    }
    
  }
  void PrintIn() {
    PrintTimeStamp("PrintIn", NULL);
    cout << "linear_solver: " << linear_solver.show() << endl;
    if(this->write_psi) {
      cout << "write_psi_filename: " << this->write_psi_filename << endl;
      VectorXd& rs = this->write_psi_rs;
      cout << format("write_psi_rs: [%5.3f, %5.3f, ..., %5.3f]\n")% rs[0] % rs[1] % rs[rs.size()-1];
    } else {
      cout << "write_psi: no" << endl;
    }
    if(this->write_basis) {
      VectorXd& rs = this->write_basis_rs;
      cout << format("write_basis_file: %s\n") % this->write_basis_filename;
      cout << format("write_basis_rs: [%5.3f, %5.3f, ..., %5.3f]\n")
	% rs[0] % rs[1] % rs[rs.size()-1];
    } else {
      cout << "write_basis: no\n";
    }
    cout << format("write_matrix: %s\n") % (write_matrix ? "yes" : "no");
    cout << "Z: " << Z << endl;
    cout << "L: " << L << endl;
    cout << "E: " << E << endl;
    cout << format("driven_term: %s\n") % driv->str();
    cout << format("basis:\n");
    cout << basis->str() << endl;
  }
  void Calc() {
    PrintTimeStamp("Calc", NULL);
    
    // -- set up basis --
    this->basis->SetUp();

    // -- build matrix --
    // D = E - T = E + 1/2 D2 - L(L+1)/2r2 +Z/r 
    MatrixXcd D, R2, R1, H, S;
    this->basis->InitMat(H);
    this->basis->InitMat(S);
    this->basis->InitMat(D);
    this->basis->InitMat(R2);
    this->basis->InitMat(R1);
    
    
    this->basis->CalcD2Mat(H);
    H *= -0.5;
    this->basis->CalcRmMat(-2, R2);
    H += 0.5*(L*(L+1))  * R2;
    this->basis->CalcRmMat(-1, R1);
    H += -Z * R1;
    this->basis->CalcRmMat(0, S);
    D = E *S - H;

    // -- print matrix --
    if(write_matrix) {
      cout << "H matrix: " << endl;
      for(int i = 0; i < H.rows(); i++)
	for(int j = 0; j < H.cols(); j++) {
	  dcomplex v(H(i,j));
	  cout << format("(%d,%d) = (%f, %f)\n") %i%j% v.real() % v.imag();
	}
      cout << "S matrix: " << endl;
      for(int i = 0; i < H.rows(); i++)
	for(int j = 0; j < H.cols(); j++) {
	  dcomplex v(S(i,j));
	  cout << format("(%d,%d) = (%f, %f)\n") %i%j% v.real() % v.imag();
	}      
    }

    // -- build vector --
    VectorXcd m; this->basis->InitVec(m);
    this->basis->CalcVec(driv, m);

    // -- solve linear problem --
    VectorXcd c;
    this->linear_solver.Solve(D, m, &c);

    // -- alpha --
    dcomplex alpha = TDot(m, c);
    cout << format("alpha: %20.15f, %20.15f\n") % alpha.real() % alpha.imag();

    // -- radial dipole --
    dcomplex k = sqrt(2.0*E);
    double rdm = sqrt(abs(alpha.imag()*k/2.0));
    
    cout << format("rad_dipole_moment: %20.15f\n") % rdm;

    // -- basis --
    if(this->write_basis) {
      ofstream f(this->write_basis_filename.c_str());
      int num_basis = this->basis->size();
      f << "r";
      for(int i = 0; i < num_basis; i++) {
	f << format(",re_%d,im_%d") % i % i;
      }
      f << endl;
      VectorXd& rs = this->write_basis_rs;
      vector<VectorXcd> yss;
      for(int i = 0; i < num_basis; i++) {
	VectorXcd c0 = VectorXcd::Zero(num_basis);
	c0[i] = 1.0;
	yss.push_back(this->basis->AtR(rs.cast<dcomplex>(), c0));
      }
      for(int j = 0; j < rs.size(); j++) {
	f << rs[j];
	for(int i = 0; i < num_basis; i++) {
	  dcomplex& y(yss[i][j]);
	  f << format(",%f,%f") % y.real() % y.imag();
	}
	f << endl;
      }
    }
    
    // -- wave function --
    if(this->write_psi) {
      VectorXd& rs = this->write_psi_rs;
      VectorXcd ys = this->basis->AtR(rs.cast<dcomplex>(), c);
      int num = ys.size();
      ofstream f(this->write_psi_filename.c_str(), ios::out);
      f << "r,re_y,im_y\n";
      for(int i = 0; i < num; i++) {
	f << format("%f,%f,%f\n") % rs[i] % ys[i].real() % ys[i].imag();
      }
    }
  }  
  void Run() {
    this->Parse();
    this->PrintIn();
    this->Calc();
  }
};
int main(int argc, char *argv[]) {
  
  cout << ">>>> r1_driv >>>>\n";

  // -- argument number --
  if(argc != 2) {
    cout << "one argument is necessary\n";
    exit(1);
  }

  // -- open and read file --  
  ifstream f(argv[1]);
  if(f.fail()) {
    throw runtime_error("open file failed");
  }
  value json; f >> json;
  if(not json.is<object>()) {
    throw runtime_error("input json file is not object");
  }
  object& obj = json.get<object>();
  cout << format("in_file: %s\n") % argv[1];
  
  // -- check function type --
  try {
    CheckObject<object>(obj, "target");
    object& target_obj = obj["target"].get<object>();
    CheckObject<object>(obj, "basis");
    object& basis_obj = obj["basis"].get<object>();
    int MS = 1;
    string target_type = ReadJson<string>(target_obj, "type");
    if(target_type == "h_pi") {
      MS = 1;
    } else if(target_type == "custom") {
      CheckObject<object>(target_obj, "driv");
      object& driv_obj = target_obj["driv"].get<object>();
      string driv_type = ReadJson<string>(driv_obj, "type");
      if(driv_type == "STO")
	MS = 1;
      else if(driv_type == "GTO")
	MS = 2;
      else
	throw runtime_error("unsupported type in driv in target");
    }
    string basis_type = ReadJson<string>(basis_obj, "type");
    int MB = 1;
    if(basis_type == "STO") {
      MB = 1;
    } else if(basis_type == "GTO") {
      MB = 2;
    } else {
      cerr << "unsupported type in basis\n";
      exit(1);
    }
    
    // -- start main --
    if(MB == 1 && MS == 1) {
      R1Driv<1,1> r1driv(obj);
      r1driv.Run();
    } else if(MB == 2 && MS == 1) {
      R1Driv<2,1> r1driv(obj);
      r1driv.Run();
    } else {
      cerr << "not supported combination of basis and driv type.\n";
      exit(1);
    }
  } catch(std::exception& e) {
    cerr << "error on calculation\n";
    cerr << e.what() << endl;
    exit(1);
  }
  
  //  r1driv.in_file = argv[1];
  //  
  //  r1driv.Parse(obj);
  //  cout << "<<<< r1_driv <<<<\n";
}

