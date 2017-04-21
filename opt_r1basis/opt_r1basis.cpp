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
#include "../math/nderiv.hpp"
#include "../r1basis/r1basis.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace picojson;
using namespace cbasis;

/*
dcomplex CalcFDP1(dcomplex z) {
  int n(1);
  
  STOs sto = Create_STOs();  
  sto->AddPrim(2, z);
  sto->SetUp();

  LC_STOs driv = Create_LC_STOs();
  driv->Add(2.0, 2, 1.0);

  MatrixXcd M00(n, n);
  MatrixXcd D00(n, n);
  sto->CalcD2Mat(D00); D00 *= 0.5;
  sto->CalcRmMat(-1, M00); D00 += M00;
  sto->CalcRmMat(-2, M00); D00 += -M00;
  sto->CalcRmMat(0,  M00); D00 += 0.5*M00;

  VectorXcd m0(n);
  sto->CalcVec(driv, m0);

  VectorXcd c(n);
  LinearSolver linear_solver("householderQr");
  linear_solver.Solve(D00, m0, &c);

  return TDot(m0, c);
}
dcomplex CalcFDP1_dzeta(dcomplex z) {
  int n(1);
  
  STOs sto = Create_STOs();  
  sto->AddPrim(2, z);
  sto->SetUp();

  STOs dsto = Create_STOs();
  sto->DerivOneZeta(dsto);
  //  cout << dsto->str() << endl;
  LC_STOs driv = Create_LC_STOs();
  driv->Add(2.0, 2, 1.0);

  MatrixXcd M00(n, n);
  MatrixXcd D00(n, n); 
  sto->CalcD2Mat(D00); D00 *= 0.5; 
  sto->CalcRmMat(-1, M00); D00 += M00;
  sto->CalcRmMat(-2, M00); D00 += -M00;
  sto->CalcRmMat(0,  M00); D00 += 0.5*M00;
  MatrixXcd D10(n, n);
  CalcD2Mat<1,1>(dsto, sto, D10);     D10 *= 0.5;
  CalcRmMat<1,1>(dsto, -1, sto, M00); D10 += M00;
  CalcRmMat<1,1>(dsto, -2, sto, M00); D10 += -M00;
  CalcRmMat<1,1>(dsto, 0, sto, M00); D10 += 0.5*M00;
  
  VectorXcd m0(n);
  sto->CalcVec(driv, m0);

  VectorXcd c(n);
  LinearSolver linear_solver("householderQr");
  linear_solver.Solve(D00, m0, &c);

  return TDot(c, D10*c);
}
*/
void ZeroVector(VectorXcd &a) {
  int n(a.size());
  for(int i = 0; i < n; i++)
    a.coeffRef(i) = 0.0;
}
void ZeroMatrix(MatrixXcd &a) {
  int n(a.rows());
  int m(a.cols());
  for(int i = 0; i < n; i++)
    for(int j = 0; j < m; j++)
      a.coeffRef(i, j) = 0.0;
}


template<int MB, int MD>
class OptR1Driv {
public:
  typedef typename _EXPs<MB>::EXPs BasisEXPs;
  typedef typename _EXPs<MD>::LC_EXPs DrivLC_EXPs;

  object& obj;

  // -- Inputs --
  string comment;
  string mode;
  string in_file;
  string out_file;
  
  LinearSolver linear_solver;
  int max_iter;
  double conv_eps;
  
  dcomplex Z;
  int L;
  dcomplex E;
  BasisEXPs basis;  
  DrivLC_EXPs driv;

  bool write_matrix;

  // -- Intermediate --
  map<int, dcomplex> coef_op;
  vector<int> opt_idx;
  BasisEXPs d_basis;
  BasisEXPs d2_basis;
  MatrixXcd D00, D10, D20, D11, M00, M10, M11;
  VectorXcd S0, S1, S2, R0, R1, R2;

  MatrixXcd G;
  MatrixXcd Di, Dj, Dij;
  VectorXcd Si, Sj, Sij, Ri, Rj, Rij;
  VectorXcd C0, Ci, Cj, Cij;

  dcomplex val;
  VectorXcd grad;
  MatrixXcd hess;
  
  // -- Methods --
  OptR1Driv(object& _obj);
  void Run();
  void Parse();
  void PrintIn();
  void Calc();
  void Calc_Pre();
  void Calc_Mat();
  void Calc_GradHess();
  void Check();  
};

template<int MB, int MD>
OptR1Driv<MB, MD>::OptR1Driv(object& _obj):
  obj(_obj), basis(new _EXPs<MB>()),driv(new _LC_EXPs<MD>())
  , d_basis(new _EXPs<MB>()), d2_basis(new _EXPs<MB>()) {}
template<int MB, int MD> void OptR1Driv<MB, MD>::Run() {
  this->Parse();
  this->PrintIn();
  if(this->mode == "check") {
    this->Check();
  }
  if(this->mode == "run") {
    this->Calc();
  }    
}
template<int MB, int MD> void OptR1Driv<MB, MD>::Parse() {

  PrintTimeStamp("Parse", NULL);

  // -- parse values --
  try {
    this->comment = ReadJson<string>(obj, "comment");
    cout << format("comment: %s\n") % this->comment;
    this->mode = ReadJson<string>(obj, "mode");
    cout << format("mode: %s\n") % this->mode;
    if(mode != "run" && mode != "check") {
      throw runtime_error("mode <- {run, check}");
    }
    this->linear_solver = ReadJsonWithDefault
      <LinearSolver>(obj, "linear_solver", LinearSolver());
    cout << format("linear_solver: %s\n") % this->linear_solver.str();
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
    write_matrix = ReadJsonWithDefault<bool>
      (opt_obj, "write_matrix", false);
    cout << format("write_matrix: %s\n") % (write_matrix ? "yes" : "no");
  }

  // -- parse optimization method --
  try {
    CheckObject<object>(obj, "method");
    object& method_obj = obj["method"].get<object>();
    this->max_iter = ReadJson<int>(method_obj, "max_iter");
    cout << "max_iter: " << this->max_iter << endl;
    this->conv_eps = ReadJson<double>(method_obj, "eps");
    cout << "eps: " << conv_eps << endl;
  } catch(std::exception& e) {
    cerr << "error on parsing method\n";
    cerr << e.what() << endl;
    exit(1);      
  }
  
  // -- parse target --
  try {      
    CheckObject<object>(obj, "target");
    object& target_obj = obj["target"].get<object>();
    string target_type = ReadJson<string>(target_obj, "type");
    if(target_type == "custom") {
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
	 normalization_const != "sqrt" ) {
	throw runtime_error("normalization_const <- {none, sqrt}" );
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
	  if(normalization_const == "sqrt") {
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

  // -- Set Up --
  try {
    for(int i = 0; i < this->basis->size(); ++i) {
      this->opt_idx.push_back(i);
    }
    this->basis->SetUp();
    this->basis->DerivOneZeta(opt_idx, INITIAL, &this->d_basis);
    this->basis->DerivTwoZeta(opt_idx, INITIAL, &this->d2_basis);
    coef_op[0] = this->E;
    coef_op[-1] = this->Z;
    coef_op[-2] = -0.5*L*(L+1);
  } catch(std::exception& e) {
    cerr << "error on SetUp\n";
    cerr << e.what() << endl;
    exit(1);
  }
  
}
template<int MB, int MD> void OptR1Driv<MB, MD>::PrintIn() {
  PrintTimeStamp("PrintIn", NULL);

  cout << format("comment: %s\n") % this->comment;
  cout << format("in_file: %s\n") % this->in_file;
  cout << format("out_file: %s\n") % this->out_file;
  cout << format("linear_solver: %s\n")
    % this->linear_solver.str();

  cout << format("Z = (%f, %f)\n")
    % this->Z.real() % this->Z.imag();
  cout << format("L = %d\n") % this->L;
  cout << format("E = (%f, %f)\n")
    % this->E.real() % this->E.imag(); 
  cout << format("basis: \n%s") % this->basis->str();
  cout << format("driv: \n%s\n") % this->driv->str();
  cout << format("write_matrix: %s\n") %
    (this->write_matrix ? "yes" : "no");
}
template<int MB, int MD> void OptR1Driv<MB, MD>::Calc() {
  PrintTimeStamp("Calc", NULL);
  this->Calc_Pre();

  int n(this->basis->size());
  int no(this->d_basis->size());
  VectorXcd dz(no);
  bool conv_q = false;
  for(int it = 0; it < this->max_iter; it++) {
    this->Calc_Mat();
    this->Calc_GradHess();
    this->linear_solver.Solve(this->hess, this->grad, &dz);
    if(grad.norm()/n < this->conv_eps) {
      conv_q = true;      
      break;
    }
    for(int i = 0; i < n; i++) {
      this->basis->basis(i)->z(0) -= dz(i);
    }    
  }

  if(conv_q) {
    cout << "Convergence!" << endl;
  } else {
    cout << "Faled!" << endl;
  }

  cout << "Basis:" << endl;
  cout << this->basis->str();
}
template<int MB, int MD> void OptR1Driv<MB, MD>::Calc_Pre() {
  int n(this->basis->size());
  int no(this->d_basis->size());
  D00 = MatrixXcd::Zero(n, n);
  D10 = MatrixXcd::Zero(no, n);
  D20 = MatrixXcd::Zero(no, n);
  D11 = MatrixXcd::Zero(no, no);
  M00 = MatrixXcd::Zero(n, n);
  M10 = MatrixXcd::Zero(no, n);
  M11 = MatrixXcd::Zero(no, no);

  G = MatrixXcd::Zero(n, n);

  S0 =VectorXcd::Zero(n);
  S1 =VectorXcd::Zero(no);
  S2 =VectorXcd::Zero(no);
  R0 =VectorXcd::Zero(n);
  R1 =VectorXcd::Zero(no);
  R2 =VectorXcd::Zero(no);
  
  
  Di = MatrixXcd::Zero(n, n);
  Dj = MatrixXcd::Zero(n, n);
  Dij= MatrixXcd::Zero(n, n);

  Si =VectorXcd::Zero(n);
  Sj =VectorXcd::Zero(n);
  Sij=VectorXcd::Zero(n);
  Ri =VectorXcd::Zero(n);
  Rj =VectorXcd::Zero(n);
  Rij=VectorXcd::Zero(n);
  
  C0 =VectorXcd::Zero(n);
  Ci =VectorXcd::Zero(n);
  Cj =VectorXcd::Zero(n);
  Cij =VectorXcd::Zero(n);

  grad = VectorXcd::Zero(no);
  hess = MatrixXcd::Zero(no, no);
}
template<int MB, int MD> void OptR1Driv<MB, MD>::Calc_Mat() {
  
  try {

    this->basis->SetUp();
    this->basis->DerivOneZeta(opt_idx, REUSE, &d_basis);
    this->basis->DerivTwoZeta(opt_idx, REUSE, &d2_basis);

    /*
    CalcD2Mat_Numeric<MB,MB>(basis,     basis, M00); D00 = 0.5 * M00;
    CalcRmMat_Numeric<MB,MB>(basis, -2, basis, M00); D00 += (-L*(L+1)*0.5) * M00;
    CalcRmMat_Numeric<MB,MB>(basis, -1, basis, M00); D00 += (+Z * M00);
    CalcRmMat_Numeric<MB,MB>(basis,  0, basis, M00); D00 += (+E * M00);
    
    CalcD2Mat_Numeric<MB,MB>(d_basis, basis, M10);     D10 = 0.5 * M10;
    CalcRmMat_Numeric<MB,MB>(d_basis, -2, basis, M10); D10 += (-L*(L+1)*0.5) * M10;
    CalcRmMat_Numeric<MB,MB>(d_basis, -1, basis, M10); D10 += (+Z* M10);
    CalcRmMat_Numeric<MB,MB>(d_basis, +0, basis, M10); D10 += (+E* M10);

    CalcD2Mat_Numeric<MB,MB>(d_basis,   d_basis, M11); D11 = 0.5 * M11;
    CalcRmMat_Numeric<MB,MB>(d_basis,-2,d_basis, M11); D11 += (-L*(L+1)*0.5) * M11;
    CalcRmMat_Numeric<MB,MB>(d_basis,-1,d_basis, M11); D11 += (+Z* M11);
    CalcRmMat_Numeric<MB,MB>(d_basis,+0,d_basis, M11); D11 += (+E* M11);    

    CalcD2Mat_Numeric<MB,MB>(d2_basis,     basis, M10); D20 = 0.5 * M10;
    CalcRmMat_Numeric<MB,MB>(d2_basis, -2, basis, M10); D20 += (-L*(L+1)*0.5) * M10;
    CalcRmMat_Numeric<MB,MB>(d2_basis, -1, basis, M10); D20 += (+Z* M10);
    CalcRmMat_Numeric<MB,MB>(d2_basis, +0, basis, M10); D20 += (+E* M10);
    */

    CalcHAtomMat<MB,MB>(basis,   0.5, coef_op, basis,   REUSE, &D00);
    CalcHAtomMat<MB,MB>(d_basis, 0.5, coef_op, basis,   REUSE, &D10);
    CalcHAtomMat<MB,MB>(d_basis, 0.5, coef_op, d_basis, REUSE, &D11);
    CalcHAtomMat<MB,MB>(d2_basis,0.5, coef_op, basis,   REUSE, &D20);
    
    CalcVec_Numeric<MB,MD>(basis,   driv, S0);
    CalcVec_Numeric<MB,MD>(d_basis, driv, S1);
    CalcVec_Numeric<MB,MD>(d2_basis,driv, S2);
    CalcVec_Numeric<MB,MD>(basis,   driv, R0);
    CalcVec_Numeric<MB,MD>(d_basis, driv, R1);
    CalcVec_Numeric<MB,MD>(d2_basis,driv, R2);

    this->linear_solver.Inv(D00, &G);
    C0 = G * S0;
  } catch(std::exception& e) {
    cerr << "Error on calculaing grad/hess\n";
    cerr << e.what() << endl;
    exit(1);
  }
}
template<int MB, int MD> void OptR1Driv<MB, MD>::Calc_GradHess() {
  int n(this->basis->size());
  int no(this->d_basis->size());
  
  try {
    val = TDot(R0, C0);
    ZeroVector(grad);
    ZeroMatrix(hess);

    for(int io = 0; io < no; io++) {
      int i = io;
      ZeroVector(Si); Si(i) = S1(io);
      ZeroVector(Ri); Ri(i) = R1(io);
      ZeroMatrix(Di);
      for(int ii = 0; ii < n; ii++) {
	Di(i, ii) += D10(io, ii);
	Di(ii, i) += D10(io, ii);
      }
      Ci = -G*Di*C0 + G*Si;
      grad(io) = TDot(Ri, C0) + TDot(R0, Ci);

      for(int jo = 0; jo < no; jo++) {
	int j = jo;
	ZeroVector(Sj); Sj(j) = S1(jo);
	ZeroVector(Rj); Rj(j) = R1(jo);
	ZeroMatrix(Dj);
	for(int jj = 0; jj < n; jj++) {
	  Dj(j, jj) += D10(jo, jj);
	  Dj(jj, j) += D10(jo, jj);
	}
	Cj = -G*Dj*C0 + G*Sj;
	
	ZeroVector(Sij);
	ZeroVector(Rij);
	ZeroMatrix(Dij);
	if(i == j) {
	  Sij(i) = S2(io);
	  Rij(i) = R2(io);
	  for(int jj = 0; jj < n; jj++) {
	    Dij(jj, j) += D20(jo, jj);
	    Dij(j, jj) += D20(jo, jj);
	  }
	}
	Dij(i, j) += D11(io, jo);
	Dij(i, j) += D11(jo, io);

	ZeroVector(Cij);
	Cij = G*Dj*G*Di*C0 - G*Dij*C0 - G*Di*Cj - G*Dj*G*Si + G*Sij;
	hess(io, jo) = (TDot(Rj,G*Si) - TDot(Rj,G*Di*G*R0) + TDot(Ri,G*Sj)
			-TDot(R0,G*Dj*G*Si) + TDot(R0,G*Dj*G*Di*G*S0)
			-TDot(Ri,G*Dj*G*S0) - TDot(R0,G*Dij*G*S0)
			+TDot(R0,G*Di*G*Dj*G*S0) - TDot(R0,G*Di*G*Sj)
			+TDot(Rij,G*S0) + TDot(R0,G*Sij));
	hess(io, jo) = (+TDot(Rij, C0)
			+TDot(Ri,  Cj)
			+TDot(Rj,  Ci)
			+TDot(R0,  Cij));
      }
    }
  } catch(std::exception& e) {
    cerr << "Error on calculaing matrix\n";
    cerr << e.what() << endl;
    exit(1);
  }
}
template<int MB, int MD> void OptR1Driv<MB, MD>::Check() {
  PrintTimeStamp("Check", NULL);
  int idx(0), jdx(0);
  double eps(0.0001);

  BasisEXPs basis_00 = this->basis->Clone();
  this->Calc_Pre(); this->Calc_Mat(); this->Calc_GradHess();
  dcomplex alpha_0 = TDot(this->R0, this->C0);
  VectorXcd g_calc = this->grad;
  MatrixXcd h_calc = this->hess;

  // -- compute infinte small difference --
  BasisEXPs basis_p = basis_00->Clone();
  basis_p->basis(idx)->z(0) += eps; basis_p->SetUp();
  this->basis.swap(basis_p);
  this->Calc_Pre(); this->Calc_Mat(); this->Calc_GradHess(); 
  dcomplex alpha_p = TDot(this->R0, this->C0);
  VectorXcd grad_p = this->grad;

  BasisEXPs basis_m = basis_00->Clone();
  basis_m->basis(idx)->z(0) -= eps; basis_m->SetUp();  
  this->basis.swap(basis_m);
  this->Calc_Pre(); this->Calc_Mat(); this->Calc_GradHess(); 
  dcomplex alpha_m = TDot(this->R0, this->C0);
  VectorXcd grad_m = this->grad;
  
  dcomplex g_nume = (alpha_p - alpha_m)/(2.0*eps);
  //  dcomplex h_nume = (grad_p(jdx) - grad_m(jdx))/(2.0*eps);

  cout << g_calc(idx)      << " vs " << g_nume << endl;
  //  cout << h_calc(idx, jdx) << " vs " << h_nume << endl;
  if(idx == jdx) {
    dcomplex h_nume = (alpha_p + alpha_m - 2.0*alpha_0)/(eps*eps);
    cout << h_calc(idx, jdx) << " vs "  << h_nume << endl;
  }
}

int main(int argc, char *argv[]) {
  
  cout << ">>>> opt_r1basis >>>>\n";

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


  // -- Check function type --
  int MS = 1;
  int MB = 1;
  try {    
    CheckObject<object>(obj, "target");
    object& target_obj = obj["target"].get<object>();
    CheckObject<object>(obj, "basis");
    object& basis_obj = obj["basis"].get<object>();
    
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
    if(basis_type == "STO") {
      MB = 1;
    } else if(basis_type == "GTO") {
      MB = 2;
    } else {
      cerr << "unsupported type in basis\n";
      exit(1);
    }
  } catch(std::exception& e) {
    cerr << "error on function type\n";
    cerr << e.what() << endl;
    exit(1);
  }
  
  // -- start main --
  try{
    if(MB == 1 && MS == 1) {
      OptR1Driv<1,1> opt(obj); opt.Run();
    } else if(MB == 2 && MS == 1) {
      OptR1Driv<2,1> opt(obj); opt.Run();
    } else {
      cerr << "not supported combination of basis and driv type.\n";
      exit(1);
    }
  } catch(std::exception& e) {
    cerr << "error on calculation\n";
    cerr << e.what() << endl;
    exit(1);
  }  

  
  cout << "<<<< opt_r1basis <<<<\n";
}

