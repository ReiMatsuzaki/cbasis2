#include <fstream>
#include <stdio.h>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <boost/format.hpp>
#include <unsupported/Eigen/FFT>
#include "../../utils/timestamp.hpp"
#include "../../utils/eigen_plus.hpp"
#include "../../utils/read_json.hpp"
#include "../../math/nderiv.hpp"
#include "../../r1basis/r1basis.hpp"

using namespace std;
using namespace boost;
using namespace Eigen;
using namespace picojson;
using namespace cbasis;

enum ETimeStepMethod {
  EULER, RK4
};

int ParseBasisType(object& obj) {

  try{
    
    CheckObject<object>(obj, "basis");
    object& basis_obj = obj["basis"].get<object>();
    string basis_type = ReadJson<string>(basis_obj, "type");
    if(basis_type == "STO") {
      return 1;
    } else if(basis_type == "GTO") {
      return 2;
    }
    throw runtime_error("unsupported type");
  } catch(std::exception& e) {
    cerr << "error on parsing basis type\n";
    cerr << e.what();    
  }
  exit(1);
}

class IWindowFunc {
  virtual double at(double t) = 0;
};
class WindowFuncATan : public IWindowFunc {
public:
  double t0_, a_;
  WindowFuncATan(double _t0, double _a): t0_(_t0), a_(_a) {}
  double at(double t) {
    return 0.5 + atan(a_*(-t + t0_)) / M_PI;
  }
};

template<int MB> class TDVar {
public:
  typedef typename _EXPs<MB>::EXPs EXPs;
  typedef boost::function<dcomplex (dcomplex)> CFunc;

  object& obj;

  // -- const --
  static const double au2ev = 27.2114;
  static const double c_light = 137.035999139;
  static const double au2mb = 5.291772 * 5.291772;
  

  // -- Inputs --
  string comment;
  string mode;
  string in_file;
  string out_file;

  bool write_matrix;
  LinearSolver linear_solver;
  string czs_out;
  FILE *czs_fp;
  string psi_out;
  FILE *psi_fp;
  string ac_out;
  FILE *ac_fp;
  FILE *cs_fp;

  bool is_hermite;
  double dt;
  double t1;
  int numt;
  ETimeStepMethod time_step;
  WindowFuncATan *window_func;

  int L;
  double Z;
  double E0;

  VectorXcd cs;  
  EXPs basis;
  vector<int> opt_index_list;

  // -- Results --
  MatrixXcd czs_t;
  vector<dcomplex> ac_t, window_ac_t;
  vector<dcomplex> ac_w;
  vector<double> cross_sec_w;

  // -- Intermediate --
  map<int, dcomplex> coef_op;
  EXPs bra_basis;
  EXPs d_basis;
  EXPs bra_d_basis;
  
  MatrixXcd S00, S10, S11, H00, H10, X;
  VectorXcd H00c_H10c;  // [H00c, H10c], RHE of update equation
  VectorXcd m0;       // temporary

  // -- Methods --
  TDVar(object& _obj);
  void Run();
  void Parse();
  void PrintIn();
  void Calc();
  void Calc_Euler();
  void Calc_RK4();
  void Calc_UpdateVec(VectorXcd *vec);
  void UpdateBasis(int reuse);
  void ShiftCoefZeta(const VectorXcd& dcz);
  void PrintOut();
  void Check();    
};

template<int MB> TDVar<MB>::TDVar(object& _obj):
  obj(_obj), basis(new _EXPs<MB>()), d_basis(new _EXPs<MB>()) {}
template<int MB> void TDVar<MB>::Run() {
  this->Parse();
  this->PrintIn();
  if(mode == "calc") {
    this->Calc();
    this->PrintOut();
  } else if(mode == "check") {
    this->Check();
  }
}
template<int MB> void TDVar<MB>::Parse() {
  
  PrintTimeStamp("Parse", NULL);

  // -- parse values --
  try {
    this->comment = ReadJson<string>(obj, "comment");
    cout << format("comment: %s\n") % this->comment;
    this->mode = ReadJson<string>(obj, "mode");
    cout << format("mode: %s\n") % this->mode;
    if(mode != "calc" && mode != "check") {
      throw runtime_error("mode <- {calc, check}");
    }
    
  } catch(std::exception& e) {
    cerr << "error on parsing json\n";
    cerr << e.what() << endl;
    exit(1);
  }

  // -- parse options --
  bool find_opt = (obj.find("opts") != obj.end());
  try {
    if(find_opt) {
      CheckObject<object>(obj, "opts");
      object& opt_obj = obj["opts"].get<object>();

      string write_matrix_str =
 	ReadJsonWithDefault<string>(opt_obj, "write_matrix", "f");
      cout << format("write_matrix: %s\n") % write_matrix_str;
      this->write_matrix = (write_matrix_str == "t");
      
      this->linear_solver = ReadJsonWithDefault
	<LinearSolver>(opt_obj, "linear_solver", LinearSolver());
      cout << format("linear_solver: %s\n") % this->linear_solver.str();    

      this->czs_out = ReadJsonWithDefault<string>(opt_obj, "czs_out", "");
      cout << format("czs_out: %s\n") % this->czs_out;
      czs_fp = NULL;
      if(czs_out != "") {
	if((czs_fp = fopen(this->czs_out.c_str(), "w")) == NULL) {
	  cerr << "failed to open czs_out file. File name is " << this->czs_out << endl;
	  exit(1);
	}
      }
      
      this->psi_out = ReadJsonWithDefault<string>(opt_obj, "psi_out", "");
      cout << format("psi_out: %s\n") % this->psi_out;
      psi_fp = NULL;
      if(this->psi_out != "") {
	if((psi_fp = fopen(this->psi_out.c_str(), "w")) == NULL) {
	  cerr << "failed to open czs_out file. File name is " << this->psi_out << endl;
	  exit(1);
	}
      }

      this->ac_out = ReadJsonWithDefault<string>(opt_obj, "ac_out", "");
      cout << format("ac_out: %s\n") % this->ac_out;
      ac_fp = NULL;
      if(this->ac_out != "") {
	if((ac_fp = fopen(this->ac_out.c_str(), "w")) == NULL) {
	  cerr << "failed to open ac_out file. File name is " << ac_out << endl;
	  exit(1);
	}
      }

      string cs_out = ReadJsonWithDefault<string>(opt_obj, "cs_out", "");
      cout << format("cs_out: %s\n") % cs_out;
      this->cs_fp = NULL;
      if(cs_out != "") {
	if((cs_fp = fopen(cs_out.c_str(), "w")) == NULL) {
	  cerr << "failed to open cs_out file. File name is " << cs_out << endl;
	  exit(1);
	}
      }      
    }
  } catch(std::exception& e) {
    cerr << "error on parsing options\n";
    cerr << e.what() << endl;
    exit(1);
  }

  // -- method --
  try {
    CheckObject<object>(obj, "method");
    object& method_obj = obj["method"].get<object>();

    string ip_str = ReadJson<string>(method_obj, "inner_product");
    cout << format("inner_product: %s\n") % ip_str;
    if(ip_str == "hermite") {
      this->is_hermite = true;
    } else if(ip_str == "complex") {
      this->is_hermite = false;
    } else {
      throw runtime_error("inner_product <- {complex, hermite}");
    }
    
    dt = ReadJson<double>(method_obj, "dt");
    cout << format("dt: %f\n") % dt;
    
    t1 = ReadJson<double>(method_obj, "t1");
    cout << format("t1: %f\n") % t1;

    numt = (int)(t1/dt);

    string time_step_name = ReadJson<string>(method_obj, "time_step");
    if(time_step_name == "Euler") {
      this->time_step = EULER;
    } else if(time_step_name == "RK4") {
      this->time_step = RK4;
    } else {
      throw runtime_error("time_step <- {Euler, RK4}");
    }
    cout << format("time_step: %s\n") % time_step_name;

    this->window_func = NULL;
    if(method_obj.find("window_func") != method_obj.end()) {
      CheckObject<object>(method_obj, "window_func");
      object& window_obj = method_obj["window_func"].get<object>();
      string window_type = ReadJson<string>(window_obj, "type");
      cout << format("window_type: %s\n") % window_type;
      if(window_type == "arctan") {
	double t0 = ReadJson<double>(window_obj, "t0");
	cout << format("t0: %f\n") % t0;
	double a = ReadJson<double>(window_obj, "a");
	cout << format("a: %f\n") % a;
	this->window_func = new WindowFuncATan(t0, a);
      }
    }
    
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
      CheckObject<object>(target_obj, "psi0");
      object& obj0 = target_obj["psi0"].get<object>();
      string target_psi0_type = ReadJson<string>(obj0, "type");
      if(target_psi0_type == "basis") {
      } else {
	throw runtime_error("unsupported type in target/psi0");
      }

      this->E0 = ReadJson<double>(target_obj, "E0");
      cout << format("E0: %f\n") % this->E0;
      
      CheckObject<object>(target_obj, "pot");
      object& obj_pot = target_obj["pot"].get<object>();
      string target_pot_type = ReadJson<string>(obj_pot, "type");
      if(target_pot_type == "hatom") {
	L = ReadJson<int>(obj_pot, "L");
	cout << format("L: %d\n") % L;
	Z = ReadJson<double>(obj_pot, "Z");
	cout << format("Z: %f\n") % Z;
      } else {
	throw runtime_error("unsupported type in target/pot");
      }
    } else {
      throw runtime_error("unsupported type in target");
    }    
  } catch(std::exception& e) {
    cerr << "error on parsing target\n";
    cerr << e.what() << endl;
    exit(1);
  } 
  
  // -- parse basis --
  try {

    bool coef_include_normal_const = false;

    // -- get object --
    CheckObject<object>(obj, "basis");
    object& basis_obj = obj["basis"].get<object>();

    // -- options --
    if(basis_obj.find("opts") != basis_obj.end()) {
      CheckObject<object>(basis_obj, "opts");
      object& opt_obj = basis_obj["opts"].get<object>();
      string normal_const_str = ReadJson<string>(opt_obj, "include_normalization_const");
      if(normal_const_str == "t") {
	coef_include_normal_const = true;
      } else if(normal_const_str == "f") {
	coef_include_normal_const = false;
      } else {
	throw runtime_error("include_normalization_const <- {t, f}" );
      }
    }
    
    // -- value --
    CheckObject<picojson::array>(basis_obj, "value");
    picojson::array& basis_vals = basis_obj["value"].get<picojson::array>();
    vector<dcomplex> cs_vector;
    for(int idx = 0; idx < (int)basis_vals.size(); idx++) {
      value& val(basis_vals[idx]);
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
	dcomplex c = ReadJson<dcomplex>(obj0, "c");
	cout << "c: " << c << endl;	
	
	dcomplex z = ReadJson<dcomplex>(obj0, "z");
	if(coef_include_normal_const) {
	  this->basis->AddNotNormalPrim(1.0, n, z);
	} else {
	  this->basis->AddPrim(n, z);
	}
	cout << "z: " << z << endl;
	string optq = ReadJson<string>(obj0, "opt");
	cout << format("opt: %d\n") % optq;
	if(optq == "t") {
	  opt_index_list.push_back(idx);
	} else {
	  throw runtime_error("unsupported basis/value/opt");
	}
	cs_vector.push_back(c);
      } else {
	throw runtime_error("key \"z\" are necessary in value in basis");
      }
    }
    this->cs = Map<VectorXcd>(&cs_vector[0], cs_vector.size());
    //this->zs = Map<VectorXcd>(&zs_vector[0], zs_vector.size());
  } catch(std::exception& e) {
    cerr << "error on parsing basis\n";
    cerr << e.what() << endl;
    exit(1);
  }

  // -- Set Up --
  try {
    // basis update
    this->UpdateBasis(INITIAL);

    // number of data
    int numt = (int)(t1/dt);    
    int num(basis->size());
    int numo(d_basis->size());

    // store matrix
    S00 = MatrixXcd::Zero(num, num);
    S10 = MatrixXcd::Zero(numo, num);
    S11 = MatrixXcd::Zero(numo, numo);
    H00 = MatrixXcd::Zero(num, num);
    H10 = MatrixXcd::Zero(numo, num);
    czs_t = MatrixXcd::Zero(num*2, numt);
    X = MatrixXcd::Zero(num*2, num*2);
    H00c_H10c = VectorXcd::Zero(num*2);
    m0 = VectorXcd::Zero(num);

    // other
    coef_op[-1] = -Z; coef_op[-2] = 0.5*L*(L+1);
    for(int i = 0; i < num; i++) {
      czs_t(i, 0)     = this->cs(i);
      czs_t(i+num, 0) = this->basis->basis(i)->z(0);
    }

  } catch(std::exception& e) {
    cerr << "error on SetUp\n";
    cerr << e.what() << endl;
    exit(1);    
  }
  
}
template<int MB> void TDVar<MB>::PrintIn() {
  PrintTimeStamp("PrintIn", NULL);

  cout << "td: " << dt << endl;
  cout << "t1: " << t1 << endl;
  cout << "numt: " << numt << endl;
  
  cout << "cs: \n" << this->cs << endl;
  cout << format("basis: \n%s\n") % basis->str();
  cout << format("d_basis: \n%s\n") % d_basis->str();
}
template<int MB> void TDVar<MB>::Calc() {
  PrintTimeStamp("Calc", NULL);

  int num(this->basis->size());
  VectorXcd c0 = this->cs;
  typename _EXPs<MB>::EXPs bra_basis0;
  this->basis->Conj(INITIAL, &bra_basis0);

  PrintTimeStamp("TimeStep", NULL);
  if(time_step == EULER) {
    this->Calc_Euler();
  } else if(time_step == RK4) {
    this->Calc_RK4();
  }

  PrintTimeStamp("ACF", NULL);
  for(int idx_t = 0; idx_t < numt; idx_t++) {
    double t = idx_t * dt;
    for(int i = 0; i < num; i++) {
      this->cs(i) = czs_t(i, idx_t);
      this->basis->basis(i)->z(0) = czs_t(i+num, idx_t);
    }
    this->UpdateBasis(REUSE);      
    CalcRmMat<MB, MB>(bra_basis0, 0, basis, REUSE, &S00);
    dcomplex ac = exp(dcomplex(0, E0*t)) * TDot(c0, S00*cs);
    ac_t.push_back(ac);
    if(this->window_func) {
      window_ac_t.push_back(ac * this->window_func->at(t));
    } else {
      window_ac_t.push_back(ac);
    }
  }

  PrintTimeStamp("FFT", NULL);
  /*
    see td_var/test/fft for learning the behavior of FFT<double>
   */
  FFT<double> fft;  
  fft.inv(ac_w, window_ac_t);
  for(int idx_w = 0; idx_w < numt; idx_w++) {
    double w = 2.0*M_PI * idx_w / t1;    
    double integrand = ac_w[idx_w].real() * 2.0 * this->t1;
    double cross_sec = 2.0*M_PI*w/c_light * integrand * au2mb / 3.0;
    this->cross_sec_w.push_back(cross_sec);
  }
}
template<int MB> void TDVar<MB>::Calc_Euler() {
  
  int num(basis->size());
  VectorXcd dcz(num*2);

  for(int idx_t = 1; idx_t < numt; idx_t++) {
    
    this->Calc_UpdateVec(&dcz);
    this->ShiftCoefZeta(dcz);
    this->UpdateBasis(REUSE);
    
    for(int i = 0; i < num; i++) {
      czs_t(i,     idx_t) = this->cs(i);
      czs_t(i+num, idx_t) = this->basis->basis(i)->z(0);
    }
  }
  
}
template<int MB> void TDVar<MB>::Calc_RK4() {
  
  int num(basis->size());
  VectorXcd c4s(num), c2s(num), c3s(num);
  VectorXcd z4s(num), z2s(num), z3s(num);
  VectorXcd dcz_1(num*2), dcz_1_half(num*2), m_dcz_1_half(num*2);
  VectorXcd dcz_2(num*2), dcz_2_half(num*2), m_dcz_2_half(num*2);
  VectorXcd dcz_3(num*2), m_dcz_3(num*2);
  VectorXcd dcz_4(num*2);
  VectorXcd dcz(num*2);

  for(int idx_t = 1; idx_t < this->numt; idx_t++) {

    this->Calc_UpdateVec(&dcz_1);
    dcz_1_half = dcz_1 * 0.5;
    m_dcz_1_half = -dcz_1_half;
    
    this->ShiftCoefZeta(dcz_1_half);
    this->UpdateBasis(REUSE);
    this->Calc_UpdateVec(&dcz_2);
    dcz_2_half = dcz_2 * 0.5;
    m_dcz_2_half = -dcz_2_half;
    
    this->ShiftCoefZeta(m_dcz_1_half);
    this->ShiftCoefZeta(dcz_2_half);
    this->UpdateBasis(REUSE);
    this->Calc_UpdateVec(&dcz_3);
    m_dcz_3 = -dcz_3;
    
    this->ShiftCoefZeta(m_dcz_2_half);
    this->ShiftCoefZeta(dcz_3);
    this->UpdateBasis(REUSE);
    this->Calc_UpdateVec(&dcz_4);
    
    this->ShiftCoefZeta(m_dcz_3);
    dcz = (dcz_1 + 2.0*dcz_2 + 2.0*dcz_3 + dcz_4) / 6.0;
    this->ShiftCoefZeta(dcz);
    this->UpdateBasis(REUSE);

    for(int i = 0; i < num; i++) {
      this->czs_t(i,     idx_t) = this->cs(i);
      this->czs_t(i+num, idx_t) = this->basis->basis(i)->z(0);
    }    
  }
  
}
template<int MB> void TDVar<MB>::Calc_UpdateVec(VectorXcd *vec) {
  
  int num(this->basis->size());
  
  // -- Update matrix --
  CalcHAtomMat<MB, MB>(bra_basis, -0.5, coef_op, basis, REUSE, &H00);
  CalcHAtomMat<MB, MB>(bra_d_basis, -0.5, coef_op, basis, REUSE, &H10);
  CalcRmMat<MB, MB>(bra_basis,   0, basis, REUSE, &S00);
  CalcRmMat<MB, MB>(bra_d_basis, 0, basis, REUSE, &S10);
  CalcRmMat<MB, MB>(bra_d_basis, 0, d_basis, REUSE, &S11);

  if(this->write_matrix) {
    cout << "H00:\n" << H00 << endl;
    cout << "H10:\n" << H10 << endl;
    cout << "S00:\n" << S00 << endl;
    cout << "S10:\n" << S10 << endl;
    cout << "S11:\n" << S11 << endl;
  }

  /*
    [S01 * diag(c)]_{ij} = sum_{k} S01_{ik} diag(c)_{kj} 
    .                    = sum_{k} S01_{ik} c_j delta_{kj} 
    .                    = S01_{ij} c_j
    .                    = S10_{ji} c_j
    [S11 * diag(c)]_{ij} = sum_{k} S11_{ik} diag(c)_{kj} 
    .                    = S11_{ij} c_j
   */
  
  m0 = H00*cs;
  for(int i = 0; i < num; i++)
    H00c_H10c(i) = m0(i);
  m0 = H10*cs;
  for(int i = 0; i < num; i++)
    H00c_H10c(i+num) = m0(i);    
  for(int i = 0; i < num; i++) {      
    for(int j = 0; j < num; j++) {
      X(i,     j)     = S00(i, j);
      X(i+num, j)     = S10(i, j);
      X(i, num+j)     = S10(j, i) * cs(j);
      X(i+num, j+num) = S11(i, j) * cs(j);
    }
  }
  linear_solver.Solve(X, H00c_H10c, vec);
}
template<int MB> void TDVar<MB>::UpdateBasis(int reuse) {
  /**
     from zs, build basis functions and its derivative.
   */
  basis->SetUp();
  this->basis->DerivOneZeta(opt_index_list, reuse, &this->d_basis);
  if(this->is_hermite) {
    this->basis->Conj(reuse, &bra_basis);
  } else {
    this->basis->Clone(reuse, &bra_basis);
  }
  this->bra_basis->DerivOneZeta(opt_index_list, reuse, &bra_d_basis);
  
}
template<int MB> void TDVar<MB>::ShiftCoefZeta(const VectorXcd& dcz) {
  int numo = this->opt_index_list.size();
  if(dcz.size() != 2*numo) {
    throw runtime_error("size mismatch");
  }  
  dcomplex ii(0, 1);
  for(int io = 0; io < numo; io++) {
    int i = this->opt_index_list[io];
    this->cs[i]                 += -ii * dt * dcz[i];
    this->basis->basis(i)->z(0) += -ii * dt * dcz[i+numo];
  }
}
template<int MB> void TDVar<MB>::PrintOut() {

  PrintTimeStamp("PrintOut", NULL);

  int num = basis->size();

  if(czs_fp != NULL) {
    cout << "write czs file.\n";
    fprintf(czs_fp, "t");
    for(int i = 0; i < num; i++) {
      fprintf(czs_fp, ",re_c%d,im_c%d,re_z%d,im_z%d", i, i, i, i);
    }
    fprintf(czs_fp, "\n");

    
    for(int it = 0; it < numt; it++) {
      double t = it * dt;
      fprintf(czs_fp, "%f", t);
      for(int i = 0; i < num; i++) {
	dcomplex c = czs_t(i,     it);
	dcomplex z = czs_t(i+num, it);
	fprintf(czs_fp, ",%f,%f,%f,%f", c.real(), c.imag(), z.real(), z.imag());
      }
      fprintf(czs_fp, "\n");      
    }
  }

  if(ac_fp != NULL) {
    cout << "write ac file.\n";
    fprintf(ac_fp, "t,re,im,re_w,im_w\n");
    for(int it = 0; it < numt; it++) {
      double t = it * dt;
      dcomplex ac = ac_t[it];
      dcomplex w_ac= window_ac_t[it];
      fprintf(ac_fp, "%f,%f,%f,%f,%f\n",
	      t, ac.real(), ac.imag(), w_ac.real(), w_ac.imag());
    }
  }

  if(cs_fp != NULL) {
    cout << "write cs file.\n";
    fprintf(cs_fp, "w,cs\n");
    for(int iw = 0; iw < numt; iw++) {
      double w = 2.0*M_PI / t1 * iw;
      double cross_sec = this->cross_sec_w[iw];
      fprintf(cs_fp, "%f,%f\n", w, cross_sec);
    }
  }
  
}
template<int MB> void TDVar<MB>::Check() {
  PrintTimeStamp("Check", NULL);
}
int main(int argc, char *argv[]) {
  cout << ">>>> fd_var >>>>\n";

  // argument check
  if(argc != 2) {
    cerr << "one argument is necessary\n";
    exit(1);
  }

  // open file
  ifstream f(argv[1]);
  if(f.fail()) {
    cerr << "failed to open file\n";
    exit(1);
  }
  cout << format("in_file: %s\n") % argv[1];

  // parse json
  value json; f >> json;
  if(not json.is<object>()) {
    cerr << "input json file is not object\n";
    exit(1);
  }
  object& obj = json.get<object>();
  
  int MB = ParseBasisType(obj);
  if(MB == 1) {
    TDVar<1> var(obj); var.Run();
  } else if(MB == 2) {
    TDVar<2> var(obj); var.Run();
  }
  
  cout << "<<<< fd_var <<<<\n";
}
