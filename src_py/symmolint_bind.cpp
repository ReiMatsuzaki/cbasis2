#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
// #include <boost/numpy.hpp>
#include <Eigen/Core>

#include "../utils/macros.hpp"
#include "../utils/eigen_plus.hpp"
#include "../src_cpp/angmoment.hpp"
#include "../src_cpp/bmatset.hpp"
#include "../src_cpp/one_int.hpp"
#include "../src_cpp/symmolint.hpp"


namespace {
  using namespace boost::python;
  using namespace cbasis;
  using namespace std;
  using namespace Eigen;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  // namespace np = boost::numpy;
}

template<class T>
boost::python::list PyListFromCppVector(vector<T>& cpp_vector) {
  boost::python::list py_list;
  typename vector<T>::const_iterator it;
  for(it = cpp_vector.begin(); it != cpp_vector.end(); ++it) {
    py_list.append(*it);
  }
  return py_list;
}
void printVectorXi(VectorXi& a) {
  std::cout << a << std::endl;
}

MatrixXcd& BMat_get(BMat* bmat, tuple args) {
  int i = extract<int>(args[0]);
  int j = extract<int>(args[1]);
  return (*bmat)[make_pair(i, j)];
}
void BMat_set(BMat* bmat, int i, int j, MatrixXcd& mat) {
  (*bmat)[make_pair(i, j)] = mat;
}
void BMatSetZero(BMat& C, BMat& mat) {

  for(int I = 0; I < (int)C.size(); I++) {
    pair<Irrep, Irrep> II(I, I);
    MatrixXcd& CII = C[II];
    mat[II] = MatrixXcd::Zero(CII.rows(), CII.cols());
  }
  
}
void BMatSetZero(SymGTOs g0, SymGTOs g1, BMat& bmat) {

  //
  // g0 and g1 are assumed to have same symmetry
  //

  if(! g0->sym_group()->IsSame(g1->sym_group())) {
    string msg; SUB_LOCATION(msg);
    msg += ": symmetry group of g0 and g1 are different";
    throw runtime_error(msg);
  }
  
  for(int irrep = 0; irrep < g0->sym_group()->order(); irrep++) {
    
    int n0(g0->size_basis_isym(irrep));
    int n1(g1->size_basis_isym(irrep));
    bmat[make_pair(irrep, irrep)] = MatrixXcd::Zero(n0, n1);
  }

}

void AddAngmoemnt() {

  def("cg_coef", cg_coef);

}

void AddSymmetryGroup() {

  register_ptr_to_python<SymmetryGroup>();
  class_<_SymmetryGroup>("SymmetryGroup", no_init)
    .add_property("order", &_SymmetryGroup::order)
    .add_property("name", &_SymmetryGroup::name)
    .add_property("irrep_s", &_SymmetryGroup::irrep_s)
    .add_property("irrep_x", &_SymmetryGroup::irrep_x)
    .add_property("irrep_y", &_SymmetryGroup::irrep_y)
    .add_property("irrep_z", &_SymmetryGroup::irrep_z)
    .def("get_irrep", &_SymmetryGroup::GetIrrep)
    .def("get_irrep_name", &_SymmetryGroup::GetIrrepName)
    .def("__str__", &_SymmetryGroup::str)
    .def("__repr__", &_SymmetryGroup::str);
  def("C1",  SymmetryGroup_C1);
  def("Cs",  SymmetryGroup_Cs);
  def("C2h", SymmetryGroup_C2h);
  def("D2h", SymmetryGroup_D2h);
  def("C4",  SymmetryGroup_C4);

}

MatrixXcd* CanonicalMatrix_py(const MatrixXcd& S, double eps) {

  MatrixXcd* X = new MatrixXcd();
  CanonicalMatrix(S, eps, X);
  return X;
}

tuple CEigenSolveCanonicalNum_py(const CM& F, const CM& S, int num0) {

  MatrixXcd C;
  VectorXcd eig;
  CEigenSolveCanonicalNum(F, S, num0, &C, &eig);
  return make_tuple(eig, C);

}

tuple generalizedComplexEigenSolve_py(const CM& F, const CM& S) {

  MatrixXcd C;
  VectorXcd eig;
  generalizedComplexEigenSolve(F, S, &C, &eig);
  return make_tuple(eig, C);
}

void AddLinearAlgebra() {

  def("canonical_matrix", CanonicalMatrix_py,
      return_value_policy<manage_new_object>());
  def("ceig", generalizedComplexEigenSolve_py);
  def("ceig_canonical", CEigenSolveCanonicalNum_py);

  class_<vector<int> >("vector_i")
    .def(vector_indexing_suite<vector<int> >());

  class_<BVec>("BVec", no_init);

  class_<BMat>("BMat", init<>())
    .def("write", BMatWrite)
    .def("read", BMatRead)
    .def("set_matrix", BMat_set)
    .def("__getitem__", BMat_get, return_internal_reference<>());


  // -- old --
  register_ptr_to_python<BMatSet>();
  class_<_BMatSet>("BMatSet", init<int>())
    .def("get_matrix",  &_BMatSet::GetMatrix, return_internal_reference<>())
    .def("__getitem__", &_BMatSet::GetBlockMatrix, return_internal_reference<>())
    .def("set_matrix", &_BMatSet::SetMatrix);

}

tuple SymGTOs_AtR_Ylm(SymGTOs gtos,
		      int L, int M,  int irrep,
		      const VectorXcd& cs_ibasis,
		      VectorXcd rs) {
  
  VectorXcd* ys = new VectorXcd();
  VectorXcd* ds = new VectorXcd();
  gtos->AtR_Ylm(L, M, irrep, cs_ibasis, rs, ys, ds);
  return make_tuple(ys, ds);
  
}
VectorXcd* SymGTOs_CorrectSign(SymGTOs gtos, int L, int M,  int irrep,
			       const VectorXcd& cs_ibasis) {
  VectorXcd* cs = new VectorXcd(cs_ibasis); // call copy constructor
  gtos->CorrectSign(L, M, irrep, *cs);
  return cs;
}
void AddSymMolInt() {

  def("print_vec", printVectorXi);

  register_ptr_to_python<Atom>();
  class_<_Atom>("_atom", init<string, dcomplex>())
    .def("add", (Atom(_Atom::*)(Vector3cd))&_Atom::Add, return_self<>())
    .def("add", (Atom(_Atom::*)(dcomplex, dcomplex, dcomplex))&_Atom::Add, return_self<>());
  def("new_Atom", (Atom(*)(string, dcomplex))NewAtom);
  //  Atom NewAtom(std::string _name, dcomplex _q);
  
  register_ptr_to_python<Molecule>();
  class_<_Molecule>("_molecule", init<SymmetryGroup>())
    .def("__str__", &_Molecule::show)
    .def("__repr__", &_Molecule::str)
    .def("sym_group", &_Molecule::sym_group)
    .def("atom", &_Molecule::atom)
    .def("add", &_Molecule::Add, return_self<>());
  
  def("new_Molecule", NewMolecule);
    
  //  def("molecule", (Molecule(*)(SymmetryGroup))NewMolecule);
  
  class_<Reduction>("Reduction", init<int, MatrixXcd>())
    .add_property("irrep",        &Reduction::irrep)
    .def("coef_iat_ipn", &Reduction::get_coef_iat_ipn,
	 return_internal_reference<>())
    .def("__str__", &Reduction::str)
    .def("__repr__", &Reduction::str);

  class_<SubSymGTOs>("SubSymGTOs", init<SymmetryGroup, Atom>())
    .def("ns",  (SubSymGTOs&(SubSymGTOs::*)(Vector3i))&SubSymGTOs::AddNs,  return_self<>())
    .def("rds", &SubSymGTOs::AddRds,  return_self<>())
    .def("zeta",&SubSymGTOs::AddZeta, return_self<>())
    .def("get_zeta", &SubSymGTOs::zeta)
    .def("x", &SubSymGTOs::x)
    .def("y", &SubSymGTOs::y)
    .def("z", &SubSymGTOs::z)
    .def("nx", &SubSymGTOs::nx)
    .def("ny", &SubSymGTOs::ny)
    .def("nz", &SubSymGTOs::nz)
    .def("__str__", &SubSymGTOs::str)
    .def("__repr__", &SubSymGTOs::str);
  def("sub_mono", Sub_Mono);
  def("sub_solid_sh", Sub_SolidSH_M);
  def("sub_solid_sh", Sub_SolidSH_Ms);


  register_ptr_to_python<SymGTOs>();
  class_<_SymGTOs>("ClassSymGTOs", no_init)
    .add_property("sym_group",
		  (SymmetryGroup(_SymGTOs::*)())&_SymGTOs::sym_group,
		  &_SymGTOs::set_sym_group)
    .add_property("molecule",
		  (Molecule(_SymGTOs::*)())&_SymGTOs::molecule,
		  &_SymGTOs::set_molecule)
    .def("sub", &_SymGTOs::AddSub, return_self<>())
    .def("clone", &_SymGTOs::Clone)
    .def("conj", &_SymGTOs::Conj)
    .def("setup", &_SymGTOs::SetUp, return_self<>())
    .def("at_r_ylm_cpp", SymGTOs_AtR_Ylm)    
    .def("correct_sign", SymGTOs_CorrectSign,
	 return_value_policy<manage_new_object>())
    .def("__str__", &_SymGTOs::show)
    .def("__repr__", &_SymGTOs::str);
  def("new_SymGTOs", NewSymGTOs);
}

void AddOneTwoInt() {

  class_<ERIMethod>("ERI_method", init<>())
    .def("use_symmetry",&ERIMethod::set_symmetry,
	 return_self<>())
    .def("coef_R_memo", &ERIMethod::set_coef_R_memo,
	 return_self<>());

  def("calc_mat_complex", CalcMat_Complex);
  def("calc_mat_hermite", CalcMat_Hermite);
  def("calc_mat", CalcMat);
  
  //  def("calc_ERI_complex", CalcERI_Complex);
  //  def("calc_ERI_hermite", CalcERI_Hermite);
  //  def("calc_ERI", CalcERI);

}


BOOST_PYTHON_MODULE(symmolint_bind) {

  Py_Initialize();
  def("CalcCs",  SymmetryGroup_Cs);
  // np::initialize();
  AddAngmoemnt();
  //  AddMO();
  AddSymmetryGroup();
  AddLinearAlgebra();
  //  AddB2EInt();
  AddSymMolInt();
  AddOneTwoInt();

}

