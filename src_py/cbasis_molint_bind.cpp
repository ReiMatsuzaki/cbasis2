#include <iostream>
#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include "../eigen_plus.hpp"
#include "../molint.hpp"
using namespace std;

namespace {
  using namespace boost::python;
  using namespace cbasis;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  namespace np = boost::numpy;
  using namespace std;

}

np::ndarray array_transform(dcomplex* xs, int ni, int nj) {  

  if(ni < 1 || nj < 1) {
    std::string msg; SUB_LOCATION(msg);
    msg += "ni and nj must be bigger than 0";
    throw std::runtime_error(msg);
  }

  if(ni > 1 && nj > 1) {
    return np::from_data(xs,
			 np::dtype::get_builtin<dcomplex>(),
			 bp::make_tuple(ni, nj),
			 bp::make_tuple(sizeof(dcomplex), sizeof(dcomplex)*ni),
			 bp::object());
  } else if(ni == 1) {
    return np::from_data(xs,
			 np::dtype::get_builtin<dcomplex>(),
			 bp::make_tuple(nj),
			 bp::make_tuple(sizeof(dcomplex)),
			 bp::object());
  } else if(nj == 1) {
    return np::from_data(xs,
			 np::dtype::get_builtin<dcomplex>(),
			 bp::make_tuple(ni),
			 bp::make_tuple(sizeof(dcomplex)),
			 bp::object());
  }

  std::string msg; SUB_LOCATION(msg);
  msg += "ni and nj must be bigger than 0";
  throw std::runtime_error(msg);
}
np::ndarray MatrixSet_get(MatrixSet* mset, std::string label) {
  
  int ni(mset->size_basis_i());
  int nj(mset->size_basis_j());
  dcomplex* ptr = mset->get(label);
  mset->set(label, NULL);
  return array_transform(ptr, ni, nj);

}
bp::dict MatrixSet_to_dict(MatrixSet* mset) {

  typedef std::map<std::string, dcomplex*> MapPtr;
  bp::dict res;

  for(MapPtr::iterator it = mset->mat_map_.begin();
      it != mset->mat_map_.end(); ++it) {
    res[it->first] = array_transform(it->second,
				     mset->size_basis_i(),
				     mset->size_basis_j());
  }

  return res;

}

np::ndarray GTOs_at_r_ylm(GTOs* gtos, int L, int M, const np::ndarray& rs,
			  const np::ndarray& cs) {

  if(rs.get_nd() != 1) {
    std::string msg; SUB_LOCATION(msg);
    msg += "rs must be 1-dimensional";
    throw std::runtime_error(msg);
  }
  
  if(cs.get_nd() != 1) {
    std::string msg; SUB_LOCATION(msg);
    msg += "cs must be 1-dimensional";
    throw std::runtime_error(msg);
  }

  if(cs.shape(0) != gtos->size_basis()) {
    std::string msg; SUB_LOCATION(msg);
    msg += "Size mismatch between cs and GTOs";
    throw std::runtime_error(msg);
  }

  dcomplex* ptr_rs = reinterpret_cast<dcomplex*>(rs.get_data());
  dcomplex* ptr_cs = reinterpret_cast<dcomplex*>(cs.get_data());

  int num_r = rs.shape(0);
  dcomplex* ptr_vs = new dcomplex[num_r];

  gtos->AtR_Ylm(L, M, ptr_rs, num_r, ptr_cs, ptr_vs);

  std::cout << ptr_rs[0] << ptr_rs[1] << std::endl;
  std::cout << ptr_cs[0] << ptr_cs[1] << std::endl;
  std::cout << ptr_vs[0] << ptr_vs[1] << std::endl;

  np::ndarray vs = np::from_data(ptr_vs,
				 np::dtype::get_builtin<dcomplex>(),
				 bp::make_tuple(num_r),
				 bp::make_tuple(sizeof(dcomplex)),
				 bp::object());
  return vs;
}

// TO BE REMOVED
bp::tuple sym_eig(np::ndarray A, np::ndarray B) {

  if(A.get_nd() != 2 || B.get_nd() != 2) {
    std::string msg; SUB_LOCATION(msg);
    msg += "A and B must be 2-dimensional";
    throw std::runtime_error(msg);
  }

  int n = A.shape(0);
  if(n != A.shape(1)) {
    std::string msg; SUB_LOCATION(msg);
    msg += "A must be square matrix";
    throw std::runtime_error(msg);
  }

  if(n != B.shape(0) || n != B.shape(1)) {
    std::string msg; SUB_LOCATION(msg);
    msg += "A and B must be same size";
    throw std::runtime_error(msg);
  }

  dcomplex* ptr_a = reinterpret_cast<dcomplex*>(A.get_data());
  dcomplex* ptr_b = reinterpret_cast<dcomplex*>(B.get_data());

  Eigen::MatrixXcd AMat = Eigen::Map<Eigen::MatrixXcd>(ptr_a, n, n);
  Eigen::MatrixXcd BMat = Eigen::Map<Eigen::MatrixXcd>(ptr_b, n, n);

  Eigen::MatrixXcd CMat(n, n);
  Eigen::VectorXcd Eig(n, 1);
  generalizedComplexEigenSolve(AMat, BMat, &CMat, &Eig);

  std::cout << "molint_bind::eig" << std::endl;
  std::cout << Eig << std::endl;

  dcomplex* ptr_c;
  dcomplex *ptr_eig;
  Eigen::Map<Eigen::MatrixXcd>(ptr_c, n, n) = CMat;
  Eigen::Map<Eigen::VectorXcd>(ptr_eig, n) = Eig;

  np::ndarray c = np::from_data(ptr_c,
				np::dtype::get_builtin<dcomplex>(),
				bp::make_tuple(n, n),
				bp::make_tuple(sizeof(dcomplex), sizeof(dcomplex)),
				bp::object());
  np::ndarray eig = np::from_data(ptr_eig,
				  np::dtype::get_builtin<dcomplex>(),
				  bp::make_tuple(n),
				  bp::make_tuple(sizeof(dcomplex)),
				  bp::object());
  return bp::make_tuple(CMat, Eig);
}

BOOST_PYTHON_MODULE(cbasis_molint_bind) {

  Py_Initialize();
  np::initialize();

  def("sym_eig", sym_eig);

  class_<MatrixSet>("MatrixSet", init<int, int>())
    .def("to_dict", MatrixSet_to_dict)
    .def("get"     , MatrixSet_get);

  class_<GTOs>("GTOs", init<>())
    .def("size_basis", &GTOs::size_basis)
    .def("calc", &GTOs::Calc)
    .def("calc_zmat_other", &GTOs::CalcZMatOther)
    .def("add_atom_cpp", &GTOs::AddAtom)    
    .def("add_gtos_cpp", &GTOs::AddSphericalGTOs)
    .def("add_one_gto_cpp", &GTOs::AddOneSphericalGTO)
    .def("at_r_ylm_cpp", GTOs_at_r_ylm)
    .def("show", &GTOs::Show);
  
}

