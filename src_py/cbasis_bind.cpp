#include <boost/python.hpp>
#include <boost/numpy.hpp>
#include "../src_cpp/typedef.hpp"
#include "../src_cpp/linspace.hpp"
#include "../src_cpp/exp_func.hpp"
#include "../src_cpp/cut_exp.hpp"
#include "../src_cpp/delta.hpp"
#include "../src_cpp/cip_exp.hpp"
#include "../src_cpp/op.hpp"

namespace {
  using namespace boost::python;
  using namespace cbasis;
  typedef std::complex<double> F;

  namespace bp = boost::python;
  namespace np = boost::numpy;
  using namespace std;

}

using namespace std;
using namespace boost;
using namespace boost::python;

np::ndarray array_transform(dcomplex* xs, int n) {  
  return np::from_data(xs,
		       np::dtype::get_builtin<dcomplex>(),
		       bp::make_tuple(n),
		       bp::make_tuple(sizeof(dcomplex)),
		       bp::object());
}

/*
np::ndarray CalcSMat(const SphericalGTOSet& a, const SphericalGTOSet& b) {
  dcomplex* vs = a.SMat(b);
  np::ndarray res = array_transform(vs, a.size() * b.size());
  return res;
}
np::ndarray CalcTMat(const SphericalGTOSet& a, const SphericalGTOSet& b) {
  dcomplex* vs = a.TMat(b);
  return array_transform(vs, a.size() * b.size());
}
np::ndarray CalcVMat(const SphericalGTOSet& a,
		     dcomplex q, dcomplex x, dcomplex y, dcomplex z,
		     const SphericalGTOSet& b) {
  dcomplex* vs = a.VMat(q, x, y, z, b);
  return array_transform(vs, a.size() * b.size());
}
np::ndarray CalcXyzMat(const SphericalGTOSet& a,
		       int nx, int ny, int nz,
		       const SphericalGTOSet& b) {
  dcomplex* vs = a.XyzMat(nx, ny, nz, b);
  return array_transform(vs, a.size() * b.size());
}
*/
BOOST_PYTHON_MODULE(cbasis_bind) {

  Py_Initialize();
  np::initialize();

  /*
  class_<array3<F> >("tuple_c3", init<F, F, F>())
    .add_property("x", &array3<F>::x, &array3<F>::set_x)
    .add_property("y", &array3<F>::y, &array3<F>::set_y)
    .add_property("z", &array3<F>::z, &array3<F>::set_z);

  class_<array3<int> >("tuple_i3", init<int, int, int>())
    .add_property("x", &array3<int>::x, &array3<int>::set_x)
    .add_property("y", &array3<int>::y, &array3<int>::set_y)
    .add_property("z", &array3<int>::z, &array3<int>::set_z);
    */
  
  typedef ExpFunc<F,1> STO;
  class_<STO>("STO", init<F, int, F>())
    .add_property("c", &STO::c, &STO::set_c)
    .add_property("n", &STO::n, &STO::set_n)
    .add_property("z", &STO::z, &STO::set_z)
    .def("__str__", &STO::str)
    .def("at", &STO::at);

  typedef ExpFunc<F,2> GTO;
  class_<GTO>("GTO", init<F, int, F>())
    .add_property("c", &GTO::c, &GTO::set_c)
    .add_property("n", &GTO::n, &GTO::set_n)
    .add_property("z", &GTO::z, &GTO::set_z)
    .def("__str__", &GTO::str)
    .def("at", &GTO::at);

  typedef CutExpFunc<F, 1> CutSTO;
  class_<CutSTO>("CutSTO", init<F, int, F, double>())
    .add_property("c", &CutSTO::c, &CutSTO::set_c)
    .add_property("n", &CutSTO::n, &CutSTO::set_n)
    .add_property("z", &CutSTO::z, &CutSTO::set_z)
    .add_property("r0", &CutSTO::r0, &CutSTO::set_r0)
    .def("__str__", &CutSTO::str)
    .def("at", &CutSTO::at);

  /*
  typedef CutExpFunc<F, 2> CutGTO;
  class_<CutGTO>("CutGTO", init<F, int, F, double>())
    .add_property("c", &CutGTO::c, &CutGTO::set_c)
    .add_property("n", &CutGTO::n, &CutGTO::set_n)
    .add_property("z", &CutGTO::z, &CutGTO::set_z)
    .add_property("r0", &CutGTO::r0, &CutGTO::set_r0)
    .def("__str__", &CutGTO::str)
    .def("at", &CutGTO::at);
  */

  def("cip_ss", CIP<STO, STO>);
  def("cip_sg", CIP<STO, GTO>);
  def("cip_gs", CIP<GTO, STO>);
  def("cip_gg", CIP<GTO, GTO>);
  def("cip_cut_ss", CIP<CutSTO, CutSTO>);
  //  def("cip_cut_gg", CIP<CutGTO, CutGTO>);
  

  class_<CRm>("Rm", init<int>())
    .add_property("m", &CRm::m);
  class_<CD1>("D1", init<>());
  class_<CD2>("D2", init<>());
  class_<Cut<F, double> >("Cut", init<double>())
    .add_property("r0", &Cut<F, double>::r0);

  def("cip_s_d2_s", CIP<STO, CD2, STO>);
  def("cip_s_d2_g", CIP<STO, CD2, GTO>);
  def("cip_g_d2_s", CIP<GTO, CD2, STO>);
  def("cip_g_d2_g", CIP<GTO, CD2, GTO>);
  def("cip_cut_s_d2_s", CIP<CutSTO, CD2, CutSTO>);

  def("cip_s_d1_s", CIP<STO, CD1, STO>);
  def("cip_s_d1_g", CIP<STO, CD1, GTO>);
  def("cip_g_d1_s", CIP<GTO, CD1, STO>);
  def("cip_g_d1_g", CIP<GTO, CD1, GTO>);
  def("cip_cut_s_d1_s", CIP<CutSTO, CD1, CutSTO>);

  def("cip_s_rm_s", CIP<STO, CRm, STO>);
  def("cip_s_rm_g", CIP<STO, CRm, GTO>);
  def("cip_g_rm_s", CIP<GTO, CRm, STO>);
  def("cip_g_rm_g", CIP<GTO, CRm, GTO>);
  def("cip_cut_s_rm_s", CIP<CutSTO, CRm, CutSTO>);


  // ==== 3D ====
  /*
  class_<SphericalGTOSet>("SphericalGTOSet", init<>())
    .add_property("size", &SphericalGTOSet::size)
    .def("add_one_basis_cpp", &SphericalGTOSet::AddOneBasis)
    .def("add_basis_cpp",     &SphericalGTOSet::AddBasis);

  def("calc_s_mat", CalcSMat);
  def("calc_t_mat", CalcTMat);
  def("calc_v_mat", CalcVMat);
  def("calc_xyz_mat", CalcXyzMat);
  */
}

