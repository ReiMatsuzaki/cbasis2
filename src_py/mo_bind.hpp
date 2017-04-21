void AddB2EInt() {

  register_ptr_to_python<B2EInt>();
  class_<IB2EInt, boost::noncopyable>("IB2EInt", no_init);
  class_<B2EIntMem, bases<IB2EInt> >("B2EIntMem", init<>())
    .def("get", &IB2EInt::Get)
    .def("set", &IB2EInt::Set)
    .def("reset", &IB2EInt::Reset)
    .def("at", &IB2EInt::At)
    .def("exist", &IB2EInt::Exist)
    .def("write", &IB2EInt::Write)
    .def("size", &IB2EInt::size)
    .def("capacity", &IB2EInt::capacity);
  def("ERI_read", ERIRead);
  
}


boost::python::list MO_get_num_occ(MO mo) {
  boost::python::list num_occ = PyListFromCppVector(mo->num_occ_irrep);
  return num_occ;
}
MO CalcRHF1(SymGTOs gtos, int nele, int max_iter, double eps, int debug_lvl) {
  bool is_conv;
  MO mo = CalcRHF(gtos, nele, max_iter, eps, &is_conv, debug_lvl);
  return mo;
}
MO CalcRHF2(pSymmetryGroup sym, BMatSet mat_set, B2EInt eri, int nele, 
	    int max_iter, double eps, int debug_lvl) {
  bool is_conv;
  MO mo = CalcRHF(sym, mat_set, eri, nele, max_iter, eps, &is_conv, debug_lvl);
  return mo;
}
BMat* py_CalcSEHamiltonian(MO mo, B2EInt eri, Irrep I0, int i0) {

  BMat *res = new BMat();
  CalcSEHamiltonian(mo, eri, I0, i0, res);
  return res;

}
BMat* py_JK(B2EInt eri, BMat& C, Irrep I0, int i0, dcomplex c_J, dcomplex c_K) {

  BMat *mat = new BMat();
  BMatSetZero(C, *mat);
  
  AddJK(eri, C, I0, i0, c_J, c_K, *mat);

  return mat;
}
BMat* py_J(B2EInt eri, VectorXcd& ca, Irrep irrep_a, 
	   SymGTOs g0, SymGTOs g1, dcomplex coef_J) {

  BMat *mat = new BMat();
  BMatSetZero(g0, g1, *mat);
  
  AddJ(eri, ca, irrep_a, coef_J, *mat);

  return mat;

}
BMat* py_K(B2EInt eri, VectorXcd& ca, Irrep ir_a,
	   SymGTOs g0, SymGTOs g1, dcomplex coef_K) {

  BMat *mat = new BMat();
  BMatSetZero(g0, g1, *mat);
  
  AddK(eri, ca, ir_a, coef_K, *mat);

  return mat;

}

void AddMO() {

  register_ptr_to_python<MO>();
  class_<_MO>("MO", no_init)
    .def_readonly("H", &_MO::H)
    .def_readonly("S", &_MO::S)
    .def_readonly("C", &_MO::C)
    .def_readonly("eigs", &_MO::eigs)
    .def("num_occ", MO_get_num_occ)
    .def_readonly("energy", &_MO::energy);


  def("calc_RHF", CalcRHF1);
  def("calc_RHF", CalcRHF2);

  def("calc_SEHamiltonian", py_CalcSEHamiltonian,
      return_value_policy<manage_new_object>());

  def("calc_JK", py_JK, return_value_policy<manage_new_object>());

  def("calc_J", py_J, return_value_policy<manage_new_object>());
  def("calc_K", py_K, return_value_policy<manage_new_object>());

  def("calc_alpha", CalcAlpha);
  def("pi_total_crosssection", PITotalCrossSection);

}
