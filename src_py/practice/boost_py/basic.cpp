#include <iostream>
#include <boost/python.hpp>

using namespace std;

int add(int a, int b) {
  return a+b;
}

BOOST_PYTHON_MODULE(basic) {
  using namespace boost::python;
  //  Py_Initialize();
  def("add", &add);
}

