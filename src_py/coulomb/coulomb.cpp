#include <gsl/gsl_sf_coulomb.h>
#include <boost/python.hpp>

using namespace boost::python;

class CoulombFunc {
public:
  double F;
  double G;
  double Fp;
  double Gp;
  double err;
  CoulombFunc(double eta, double x, double L_F, int k) {
  
    gsl_sf_result f, fp, g, gp;
    double exp_f, exp_g;
    gsl_sf_coulomb_wave_FG_e(eta, x, L_F, k, &f, &fp, &g, &gp, &exp_f, &exp_g);
  
    F = f.val;
    Fp = fp.val;
    G = g.val;
    Gp = gp.val;

    err = f.err;
    if(f.err < g.err)
      err = g.err;
    if(g.err < fp.err)
      err = fp.err;
    if(fp.err < gp.err)
      err = gp.err;
  }
  double get_F() { return F; }
  double get_G() { return G; }
  double get_Fp() { return Fp; }
  double get_Gp() { return Gp; }
  double get_err() { return err; }
};

BOOST_PYTHON_MODULE(coulomb_bind) {

  class_<CoulombFunc>("coulomb", init<double,double,double,int>())
    .add_property("f", &CoulombFunc::get_F)
    .add_property("g", &CoulombFunc::get_G)
    .add_property("fp", &CoulombFunc::get_Fp)
    .add_property("gp", &CoulombFunc::get_Gp)
    .add_property("err", &CoulombFunc::get_err);
}
