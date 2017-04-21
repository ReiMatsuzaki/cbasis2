#include "r1gtoint.hpp"
#include "op_driv.hpp"

#include <Eigen/Core>

using namespace Eigen;

namespace cbasis {
  
  // ==== Driven term ====
  IDriv::~IDriv() {}
  
  DrivSTO::DrivSTO(const R1STOs& _stos): IDriv(), stos(new R1STOs(_stos)) {}
  DrivSTO::~DrivSTO() {
    delete stos;
  }
  void DrivSTO::Calc(const R1GTOs& gs, VectorXcd& m) {
    gs.CalcVec(*stos, m);
  }

  
  // ==== Operator ====
  IOp::~IOp() {}
  void IOp::Calc(const R1GTOs& gi, MatrixXcd& m) {
    this->Calc(gi, gi, m);
  }
  OpCoulomb::OpCoulomb(int _L, dcomplex _ene): L(_L), energy(_ene) {}
  OpCoulomb::~OpCoulomb() {}
  void OpCoulomb::Calc(const R1GTOs& gi, const R1GTOs& gj,MatrixXcd& m) {
    static MatrixXcd s, t, v;
    gi.CalcMatSTV(gj, L, s, t, v);
    m = t + v - energy * s;
  }

  OpCoulombShort::OpCoulombShort(int _L, dcomplex _ene, const R1STOs& _short_pot) {
    L = _L;
    energy = _ene;
    short_pot = new R1STOs(_short_pot);
  }
  OpCoulombShort::~OpCoulombShort() {
    delete short_pot;
  }
  void OpCoulombShort::Calc(const R1GTOs& gi, const R1GTOs& gj, MatrixXcd& m) {
    static MatrixXcd s, t, v, v2;
    gi.CalcMatSTV(gj, L, s, t, v);
    gi.CalcMatSTO(gj, *short_pot, v2);
    m = t + v + v2 - energy * s;
  }  
  

}
