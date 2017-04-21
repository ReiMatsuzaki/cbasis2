#ifndef OP_FUNC_H
#define OP_FUNC_H

#include <Eigen/Core>
#include "typedef.hpp"



namespace cbasis {

  class R1GTOs;
  class R1STOs;

  // ==== Driven term ====
  class IDriv {
  public:
    virtual ~IDriv();
    virtual void Calc(const R1GTOs&, Eigen::VectorXcd&)=0;
  };
  class DrivSTO :public IDriv {
  public:
    R1STOs* stos;
    DrivSTO(const R1STOs& stos);
    ~DrivSTO();
    void Calc(const R1GTOs&, Eigen::VectorXcd&);
  };

  // ==== Operator ====
  class IOp {
  public:
    virtual void Calc(const R1GTOs&, Eigen::MatrixXcd&);
    virtual void Calc(const R1GTOs&, const R1GTOs&, Eigen::MatrixXcd&)=0;
    virtual ~IOp();
  };
  class OpCoulomb : public IOp{
  public:
    int L;
    dcomplex energy;
    OpCoulomb(int, dcomplex);
    ~OpCoulomb();
    void Calc(const R1GTOs&, const R1GTOs&, Eigen::MatrixXcd&);
  };
  class OpCoulombShort: public IOp {
  public:
    int L;
    dcomplex energy;
    R1STOs* short_pot;
    OpCoulombShort(int _L, dcomplex _ene, const R1STOs& _short_pot);
    ~OpCoulombShort();
    void Calc(const R1GTOs&, const R1GTOs&, Eigen::MatrixXcd&);
  };


}

#endif
