#ifndef OPT_ALPHA_H
#define OPT_ALPHA_H

#include <Eigen/Core>
#include "typedef.hpp"

namespace cbasis {

  /*
  void AlphaGrad(const Eigen::VectorXcd& d, const Eigen::VectorXcd& c,
   		 const Eigen::MatrixXcd& D10,
		 const Eigen::VectorXcd& r1,
		 const Eigen::VectorXcd& s0, const Eigen::VectorXcd& s1,
		 dcomplex *a, Eigen::VectorXcd *g);
  
  void AlphaGradHess(const Eigen::VectorXcd& c,
		     const Eigen::VectorXcd& d,
		     const Eigen::MatrixXcd& U,
		     const Eigen::MatrixXcd& D00,
		     const Eigen::MatrixXcd& D10,
		     const Eigen::MatrixXcd& D20,
		     const Eigen::MatrixXcd& D11,
		     const Eigen::VectorXcd& r0,
		     const Eigen::VectorXcd& r1,
		     const Eigen::VectorXcd& r2,
		     const Eigen::VectorXcd& s0,
		     const Eigen::VectorXcd& s1,
		     const Eigen::VectorXcd& s2,
		     dcomplex* a, Eigen::VectorXcd* g, Eigen::MatrixXcd* h);
  */

  /*
    Optimize orbital exponents of cGTO for alpha.    
   */
  /*
  class R1STOs;
  class R1GTOs;

  void OptAlphaShift(const R1STOs& driv,
		     const Eigen::VectorXi& opt_idx,
		     int L,
		     dcomplex ene,
		     double h,
		     int max_iter,
		     double eps,
		     double eps_canonical,
		     R1GTOs& gtos,
		     dcomplex& z_shift,
		     bool* convq,
		     dcomplex* alpha);
		     */

  // ==== calculation for driv eq ====
  class R1STOs;
  class R1GTOs;
  class IDriv;
  class IOp;
  dcomplex CalcAlpha(IDriv* driv, IOp* op, const R1GTOs& gs);
  void SolveAlpha(IDriv* driv, IOp* op, const R1GTOs& gs, Eigen::VectorXcd&);

  // ==== optimize target ====
  class IOptTarget {
  public:
    //    IOptTarget();
    virtual ~IOptTarget();
    virtual int dim() = 0;
    virtual dcomplex Val(const Eigen::VectorXcd& zs);
    virtual void ValGradHess(const Eigen::VectorXcd&, dcomplex*,
			     Eigen::VectorXcd&, Eigen::MatrixXcd&) = 0;
  };
  bool CheckOptTarget(IOptTarget* target, const Eigen::VectorXcd& zs,
		      double h, double eps);
  bool CheckOptTarget(IOptTarget* target, dcomplex z, double h, double eps);
		      

  // ==== Optimization result ====
  class OptResult {
  public:
    bool conv_q;
    Eigen::VectorXcd zs;
    dcomplex val;
    Eigen::VectorXcd dz;
    Eigen::VectorXcd grad;
    Eigen::MatrixXcd hess;
    OptResult(int n);
  };
  std::ostream& operator<<(std::ostream& out, const OptResult&);

  // ==== Optimizer ====
  // 
  // max_iter : maximum iteration for optimization. (max_iter > 0)
  // eps      :convergence epsilon. (eps > 0)
  // target   : optimization target. (not null)
  // debug_lvl: debug level.
  class IOptimizer {
  public:
    // -- input --
    int max_iter;   
    double eps;
    IOptTarget* target;
    int debug_lvl;
    // -- results --
    /*
    bool conv_q;
    Eigen::VectorXcd zs;
    dcomplex val;
    Eigen::VectorXcd dz;
    Eigen::VectorXcd grad;
    Eigen::MatrixXcd hess;
    */
  public:
    IOptimizer(int _max_iter, double _eps, IOptTarget* _target, int lvl);
    virtual ~IOptimizer();
    OptResult* Optimize(const Eigen::VectorXcd& z0s);
    OptResult* Optimize(dcomplex);
    virtual void OneStep(const Eigen::VectorXcd& _grad, const Eigen::MatrixXcd& _hess,
			 Eigen::VectorXcd& _dz) = 0;
  };
  class OptNewton :public IOptimizer {
  public:
    OptNewton(int _max_iter, double _eps, IOptTarget* _target, int lvl=0);
    ~OptNewton();
    void OneStep(const Eigen::VectorXcd& _grad, const Eigen::MatrixXcd& _hess,
		 Eigen::VectorXcd& _dz);
  };

  // ==== alpha ====
  class OptAlpha :public IOptTarget {
  public:
    IDriv* driv;
    IOp*   op;
    R1GTOs& ref_g0s;
    R1GTOs* g1s;
    R1GTOs* g2s;
    OptAlpha(IDriv* _driv, IOp* _op, R1GTOs& gtos);
    ~OptAlpha();
    dcomplex Val(const Eigen::VectorXcd& zs);
    int dim();
    void Solve(Eigen::VectorXcd&);
    void ValGradHess(const Eigen::VectorXcd& zs, dcomplex *,
		     Eigen::VectorXcd&, Eigen::MatrixXcd&);
    void Update(const Eigen::VectorXcd&);
  };
  class OptAlphaShift : public IOptTarget {
  public:
    OptAlpha* opt_alpha;
    Eigen::VectorXi index;
    Eigen::VectorXcd z0s;
    OptAlphaShift(IDriv* _driv, IOp* _op, R1GTOs& gtos,
		  const Eigen::VectorXi& idx);
    ~OptAlphaShift();
    int dim();
    //dcomplex Val(const Eigen::VectorXcd& zs);
    void ValGradHess(const Eigen::VectorXcd& zs, dcomplex *,
		     Eigen::VectorXcd&, Eigen::MatrixXcd&);    
  };
  class OptAlphaPartial : public IOptTarget {
  public:
    OptAlpha* opt_alpha;
    Eigen::VectorXi index;
    Eigen::VectorXcd z0s;
    OptAlphaPartial(IDriv* _driv, IOp* _op, R1GTOs& gtos,
		    const Eigen::VectorXi& idx);
    ~OptAlphaPartial();
    void Update(const Eigen::VectorXcd& zs);
    int dim();
    //dcomplex Val(const Eigen::VectorXcd&);
    void ValGradHess(const Eigen::VectorXcd&, dcomplex *,
		     Eigen::VectorXcd&, Eigen::MatrixXcd&);    
  };

  /*
    Compute coefficient for driven type equation.
    Notice that calculation speed is not fast.
   */
  /*
  Eigen::VectorXcd SolveAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene);
  dcomplex         CalcAlpha( const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene);
  Eigen::VectorXcd SolveAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene, const R1STOs& sto_pot);
  dcomplex         CalcAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene, const R1STOs& sto_pot);
  */
  //void SolveAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene,
  //		  Eigen::VectorXcd& c);
  //  void SolveAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene,
  //		  const R1STOs& sto_pot, Eigen::VectorXcd& c, dcomplex* a);

  /*
  dcomplex CalcAlphaFull(const R1STOs& driv,R1GTOs& gtos, int L, dcomplex ene);
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos, int L,  dcomplex ene,
		     double eps);
  dcomplex CalcAlpha(const R1STOs& driv, R1GTOs& gtos, int L, dcomplex ene, const R1STOs&);
  */
		     
}
#endif
