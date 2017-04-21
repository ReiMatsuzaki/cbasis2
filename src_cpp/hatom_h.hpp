#ifndef HATOM_H_TEMPLATE_H
#define HATOM_H_TEMPLATE_H

#include "exp_func.hpp"
#include "op.hpp"

namespace cbasis { namespace h_atom_rad {

  // ==== radial hamiltonian ====
  template<class Field, int l>
  struct HOp {
    typedef OpD<Field, double, 2> D2;
    typedef OpRm<Field, double> Rm;
    typedef OpAdd<ScalarOpMult<D2>,
		  OpAdd<ScalarOpMult<Rm>,
			ScalarOpMult<Rm> > > type; 
    type operator()() {
      return op_add_op(scalar_mult_op(-0.5, D2()),
		       op_add_op(scalar_mult_op((l*l+l)/Field(2), Rm(-2)),
				 scalar_mult_op(-Field(1), Rm(-1))));
    }
  };
  template<class Field>
  struct HOp<Field, 0> {
    typedef OpD<Field, double, 2> D2;
    typedef OpRm<Field, double> Rm;
    typedef OpAdd<ScalarOpMult<D2>,
		  ScalarOpMult<Rm> > type; 
    type operator()() {
      return op_add_op(scalar_mult_op(-0.5, D2()),
		       scalar_mult_op(-Field(1), Rm(-1)));    
      
    }
  };

  // ==== radial H-E ====
  template<class Field, int l>
  struct HminusE {
    typedef OpAdd<typename HOp<Field, l>::type,
		  ScalarOpMult<One> > type;
    type operator() (Field ene){
      return op_add_op(HOp<Field, l>()(),
		       scalar_mult_op(ene, One()));
    }
  };

  // ==== eigen function ====
  template<class F, int n, int l> struct EigenFunc;
  template<class F> struct EigenFunc<F, 1, 0> {
    typedef ExpFunc<F, 1> type;
    type operator ()(){ return type(2, 1, 1); }
  };
  template<class F> struct EigenFunc<F, 2, 0> {
    typedef ExpFunc<F, 1> S;
    typedef FuncAdd<S, S> type;
    type operator ()(){ 
      F z = F(1)/F(2);
      return func_add_func(S(F(1)/sqrt(F(2)), 1, z),
			   S(-F(1)/(F(2)*sqrt(F(2))), 2, z));
    }
  };
  template<class F> struct EigenFunc<F, 2, 1> {
    typedef ExpFunc<F, 1> type;
    type operator ()(){ return type(F(1)/(F(2)*sqrt(F(6))), 2, F(1)/F(2)); }
  };  

  // ==== length ====
  template<class F, int n0, int l0, int l1> struct Length;
  template<class F> struct Length<F, 1, 0, 1> {
    typedef ExpFunc<F, 1> type;
    type operator()() {
      return type(2, 2, 1);
    }
  };

  // ==== velocity ====
  template<class F, int n0, int l0, int l1> struct Velocity;
  template<class F> struct Velocity<F, 1, 0, 1> {
    typedef ExpFunc<F, 1> type;
    type operator()() {
      return type(-2, 1, 1);
    }
  };

  }}
#endif
