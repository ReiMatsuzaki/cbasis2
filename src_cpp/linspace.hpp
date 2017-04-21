#ifndef LINSPACE_TEMPLATE_H
#define LINSPACE_TEMPLATE_H

#include <complex>
#include <list>
#include <utility>
#include <boost/array.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/or.hpp>
#include <boost/static_assert.hpp>

namespace cbasis {

  // ==== L2 function space ====  
  // ---- base class ----
  template<class _Field, class _Coord>
  struct L2funcComponent {
    typedef _Field Field;
    typedef _Coord Coord;
  };

  // ---- (only declation) func and op ----
  template<class _Field, class _Coord>
  struct Func;// : public L2funcComponent<_Field, _Coord>;

  template<class _Field, class _Coord>
  struct Op;// : public L2funcComponent<_Field, _Coord>;

  // ---- meta function ----
  template<class A>
  struct is_lin_separable: boost::mpl::false_ {};

  template<class A, class B>
  struct in_same_space : boost::mpl::and_<
    boost::is_same<typename A::Field, typename B::Field>,
    boost::is_same<typename A::Coord, typename B::Coord> > {};

  template<class A> struct is_func: boost::mpl::false_ {};

  // ==== Linear space operation ====
  // ---- zero ----
  template<class _Field, class _Coord> 
  class ZeroFunc : public Func<_Field, _Coord> {
    _Field at(_Coord ) { return _Field(0); }
  };

  // ---- one ----
  struct One {};
  template<class Field, class Coord>
  class OneOp : public Op<Field, Coord>, One {};

  // ---- func_add_func ----
  template<class A, class B>
  struct FuncAdd : public Func<typename A::Field, typename A::Coord> {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    typedef Func<typename A::Field, typename A::Coord> base;
    A left;
    B right;
    FuncAdd(A a, B b): left(a), right(b) {}
    typename base::Field at(typename base::Coord x) {
      return left.at(x) + right.at(x);
    }
  };

  template<class A, class B>
  struct is_lin_separable<FuncAdd<A, B> > : boost::mpl::true_ {};

  template<class A, class B>
  FuncAdd<A, B> func_add_func(const A& a, const B& b) {
    return FuncAdd<A, B>(a, b);
  }
  template<class B>
  B func_add_func(const ZeroFunc<typename B::Field, typename B::Coord>&, const B& b) {
    return b; }
  template<class A>
  A func_add_func(const A& a, const ZeroFunc<typename A::Field, typename A::Coord>&) {
    return a; }


  // ---- scalar_mult_func ----
  template<class A>
  struct ScalarFuncMult : public Func<typename A::Field, typename A::Coord> {
    typedef Func<typename A::Field, typename A::Coord> base;
    typename base::Field scalar;
    A func;
    ScalarFuncMult(typename A::Field c, A a): scalar(c), func(a) {}
    typename base::Field at(typename base::Coord x) {
      return scalar * func.at(x);
    }
  };

  template<class A>
  struct is_lin_separable<ScalarFuncMult<A> > : boost::mpl::true_ {};

  template<class A>
  ScalarFuncMult<A> scalar_mult_func(typename A::Field c, const A& a) {
    return ScalarFuncMult<A>(c, a);
  }


  // ---- linear combination of func ----
  template<class A>
  struct LinFunc: public Func<typename A::Field, typename A::Coord> {
    typedef Func<typename A::Field, typename A::Coord> base;
    typedef typename A::Field Field;
    typedef typename A::Coord Coord;
    typedef std::list<std::pair<typename A::Field, A> > CFList;
    typedef typename CFList::const_iterator const_iterator;
    // typedef typename CFList::iterator iterator;
    CFList c_f_list_;    
    void Add(Field c, const A& a)  {
      c_f_list_.push_back(std::make_pair(c, a));
    }
    Field at(Coord x) const { 

      typename A::Field cumsum(0);
      for(const_iterator it = this->begin(); it != this->end(); ++it) {
	cumsum += it->first * it->second.at(x);
      }
      return cumsum;
    }
    void SetScalarProd(Field c) {
      for(typename CFList::iterator it = c_f_list_.begin();
	  it != c_f_list_.end(); ++it) {
	it->first *= c;
      }
    }
    void SetNormalize() {
      Field norm2 = CIP(*this, *this);
      this->SetScalarProd(Field(1)/sqrt(norm2));
    }
    typename CFList::const_iterator begin() const { return c_f_list_.begin(); }
    typename CFList::const_iterator end() const { return c_f_list_.end(); }
  };

  template<class A>
  struct is_lin_separable<LinFunc<A> > : boost::mpl::true_ {};

  // ---- op_add_op ----
  template<class A, class B>
  struct OpAdd :public Op<typename A::Field, typename A::Coord> {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    A left;
    B right;
    OpAdd(A a, B b): left(a), right(b) {}
  };
  
  template<class A, class B>
  struct is_lin_separable<OpAdd<A, B> > : boost::mpl::true_ {};

  template<class A, class B>
  OpAdd<A, B> op_add_op(const A& a, const B& b) { 
    return OpAdd<A, B>(a, b);
  }

  template<class A, class B>
  std::ostream& operator << (std::ostream& os, const OpAdd<A, B>& a) {
    os << a.left << " + " << a.right;
    return os;
  }

  // ---- scalar_mult_op ----
  template<class A>
  struct ScalarOpMult :public Op<typename A::Field, typename A::Coord> {
    typename A::Field scalar;
    A op;
    ScalarOpMult(typename A::Field c, A a): scalar(c), op(a) {}    
  };

  template<class A>
  struct is_lin_separable<ScalarOpMult<A> > : boost::mpl::true_{};

  template<class A>
  ScalarOpMult<A> scalar_mult_op(typename A::Field c, const A& a) {
    return ScalarOpMult<A>(c, a);
  }

  // ==== impl of func/op ====
  template<class _Field, class _Coord>
  struct Func : public L2funcComponent<_Field, _Coord> {
    typedef Func<_Field, _Coord> type;

    template<class Other>
    FuncAdd<type, Other> operator+(const Other& a) {
      return func_add_func(*this, a);
    }
  };

  template<class _Field, class _Coord>
  struct Op : public L2funcComponent<_Field, _Coord> {

  };


  // ==== operator apply ====
  /*
  // ---- primitive ----
  template<class OpT, class A>
  struct op_func_res;

  template<class OpT, class A>
  typename op_func_res<OpT, A>::type OP_prim(const OpT& o, const A& a);

  // ---- linearly separable ----
  template<class OpT, class A>
  typename op_func_res<OpT, A>::type
  OP_lin(const OpT& o, const A& a,
	      typename boost::disable_if<is_lin_separable<OpT> >::type* =0,
	      typename boost::disable_if<is_lin_separable<A> >::type* =0
	      ) {
    return OP_prim(o, a);
  } 

  template<class OpT, class A, class B>
  FuncAdd<typename op_func_res<OpT, A>::type, 
	  typename op_func_res<OpT, B>::type > 
  OP_lin(const OpT& o, const FuncAdd<A,B>& ab,
	      typename boost::disable_if<is_lin_separable<OpT> >::type* =0
	      ) {
    return func_add_func(OP_lin(o, ab.left), OP_lin(o, ab.left));
  } 

  template<class OpT, class A>
  ScalarFuncMult<typename op_func_res<OpT, A>::type > 
  OP_lin(const OpT& o, const ScalarFuncMult<A>& ca,
	      typename boost::disable_if<is_lin_separable<OpT> >::type* =0
	      ) {
    return scalar_mult_func(ca.scalar, OP_lin(o, ca.func));
  } 


  //  template<class OpT, class A>
  //  LinFunc<typename op_func_res<OpT, A>::type >
  //  OP_lin(const OpT& o, const LinFunc<A>& as) {  }


  template<class OpA, class OpB, class C>
  FuncAdd<typename op_func_res<OpA, C>::type, 
	  typename op_func_res<OpB, C>::type > 
  OP_lin(const OpAdd<OpA, OpB>& o, const C& c) {
    return func_add_func(OP_lin(o.left, c), OP_lin(o.right, c));
  } 

  template<class OpT, class A>
  ScalarFuncMult<typename op_func_res<OpT, A>::type >
  OP_lin(const ScalarOpMult<OpT>& o, const A& a) {
    return scalar_mult_func(o.scalar, OP_lin(o.op, a));
  } 

  // ---- interface ----
  //  template<class OpT, class A>
  //  OP(const OpT& o, const A& a) {
  //    return OP_lin(o, a);
  //  }
	      
  */
  // ==== inner product ====
  // ---- primitive ----
  template<class A, class OpT, class B>
  typename A::Field CIP_impl_prim(const A&, const OpT&, const B&);// {
    //    return typename A::Field(777);
    //  }

  // ---- linearly separate ----
  template<class A, class OpT, class B>
  typename A::Field 
  CIP_impl_lin(const A& a, const OpT& o, const B& b, 
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0,
	       typename boost::disable_if<is_lin_separable<B> >::type* =0) {
    return CIP_impl_prim(a, o, b);
  }

  template<class A,class OpT, class B, class C>
  typename A::Field 
  CIP_impl_lin(const A& a, const OpT& o, const FuncAdd<B, C>& bc,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0) {
    return CIP_impl_lin(a, o, bc.left) + CIP_impl_lin(a, o, bc.right);
  }

  template<class A, class OpT, class B>
  typename A::Field 
  CIP_impl_lin(const A& a, const OpT& o, const ScalarFuncMult<B>& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0) {
    return b.scalar * CIP_impl_lin(a, o, b.func);
  }

  template<class A, class OpT, class B>
  typename A::Field 
  CIP_impl_lin(const A& a, const OpT& o, const LinFunc<B>& bs,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0,
	       typename boost::disable_if<is_lin_separable<OpT> >::type* =0) {

    typename A::Field cumsum(0);
    for(typename LinFunc<B>::const_iterator it = bs.begin(); it != bs.end(); ++it) {

      cumsum += it->first * CIP_impl_lin(a, o, it->second);
      
    }
    return cumsum;
  }


  template<class A, class OpT, class OpS, class B>
  typename A::Field
  CIP_impl_lin(const A& a, const OpAdd<OpT, OpS>& o, const B& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0) {
    return CIP_impl_lin(a, o.left, b) + CIP_impl_lin(a, o.right, b);
  }
  template<class A, class OpT, class B>
  typename A::Field
  CIP_impl_lin(const A& a, const ScalarOpMult<OpT>& o, const B& b,
	       typename boost::disable_if<is_lin_separable<A> >::type* =0) {
    return o.scalar * CIP_impl_lin(a, o.op, b);
  }

  template<class A, class B, class OpT, class C>
  typename A::Field 
  CIP_impl_lin(const FuncAdd<A, B>& ab, const OpT& o, const C& c) {
    return CIP_impl_lin(ab.left, o, c) + CIP_impl_lin(ab.right, o, c);
  }

  template<class A, class OpT, class B>
  typename A::Field 
  CIP_impl_lin(const ScalarFuncMult<A>& a, const OpT& o, const B& b) {
    return a.scalar * CIP_impl_lin(a.func, o, b);
  }

  template<class A, class OpT, class B>
  typename A::Field
  CIP_impl_lin(const LinFunc<A>& as, const OpT& o, const B& b) {

    typename A::Field cumsum(0);
    for(typename LinFunc<A>::const_iterator it = as.begin(); it != as.end(); ++it) {

      cumsum += it->first * CIP_impl_lin(it->second, o, b);
      
    }
    return cumsum;    
  }

  // ---- general interface ----
  template<class A, class B>
  typename A::Field CIP(const A& a, const B& b) {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    return CIP_impl_lin(a, One(), b);
  }

  template<class A, class O, class B>
  typename A::Field CIP(const A& a, const O& o, const B& b) {
    BOOST_STATIC_ASSERT((in_same_space<A, B>::value));
    BOOST_STATIC_ASSERT((in_same_space<A, O>::value));
    return CIP_impl_lin(a, o, b);
  }


  // ==== Utilities ====
  // ---- n term ----
  template<int N, class FuncT>
  struct NTermFunc {
    typedef FuncAdd<FuncT, typename NTermFunc<N-1, FuncT>::type> type;
  };
  template<class FuncT>
  struct NTermFunc<1, FuncT> {
    typedef FuncT type;
  };

  // ---- normalization ----
  template<class FuncT>
  typename FuncT::Field CNorm2(const FuncT& a) { return CIP(a, a); }
  template<class FuncT>
  typename FuncT::Field CNorm(const FuncT& a) { return sqrt(CNorm2(a)); }
  template<class FuncT>
  ScalarFuncMult<FuncT> CNormalize(const FuncT& a) {
    typename FuncT::Field one(1);
    return scalar_mult_func(one/CNorm(a), a);
  }
  
  // ==== type ====
  typedef double R1d;
  typedef boost::array<double, 3> R3d;

  // ==== For practice ====
  /*
  class RealBasisOn1D : public Func<double, double> { 
  public:
    double param; 
    RealBasisOn1D(double p) : param(p) {}
  };
  class RealBasisOn3D : public Func<double, boost::array<double, 3> > {};
  class ComplexBasisOn1D : public Func< std::complex<double>, double> {};
  template<>double CIP_impl_prim<RealBasisOn1D, One, RealBasisOn1D>
  (const RealBasisOn1D& a, const One&, const RealBasisOn1D& b) {
    return a.param * b.param;
  }
  */
}

#endif
