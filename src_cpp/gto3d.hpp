#ifndef GTO3D_HPP
#define GTO3D_HPP

#include <iostream>

#include "math_utils.hpp"
#include "macros.hpp"
#include "linspace.hpp"
#include "op.hpp"

#include "angmoment.hpp"
#include "cints.hpp"
#include "molint.hpp"

namespace cbasis {

  typedef cbasis::array3<double> d3;
  typedef cbasis::array3<int>    i3;
  typedef cbasis::array3<dcomplex> c3;

  // ==== Operators ====
  template<class Field, class Coord>
  class OpKE : public Op<Field, Coord> {};
  template<class Field, class Coord>
  class OpNA :public Op<Field, Coord>  {
  private:
    Field q_;
    Coord xyz_;
  public:
    OpNA(Field q, const Coord& xyz): q_(q), xyz_(xyz) {}
    Field q() const { return q_;}
    const Coord& xyz() const { return xyz_;}
  };
  template<class Field, class Coord>
  class OpXyz :public Op<Field, Coord>  {
  private:
    i3 nml_;
  public:
    OpXyz(int n, int m, int l) : nml_(n, m, l) {}
    const i3& nml() const { return nml_; }
    int n() const { return nml_[0]; }
    int m() const { return nml_[1]; }
    int l() const { return nml_[2]; }
    void show() const {
      std::cout << nml_[0] << nml_[1] << nml_[2] << std::endl;
    }
  };


  // ==== Cartesian GTO ====
  // F : Field
  // FC : Field of coordinate
  template<class F, class FC>
  class CartGTO :public Func<F, FC> {
  public:
    // ---- type ----
    typedef F Field;
    typedef FC FieldCoord;
    typedef typename FC::Field CoordRad;
    //typedef cbasis::array3<FC> FC3;
    
  private:
    // ---- Member Field ----
    F c_;
    i3 nml_;
    FC xyz_;
    F zeta_;

  public:
    // ---- Constructors ----
    CartGTO(F c, const i3& nml, const FC& xyz, F z) : c_(c), nml_(nml), xyz_(xyz), zeta_(z) {}
    template<class F2, class FC2>
    CartGTO(const CartGTO<F2, FC2>& o): c_(o.c()), nml_(o.nml()), xyz_(o.xyz()), zeta_(o.z()) {}
    
    // ---- Accessors ----
    F c() const { return this->c_; }
    const i3& nml() const { return this->nml_; }
    int n() const {return this->nml_[0];}
    int m() const {return this->nml_[1];}
    int l() const {return this->nml_[2];}
    const FC& xyz() const { return this->xyz_; }
    FC x() const {return this->xyz_[0];}
    FC y() const {return this->xyz_[1];}
    FC z() const {return this->xyz_[2];}
    F zeta() const { return this->zeta_; }

    void set_c(F c) {this->c_ = c; } 
    void set_nml(i3 nml) {this->nml_ = nml; } 
    void set_xyz(FC xyz) {this->xyz_ = xyz; } 
    void set_z(F z) {this->z_ = z; } 

    // ---- Method ----
    F at(FC x) const {
      return c_ *
	pow(x[0]-xyz_[0], nml_[0]) *
	pow(x[1]-xyz_[1], nml_[1]) *
	pow(x[2]-xyz_[2], nml_[2]) *
	exp(-zeta_ * pow(x[0]-xyz_[0], 2) * pow(x[1]-xyz_[1], 2) * pow(x[2]-xyz_[2], 2));
    }
    std::string str() const;
    
    void SetComplexConjugate() {
      this->c_ = conjug(this->c_);
      this->z_ = conjug(this->z_);
    }
    void SetScalarProd(F c) {
      this->c_ *= c;
    }
    void SetNormalize() {
      F norm2 = CIP(*this, *this);
      this->SetScalarProd(F(1)/sqrt(norm2));
    }

  };

  // ---- External ----
  template<class F, class FC>
  std::ostream& operator << (std::ostream& os, const CartGTO<F, FC>& a);

  // ---- inner product ----
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const One&, const CartGTO<F, FC>& b) {
    /*
#ifndef USE_MOLINT
    return a.c() * b.c() * overlap(a.zeta(),
				   a.nml()[0], a.nml()[1], a.nml()[2],
				   a.xyz()[0], a.xyz()[1], a.xyz()[2],
				   b.zeta(),
				   b.nml()[0], b.nml()[1], b.nml()[2],
				   b.xyz()[0], b.xyz()[1], b.xyz()[2]);
				       
#endif
    cout << "molint" << endl;
    */
    return a.c() * b.c() * gto_overlap(a.nml()[0], a.nml()[1], a.nml()[2],
				       a.xyz()[0], a.xyz()[1], a.xyz()[2],
				       a.zeta(),
				       b.nml()[0], b.nml()[1], b.nml()[2],
				       b.xyz()[0], b.xyz()[1], b.xyz()[2],
				       b.zeta());
	

  }
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const OpKE<F, FC>&, const CartGTO<F, FC>& b) {
    /*
#ifndef USE_MOLINT
    return a.c() * b.c() * kinetic(a.zeta(),
				   a.nml()[0], a.nml()[1], a.nml()[2],
				   a.xyz()[0], a.xyz()[1], a.xyz()[2],
				   b.zeta(),
				   b.nml()[0], b.nml()[1], b.nml()[2],
				   b.xyz()[0], b.xyz()[1], b.xyz()[2]);
#endif
*/
    return a.c() * b.c() * gto_kinetic(a.nml()[0], a.nml()[1], a.nml()[2],
				       a.xyz()[0], a.xyz()[1], a.xyz()[2],
				       a.zeta(),
				       b.nml()[0], b.nml()[1], b.nml()[2],
				       b.xyz()[0], b.xyz()[1], b.xyz()[2],
				       b.zeta());

  }
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const OpNA<F, FC>& v, const CartGTO<F, FC>& b) {
/*
#ifndef USE_MOLINT
    return v.q() *  nuclear_attraction(a.xyz()[0], a.xyz()[1], a.xyz()[2],
				       a.c(),
				       a.nml()[0], a.nml()[1], a.nml()[2],
				       a.zeta(),
				       b.xyz()[0], b.xyz()[1], b.xyz()[2],
				       b.c(),
				       b.nml()[0], b.nml()[1], b.nml()[2],
				       b.zeta(),
				       v.xyz()[0], v.xyz()[1], v.xyz()[2]);
#endif 
    */
    return v.q() * a.c() *b.c() *
      gto_nuclear_attraction(a.nml()[0], a.nml()[1], a.nml()[2],
			     a.xyz()[0], a.xyz()[1], a.xyz()[2],
			     a.zeta(),
			     b.nml()[0], b.nml()[1], b.nml()[2],
			     b.xyz()[0], b.xyz()[1], b.xyz()[2],
			     b.zeta(),
			     v.xyz()[0], v.xyz()[1], v.xyz()[2]);

  }
  template<class F, class FC>
  F CIP_impl_prim(const CartGTO<F, FC>& a, const OpXyz<F, FC> & op, const CartGTO<F, FC>& b) {

    if(op.nml()[0] == 0 && op.nml()[1] == 0 && op.nml()[2] == 0) {
      return CIP_impl_prim(a, One(), b);
    }

    if(op.nml()[0] != 0 || op.nml()[1] != 0 || op.nml()[2] != 1) {
      std::string msg;
      SUB_LOCATION(msg);
      msg += "only OpXyz(0,0,1) is supported";
      std::cout << op.nml()[0] << op.nml()[1] << op.nml()[2] << std::endl;
      throw std::runtime_error(msg);
    }

    return a.c() * b.c() * (overlap(a.zeta(),
				    a.nml()[0],
				    a.nml()[1],
				    a.nml()[2] + 1,
				    a.xyz()[0], a.xyz()[1], a.xyz()[2],
				    b.zeta(),
				    b.nml()[0], b.nml()[1], b.nml()[2],
				    b.xyz()[0], b.xyz()[1], b.xyz()[2])
			    + a.xyz()[2] *
			    overlap(a.zeta(),
				    a.nml()[0],
				    a.nml()[1],
				    a.nml()[2],
				    a.xyz()[0], a.xyz()[1], a.xyz()[2],
				    b.zeta(),
				    b.nml()[0], b.nml()[1], b.nml()[2],
				    b.xyz()[0], b.xyz()[1], b.xyz()[2]));    
  }


  // ==== Spherical GTO ====
  // ---- external ----
  template<class Field, class FieldCoordOne>
  void SetSphericalGTO(int L, int M, array3<FieldCoordOne> xyz, Field zeta,
		       LinFunc<CartGTO<Field, array3<FieldCoordOne> > >* target) {

    typedef array3<FieldCoordOne> FC3;
    if(L == 0 && M == 0) {
      CartGTO<Field, FC3> gto(Field(1), i3(0, 0, 0), xyz, zeta);
      gto.SetNormalize();
      target->Add(Field(1), gto);
    } else if(L == 1 && M == 0) {
      CartGTO<Field, FC3> g3(Field(1), i3(0, 0, 1), xyz, zeta);
      g3.SetNormalize(); target->Add(1, g3);
    } else if(L == 1 && M == 1) {
      CartGTO<Field, FC3> g1(1, i3(1, 0, 0), xyz, zeta);
      g1.SetNormalize(); target->Add(1, g1);
    } else if(L == 1 && M == -1) {
      CartGTO<Field, FC3> g1(Field(1), i3(0, 1, 0), xyz, zeta);
      g1.SetNormalize(); target->Add(1, g1);
    } else if(L == 2 && M == 1) {
      CartGTO<Field, FC3> g(Field(1), i3(1, 0, 1), xyz, zeta); 
      g.SetNormalize(); target->Add(Field(1), g);
    } else if(L == 2 && M == -1) {
      CartGTO<Field, FC3> g(Field(1), i3(0, 1, 1), xyz, zeta); 
      g.SetNormalize(); target->Add(Field(1), g);
    } else if(L == 2 && M == 2) {
      CartGTO<Field, FC3> g1(Field(1), i3(2, 0, 0), xyz, zeta); 
      g1.SetNormalize(); target->Add(Field(1), g1);
      CartGTO<Field, FC3> g2(Field(1), i3(0, 2, 0), xyz, zeta); 
      g2.SetNormalize(); target->Add(-Field(1), g2);	
    } else if(L == 2 && M == -2) {
      CartGTO<Field, FC3> g1(Field(1), i3(1, 1, 0), xyz, zeta); 
      g1.SetNormalize(); target->Add(1, g1);
    } else if(L == 2 && M == 0) {
      CartGTO<Field, FC3> g1(Field(1), i3(2, 0, 0), xyz, zeta); 
      g1.SetNormalize(); target->Add(-Field(1), g1);
      CartGTO<Field, FC3> g2(Field(1), i3(0, 2, 0), xyz, zeta); 
      g2.SetNormalize(); target->Add(-Field(1), g2);
      CartGTO<Field, FC3> g3(Field(1), i3(0, 0, 2), xyz, zeta); 
      g3.SetNormalize(); target->Add(Field(2), g3);
    } else {
      std::string msg;
      SUB_LOCATION(msg);
      throw(ExceptionBadYlm(L, M, msg));
    }
    target->SetNormalize();
  }

  // ---- class ----
  template<class F, class FC1>
  class SphericalGTO :public Func<F, cbasis::array3<FC1> > {
  public:
    // ---- type ----
    typedef F Field;
    typedef FC1 FieldCoord;
    typedef cbasis::array3<FC1> FC3;

  private:
    // ---- Member Field ----
    int L_;
    int M_;
    LinFunc<CartGTO<F, FC3> > funcs_;
  public:
    SphericalGTO(int L, int M, FC3 xyz, Field zeta): L_(L), M_(M) {
      try {
	SetSphericalGTO(L, M, xyz, zeta, &funcs_);
      } catch(const ExceptionBadYlm& e) {
	throw(e);
      }
    }
    ~SphericalGTO() {}
    const LinFunc<CartGTO<F, FC3> >& GetLinFunc() const { return funcs_;}
    F at(FC3 x) const { return funcs_.at(x); }
    int L() const { return L_; }
    int M() const { return M_; }
    F zeta() const { return funcs_.begin()->second.zeta(); }
    FC3 xyz() const { return funcs_.begin()->second.xyz(); }
  };
  
  // ---- inner product ----
  template<class F, class FC>
  F CIP_impl_prim(const SphericalGTO<F, FC>& a,
		  const One&,
		  const SphericalGTO<F, FC>& b) {
    return CIP(a.GetLinFunc(), b.GetLinFunc());
  }
  template<class F, class FC>
   F CIP_impl_prim(const SphericalGTO<F, FC>& a,
		   const OpKE<F,array3<FC> >& op,
		   const SphericalGTO<F, FC>& b) {
    return CIP(a.GetLinFunc(), op, b.GetLinFunc());
  }
  template<class F, class FC>
  F CIP_impl_prim(const SphericalGTO<F, FC>& a,
		  const OpNA<F,array3<FC> >& op,
		  const SphericalGTO<F, FC>& b) {
    return CIP(a.GetLinFunc(), op, b.GetLinFunc());
  }
  template<class F, class FC>
  F CIP_impl_prim(const SphericalGTO<F, FC>& a,
		  const OpXyz<F,array3<FC> >& op,
		  const SphericalGTO<F, FC>& b) {
    return CIP(a.GetLinFunc(), op, b.GetLinFunc());
  }
}

#endif
