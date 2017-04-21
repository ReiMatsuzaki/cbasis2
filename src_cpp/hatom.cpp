#include "hatom.hpp"
#include <boost/lexical_cast.hpp>
#include "macros.hpp"
#include "op.hpp"
#include "op_func.hpp"


namespace cbasis {

  template<class F>
  HLikeAtom<F>::HLikeAtom(): n_(1), l_(0), z_() {}
  template<class F>
  HLikeAtom<F>::HLikeAtom(int _n, int _l): n_(_n), l_(_l), z_() {}
  template<class F>
  HLikeAtom<F>::HLikeAtom(int _n, int _l, F _z): n_(_n), l_(_l), z_(_z) {}
  
  template<class F>
  typename HLikeAtom<F>::STOs HLikeAtom<F>::EigenState() const {

      STOs lc;

      if( n_ == 1 && l_ ==0) {
      
	lc.Add(F(1), STO(F(2), 1, F(1)));
      
      } else if ( n_ == 2 && l_ == 0) {

	lc.Add(F(1), STO(F(1)/sqrt(F(2)), 1, F(1)/F(2)));
	lc.Add(F(1), STO(-F(1)/(F(2)*sqrt(F(2))), 2, F(1)/F(2)));
      
      } else if ( n_ == 2 && l_ == 1) {
	
	lc.Add(F(1), STO(F(1) / (F(2)*sqrt(F(6))), 2, F(1)/F(2)));
	
      } else if (n_ == 3 && l_ == 0) {

	lc.Add(F(1), STO(F(4)/(F(81)*sqrt(F(30))), 3, F(1)/F(3)));
	
      } else if (n_ == 3 && l_ == 1) {
	F c = F(8) / (F(27) * sqrt(F(6)));
	F z = F(1) / F(3);
	lc.Add(F(1), STO(c, 2, z));
	lc.Add(F(1), STO(-c/F(6), 3, z));
	
      } else if (n_ == 3 && l_ == 2) {
	
	F c = F(1) / F(81) * sqrt(F(8)/F(15));
	lc.Add(F(1), STO(c, 3, F(1)/F(3)));
	
      } else {
	
	string msg;
	SUB_LOCATION(msg);
	msg += "\ninputted n and l is not supported for HLikeAtom::EigenState\n";
	msg += "n: ";
	msg += boost::lexical_cast<string>(n_);
	msg += "\nl: ";
	msg += boost::lexical_cast<string>(l_);
	throw std::invalid_argument(msg);
      
      }
      
      return lc;
    
    }
    
  template<class F>
  LinFunc<typename  HLikeAtom<F>::STOs> HLikeAtom<F>::DipoleInitLength(int l1) const {

    if(l_ != l1 + 1 && l_ != l1 - 1) {
	
      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitL";
      throw std::invalid_argument(msg);
      
    }
      
    STOs psi_n = this->EigenState();
    return OP(OpRm(1), psi_n);
    
  }

  template<class F>
  LinFunc<typename  HLikeAtom<F>::STOs>  HLikeAtom<F>::DipoleInitVelocity(int l1) const {

      if(l_ != l1 + 1 && l_ != l1 - 1) {
	string msg;
	msg = "|l0 - l1| != 1 in HAtomDipoleInitV";
	throw std::invalid_argument(msg);
      }
    
      STOs psi_n = this->EigenState();
      LinFunc<STOs> res = OP(OpD1(), psi_n);
      
      if( l_ > 0) {
	F coef = F(l_ < l1 ? - (l_ + 1) : l_);
	STOs b = Expand(OP(OpRm(-1), psi_n));
	res.Add(coef, b);
      }
      
      return res;    	
    }
  
  template<class F>
  F HLikeAtom<F>::EigenEnergy() const {
      return -F(1) / F(2 * n_ * n_); 
    }

  template<class F>
  typename HLikeAtom<F>::HamiltonianOp HLikeAtom<F>::Hamiltonian() const {

      return AddOp(ProdOp(       -F(1)/F(2),    OpD2()),
		   AddOp(ProdOp(  -z_,          OpRm(-1)),
			 ProdOp( F(l_*l_-l_)/F(2), OpRm(-2))));

    }

  template<class F>
  typename HLikeAtom<F>::HMinusEnergyOp HLikeAtom<F>::HMinusEnergy(F ene) const {

    return AddOp(Hamiltonian(), 
		 ProdOp(-ene, OpRm(0)));

  }


  // ----------- explicit instance ---------
  typedef std::complex<double> CD;
  template class HLikeAtom<double>;
  template class HLikeAtom<CD>;

  /*
  template Op<RSTO> HLikeAtom<double>::Hamiltonian<RSTO>();
  template Op<CSTO> HLikeAtom<CD>::Hamiltonian<CSTO>();
  template Op<RGTO> HLikeAtom<double>::Hamiltonian<RGTO>();
  template Op<CGTO> HLikeAtom<CD>::Hamiltonian<CGTO>();

  template Op<RSTO> HLikeAtom<double>::HMinusEnergy<RSTO>(double e);
  template Op<CSTO> HLikeAtom<CD>::HMinusEnergy<CSTO>(CD);
  template Op<RGTO> HLikeAtom<double>::HMinusEnergy<RGTO>(double);
  template Op<CGTO> HLikeAtom<CD>::HMinusEnergy<CGTO>(CD);
  */
  
}
