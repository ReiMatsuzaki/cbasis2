#ifndef BASIS_SET_TEMPLATE_H
#define BASIS_SET_TEMPLATE_H

#include <vector>
#include <Eigen/Core>

namespace cbasis {

  template<class FuncT>
  class BasisSet {
  public:
    // ---- type ----
    typedef typename FuncT::Field Field;
    typedef Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic> Mat;
    typedef Eigen::Matrix<Field, Eigen::Dynamic, 1> Vec;
    typedef typename std::vector<FuncT>::const_iterator const_iterator;

  private:
    std::vector<FuncT> us_; // basis set
    //    std::vector<FuncT> vs_; // complex conjugate of adjoint functions

  public:
    int size() const { return us_.size(); }
    const_iterator begin() const { return us_.begin(); }
    const_iterator end() const { return us_.end(); }
    template<class OpT>
    Mat BuildMatrix(OpT op) {

      int n(size());
      Mat m = Mat::Zero(n, n);
	
      int i = 0; 
      for(const_iterator it_i = begin(); it_i != end(); ++it_i, ++i) {
	int j = 0;
	for(const_iterator it_j = begin(); it_j != end(); ++it_j, ++j) {
	  m(i, j) = CIP(*it_i, op, *it_j);
	}
      }

      return m;
    }
    template<class OpT, class FuncOther>
    Mat BuildMatrix(OpT op, const BasisSet<FuncOther>& other) {

      Mat m(this->size(), other.size());
      int i = 0; 
      for(const_iterator it_i = begin(); it_i != end(); ++it_i, ++i) {
	int j = 0;
	for(const_iterator it_j = other.begin(); it_j != other.end(); ++it_j, ++j) {
	  m(i, j) = CIP(*it_i, op, *it_j);
	}
      }      

      return m;
    }
    template<class FuncOther>
    Vec BuildVector(const FuncOther& other) {

      Vec v(size());

      int i(0);
      for(const_iterator it = begin(); it != end(); ++it, ++i) {
	v(i)  = CIP(*it, other);
      }

    }
    
  };


}

#endif
