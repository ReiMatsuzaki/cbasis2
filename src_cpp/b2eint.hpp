#ifndef B2EINT_H
#define B2EINT_H

//#include "symmolint.hpp"
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include "../utils/macros.hpp"
#include "../utils/typedef.hpp"

using namespace std;

namespace cbasis {

  /**
    Interface for store of two electron integrals.
   */
  class IB2EInt {
  public:
    IB2EInt(){};
    virtual ~IB2EInt();

    /*
      Obtain next index list and value.
     */
    virtual bool Get(int *ib, int *jb, int *kb, int *lb,
		     int *i, int *j, int *k, int *l, int *type, dcomplex *val) = 0;
    /*
      Set value at given index list.
     */
    virtual bool Set(int ib, int jb, int kb, int lb,
	     int i, int j, int k, int l, dcomplex val) = 0;
    /*
      Reset internal index counter "idx_".
     */
    virtual void Reset() = 0;
    virtual void Init(int num) = 0;
    /*
      Obtain value at given index list. This function is for debug and
      very slow for numerical calcualtion. 
      Additionaly, in this function Reset method is called.
     */    
    dcomplex At(int ib, int jb, int kb, int lb,
		int i, int j, int k, int l);
    bool Exist(int ib, int jb, int kb, int lb,
	       int i, int j, int k, int l);
    /*
      Write to file
     */
    virtual void Write(std::string fn) = 0;

    /*
      Accessor for size and capacity
     */
    //    int basis_size_irrep_i(int irrep) const;
    //    int basis_size_irrep_j(int irrep) const;
    //    int basis_size_irrep_k(int irrep) const;
    //    int basis_size_irrep_l(int irrep) const;
    //    int num_irrep() const;
    virtual int size() const = 0;
    virtual int capacity() const = 0;
  };

  class B2EIntMem :public IB2EInt {
  private:
    int capacity_; // capacity of each array.
    int size_;     // size of data
    int idx_;      // used for Get function.
    std::vector<int> ibs, jbs, kbs, lbs;
    std::vector<int> is, js, ks, ls;
    std::vector<int> ts;
    std::vector<dcomplex> vs;
  public:
    B2EIntMem();
    B2EIntMem(int num);
    ~B2EIntMem();

    void Init(int num);
    bool Get(int *ib, int *jb, int *kb, int *lb,
	     int *i, int *j, int *k, int *l, int *type, dcomplex *val);
    bool Set(int ib, int jb, int kb, int lb,
	     int i, int j, int k, int l, dcomplex val);
    dcomplex At(int ib, int jb, int kb, int lb,
		int i, int j, int k, int l);    
    void Reset();
    void Write(std::string fn);
    int size() const;
    int capacity() const;
    
  };

  typedef boost::shared_ptr<IB2EInt> B2EInt;

  B2EInt ERIRead(std::string fn);
}
#endif
