#ifndef BMATSET_H
#define  BMATSET_H

#include <map>
#include <Eigen/Core>
#include <boost/shared_ptr.hpp>
#include "../utils/typedef.hpp"

/**
   BMatSet is Block Matrix set.
   if ARG_NO_CHECK is defined, given argument for map and vector is not checked.
*/

namespace cbasis {

  // ==== Block Vector ====
  class BVec {
  public:
    typedef std::map<int, Eigen::VectorXcd>::iterator iterator;
    typedef std::map<int, Eigen::VectorXcd>::const_iterator const_iterator;
  private:
    std::string name_;
    std::map<int, Eigen::VectorXcd> map_;
  public:
    BVec() {}
    BVec(std::string _name) :name_(_name) {}
    iterator begin() { return map_.begin(); }
    const_iterator begin() const { return map_.begin(); }
    iterator end() { return map_.end(); }
    const_iterator end() const { return map_.end(); }
    void set_name(std::string _name) { name_ = _name; }
    std::string get_name() const { return name_; }
    bool has_block(int irrep) { return map_.find(irrep) != map_.end(); }
    int size() const { return map_.size(); }
    Eigen::VectorXcd& at(int irrep) {return map_[irrep];}
    const Eigen::VectorXcd& at(int irrep) const {return map_.find(irrep)->second;}    
    Eigen::VectorXcd& operator()(int irrep) {return at(irrep); }
    const Eigen::VectorXcd& operator()(int irrep) const {return at(irrep); }
    Eigen::VectorXcd& operator[] (int irrep) { return at(irrep); }
    const Eigen::VectorXcd& operator[](int irrep) const {return at(irrep); }
    void Write(std::string filename) const;
    void Read(std::string filename);
    void Shift(dcomplex d);
    void Add(dcomplex a, const BVec& o);
  };
  std::ostream& operator << (std::ostream& os, const BVec& a);
  void BVecSqrt(const BVec& a, BVec *b);  

  // ==== Block Matrix ====
  class BMat {
  public:
    typedef std::pair<int, int> Key;
    typedef Eigen::MatrixXcd Value;
    typedef std::map<Key, Value> Map;
    typedef Map::iterator iterator;
    typedef Map::const_iterator const_iterator;
    
  private:
    std::string name_;
    Map map_;
    bool simple_print_;
  public:
    BMat(): simple_print_(false) {}
    BMat(std::string _name): name_(_name), simple_print_(false) {}
    void Clear();
    iterator begin() { return map_.begin(); }
    const_iterator begin() const { return map_.begin(); }
    iterator end() { return map_.end(); }
    const_iterator end() const { return map_.end(); }
    iterator find(Key& k) { return map_.find(k); }
    const_iterator find(const Key& k) const { return map_.find(k); }
    void set_name(std::string _name) { name_ = _name; }
    BMat& simple_print() { simple_print_=true; return *this; }
    BMat& full_print() { simple_print_=false; return *this; }
    bool get_simple_print() const { return simple_print_;}
    std::string get_name() const { return name_; }
    bool has_block(Key ijrrep) const {
      return (map_.find(ijrrep) != map_.end()) &&
	(map_.find(ijrrep)->second.rows() != 0);
    }    
    bool has_block(int irrep, int jrrep) const {
      return this->has_block(std::make_pair(irrep, jrrep));
    }
    bool is_same_structure(const BMat& o) const;
    bool is_block_diagonal() const;
    int size() const { return map_.size(); }
    Value& operator()(int irrep, int jrrep) {
      return map_[std::make_pair(irrep, jrrep)];
    }
    const Value& operator()(int irrep, int jrrep) const {
      return map_.find(std::make_pair(irrep, jrrep))->second;
    }
    Value& operator[] (Key k) { return map_[k]; }
    const Value& operator[](Key k) const {return map_.find(k)->second; }
    void SetZero();
    void SetId();
    void Scale(dcomplex c);
    void Add(dcomplex c, const BMat&);
    void Shift(dcomplex c);    
    void swap(BMat& o);
    void Write(std::string filename) const;
    void Read(std::string filename);
    void SimplePrint() const;
  };
  std::ostream& operator << (std::ostream& os, const BMat& a);
  void Copy(const BMat& a, BMat& b);
  void BMatDiag(const BVec& a, BMat *m);
  void Multi(const BMat& a, const BMat& b, BMat& c);
  void Multi3(const BMat& a, const BMat& b, const BMat& c, BMat& res);
  void BMatInvSqrt(const BMat& a, BMat& b);
  void BMatSqrt(const BMat& a, BMat& b);
  void BMatCtAC(const BMat& C, const BMat& A, BMat *res);
  void BMatCtAD(const BMat& C, const BMat& D, const BMat& A, BMat *res);
  void BMatCtA(const BMat& C, const BMat& A, BMat *res);
  void BVecCtx(const BMat& C, const BVec& x, BVec *y);
  void BMatEigenSolve(const BMat& H, const BMat& S, BMat *C, BVec *E);
  void BMatEigenSolve(const BMat& H, BMat *C, BVec *E);
  
  // ==== Old ====
  void BMatRead(BMat::Map& bmat, std::string fn);
  void BMatWrite(BMat::Map& bmat, std::string fn);
  typedef std::map<std::string, BMat> BMatMap;
  class _BMatSet {
  private:
    int     block_num_;
    BMatMap mat_map_;
  public:
    _BMatSet();
    _BMatSet(int _block_num);
    int block_num() const { return block_num_; }
    void SetMatrix(std::string name, int i, int j, Eigen::MatrixXcd& a);
    bool Exist(std::string, int i, int j);
    const Eigen::MatrixXcd& GetMatrix(std::string name, int i, int j);
    const BMat& GetBlockMatrix(std::string name);
    void SelfAdd(std::string name, int ib, int jb, int i, int j, dcomplex v);
    dcomplex GetValue(std::string name, int ib, int jb, int i, int j);
    void swap(_BMatSet& o);
    std::string str() const;
  };
  typedef boost::shared_ptr<_BMatSet> BMatSet;
}

#endif
