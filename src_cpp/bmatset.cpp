#include <string>
#include <ostream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <boost/format.hpp>
#include "../utils/macros.hpp"
#include "../utils/eigen_plus.hpp"
#include "bmatset.hpp"

using namespace std;
using namespace Eigen;
using boost::format;

namespace cbasis {

  // ==== BVec ====  
  void BVec::Write(string filename) const {

    ofstream f;
    f.open(filename.c_str(), ios::out|ios::binary|ios::trunc);

    if(!f) {
      string msg; SUB_LOCATION(msg); msg = "\n" + msg + "file not found";
      throw runtime_error(msg);
    }

    int id(6677);
    f.write((char*)&id, sizeof(int));
    
    int num = this->size();
    f.write((char*)&num, sizeof(int));

    for(const_iterator it = this->begin(); it != this->end(); ++it) {
      int irrep = it->first;
      f.write((char*)&irrep, sizeof(int));
      const VectorXcd& V = it->second;
      int n = V.size();
      f.write((char*)&n, sizeof(int));
      for(int i = 0; i < n; i++) {
	f.write((char*)&V.coeff(i), sizeof(dcomplex));
      }
    }
  }
  void BVec::Read(string filename) {
    ifstream f(filename.c_str(), ios::in|ios::binary);
    
    if(!f) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "file not found. filename: " + filename;
      throw runtime_error(msg);
    }

    int id;
    f.read((char*)&id, sizeof(int));
    if(id != 6677) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + ": invalid format. filename: " + filename;
      throw runtime_error(msg);
    }

    int num;
    f.read((char*)&num, sizeof(int));
    for(int i = 0; i < num; i++) {
      
      int irrep, n;
      f.read((char*)&irrep, sizeof(int));
      f.read((char*)&n, sizeof(int));
      this->at(irrep) = VectorXcd::Zero(n);
      VectorXcd& V = (*this)(irrep);
      for(int i = 0; i < n; i++) {
	dcomplex v;
	f.read((char*)&v, sizeof(dcomplex));
	V(i) = v;
      }
    }

  }
  void BVec::Shift(dcomplex d) {
    for(BVec::iterator it = this->begin(); it != this->end(); ++it) {	
      VectorXcd& vec = it->second;
      for(int i = 0; i < vec.size(); i++) {
	vec(i) += d;
      }
    }
  }
  void BVec::Add(dcomplex a, const BVec& o) {
    for(BVec::iterator it = this->begin(); it != this->end(); ++it) {	
      int k = it->first;
      VectorXcd& vec = it->second;
      const VectorXcd& ovec = o[k];
      for(int i = 0; i < vec.size(); i++) {
	vec(i) += a * ovec[i];
      }
    }
  }
  ostream& operator << (ostream& os, const BVec& a) {
    os << "==== BVec ====" << endl;
    os << "Block vector object" << endl;
    os << "name: " << a.get_name() << endl;
    for(BVec::const_iterator it = a.begin(); it != a.end(); ++it) {
      int irrep = it->first;
      const VectorXcd& vec = it->second;
      os << "irrep = " << irrep << endl;
      os << vec << endl;
    }
    os << "==============" << endl;
    return os;
  }
  void BVecSqrt(const BVec& a, BVec *b) {
    for(BVec::const_iterator it = a.begin();
	it != a.end(); ++it) {
      int k = it->first;
      const VectorXcd& vec = it->second;      
      if(not b->has_block(it->first)) {
	(*b)(k) = vec;
      }
      for(int i = 0; i < vec.size(); i++) {
	(*b)(k)(i) = sqrt(vec(i));
      }
    }
  }  
  // ==== BMat ====
  void BMat::Clear() {
    this->map_.clear();
  }
  void BMat::swap(BMat& o) {
    this->name_.swap(o.name_);
    this->map_.swap(o.map_);
  }
  bool BMat::is_block_diagonal() const {
    typedef const_iterator It;
    bool acc(true);
    for(It it = begin(); it != end(); ++it) {
      Key k = it->first;
      acc = acc && (k.first == k.second);
    }
    return acc;
  }
  void BMat::Write(string filename) const {
    
    ofstream f;
    f.open(filename.c_str(), ios::out|ios::binary|ios::trunc);

    if(!f) {
      string msg; SUB_LOCATION(msg); msg = "\n" + msg + "file not found";
      throw runtime_error(msg);
    }

    int id(668778);
    f.write((char*)&id, sizeof(int));
    int num = this->size();
    f.write((char*)&num, sizeof(int));

    for(const_iterator it = this->begin(); it != this->end(); ++it) {
      int irrep = it->first.first;
      int jrrep = it->first.second;
      f.write((char*)&irrep, sizeof(int));
      f.write((char*)&jrrep, sizeof(int));
      const MatrixXcd& M = it->second;
      int n = M.rows();
      int m = M.cols();
      f.write((char*)&n, sizeof(int));
      f.write((char*)&m, sizeof(int));
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  f.write((char*)&M.coeff(i, j), sizeof(dcomplex));
	}
      }
    }
  }
  void BMat::Read(string filename) {
    ifstream f(filename.c_str(), ios::in|ios::binary);
    
    if(!f) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "file not found. filename: " + filename;
      throw runtime_error(msg);
    }

    int id;
    f.read((char*)&id, sizeof(int));
    if(id != 668778) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + "invalid format. filename: "+filename;
      throw runtime_error(msg);
    }
    
    int num;
    f.read((char*)&num, sizeof(int));
    for(int i = 0; i < num; i++) {
      
      int irrep, jrrep, n, m;
      f.read((char*)&irrep, sizeof(int));
      f.read((char*)&jrrep, sizeof(int));
      f.read((char*)&n, sizeof(int));
      f.read((char*)&m, sizeof(int));
      (*this)(irrep, jrrep) = MatrixXcd::Zero(n, m);
      MatrixXcd& M = (*this)(irrep, jrrep);
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  dcomplex v;
	  f.read((char*)&v, sizeof(dcomplex));
	  M(i, j) = v;
	}
      }
    }

  }
  void BMat::SetZero() {
    for(iterator it = this->begin(); it != this->end(); ++it) {
      it->second.setZero();
    }
  }
  void BMat::SetId() {
    for(iterator it = this->begin(); it != this->end(); ++it) {
      Key k = it->first;
      MatrixXcd& mat = it->second;
      int n = mat.rows();
      int m = mat.cols();
      if(k.first == k.second) {
	it->second = MatrixXcd::Identity(n, m);
      } else {
	it->second = MatrixXcd::Zero(n, m);
      }
    }
    void SetId();
  }
  void BMat::Scale(dcomplex c) {
    
    for(iterator it = this->begin(); it != this->end(); ++it) {	
      Key k = it->first;
      MatrixXcd& m = it->second;
      m = c * m;
    }
  }
  void BMat::Add(dcomplex c, const BMat& o) {
    
    for(iterator it = this->begin(); it != this->end(); ++it) {	
      if(not o.has_block(it->first)) {
	stringstream ss;
	Key kthis = it->first;
	ss << "block mismatch\n";
	ss << format("key(this): (%d, %d)\n") % kthis.first % kthis.second;
	ss << "this:\n";
	ss << *this;
	ss << "other:\n";
	ss << o;
	THROW_ERROR(ss.str());
      }
      Key k = it->first;
      if(it->second.rows() != o[k].rows()) {
	THROW_ERROR("size mismatch");
      }
      if(it->second.cols() != o[k].cols()) {
	THROW_ERROR("size mismatch");
      }
      (*this)[k] += c * o[k];
    }
  }
  void BMat::Shift(dcomplex c) {
    for(iterator it = this->begin(); it != this->end(); ++it) {

      if(it->first.first == it->first.second) {
	MatrixXcd& M = it->second;
	if(M.rows() != M.cols()) {
	  THROW_ERROR("Shift works only square matrix");
	}
	int n = M.rows();
	for(int i = 0; i < n; i++) {
	  M(i, i) += c;
	}
      }                      
    }
  }
  bool BMat::is_same_structure(const BMat& o) const {
    
    for(const_iterator k = this->begin(); k != this->end(); ++k) {
      if(not o.has_block(k->first))
	return false;
    }

    for(const_iterator k = o.begin(); k != o.end(); ++k) {
      if(not o.has_block(k->first))
	return false;
    }

    for(const_iterator k = o.begin(); k != o.end(); ++k) {
      const MatrixXcd& a = (*this)[k->first];
      const MatrixXcd& b = k->second;
      if(a.rows() != b.rows() or a.cols() != b.cols())
	return false;
    }

    return true;
  }
  void BMatDiag(const BVec& a, BMat *m) {

    for(BVec::const_iterator it_a = a.begin();
	it_a != a.end(); ++it_a) {
      int irr = it_a->first;
      const VectorXcd& vec = it_a->second;
      int n = vec.size();
      (*m)(irr, irr) = MatrixXcd::Zero(n, n);
      for(int i = 0; i < n; ++i) {
	(*m)(irr, irr)(i, i) = vec(i);
      }
    }
    
  }
  void Multi(const BMat& a, const BMat& b, BMat& c) {
    typedef BMat::const_iterator It;
    typedef BMat::Key Key;
    for(It ia = a.begin(); ia != a.end(); ++ia) {
      for(It ib = b.begin(); ib != b.end(); ++ib) {	
	Key irr_a = ia->first;
	Key irr_b = ib->first;
	if(irr_a.second != irr_b.first) {
	  continue;
	}

	const MatrixXcd& am = ia->second;
	const MatrixXcd& bm = ib->second;

	if(am.cols() != bm.rows()) {
	    stringstream ss;
	    ss << "size mismatch.\n";
	    ss << format("a = (%d, %d)\n") % am.rows() % am.cols();
	    ss << format("b = (%d, %d)\n") % bm.rows() % bm.cols();
	    THROW_ERROR(ss.str());
	}
	
	Key irr_c(irr_a.first, irr_b.second);
	if(not c.has_block(irr_c) ||
	   c[irr_c].rows() != am.rows() ||
	   c[irr_c].cols() != bm.cols()) {
	  c[irr_c] = MatrixXcd::Zero(am.rows(), bm.cols());
	}
      }
    }

    c.SetZero();

    for(It ia = a.begin(); ia != a.end(); ++ia) {
      for(It ib = b.begin(); ib != b.end(); ++ib) {
	Key irr_a = ia->first;
	Key irr_b = ib->first;
	if(irr_a.second != irr_b.first) {
	  continue;
	}
	Key irr_c(irr_a.first, irr_b.second);
	c[irr_c] += a[irr_a] * b[irr_b];
      }
    }
  }
  void Multi3(const BMat& a, const BMat& b, const BMat& c, BMat& res) {
    
    typedef BMat::const_iterator It;
    typedef BMat::Key Key;

    // -- build block matrix if necessary --
    for(It ia = a.begin(); ia != a.end(); ++ia) {
      for(It ib = b.begin(); ib != b.end(); ++ib) {
	for(It ic = c.begin(); ic != c.end(); ++ic) {	
	  Key irr_a = ia->first;
	  Key irr_b = ib->first;
	  Key irr_c = ic->first;
	  if(irr_a.second != irr_b.first || irr_b.second != irr_c.first) {	     
	    continue;
	  }
	  const MatrixXcd& am = ia->second;
	  const MatrixXcd& bm = ib->second;
	  const MatrixXcd& cm = ic->second;

	  if(am.cols() != bm.rows() ||
	     bm.cols() != cm.rows() ) {
	    stringstream ss;
	    ss << "size mismatch.\n";
	    ss << format("a = (%d, %d)\n") % am.rows() % am.cols();
	    ss << format("b = (%d, %d)\n") % bm.rows() % bm.cols();
	    ss << format("c = (%d, %d)\n") % cm.rows() % cm.cols();
	    THROW_ERROR(ss.str());
	  } 
	  
	  Key irr_res(irr_a.first, irr_c.second);
	  if(not res.has_block(irr_res)) {
	    int n(ia->second.rows());
	    int m(ic->second.cols());
	    res[irr_res] = MatrixXcd::Zero(n, m);
	  }
	}
      }
    }

    // -- set all values to zero --
    res.SetZero();

    // -- calculation --
    for(It ia = a.begin(); ia != a.end(); ++ia) {
      for(It ib = b.begin(); ib != b.end(); ++ib) {
	for(It ic = c.begin(); ic != c.end(); ++ic) {	
	  Key irr_a = ia->first;
	  Key irr_b = ib->first;
	  Key irr_c = ic->first;
	  if(irr_a.second != irr_b.first || irr_b.second != irr_c.first) {	     
	    continue;
	  }
	  Key irr_res(irr_a.first, irr_c.second);
	  res[irr_res] += a[irr_a] * b[irr_b] * c[irr_c];	    
	}
      }
    }
  }
  ostream& operator << (ostream& os, const BMat& a) {
    os << "==== BMat ====" << endl;
    os << "Block matrix object" << endl;
    os << "name: " << a.get_name() << endl;
    os << "non0_block: " << a.size() << endl;
    for(BMat::const_iterator it = a.begin(); it != a.end(); ++it) {
      BMat::Key key = it->first;
      int irrep = key.first;
      int jrrep = key.second;
      const MatrixXcd& mat = it->second;
      os << format("(irrep, jrrep; inum, jnum) = (%d, %d, %d, %d)\n")
	% irrep % jrrep % mat.rows() % mat.cols();
      if(not a.get_simple_print())
	os << mat << endl;
    }
    os << "==============" << endl;
    return os;
  }
  void Copy(const BMat& a, BMat& b) {
    b.Clear();
    for(BMat::const_iterator it = a.begin(); it != a.end(); ++it) {      
      b[it->first]  = it->second;
    }	
  }  
  void BMatSqrt(const BMat& a, BMat& b) {
    
    if(not a.is_block_diagonal()) {
      THROW_ERROR("a must be block diagonal");
    }

    typedef BMat::const_iterator It;
    typedef BMat::Key Key;
    
    for(It it = a.begin(); it != a.end(); ++it) {
      Key k = it->first;
      if(not b.has_block(k)) {
	b[k] = a[k];
      }
    }

    for(It it = a.begin(); it != a.end(); ++it) {
      Key k = it->first;
      matrix_sqrt(a[k], b[k]);
    }
    
  }
  void BMatInvSqrt(const BMat& a, BMat& b) {
    
    if(not a.is_block_diagonal()) {
      THROW_ERROR("a must be block diagonal");
    }

    typedef BMat::const_iterator It;
    typedef BMat::Key Key;
    
    for(It it = a.begin(); it != a.end(); ++it) {
      Key k = it->first;
      if(not b.has_block(k)) {
	b[k] = a[k];
      }
    }

    for(It it = a.begin(); it != a.end(); ++it) {
      Key k = it->first;
      matrix_inv_sqrt(a[k], &b[k]);
    }
    
  }
  void BMatCtAC(const BMat& C, const BMat& A, BMat *res) {
    if(not C.is_block_diagonal()) {
      THROW_ERROR("C must be block diagonal");      
    }
    if(not C.is_same_structure(A)) {
      THROW_ERROR("C and A must be same structure");
    }

    typedef BMat::const_iterator It;
    typedef BMat::Key Key;
    
    for(It it = A.begin(); it != A.end(); ++it) {
      Key k = it->first;
      if(not res->has_block(k)) {
	(*res)[k] = A[k];
      }
    }

    for(It it = A.begin(); it != A.end(); ++it) {
      Key k = it->first;
      CtAC(C[k], A[k], &(*res)[k]);
    }
    
  }
  void BMatCtAD(const BMat& C, const BMat& D, const BMat& A, BMat *res) {
    
    if(not C.is_block_diagonal()) {
      THROW_ERROR("C must be block diagonal");      
    }
    if(not D.is_block_diagonal()) {
      THROW_ERROR("D must be block diagonal");      
    }
    
    typedef BMat::const_iterator It;
    typedef BMat::Key Key;
    
    for(It it = A.begin(); it != A.end(); ++it) {
      Key k = it->first;
      if(not res->has_block(k)) {
	(*res)[k] = MatrixXcd::Zero(C[k].rows(), D[k].rows());
      }
    }

    for(It itc = C.begin(); itc != C.end(); ++itc) {
      for(It ita = A.begin(); ita != A.end(); ++ita) {
	for(It itd = D.begin(); itd != D.end(); ++itd) {
	  Key kc = itc->first;
	  Key ka = ita->first;
	  Key kd = itd->first;
	  const MatrixXcd& mc = itc->second;
	  const MatrixXcd& ma = ita->second;
	  const MatrixXcd& md = itd->second;
	  if(kc.first == ka.first and ka.second == kd.first) {
	    Key kres(kc.second, kd.second);
	    CtAD(mc, ma, md, &(*res)[kres]);
	  }
	}
      }
    }
  }
  void BMatRead(BMat::Map& bmat, string fn) {
    ifstream f(fn.c_str(), ios::in|ios::binary);
    
    if(!f) {
      string msg; SUB_LOCATION(msg);
      msg += ": file not found";
      throw runtime_error(msg);
    }
    
    int id;
    f.read((char*)&id, sizeof(int));
    if(id != 8371) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + ": invalid format";
      throw runtime_error(msg);
    }
    int num;
    f.read((char*)&num, sizeof(int));
    for(int i = 0; i < num; i++) {
      int irrep, jrrep, n, m;
      f.read((char*)&irrep, sizeof(int));
      f.read((char*)&jrrep, sizeof(int));
      f.read((char*)&n, sizeof(int));
      f.read((char*)&m, sizeof(int));
      
      MatrixXcd M(n, m);
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  dcomplex v;
	  f.read((char*)&v, sizeof(dcomplex));
	  M(i, j) = v;
	}
      }
      bmat[make_pair(irrep, jrrep)] = MatrixXcd::Zero(1, 1);
      bmat[make_pair(irrep, jrrep)].swap(M);
    }
  }
  void BMatWrite(BMat::Map& bmat, string fn) {

    ofstream f;
    f.open(fn.c_str(), ios::out|ios::binary|ios::trunc);
    
    if(!f) {
      string msg; SUB_LOCATION(msg); msg+=": file not found";
      throw runtime_error(msg);
    }
    int id(8371);
    f.write((char*)&id, sizeof(int));
    int num = bmat.size();
    f.write((char*)&num, sizeof(int));

    for(BMat::iterator it = bmat.begin(); it != bmat.end(); ++it) {
      int irrep = it->first.first;
      int jrrep = it->first.second;
      f.write((char*)&irrep, sizeof(int));
      f.write((char*)&jrrep, sizeof(int));
      MatrixXcd& M = it->second;
      int n = M.rows();
      int m = M.cols();
      f.write((char*)&n, sizeof(int));
      f.write((char*)&m, sizeof(int));
      for(int i = 0; i < n; i++) {
	for(int j = 0; j < m; j++) {
	  f.write((char*)&M.coeff(i, j), sizeof(dcomplex));
	}
      }
    }
  }
  void BMatEigenSolve(const BMat& H, const BMat& S, BMat *C, BVec *E) {

    if(not H.is_block_diagonal()) {
      THROW_ERROR("H must be block diagonal");
    }
    if(not S.is_same_structure(H)) {
      THROW_ERROR("S must be same structure with H");
    }
    
    for(BMat::const_iterator it = H.begin(); it != H.end(); ++it) {
      BMat::Key k = it->first;
      int irr = k.first;
      SymGenComplexEigenSolver solver(H[k], S[k]);
      (*E)(irr) = solver.eigenvalues();
      (*C)[k]   = solver.eigenvectors();
    }
  }
  void BMatEigenSolve(const BMat& H, BMat *C, BVec *E) {
    BMat S;
    Copy(H, S);
    S.SetId();
    BMatEigenSolve(H, S, C, E);
  }
  // ==== BlockMatrixSets ====
  _BMatSet::_BMatSet(): block_num_(1) {}
  _BMatSet::_BMatSet(int _block_num): block_num_(_block_num) {}
  void _BMatSet::SetMatrix(string name, int i, int j, MatrixXcd& mat) {
#ifndef ARG_NO_CHECK
#endif

    MatrixXcd tmp;
    mat_map_[name][make_pair(i, j)] = tmp;
    mat_map_[name][make_pair(i, j)].swap(mat);
    
    if(block_num_ <= i)
      block_num_ = i+1;
    if(block_num_ <= j)
      block_num_ = j+1;    

  }
  bool _BMatSet::Exist(string name, int i, int j) {

    if(i < 0 || block_num_ <= i ||
       j < 0 || block_num_ <= j) {
      string msg; SUB_LOCATION(msg);
      stringstream ss;
      ss << msg;
      ss << ": Error for int i or j" << endl;
      ss << "i: " << i << endl;
      ss << "j: " << j << endl;
      ss << "block_num: " << block_num_ << endl;
      throw runtime_error(ss.str());
    }

    if(mat_map_.find(name) == mat_map_.end()) 
      return false;

    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) 
      return false;

    MatrixXcd& mat = mat_map_[name][make_pair(i, j)];
    if(mat.rows() == 0 || mat.cols() == 0)
      return false;
    
    return true;
  }
  const MatrixXcd& _BMatSet::GetMatrix(string name, int i, int j) {

    if(mat_map_.find(name) == mat_map_.end()) {
      string msg; SUB_LOCATION(msg);
      stringstream ss;
      ss << ": failed to find name." << endl
	 << "name : " << name <<endl;
      msg += ss.str();
      throw runtime_error(msg);
    }
    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) {
      string msg; SUB_LOCATION(msg);
      stringstream ss;
      ss << ": failed to find index.\n"
	 << "i : " << i << endl
	 << "j : " << j << endl;
      
      msg += ": failed to find index.\n";
      msg += ss.str();
      throw runtime_error(msg);
    }

    return mat_map_[name][make_pair(i, j)];
  }
  const BMat& _BMatSet::GetBlockMatrix(std::string name) {
    return mat_map_[name];
  }
  void _BMatSet::SelfAdd(string name, int i, int j, int a, int b, dcomplex v) {

#ifndef ARG_NO_CHECK
    if(i < 0 || block_num_ <= i ||
       j < 0 || block_num_ <= j) {
      string msg; SUB_LOCATION(msg);
      msg += "Error for int i or j.";
      throw runtime_error(msg);
    }

    if(mat_map_.find(name) == mat_map_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += " : key not found.\n" ;
      msg += "name : " + name;
      throw runtime_error(msg);
    }

    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) {
      string msg; SUB_LOCATION(msg);
      msg += " : matrix (i,j) is not found.\n";
      throw runtime_error(msg);
    }
#endif

    mat_map_[name][make_pair(i, j)](a, b) += v;
  }
  dcomplex _BMatSet::GetValue(string name, int i, int j, int a, int b) {

#ifndef ARG_NO_CHECK
    if(mat_map_.find(name) == mat_map_.end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": key not found.\n" ;
      msg += "name : " + name;
      throw runtime_error(msg);
    }

    if(i < 0 || block_num_ <= i ||
       j < 0 || block_num_ <= j) {
      string msg; SUB_LOCATION(msg);
      stringstream oss;
      oss << msg << endl << "Error for int i or j" << endl;
      oss << "i = " << i << endl;
      oss << "j = " << j << endl;
      oss << "block_num = " << block_num_ << endl;
      throw runtime_error(oss.str());
    }

    if(mat_map_[name].find(make_pair(i, j)) == mat_map_[name].end()) {
      string msg; SUB_LOCATION(msg);
      msg += ": matrix (i,j) is not found.\n";
      throw runtime_error(msg);
    }
#endif

    return mat_map_[name][make_pair(i, j)](a, b);
  } 
  void _BMatSet::swap(_BMatSet& o) {
    ::swap(this->block_num_, o.block_num_);
    this->mat_map_.swap(o.mat_map_);
    
    /*
    typedef BMatMap::iterator It;
    BMatMap tmp_o;
    for(It it_o = o.mat_map_.begin(); it_o != o.mat_map_.end();) {
      tmp_o[it_o->first] = BMat::Map();
      ::swap(tmp_o[it_o->first], it_o->second);
      o.mat_map_.erase(it_o++);
    }
    
    for(It it_this = this->mat_map_.begin(); it_this != this->mat_map_.end();) {
      o.mat_map_[it_this->first] = BMat::Map();
      o.mat_map_[it_this->first].swap(it_this->second);
      this->mat_map_.erase(it_this++);
    }

    for(It it_tmp = tmp_o.begin(); it_tmp != tmp_o.end();) {
      this->mat_map_[it_tmp->first] = BMat::Map();
      this->mat_map_[it_tmp->first].swap(it_tmp->second);
      ++it_tmp;
    }

    int tmp = o.block_num_;
    o.block_num_ = this->block_num_;
    this->block_num_ = tmp;
    */
  }
  string _BMatSet::str() const {
    ostringstream oss;
    for(BMatMap::const_iterator it = mat_map_.begin(); it != mat_map_.end(); ++it) {
      oss << it->first << endl;
      for(BMat::const_iterator it_mat = it->second.begin(); it_mat != it->second.end();
	  ++it_mat) {
	oss << it_mat->first.first << ", " << it_mat->first.second<< endl
	    << it_mat->second << endl;
      }
    }
    return oss.str();
  }
	      
  // ==== External ====
  void swap(BMat::Map& a, BMat::Map& b) {

    typedef map<pair<int, int>, MatrixXcd>::iterator It;
    map<pair<int, int>, MatrixXcd> tmp_a;
    for(It it_a = a.begin(); it_a != a.end();) {
      tmp_a[it_a->first] = MatrixXcd();
      tmp_a[it_a->first].swap(it_a->second);
      a.erase(it_a++);
    }

    for(It it_b = b.begin(); it_b != b.end();) {
      a[it_b->first] = MatrixXcd();
      a[it_b->first].swap(it_b->second);
      b.erase(it_b++);
    }

    for(It it_tmp = tmp_a.begin(); it_tmp != tmp_a.end();) {
      b[it_tmp->first] = MatrixXcd();
      b[it_tmp->first].swap(it_tmp->second);
      it_tmp++;
    }
  }

}
