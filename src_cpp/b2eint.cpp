#include <fstream>
#include <iostream>
#include <stdexcept>
#include "b2eint.hpp"

namespace cbasis {

  // ==== Interface ====
  IB2EInt::~IB2EInt() {}
  dcomplex IB2EInt::At(int ib, int jb, int kb, int lb,
		       int i, int j, int k, int l) {
    
    int xib,xjb,xkb,xlb,xi,xj,xk,xl,xt;
    dcomplex xv;
    this->Reset();
    while(this->Get(&xib,&xjb,&xkb,&xlb,&xi,&xj,&xk,&xl,&xt, &xv)) {
      if(ib == xib && jb == xjb && kb == xkb && lb == xlb &&
	 i == xi   && j == xj   && k == xk   && l == xl) {
	return xv;
      }
    }
    string msg; SUB_LOCATION(msg); 
    msg += ": failed to find given index list.";
    throw runtime_error(msg);
  }  
  bool IB2EInt::Exist(int ib, int jb, int kb, int lb,
		      int i, int j, int k, int l) {
    try {
      this->At(ib, jb, kb, lb, i, j, k, l);
    } catch(exception& e) {
      return false;
    }
    return true;
  }

  // ==== Mem version ====
  // ---- Constructors ----
  B2EIntMem::B2EIntMem() {
    this->Init(1);
  }
  B2EIntMem::B2EIntMem(int n) {
    this->Init(n);
  }
  
  B2EIntMem::~B2EIntMem() {}

  // ---- Main ----
  void B2EIntMem::Init(int num) {
    capacity_ = num;
    size_ = 0;
    idx_ = 0;
    ibs.reserve(num);
    jbs.reserve(num);
    kbs.reserve(num);
    lbs.reserve(num);

    is.reserve(num);
    js.reserve(num);
    ks.reserve(num);
    ls.reserve(num);

    ts.reserve(num);
    vs.reserve(num);
  }
  bool B2EIntMem::Get(int *ib, int *jb, int *kb, int *lb,
		      int *i, int *j, int *k, int *l,
		      int *type, dcomplex *val) {
    if(this->idx_ >= this->size_ ) {
      return false;
    }
    *ib = this->ibs[this->idx_];
    *jb = this->jbs[this->idx_];
    *kb = this->kbs[this->idx_];
    *lb = this->lbs[this->idx_];
    *i  = this->is[this->idx_];
    *j  = this->js[this->idx_];
    *k  = this->ks[this->idx_];
    *l  = this->ls[this->idx_];
    *type = this->ts[this->idx_];
    *val  = this->vs[this->idx_];
    this->idx_++;
    return true;
  }
  bool B2EIntMem::Set(int ib, int jb, int kb, int lb,
		      int i, int j, int k, int l,
		      dcomplex val) {
    this->ibs[this->size_] = ib;
    this->jbs[this->size_] = jb;
    this->kbs[this->size_] = kb;
    this->lbs[this->size_] = lb;
    this->is[this->size_]  = i;
    this->js[this->size_]  = j;
    this->ks[this->size_]  = k;
    this->ls[this->size_]  = l;
    this->ts[this->size_]  = 0;
    this->vs[this->size_]  = val;
    this->size_++;
    return true;
  }
  void B2EIntMem::Reset() {
    idx_ = 0;
  }
  void B2EIntMem::Write(string fn) {

    ofstream f;
    f.open(fn.c_str(), ios::out|ios::binary|ios::trunc);
    
    if(!f) {
      string msg; SUB_LOCATION(msg); msg+=": file not found";
      throw runtime_error(msg);
    }

    int num(this->size());
    f.write((char*)&num, sizeof(int));

    this->Reset();
    int ib,jb,kb,lb,i,j,k,l,t;
    dcomplex v;
    this->Reset();
    while(this->Get(&ib,&jb,&kb,&lb,&i,&j,&k,&l, &t, &v)) {
      f.write((char*)&ib, sizeof(int));
      f.write((char*)&jb, sizeof(int));
      f.write((char*)&kb, sizeof(int));
      f.write((char*)&lb, sizeof(int));
      f.write((char*)&i,  sizeof(int));
      f.write((char*)&j,  sizeof(int));
      f.write((char*)&k,  sizeof(int));
      f.write((char*)&l,  sizeof(int));
      f.write((char*)&t,  sizeof(int));
      f.write((char*)&v,  sizeof(dcomplex));
    }
    f.close();

  }
  int B2EIntMem::size() const {
    return this->size_;
  }
  int B2EIntMem::capacity() const {
    return this->capacity_;
  }

  // ==== ERI read ====
  B2EInt ERIRead(string fn) {

    B2EInt eri(new B2EIntMem);

    ifstream f(fn.c_str(), ios::in|ios::binary);
    
    if(!f) {
      string msg; SUB_LOCATION(msg);
      msg += ": file not found";
      throw runtime_error(msg);
    }

    int num;
    f.read((char*)&num, sizeof(int));

    //this->Init(num);
    eri->Init(num);

    for(int idx = 0; idx < num; idx++) {
      int ib,jb,kb,lb,i,j,k,l,t;
      dcomplex v;      
      f.read((char*)&ib, sizeof(int));
      f.read((char*)&jb, sizeof(int));
      f.read((char*)&kb, sizeof(int));
      f.read((char*)&lb, sizeof(int));

      f.read((char*)&i, sizeof(int));
      f.read((char*)&j, sizeof(int));
      f.read((char*)&k, sizeof(int));
      f.read((char*)&l, sizeof(int));

      f.read((char*)&t, sizeof(int));
      f.read((char*)&v, sizeof(dcomplex));
      eri->Set(ib, jb, kb, lb, i, j, k, l, v);
    }

    return eri;
    
  }

}
