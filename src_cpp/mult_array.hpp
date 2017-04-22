#ifndef MULT_ARRAY_H
#define MULT_ARRAY_H

#include <string>
#include <ostream>
#include <iostream>
#include <stdexcept>
#include "../utils/macros.hpp"

/**
   MultArray gives multi indexed array. This file is header library for inline optimization.
*/

namespace cbasis {

  template<class F, int N>
  class MultArray {};

  template<class F>
  class MultArray<F, 1> {
  public:
    static const int N = 1;
    F* data_;
    int data_num_;
    int num_;
    int n0_[N];
    int n1_[N];
    std::string name_;
  public:
    MultArray(int _num) {
      data_ = new F[_num];
      data_num_ = _num;
      num_ = _num;
    }
    ~MultArray() {
      delete[] data_;
    }
    void set_name(std::string _name) { name_ = _name; }
    void get_name() { return name_; }    
    int size() const { return num_; }
    void SetRange(int n0, int n1) {
      n0_[0] = n0; n1_[0] = n1;
      num_ = n1-n0+1;
      if(data_num_ < num_) {
	std::string msg; SUB_LOCATION(msg);
	msg = "\n" + msg + " : not enough capacity";
	throw std::runtime_error(msg);
      }
    }
    F& operator()(int nx) {
      int index = nx - n0_[0];

#ifndef ARG_NO_CHECK

      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx) {
	std::string msg; SUB_LOCATION(msg);
	std::ostringstream oss;
	oss << std::endl << msg;
	oss << "index: (" << nx << ") "
	   << index << std::endl;
	throw std::runtime_error(oss.str());
      }
#endif
    return data_[index];
    }
    
  };

  template<class F>
  class MultArray<F, 2> {
  public:
    static const int N = 2;
    F* data_;
    int data_num_;  // size of data_ (Capacity)
    int num_;       // Prod_i(n1_[i] - n0_[i] + 1)
    int n0_[N];
    int n1_[N];
    std::string name_;
  public:
    MultArray(int _num0) {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
    }
    MultArray(const MultArray<F, 2>& o) {
      this->data_ = new F[o.data_num_];
      for(int i = 0; i < o.num_; i++) {
	this->data_[i] = o.data_[i];
      }
      this->data_num = o.data_num_;
      this->num_ = o.num_;
      for(int n = 0; n < N; n++) {
	int n0_[n] = o.n0_[n];
	int n1_[n] = o.n1_[n];
      }
      this->name_ = o.name_;
    }
    ~MultArray() {
      delete[] data_;
    }
    void set_name(std::string _name) { name_ = _name; }
    std::string get_name() { return name_; }    
    int size() const { return num_; }
    void SetRange(int nx0, int nx1, int ny0, int ny1) {
      n0_[0] = nx0; n0_[1] = ny0;
      n1_[0] = nx1; n1_[1] = ny1;
      num_ = (nx1-nx0+1)*(ny1-ny0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
      }
    }
    F& operator()(int nx, int ny) {
      int index = ((nx - n0_[0]) * (n1_[1] - n0_[1] + 1) + 
		   (ny - n0_[1]));
#ifndef ARG_NO_CHECK

      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny) { 
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << "index: (" << nx << ", " << ny << ") "
	   << index << std::endl;
	msg += ss.str();
	throw std::runtime_error(msg);
      }
#endif
      return data_[index];
    }
    std::string str() {
      std::ostringstream oss;
      oss << "MultArray<2> object:" << std::endl;
      oss << "  name: " << this->get_name() << std::endl;
      oss << "  size: " << this->size() << std::endl;
      for(int n0 = n0_[0]; n0 < n1_[0]+1; n0++)
	for(int n1 = n0_[1]; n1 < n1_[1]+1; n1++){
	  oss << n0 << n1 << " : "
	      << (*this)(n0, n1)
	      << std::endl;
	}
      return oss.str();
    }    
  };

  template<class F>
  class MultArray<F, 3> {
  public:
    static const int N = 3;
    F* data_;
    int data_num_;
    int num_;
    int n0_[N];
    int n1_[N];
    std::string name_;
  public:
    MultArray(int _num0) {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
    }
    MultArray(int _num0, std::string _name) {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
      name_ = _name;
    }    
    ~MultArray() {
      delete[] data_;
    }
    void set_name(std::string _name) { name_ = _name; }
    std::string get_name() { return name_; }
    void SetRange(int nx0, int nx1, int ny0, int ny1, int nz0, int nz1) {
      n0_[0] = nx0; n0_[1] = ny0; n0_[2] = nz0; 
      n1_[0] = nx1; n1_[1] = ny1; n1_[2] = nz1; 
      this->num_ = (nx1-nx0+1)*(ny1-ny0+1)*(nz1-nz0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
	/*
	std::string msg; SUB_LOCATION(msg);
	std::stringstream ss;
	ss << std::endl;
	ss << msg << " : not enough capcaity" << std::endl;
	ss << "name: " << this->get_name() << std::endl;
	for(int i = 0; i < N; i++) {
	  ss << "n" << i << "s = (" << n0_[0] << ", " << n1_[0] << ")" << std::endl;
	}
	throw std::runtime_error(ss.str());
	*/
      }

    }
    int size() const { return num_; }
    F& operator()(int nx, int ny, int nz) {
      int index = ((nx - n0_[0]) * (n1_[2] - n0_[2] + 1) * (n1_[1] - n0_[1] + 1) + 
		   (ny - n0_[1]) * (n1_[2] - n0_[2] + 1) + 
		   (nz - n0_[2]));
#ifndef ARG_NO_CHECK

      if(index < 0   || num_-1 < index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny ||
	 nz < n0_[2] || n1_[2] < nz) {
	std::string msg; SUB_LOCATION(msg);
	std::stringstream ss;
	ss << std::endl;
	ss << msg;
	ss << "index out of range" << std::endl;
	ss << "name: " << this->get_name() << std::endl;
	ss << "(n0,n1,n2) = (" << nx << ", " << ny << ", " << nz << ") " << std::endl;
	ss << "index = " << index << std::endl;
	for(int i = 0; i < N; i++) {
	  ss << "n" << i << "s = (" << n0_[0] << ", " << n1_[0] << ")" << std::endl;
	}
	   
	throw std::runtime_error(ss.str());
      }
#endif

      return data_[index];
    }
    void Show() const {
      std::cout << "MultArray<3>" << std::endl;
      for(int n = 0; n < N; n++) {
	std::cout << n << ": " << n0_[n] << ", " << n1_[n] << std::endl;
      }
      
    }
  };

  template<class F>
  class MultArray<F, 4> {
  public:
    static const int N = 4;
    F* data_;
    int data_num_;
    int num_;
    int n0_[N];
    int n1_[N];
    std::string name_;
  public:
    MultArray(int _num0, std::string _name="") {
      data_ = new F[_num0];
      data_num_ = _num0;
      num_ = _num0;
      n0_[0] = 0; n1_[0] = _num0;
      for(int i = 1; i < N; i++) {
	n0_[i] = 0; n0_[i] = _num0;
      }
      name_ = _name;
    }
    ~MultArray() {
      delete[] data_;
    }
    int size() const { return num_; }
    void set_name(std::string _name) { name_ = _name; }
    std::string get_name() { return name_; }
    void SetRange(int nx0, int nx1, int ny0, int ny1,
		  int nz0, int nz1, int nw0, int nw1) {
      n0_[0] = nx0; n0_[1] = ny0; n0_[2] = nz0; n0_[3] = nw0; 
      n1_[0] = nx1; n1_[1] = ny1; n1_[2] = nz1; n1_[3] = nw1; 
      num_ = (nx1-nx0+1)*(ny1-ny0+1)*(nz1-nz0+1)*(nw1-nw0+1);
      if(data_num_ < num_) {
	delete[] data_;
	data_num_ = num_;
	data_ = new F[num_];
	/*
	std::string msg; SUB_LOCATION(msg);
	msg = "\n" + msg + " : not enough capacity";
	throw std::runtime_error(msg);
	*/
      }
    }
    void SetValue(F val) {
      for(int i = 0; i < data_num_; i++)
	data_[i] = val;
    }
    int index(int nx, int ny, int nz, int nw) {
      int num1 = n1_[1] - n0_[1] + 1;
      int num2 = n1_[2] - n0_[2] + 1;
      int num3 = n1_[3] - n0_[3] + 1;
      int _index = ((nx - n0_[0]) * num3 * num2 * num1 + 
		   (ny - n0_[1]) * num3 * num2 + 
		   (nz - n0_[2]) * num3 +
		   (nw - n0_[3]));

#ifndef ARG_NO_CHECK
      if(_index < 0   || num_-1 < _index ||
	 nx < n0_[0] || n1_[0] < nx ||
	 ny < n0_[1] || n1_[1] < ny ||
	 nz < n0_[2] || n1_[2] < nz ||
	 nw < n0_[3] || n1_[3] < nw) {
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << ": index: (" << nx << ", " << ny << ", " << nz << ", " << nw << ") "
	   << _index << std::endl;
	for(int i = 0; i < N; i++)
	  ss << n0_[i] << ", " << n1_[i] << " ";
	msg += ss.str();
	throw std::out_of_range(msg);
      }
#endif      
      
      return _index;
    }
    F& operator()(int nx, int ny, int nz, int nw) {
      int _index = this->index(nx, ny, nz, nw);
      return data_[_index];
    }
    std::string str() {
      std::ostringstream oss;
      oss << "MultArray<4> object:" << std::endl;
      oss << "  name: " << this->get_name() << std::endl;
      oss << "  size: " << this->size() << std::endl;
      for(int n0 = n0_[0]; n0 < n1_[0]+1; n0++)
	for(int n1 = n0_[1]; n1 < n1_[1]+1; n1++)
	  for(int n2 = n0_[2]; n2 < n1_[2]+1; n2++)
	    for(int n3 = n0_[3]; n3 < n1_[3]+1; n3++) {
	      oss << n0 << n1 << n2 << n3 << " : "
		  << this->index(n0, n1, n2, n3) << " : "
		  << (*this)(n0, n1, n2, n3)
		  << std::endl;
	    }
      return oss.str();
    }
  };

  template<class F, int N>
  F MultArrayTDot(MultArray<F, N>& a, MultArray<F, N>& b) {
    for(int n = 0; n < N; n++) {
      if(a.n0_[n] != b.n0_[n] || a.n1_[n] != b.n1_[n]) {
	std::string msg; SUB_LOCATION(msg);
	msg = "\n" + msg + " : size mismatch";
	throw std::runtime_error(msg);	
      }
    }      
    int n(a.size());
    F acc(0);
    for(int i = 0; i < n; i++) {
      acc += a.data_[i] * b.data_[i];
    }
    return acc;
  }
}
/*
namespace future {
  
  template<class F, int N>
  class MultArray {

    // ---- Typedef ----
  public:
    typedef MultArray<F, N-1> Element;

    // ---- Field member ----
    Element* data;
    int step;
    int n0;
    int n1;

    // ---- Constructors ----
    MultArray(int num): elements(num), n0(0), n1(num-1) {}
    static Create
    ~MultArray() {}
    int size() const { return elements.size(); }
    Element& SetRange(int _n0, int _n1) { n0 = _n0; n1 = _n1; }
    void CheckArg(int index, int n) {
      int num(this->size());
      if(index < 0   || num-1 < index) {
	std::string msg;
	std::stringstream ss;
	SUB_LOCATION(msg);
	ss << "index: " << n << ", "
	   << index << std::endl;
	msg += ss.str();
	throw std::runtime_error(msg);
      }
    }
    Element& operator[](int n) {
      int index = (n-n0);
#ifndef ARG_NO_CHECK
      CheckArg(index, n);
#endif
      return elements[index];
    }
    const Element& operator[](int n) const {
      int index = n-n0;
#ifndef ARG_NO_CHECK
      CheckArg(index, n);
#endif
      return elements[index];      
    }
  };
  
}
*/
#endif
