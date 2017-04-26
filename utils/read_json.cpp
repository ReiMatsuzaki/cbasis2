#include <fstream>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include "../utils/macros.hpp"
#include "../utils/eigen_plus.hpp"
#include "read_json.hpp"


using namespace std;
using namespace picojson;
using namespace Eigen;
using boost::format;

namespace cbasis {

  template<> void CheckValue<object>(picojson::value& val, int n, int m) {
    
    if(not val.is<object>()) {
      throw(runtime_error("value is not object. "));
    }
    if(n > 0) {
      if((int)val.get<object>().size() != n)
	throw(runtime_error("invalid size object."));
    }
  }
  template<> void CheckValue<string>(picojson::value& val, int n, int m) {
    if(not val.is<string>()) {
      throw(runtime_error("value is not string. "));
    }    
  }
  template<> void CheckValue<double>(picojson::value& val, int n, int m) {
    
    if(not val.is<double>()) {
      throw(runtime_error("value is not double."));
    }    
  }
  template<> void CheckValue<dcomplex>(picojson::value& val, int n, int m) {
    string msg = "value is not compelx (double or two element array of double).";
    if(val.is<array>()) {
      array& ary = val.get<array>();
      if(ary.size() != 2) {
	throw(runtime_error(msg));
      }
      if(not ary[0].is<double>())
	throw(runtime_error(msg));
      if(not ary[1].is<double>())
	throw(runtime_error(msg));
    } else if(not val.is<double>()) {
      throw runtime_error(msg);
    }
  }
  template<> void CheckValue<int>(picojson::value& val, int n, int m) {
    string msg = "value is not int. ";
    if(not val.is<double>()) {
      throw runtime_error(msg);
    }
    double x(val.get<double>());
    double y((double)(int)x);
    double eps(0.0000000001);
    if(abs(x-y) > eps) {
      throw runtime_error(msg);
    }
  }
  template<> void CheckValue<array>(picojson::value& val, int n, int m) {
    if(not val.is<array>()) {
      throw(runtime_error("value  is not array. "));
    }

    if(n > 0) {
      if((int)val.get<array>().size() != n) {
	throw(runtime_error("invalid size of array. "));
      }
    }
  }
  template<> void CheckValue<VectorXcd>(picojson::value& val, int n, int m) {

    if(val.is<array>()) {
      try {
	CheckValue<array>(val, n);
      } catch(exception& e) {
	throw(runtime_error("value is not array for VectorXcd."));
      }
      array& ary = val.get<array>();
      for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
	try {
	  CheckValue<dcomplex>(*it);
	} catch(...) {
	  throw(runtime_error("element of value is not dcomplex for VectorXcd."));
	}
      }
      
    } else {
      throw runtime_error("value must be array");
    }
  }
  template<> void CheckValue<VectorXi>(picojson::value& val, int n, int m) {

    if(val.is<array>()) {
      array& ary = val.get<array>();
      for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
	try {
	  CheckValue<int>(*it);
	} catch(...) {
	  throw(runtime_error("element is not int for VectorXi."));
	}
      }
    } else {
      throw(runtime_error("value is not array for VectorXi."));
    }
    
  }
  template<> void CheckValue<MatrixXi>(picojson::value& val, int n, int m) {
    
    try {
      CheckValue<array>(val, -1);
    } catch(exception& e) {
      throw(runtime_error("value is not array for MatrixXi"));
    }
    array& ary = val.get<array>();
    
    if(n > 0) 
      if(m != (int)ary.size()) 
	throw(runtime_error("invalid rows size for MatrixXi"));
	
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<array>(*it, n);
      } catch(...) {
	throw(runtime_error("element of value is not array for MatrixXi."));
      }
      array& aryary = it->get<array>();
      if(m > 0) {
	if(m != (int)aryary.size()) {
	  throw(runtime_error("invalid cols size of element for MatrixXi."));
	}
      }
      for(array::iterator jt = aryary.begin(); jt != aryary.end(); ++jt) {
	try {
	  CheckValue<int>(*jt);
	} catch(...) {
	  throw(runtime_error("element of element of value is not int for MatrixXi"));
	}
      }
    }
        
  }
  template<> void CheckValue<MatrixXcd>(picojson::value& val, int n, int m) {
    
    try {
      CheckValue<array>(val, -1);
    } catch(exception& e) {
      throw(runtime_error("value is not array for MatrixXcd"));
    }
    array& ary = val.get<array>();
    for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
      try {
	CheckValue<array>(*it, n);
      } catch(...) {
	string msg = (
		      format("%d th element is not array")
		      % distance(ary.begin(), it)
		      ).str();
	throw runtime_error(msg);
	//throw(runtime_error("element of value is not array for MatrixXcd."));
      }
      array& aryary = it->get<array>();
      if(n > 0) {
	if(n != (int)aryary.size()) {
	  throw(runtime_error("invalid size of element for MatrixXcd."));
	}
      }
      for(array::iterator jt = aryary.begin(); jt != aryary.end(); ++jt) {
	try {
	  CheckValue<dcomplex>(*jt);
	} catch(exception& e) {
	  string msg = (
			format("element (%d, %d) cannot convert to complex\n")
			% distance(ary.begin(), it)
			% distance(aryary.begin(), jt)
			).str();
	  msg += e.what();
	  throw runtime_error(msg);
	}
      }
    }  
  }

  template void CheckObject<object>(object&, string, int , int);

  template<> string ReadJson<string>(value& json, int n, int m) {
    CheckValue<string>(json);
    return json.get<string>();
  }
  template<> double ReadJson<double>(value& json, int n, int m) {
    CheckValue<double>(json);
    return json.get<double>();
  }  
  template<> int ReadJson<int>(value& json, int n, int m) {
    CheckValue<int>(json, n, m);
    return (int)json.get<double>();;
  }
  template<> bool ReadJson<bool>(value& json, int n, int m) {    
    if(not json.is<bool>()) {
      throw runtime_error("value is not bool");
    }
    return json.get<bool>();
  }  
  template<> dcomplex ReadJson<dcomplex>(value& json, int n, int m) {
    CheckValue<dcomplex>(json);
    dcomplex x(0);
    if(json.is<double>()) {
      x = json.get<double>();
    } else if(json.is<array>()) {
      array& ary = json.get<array>();
      x = dcomplex(ary[0].get<double>(), ary[1].get<double>());
    }   
    return x;
  }  
  template<> VectorXcd ReadJson<VectorXcd>(value& json, int n, int m) {

    if(json.is<array>()) {
      //      CheckValue<VectorXcd>(json, n);
      array& ary = json.get<array>();
      VectorXcd vec = VectorXcd::Zero(ary.size());
      
      for(int i = 0; i <(int) ary.size(); i++) {
	vec.coeffRef(i) = ReadJson<dcomplex>(ary[i]);
      }
      return vec;

      
    } else if(json.is<object>()) {
      object& obj = json.get<object>();
      string type = ReadJson<string>(obj, "type");
      if(type == "lin") {
	dcomplex x0, dx;
	x0 = ReadJson<dcomplex>(obj, "x0");
	int N = ReadJson<int>(obj, "N");
	if(obj.find("dx") != obj.end()) {
	  dx = ReadJson<dcomplex>(obj, "dx");
	} else if(obj.find("x1") != obj.end()) {
	  dcomplex xN = ReadJson<dcomplex>(obj, "xN");
	  dx = (xN-x0)/dcomplex(N);
	} else {
	  throw runtime_error("dx or xN is necessary in ReadJson<VectorXcd>");
	}
	
	VectorXcd vec(N+1);
	//	dcomplex dx = (xN-x0)/(N-1);
	for(int i =0; i<=N; i++) {
	  vec(i) = x0 + dx * dcomplex(i);
	}
	return vec;
      } else if(type == "geo") {
	dcomplex x0 = ReadJson<dcomplex>(obj, "x0");
	dcomplex r = ReadJson<dcomplex>(obj, "r");
	int N = ReadJson<int>(obj, "N");
	VectorXcd vec(N);
	//	dcomplex r = (xN-x0)/(N-1);
	for(int i =0; i<N; i++) {
	  vec(i) = x0 * pow(r, i);
	}
	return vec;
      } else {
	throw runtime_error("type must be [lin or geo]");
      }
    } else {
      throw runtime_error("value must be array or object for VectorXcd");
    }
    
  }
  template<> VectorXd ReadJson<VectorXd>(value& json, int n, int m) {

    if(json.is<array>()) {

      array& ary = json.get<array>();
      VectorXd vec = VectorXd::Zero(ary.size());
      for(int i = 0; i <(int) ary.size(); i++) {
	vec.coeffRef(i) = ReadJson<double>(ary[i]);
      }
      return vec;
      
    } else if(json.is<object>()) {

      object& obj = json.get<object>();
      string type = ReadJson<string>(obj, "type");
      double x0, dx;
      int N;
      VectorXd vec;
      if(type == "lin") {
	if(obj.find("x0") != obj.end() &&
	   obj.find("dx") != obj.end() &&
	   obj.find("N") != obj.end() ) {
	  x0 = ReadJson<double>(obj, "x0");
	  dx = ReadJson<double>(obj, "dx");
	  N = ReadJson<int>(obj, "N");
	  
	} else if(obj.find("x0") != obj.end() &&
		  obj.find("x1") != obj.end() &&
		  obj.find("N") != obj.end() ) {
	  x0 = ReadJson<double>(obj, "x0");
	  double x1 = ReadJson<double>(obj, "x1");
	  N = ReadJson<int>(obj, "N");
	  dx = (x1-x0) / (N-1);	  
	} else {
	  throw runtime_error("set of options (x0,dx,N) or (x0,x1,N) is necessary");
	}
	vec = VectorXd::Zero(N);
	for(int i =0; i<N; i++) {
	  vec(i) = x0 + dx * i;
	}
	return vec;
	
      } else if(type == "geo") {
	double x0 = ReadJson<double>(obj, "x0");
	double r  = ReadJson<double>(obj, "r");
	int N = ReadJson<int>(obj, "N");
	VectorXd vec(N);
	for(int i =0; i<N; i++) {
	  vec(i) = x0 * pow(r, i);
	}
	return vec;
      
      } else {
	throw runtime_error("value must be array or object for VectorXd");
      }
    }
    throw runtime_error("un resolved error");
  }
  template<> MatrixXcd ReadJson<MatrixXcd>(value& json, int _n, int _m) {
    CheckValue<MatrixXcd>(json, _n, _m);
    array& ary = json.get<array>();
    int n(0), m(0);
    n = ary.size();
    m = ary[0].get<array>().size();

    MatrixXcd mat = MatrixXcd::Zero(n, m);

    for(int i = 0; i < n; i++) {
      array& aryary = ary[i].get<array>();
      for(int j = 0; j < m; j++) {
	mat.coeffRef(i, j) = ReadJson<dcomplex>(aryary[j]);
      }
    }
    
    return mat;
  }
  template<> MatrixXd ReadJson<MatrixXd>(value& json, int _n, int _m) {
    CheckValue<array>(json, _n);
    array& ary = json.get<array>();
    int n = ary.size();

    if(n == 0) {
      throw runtime_error("no element for MatrixXd");
    }

    CheckValue<array>(ary[0], _m);
    array& aryary0 = ary[0].get<array>();
    int m = aryary0.size();

    MatrixXd mat = MatrixXd::Zero(n, m);    
    for(int i = 0; i < n; i++) {
      CheckValue<array>(ary[i], m);
      array& aryary = ary[i].get<array>();
      for(int j = 0; j < m; j++) {
	mat.coeffRef(i, j) = ReadJson<double>(aryary[j]);
      }
    }
    
    return mat;
  }  
  template<> VectorXi ReadJson<VectorXi>(value& json, int n, int m) {
    
    CheckValue<VectorXi>(json, n);

    if(json.is<array>()) {
      array& ary = json.get<array>();
      VectorXi vec = VectorXi::Zero(ary.size());
      
      for(int i = 0; i < (int)ary.size(); i++) {
	int x(ary[i].get<double>());
	vec.coeffRef(i) = x;
      }
      return vec;
    } else {
      throw runtime_error("value is not array for VectorXi");
    }
    
  }  
  template<> MatrixXi ReadJson<MatrixXi>(value& json, int _n, int _m) {
    CheckValue<MatrixXi>(json, _n, _m);
    array& ary = json.get<array>();
    int n(0), m(0);
    n = ary.size();
    m = ary[0].get<array>().size();
    MatrixXi mat = MatrixXi::Zero(n, m);

    for(int i = 0; i < n; i++) {
      array& aryary = ary[i].get<array>();
      for(int j = 0; j < m; j++) {
	mat.coeffRef(i, j) = (int)aryary[j].get<double>();
      }
    }
    
    return mat;
  }
  template<> LinearSolver ReadJson<LinearSolver>(value& json, int n, int m) {
    CheckValue<object>(json);
    object& obj = json.get<object>();

    string type = ReadJson<string>(obj, "type");
    return LinearSolver(type);
  }
  
  template int ReadJson<int>(object& o, string k, int n, int m);
  template VectorXd ReadJson<VectorXd>(object& o, string k, int n, int m);
  template double ReadJson<double>(object& o, string k, int n, int m);
  template bool ReadJson<bool>(object& o, string k, int n, int m);
  template LinearSolver ReadJson<LinearSolver>(object& o, string k, int n, int m);

  template<class T>
  T ReadJsonWithDefault(object& json, string key, T t, int n, int m) {
    if(json.find(key) == json.end()) {
      return t;
    } else {
      return ReadJson<T>(json[key], n, m);
    }
  }
  template LinearSolver
  ReadJsonWithDefault<LinearSolver>(object&, string, LinearSolver, int, int);
  template bool
  ReadJsonWithDefault<bool>(object&, string, bool, int, int);    
  template string ReadJsonWithDefault<string>(object&, string, string, int, int);
  template int ReadJsonWithDefault<int>(object&, string, int, int, int);
  template double ReadJsonWithDefault<double>(object&, string, double, int, int);
  
  template<>
  value ToJson<dcomplex>(dcomplex& x) {
    double eps(0.000000000000001);
    if(abs(x.imag()) < eps ) {
      double re_x(x.real());
      return value(re_x);
    } else {
      array xs(2);
      xs[0] = value(x.real());
      xs[1] = value(x.imag());
      return value(xs);
    }
  }
  template<>
  value ToJson<int>(int& x) {
    double y(x);
    return value(y);
  }  

}
