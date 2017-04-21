#ifndef READ_JSON_HPP
#define READ_JSON_HPP

#include "../external/picojson/picojson.h"

namespace cbasis {

  template<class T>
  void CheckValue(picojson::value& val, int n=-1, int m=-1);
  
  template<class T>
  void CheckObject(picojson::object& obj, std::string k, int n=-1, int m=-1);
  
  void CheckJson_string(picojson::object& obj, std::string k);
  void CheckJson_double(picojson::object& obj, std::string k);
  void CheckJson_complex(picojson::object& obj, std::string k);
  void CheckJson_array(picojson::object& obj, std::string k);
  void CheckJson_VectorXcd(picojson::value& json);

  template<class T>
  void CheckObject(picojson::object& obj, std::string k, int n=-1, int m=-1) {
    
    std::string key = "key \"" + k + "\"";
    if(obj.find(k) == obj.end()) {
      throw(std::runtime_error(key + " not found."));
    }
    
    try {
      CheckValue<T>(obj[k], n, m);
    } catch(std::exception& e) {
      std::string msg = "error on parsing value of " + key + "\n";
      msg += e.what();
      throw(std::runtime_error(msg));
    }
  }  

  template<class T>
  T ReadJson(picojson::value& json, int n=-1, int m=-1);
  template<class T>
  T ReadJson(picojson::object& json, std::string key, int n=-1, int m=-1) {
    if(json.find(key) == json.end()) {
      throw std::runtime_error("key \"" + key + "\" not found");
    }
    T t;
    try {
      t = ReadJson<T>(json[key], n, m);
    } catch(std::exception& e) {
      std::ostringstream oss;
      oss << "error on parsing value of key \"" << key << "\"\n";
      oss << e.what() << "\n";
      throw std::runtime_error(oss.str());
    }
    return t;
  }
  template<class T>
  T ReadJsonWithDefault(picojson::object& obj, std::string,
			T t, int n=-1, int m=-1);  

  
  template<class T>
  picojson::value ToJson(T& t);

}

#endif
