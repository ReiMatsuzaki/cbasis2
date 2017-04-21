#ifndef SYMMOL_READ_JSON_HPP
#define SYMMOL_READ_JSON_HPP

#include "symmolint.hpp"
#include "../utils/read_json.hpp"

namespace cbasis {

  Reduction ReadJson_Reduction(picojson::value& json, SymmetryGroup sym);
  void ReadJson_Molecule(picojson::value& json, Molecule mole);
  void ReadJson_Molecule(picojson::object& obj, std::string k, Molecule mole);
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos);
  void ReadJson_SymGTOs_Subs(picojson::object& obj, std::string k, SymGTOs gtos);

  void ReadJson_Orbital_file(picojson::object& initstate,
			     SymmetryGroup _sym, dcomplex *_E0,
			     Eigen::VectorXcd *_c0, Irrep *_irrep0, int *_i0);
  void ReadJson_Orbital(picojson::object& obj, string k, SymmetryGroup _sym,
			dcomplex *_E0, Eigen::VectorXcd *_c0,
			Irrep *_irrep0, int *_i0);  
  
}

#endif
