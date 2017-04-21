#include "symmol_read_json.hpp"

using namespace std;
using namespace picojson;
using namespace Eigen;

namespace cbasis {

  template<> ERIMethod ReadJson<ERIMethod>(value& json, int _n, int _m) {
    ERIMethod method;
    CheckValue<object>(json);
    object& obj = json.get<object>();
    if(obj.find("symmetry") != obj.end()) {
      method.set_symmetry(ReadJson<int>(obj, "symmetry"));
    }
    if(obj.find("memo") != obj.end()) {
      method.set_coef_R_memo(ReadJson<int>(obj, "memo"));
    }
    if(obj.find("perm") != obj.end()) {
      method.set_perm(ReadJson<int>(obj, "perm"));
    }
    return method;
  }  

  template<> SymmetryGroup ReadJson<SymmetryGroup>( value& v, int n, int m) {
    string sym;
    try {
      sym = ReadJson<string>(v);
    } catch(exception& e) {
      string msg = "error on parsing SymmetryGroup.\n";
      msg += e.what();
      throw(runtime_error(msg));
    }
    
    if(sym == "C1")
      return SymmetryGroup_C1();
    else if(sym == "Cs")
      return SymmetryGroup_Cs();
    else if(sym == "C2h")
      return SymmetryGroup_Cs();
    else if(sym == "D2h")
      return SymmetryGroup_D2h();
    else
      throw(runtime_error("error on parsing SymmetryGroup.\nsym must be C1,Cs,C2h or D2h"));
  }
  template<> SymmetryGroup ReadJson<SymmetryGroup>(object& o, string k, int n, int m){
    
    if(o.find(k) == o.end()) {
      string msg = "key \"" + k + "\" not found in parsing SymmetryGroup";
      throw runtime_error(msg);
    }
    return ReadJson<SymmetryGroup>(o[k]);
  }
  //  template ERIMethod ReadJson<ERIMethod>(object& o, string k, int n, int m);
  
  void ReadJson_Molecule(value& v, Molecule mole) {

    try {
      CheckValue<array>(v);
    }catch(exception& e) {
      string msg = "error on parsing Molecule\n";
      msg += e.what();
      throw(runtime_error(msg));
    }

    string key;
    array& atoms = v.get<array>();
    for(array::iterator it = atoms.begin(); it != atoms.end(); ++it) {

      try {
	// -- build --
	CheckValue<object>(*it);
	object& obj_atom = it->get<object>();
	key = "atom"; string name = ReadJson<string>(obj_atom, "atom");
	key = "q";    dcomplex q = ReadJson<dcomplex>(obj_atom, "q");
	key = "xyz";  MatrixXcd xyz_mat = ReadJson<MatrixXcd>(obj_atom, "xyz", -1, 3);

	Atom atom = NewAtom(name, q);
	for(int i = 0 ; i < xyz_mat.rows(); i++) {
	  atom->Add(Vector3cd(xyz_mat(i, 0), xyz_mat(i, 1), xyz_mat(i, 2)));
	}
	mole->Add(atom);
      } catch(exception& e) {
	int i = distance(atoms.begin(), it);
	ostringstream oss;
	oss << "error on parsing value of "
	    << i << "th element of molecule" << endl;
	oss << e.what();
	throw(runtime_error(oss.str()));
      }
    }
    mole->SetSymPos();

  }
  void ReadJson_Molecule(picojson::object& obj, std::string k, Molecule mole) {
    try {
      CheckObject<array>(obj, k);
      ReadJson_Molecule(obj[k], mole);
    } catch(exception& e) {
      string msg = "error on parsing Molecule\n";
      msg += e.what();
      throw runtime_error(msg);
    }
    
  }
  Reduction ReadJson_Reduction(value& json, SymmetryGroup sym){

    try {
      CheckValue<object>(json);
      object& obj = json.get<object>();
      string str_irrep = ReadJson<string>(obj, "irrep");
      MatrixXcd coef = ReadJson<MatrixXcd>(obj, "coef");
      int irrep = sym->GetIrrep(str_irrep);
      Reduction rds(irrep, coef);
      return rds;

    } catch(exception& e) {
      string msg = "error on parsing Reduction.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }
    
  }
  void ReadJson_SymGTOs_Subs_cart( object& basis, SymGTOs gtos) {

    try {
      
      // -- pn --
      Vector3i ns(ReadJson<VectorXi>(basis, "ns", -1, 3));

      // -- atom --
      string atom_name = ReadJson<string>(basis["atom"]);
      Molecule mole = gtos->molecule();
      if(not mole->exist_atom(atom_name))  {
	throw runtime_error("atom name \"" + atom_name + "\" not found");
      }

      // -- zeta --
      VectorXcd zeta =  ReadJson<VectorXcd>(basis, "zeta");

      // -- build --
      SymmetryGroup sym = gtos->sym_group();
      gtos->NewSub(atom_name).Mono(0, ns, zeta);

    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs_cart.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }    

  }
  void ReadJson_SymGTOs_Subs_full( object& basis, SymGTOs gtos) {
    
    try {

      // -- atom --
      string atom_name = ReadJson<string>(basis, "atom");
      Molecule mole = gtos->molecule();
      Atom atom = mole->atom(atom_name);

      // -- sub --
      SymmetryGroup sym = gtos->sym_group();
      SubSymGTOs sub(sym, atom);

      // -- ns --
      MatrixXi ns = ReadJson<MatrixXi>(basis, "ns", -1, 3);
      for(int i = 0; i < ns.rows(); i++) {
	sub.AddNs(ns(i, 0), ns(i, 1), ns(i, 2));
      }
      
      // -- reduction --
      CheckObject<array>(basis, "rds");
      array& rds_list = basis["rds"].get<array>();
      for(array::iterator it = rds_list.begin(); it != rds_list.end(); ++it) {
	Reduction rds = ReadJson_Reduction(*it, sym);
	sub.AddRds(rds);
      }

      // -- zeta --
      VectorXcd zeta = ReadJson<VectorXcd>(basis, "zeta");
      sub.AddZeta(zeta);

      gtos->AddSub(sub);
      
    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs_full.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }
    
  }
  void ReadJson_SymGTOs_Subs_SolidSH( object& basis, SymGTOs gtos) {

    try {

      // -- atom --
      string atom_name = ReadJson<string>(basis, "atom");
      Molecule mole = gtos->molecule();
      if(not mole->exist_atom(atom_name))  {
	throw runtime_error("atom name \"" + atom_name + "\" not found");
      }

      // -- L --
      int L = ReadJson<int>(basis, "L");

      // -- M --
      VectorXi Ms = ReadJson<VectorXi>(basis, "Ms");
      
      // -- zeta --
      VectorXcd zeta =  ReadJson<VectorXcd>(basis, "zeta");

      // -- build --
      gtos->NewSub(atom_name).SolidSH_Ms(L, Ms, zeta);

    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs_SolidSH.\n";
      msg += e.what();
      throw(runtime_error(msg));      
    }    

  }
  void ReadJson_SymGTOs_Subs(picojson::value& json, SymGTOs gtos) {

    try {
      /*
      if(json.is<string>()) {
	string path = ReadJson<string>(json);
	fstream f; f.open(path.c_str(), ios::in);
	value new_json; f >> new_json;
	cout << "new_json: " << new_json.to_str() << endl;
	ReadJson_SymGTOs_Subs(new_json, gtos);
	return;
      } 
      */
      CheckValue<array>(json);
      array& ary = json.get<array>();
      
      for(array::iterator it = ary.begin(); it != ary.end(); ++it) {
	
	CheckValue<object>(*it);
	object& basis = it->get<object>();
	
	CheckObject<string>(basis, "type");
	string type = ReadJson<string>(basis, "type");
	if(type == "cart")
	  ReadJson_SymGTOs_Subs_cart(basis, gtos);
	else if(type == "full")
	  ReadJson_SymGTOs_Subs_full(basis, gtos);
	else if(type == "solid_sh")
	  ReadJson_SymGTOs_Subs_SolidSH(basis, gtos);
	else 
	  throw(runtime_error("unsupported type: " + type));
      }
    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs\n";
      msg += e.what();
      throw runtime_error(msg);
    }
  }
  void ReadJson_SymGTOs_Subs(picojson::object& obj, std::string k, SymGTOs gtos) {
    if(obj.find(k) == obj.end()) {
      throw runtime_error("key \"" + k + "\" not found.");
    }
    try {
      ReadJson_SymGTOs_Subs(obj[k], gtos);
    } catch(exception& e) {
      string msg = "error on parsing SymGTOs_Subs. keys = \"" + k + "\"\n";
      msg += e.what();
      throw runtime_error(msg);
    }
  }
  

  
  void ReadJson_Orbital_file(object& initstate, SymmetryGroup _sym, dcomplex *_E0,
			     VectorXcd *_c0, Irrep *_irrep0, int *_i0) {

    try {
      string irrep_str = ReadJson<string>(initstate, "irrep");
      *_irrep0 = _sym->GetIrrep(irrep_str);
      *_i0 = ReadJson<int>(initstate, "i0");
    } catch(exception& e) {
      string msg = "error on parsing Orbital_file.\n";
      msg += e.what();
      throw runtime_error(msg);
    }
    
    try {
      string eigvecs_file = ReadJson<string>(initstate, "eigvecs");
      BMat eigvecs;
      eigvecs.Read(eigvecs_file);
      *_c0 = eigvecs(*_irrep0, *_irrep0).col(*_i0);    
    } catch(exception& e) {
      string msg = "error on reading Orbital_file.\n";
      msg += e.what();
      throw runtime_error(msg);
    }

    try {
      string eigvals_file = ReadJson<string>(initstate, "eigvals");
      BVec eigvals;
      eigvals.Read(eigvals_file);
      *_E0 = eigvals(*_irrep0)(*_i0);
    } catch(exception& e) {
      string msg = "error on parsing Orbital_file.\n";
      msg += e.what();
      throw runtime_error(msg);
    }
    
  }
  void ReadJson_Orbital(object& obj, string k, SymmetryGroup _sym,
			dcomplex *_E0, VectorXcd *_c0, Irrep *_irrep0, int *_i0) {
  
    try {
      CheckObject<object>(obj, k);
      object& initstate = obj[k].get<object>();
      string type = ReadJson<string>(initstate, "type");
      
      if(type == "file")
	ReadJson_Orbital_file(initstate, _sym, _E0, _c0, _irrep0, _i0);
      else
	throw runtime_error("type <- [file]");
      
      
    } catch(exception& e) {
      string msg = "error on parsing Orbital. key = \"" + k + "\"\n";
      msg += e.what();
      throw runtime_error(msg);
    }
  }  
}
