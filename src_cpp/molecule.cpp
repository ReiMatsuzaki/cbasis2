#include <stdexcept>
#include <iomanip>
#include <boost/foreach.hpp>
#include "molecule.hpp"
#include "../utils/macros.hpp"

using namespace std;
using namespace Eigen;

namespace cbasis {  
  Atom _Atom::Add(Eigen::Vector3cd xyz) {
    xyz_list_.push_back(xyz);
    return shared_from_this();
  }
  Atom _Atom::Add(dcomplex x, dcomplex y, dcomplex z) {
    return this->Add(Vector3cd(x, y, z));
  }
  Atom _Atom::SetSymPos(SymmetryGroup sym) {

     vector<Vector3cd> new_xyz_list;
     sym->CalcSymPosList(this->xyz_list(), &new_xyz_list);
     this->xyz_list_.swap(new_xyz_list);   
     return shared_from_this();
  }
  Atom NewAtom(std::string _name, dcomplex _q) {
    return Atom(new _Atom(_name, _q));
  }
  Atom NewAtom(std::string _name, dcomplex _q, Eigen::Vector3cd _xyz) {
    Atom atom = NewAtom(_name, _q);
    atom->Add(_xyz);
    return atom;
  }
  
  _Molecule* _Molecule::Add(Atom atom) {
    
    if(atom_map_.find(atom->name()) != atom_map_.end()) {
      string msg; SUB_LOCATION(msg);
      msg = "\n" + msg + " : atom " + atom->name() + " already exist." ;
      throw runtime_error(msg);
    }
    
    atom_map_[atom->name()] = atom;
    typedef vector<Vector3cd>::const_iterator It;
    const vector<Vector3cd>& xyz_list = atom->xyz_list();
    for(It it = xyz_list.begin(); it != xyz_list.end(); ++it) {
      this->names.push_back(atom->name());
      this->qs.push_back(atom->q());
      this->xyzs.push_back(*it);
    }    
    return this;
  }
  dcomplex _Molecule::NucEnergy() {

    typedef vector<Vector3cd>::const_iterator cIt;
    
    dcomplex acc(0);
    for(Map::iterator it = atom_map_.begin(); it != atom_map_.end(); ++it) {
      Atom atomi = it->second;
      const vector<Vector3cd>& xyzi = atomi->xyz_list();
      for(cIt it_pos = xyzi.begin(); it_pos != xyzi.end(); ++it_pos) {
	cIt jt_pos = it_pos+1;
	for(; jt_pos != xyzi.end(); ++jt_pos) {
	  acc += atomi->q() * atomi->q() / (*it_pos-*jt_pos).norm();
	}
      }
      for(Map::iterator jt = it; jt != atom_map_.end(); ++jt) {
	if(jt != it) {
	  Atom atomj = jt->second;
	  const vector<Vector3cd>& xyzj = atomj->xyz_list();
	  for(cIt it_pos = xyzi.begin(); it_pos != xyzi.end(); ++it_pos)
	    for(cIt jt_pos = xyzj.begin(); jt_pos != xyzj.end(); ++jt_pos) {
	      acc += atomi->q() * atomj->q() / (*it_pos-*jt_pos).norm();
	    }
	}
      }
    }

    return acc;
  }
  _Molecule* _Molecule::SetSymPos() {

    vector<Atom> atoms;
    for(Map::iterator it = atom_map_.begin(); it != atom_map_.end(); ++it) {
      Atom atom = it->second;
      atom->SetSymPos(this->sym_group());
      atoms.push_back(atom);
    }
    this->atom_map_.clear();
    this->names.clear();
    this->qs.clear();
    this->xyzs.clear();

    BOOST_FOREACH(Atom atom, atoms) {
      this->Add(atom);
    }
    return this;
  }
  string _Molecule::str() const {
    ostringstream oss;
    oss << "==== Molecule ====" << endl;

    typedef Map::const_iterator It;
    for(It it = atom_map_.begin(); it != atom_map_.end(); ++it) {
      int i = distance(atom_map_.begin(), it);
      Atom atom = it->second;
      oss << "atom" << i << ": {";
      oss << atom->name() << ", " << atom->q() << ", [";
      const vector<Vector3cd>& xyz_list = atom->xyz_list();
      for(vector<Vector3cd>::const_iterator jt = xyz_list.begin();
	  jt != xyz_list.end(); ++jt) {
	const Vector3cd& xyz = *jt;
	oss << "[" << xyz[0] << xyz[1] << xyz[2] << "] ";
      }
      oss << "]}" << endl;
    }
    return oss.str();
  }
  string _Molecule::show() {
    string sep = " | ";
    string line = "-|-";
    ostringstream oss;

    ostringstream lines;
    lines << setfill('-')
	  << line << setw(4) << right << ""
	  << line << setw(6) << right << ""
	  << line << setw(10) << right << ""
	  << line << setw(10) << right << ""
	  << line << setw(10) << right << ""
	  << line << endl;
    oss << lines.str();
    
    oss << setfill(' ')
	<< sep << setw(4) << right << "name"
	<< sep << setw(6) << right << "q"
	<< sep << setw(10) << right << "x"
	<< sep << setw(10) << right << "y"
	<< sep << setw(10) << right << "z"
	<< sep << endl;
    
    oss << lines.str();

    BOOST_FOREACH(Map::value_type key_atom, atom_map_) {
      Atom atom = key_atom.second;
      oss << setfill(' ')
	  << sep << setw(4) << right << atom->name()	
	  << sep << setw(6) << right << atom->q()
	  << sep << setw(10) << right << ""
	  << sep << setw(10) << right << ""
	  << sep << setw(10) << right << ""
	  << sep<< endl;

      BOOST_FOREACH(Vector3cd xyz, atom->xyz_list()) {
	oss << setfill(' ')
	    << sep << setw(4) << ""
	    << sep << setw(6) << ""
	    << sep << setw(10) << right << xyz[0]
	    << sep << setw(10) << right << xyz[1]
	    << sep << setw(10) << right << xyz[2]
	    << sep << endl;
      }
    }
    oss << lines.str();
    return oss.str();
  }
  Molecule NewMolecule(SymmetryGroup _sym) {
    Molecule ptr(new _Molecule(_sym));
    return ptr;
  }
}
