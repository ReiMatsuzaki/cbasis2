#ifndef MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Eigen/Core>
#include "../utils/typedef.hpp"
#include "symgroup.hpp"

namespace cbasis {

  class _Atom;
  typedef boost::shared_ptr<_Atom> Atom;
  Atom NewAtom(std::string _name, dcomplex _q);
  Atom NewAtom(std::string _name, dcomplex _q, Eigen::Vector3cd _xyz);
  
  class _Atom : public boost::enable_shared_from_this<_Atom> {
  private:
    std::string name_;
    dcomplex q_;
    std::vector<Eigen::Vector3cd> xyz_list_;
  public:
    _Atom(std::string _name, dcomplex _q): name_(_name), q_(_q) {}
    int size() const { return xyz_list_.size(); }
    std::string& name() { return name_; }
    const std::string& name() const { return name_; }
    dcomplex& q() { return q_; }
    const dcomplex& q() const { return q_; }
    Atom Add(Eigen::Vector3cd xyz);
    Atom Add(dcomplex x, dcomplex y, dcomplex z);
    Atom SetSymPos(SymmetryGroup sym);
    const std::vector<Eigen::Vector3cd>& xyz_list() const { return xyz_list_; }    
  };    

  
  class _Molecule;
  typedef boost::shared_ptr<_Molecule> Molecule;
  Molecule NewMolecule(SymmetryGroup _sym);
  
  class _Molecule {
  public:
    typedef std::map<std::string, Atom> Map;
    SymmetryGroup sym_;
    Map atom_map_;

    // -- data for calculation --
    std::vector<std::string> names;
    std::vector<dcomplex> qs;
    std::vector<Eigen::Vector3cd> xyzs;
  public:
    _Molecule(SymmetryGroup _sym): sym_(_sym) {}
    Eigen::Vector3cd& At(int i) { return xyzs[i]; }
    const Eigen::Vector3cd& At(int i) const { return xyzs[i]; }
    const std::string& name(int i) const { return names[i]; }
    dcomplex x(int i) const { return this->At(i)[0]; }
    dcomplex y(int i) const { return this->At(i)[1]; }
    dcomplex z(int i) const { return this->At(i)[2]; }
    dcomplex& q(int i) { return this->qs[i]; }
    const dcomplex& q(int i) const { return this->qs[i]; }
    int size() const { return names.size(); }
    SymmetryGroup sym_group() { return sym_; }
    _Molecule* Add(Atom atom);
    Atom atom(std::string name) { return atom_map_[name]; }
    bool exist_atom(std::string name) {return atom_map_.find(name)!=atom_map_.end(); }
    dcomplex NucEnergy();
    _Molecule *SetSymPos();
    std::string str() const;
    std::string show();
  }; 
}

#endif
