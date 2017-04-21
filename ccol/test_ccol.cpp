#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/assign.hpp>
#include <gtest/gtest.h>
#include "../src_cpp/bmatset.hpp"
#include "../src_cpp/symmolint.hpp"
#include "../src_cpp/symgroup.hpp"
#include "read_aoint_wrapper.hpp"

using namespace std;
using boost::format;
using namespace boost::assign;
using namespace cbasis;
using namespace Eigen;


void WriteMolint(SymGTOs gtos, string label, string filename) {
  ofstream f(filename.c_str(), ios::out);
  f << label << endl;
  f << "1.0000000000000 0.0000000000 0 0\n";

  int num_atoms = 2;
  int num_basis = 20;
  int num_sub = 5;
  f << format("  1  %d  %d %d  %d 65 65\n") % num_atoms % num_sub % num_basis %num_sub;

  SymmetryGroup sym = gtos->sym_group();
  int num_irrep = sym->order();
  f << format("  %d") % num_irrep;
  for(int i = 0; i < num_irrep; i++) {
    f << format("%6s") % sym->irrep_name_[i];
  }
  f << endl;

  int num_prod_table(0);
  for(int i = 0; i < num_irrep; i++)
    for(int j = 0; j < num_irrep; j++)
      for(int k = 0; k < num_irrep; k++)
	if(sym->prod_table_(i, j, k) && i < j)
	  num_prod_table++;
  f << format(" %d\n") % num_prod_table;
  for(int i = 0; i < num_irrep; i++)
    for(int j = 0; j < num_irrep; j++)
      for(int k = 0; k < num_irrep; k++)
	if(sym->prod_table_(i, j, k) && i<j && j<k && i!=0 && j!=0)
	  f << format(" %d  %d  %d\n") % (i+1) % (j+1) % (k+1);
  
  BOOST_FOREACH(SubSymGTOs& sub, gtos->subs()) {
    f << format(" %d") % sub.rds.size();
    BOOST_FOREACH(Reduction& rds, sub.rds) {
      f << format("  %d") % (rds.irrep + 1);
    }
    f << endl;
  }

  int idx(1);
  BOOST_FOREACH(SubSymGTOs& sub, gtos->subs()) {
    f << format(" %d  %d  %d\n") % (sub.rds.size())  % (sub.maxn*2+1) % idx;
    BOOST_FOREACH(Reduction& rds, sub.rds) {
      for(int ipn = 0; ipn < rds.size_pn(); ipn++) {
	for(int iat = 0; iat < rds.size_at(); iat++) {
	  dcomplex x = rds.coef_iat_ipn(iat, ipn);
	  f << format(" %d ") % (x.real() * abs(x));
	}
      }
      f << endl;
    }
  }

  BOOST_FOREACH(SubSymGTOs& sub, gtos->subs()) {
    for(int iz = 0; iz < sub.size_zeta(); iz++) {
      f << format("1, %d\n") % (sub.maxn+1);
      f << format("%10.6f, %10.6f, 1.000000, 1.000000\n")
	% sub.zeta(iz).real() % sub.zeta(iz).imag();
    } 
  }

  map<Atom, int> atom_num;
  BOOST_FOREACH(SubSymGTOs& sub, gtos->subs()) {
    Atom atom = sub.atom();
    if(atom_num.find(atom) == atom_num.end()) {
      atom_num[atom] = 0;
    }
    atom_num[atom] += sub.size_zeta();
  }
  
  Molecule mole = gtos->molecule();
  typedef pair<string, Atom> StrAtom;
  BOOST_FOREACH(StrAtom str_atom, mole->atom_map_) {
    string name = str_atom.first;
    Atom atom = str_atom.second;
    f << format("%s %d  %d  %2.1f\n")
      % name % (atom_num[atom]) % atom->size() % atom->q().real();
    for(int i = 0; i < atom->size(); i++) {
      Vector3cd xyz = atom->xyz_list()[i];
      f << format("%f   %f   %f\n") % xyz[0].real() % xyz[1].real() % xyz[2].real();
    }
    if(atom->size() != 1) {
      f << "  2  1\n";
    }
    int idx(1), idx_sub(0);
    BOOST_FOREACH(SubSymGTOs& sub, gtos->subs()) {
      idx_sub++;
      if(sub.atom() == atom) {
	for(int i = 0; i < sub.size_zeta(); i++) {
	  f << format("%d %d\n") % idx % idx_sub;
	  idx++;
	}
      }
    }
  }

}
TEST(WriteMolint, First) {

  SymmetryGroup sym = SymmetryGroup_D2h();

  Molecule mole = NewMolecule(sym);
  mole->Add(NewAtom("H", 1.0)->Add(0,0,+1.0)->Add(0,0,-1.0));
  mole->Add(NewAtom("Cen", 0.0)->Add(0,0,0));

  int num_zeta(10);
  VectorXcd zeta(num_zeta);
  for(int n = 0; n < num_zeta; n++)
    zeta(n) = pow(2.0, n-5);

  VectorXcd zetap(3);
  zetap[0] = dcomplex(1.1, -0.3);
  zetap[1] = dcomplex(0.3, -0.1);
  zetap[2] = dcomplex(0.1, -0.05);

  VectorXi Ms(3); Ms << -1,0,+1;
  SymGTOs gtos = NewSymGTOs(mole);
  gtos->NewSub("H")
    .AddNs(  Vector3i(0, 0, 0))
    .AddZeta(zeta)
    .AddRds( Reduction(sym->irrep_s(), MatrixXcd::Ones(2, 1)));
  gtos->NewSub("Cen")
    .SolidSH_Ms(1, Ms, zetap);
  gtos->NewSub("Cen")
    .SolidSH_Ms(3, Ms, zetap);
  gtos->SetUp();
  
  WriteMolint(gtos, "produced by cpp", "out/molint_h2plus_cpp.in");
  system("cd out/; ~/src/prog/cCOLUMBUS/molint< molint_h2plus_cpp.in");
  //system("cd out/; ~/src/prog/cCOLUMBUS/molint< molint_h2plus_cpp.in > molint_h2plus_cpp.out");
}

TEST(ReadAOINTS, First) {
  /*
  int aoint_ifile;
  char aoint_filename[20] = "out/AOINTS";
  //long aoint_filename_length = 20;
  long aoint_filename_length = strlen(aoint_filename)+1;
  bool succ;

  aoint_filename[aoint_filename_length-1] = ' ';
  open_file_binary_read_(&aoint_ifile, &succ, aoint_filename, aoint_filename_length);
  
  char blabel[80];
  double repfunc;
  int    nst;
  short int nd[8];
  
  read_header_(&aoint_ifile, blabel, &repfunc, &nst, nd);
  cout << "blabel: " << blabel << endl;
  cout << "repfunc: " << repfunc << endl;
  cout << "nst: " << nst << endl;
  cout << "nd: " << nd[0] << nd[1] << endl;
  close_file_(&aoint_ifile);
  EXPECT_TRUE(succ);
  EXPECT_EQ(2, 1+1);
  */
}
TEST(ReadAOINTS, third) {

  AoIntsHeader ao;
  BMat s, t, v;
  char filename[20] = "out/AOINTS";
  ReadAOINTS(filename, &ao, &s, &t, &v);
  cout << s(0, 0)(0, 0) << endl;
  cout << s(0, 0)(0, 1) << endl;
  cout << s(0, 0)(1, 0) << endl;
}
TEST(ReadAOINTS, Second) {
  /*
  int aoint_ifile;
  char aoint_filename[20] = "out/AOINTS";
  //long aoint_filename_length = 20;
  long aoint_filename_length = strlen(aoint_filename)+1;
  bool succ;

  aoint_filename[aoint_filename_length-1] = ' ';
  open_file_binary_read_(&aoint_ifile, &succ, aoint_filename, aoint_filename_length);

  AoIntsHeader ao;
  aoints_read_header_(&aoint_ifile, &ao);

  cout << ao.blabel << endl;
  cout << ao.repfunc << endl;
  cout << ao.nst << endl;
  cout << ao.nd[0] << endl;
  cout << ao.ityp[0] << endl;
  cout << "ns: " << ao.ns << endl;
  cout << "zscale: " << ao.zscale.re << "+" << ao.zscale.im << endl;

  // -- BMat --
  BMat smat, tmat, vmat;
  int num_isym[10];
  aoints_read_mat_structure_(&aoint_ifile, num_isym);
  cout << format("num_isym = ");
  for(int isym = 0; isym < 8; isym++) {
    int ni = num_isym[isym];
    cout << format(" %d") % ni;
    smat(isym, isym) = MatrixXcd::Zero(ni, ni);
    tmat(isym, isym) = MatrixXcd::Zero(ni, ni);
    vmat(isym, isym) = MatrixXcd::Zero(ni, ni);
  }
  cout << endl;
  close_file_(&aoint_ifile);

  // re read
  open_file_binary_read_(&aoint_ifile, &succ, aoint_filename, aoint_filename_length);
  aoints_read_header_(&aoint_ifile, &ao);
  int is[1080];
  int js[1080];
  int isyms[1080];
  for_complex vs[1080];
  int num;
  bool is_end_block = false;
  for(bool is_end_block = false; not is_end_block; ) {
    aoints_read_mat_value_block_(&aoint_ifile, &num, &is_end_block, vs, is, js, isyms);
    cout << "num : " << num << endl;
    cout << "is end : " << (is_end_block ? "yes" : "no") << endl;
    for(int i = 0; i < num; i++) {
      smat(isyms[i]-1, isyms[i]-1)(is[i]-1, js[i]-1) = dcomplex(vs[i].re, vs[i].im);
    }
  }
  for(bool is_end_block = false; not is_end_block; ) {
    aoints_read_mat_value_block_(&aoint_ifile, &num, &is_end_block, vs, is, js, isyms);
    for(int i = 0; i < num; i++) {
      tmat(isyms[i]-1, isyms[i]-1)(is[i]-1, js[i]-1) = dcomplex(vs[i].re, vs[i].im);
    }
  }
  for(bool is_end_block = false; not is_end_block; ) {
    aoints_read_mat_value_block_(&aoint_ifile, &num, &is_end_block, vs, is, js, isyms);
    for(int i = 0; i < num; i++) {
      vmat(isyms[i]-1, isyms[i]-1)(is[i]-1, js[i]-1) = dcomplex(vs[i].re, vs[i].im);
    }
  }    
  close_file_(&aoint_ifile);
  */

  
}



int main(int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
