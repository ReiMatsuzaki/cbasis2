#include <iostream>
#include <boost/format.hpp>
#include "../utils/macros.hpp"
#include "read_aoint_wrapper.hpp"

using boost::format;
using namespace std;
using namespace Eigen;

namespace cbasis {

  void ReadAOINTS(char *filename, AoIntsHeader *header,
		  BMat *smat, BMat *tmat, BMat *vmat) {

    // -- file name --
    long filename_length = strlen(filename) + 1;    
    filename[filename_length-1] = ' ';

    // -- open file 1st time --
    int ifile;
    bool succ;
    open_file_binary_read_(&ifile, &succ, filename, filename_length);
    if(not succ) {
      cerr << "open failed\n";
      exit(1);
    }

    // -- header 1st time --
    aoints_read_header_(&ifile, header);

    // -- structure --
    int num_isym[10];
    aoints_read_mat_structure_(&ifile, num_isym);
    for(int isym = 0; isym < 8; isym++) {
      int ni = num_isym[isym];
      //      cout << ni << endl;
      (*smat)(isym, isym) = MatrixXcd::Zero(ni, ni);
      (*tmat)(isym, isym) = MatrixXcd::Zero(ni, ni);
      (*vmat)(isym, isym) = MatrixXcd::Zero(ni, ni);
    }    
    
    // -- close --
    close_file_(&ifile);

    // -- 2nd time --
    open_file_binary_read_(&ifile, &succ, filename, filename_length);
    aoints_read_header_(&ifile, header);
    //    cout << "zscale: " << header->zscale.re <<endl;

    // -- matrix --
    int is[1080], js[1080], isyms[1080];
    for_complex vs[1080];
    int num;
    for(bool is_end = false; not is_end; ) {
      aoints_read_mat_value_block_(&ifile, &num, &is_end, vs, is, js, isyms);
      for(int k = 0; k < num; k++) {
	dcomplex v(vs[k].re, vs[k].im);
	int isym = isyms[k]-1;
	MatrixXcd& S = (*smat)(isym, isym);
	int i = is[k]-1;
	int j = js[k]-1;
	S(i, j) = v;
	if(i!=j)
	  S(j, i) = v;
      }
    }
    for(bool is_end = false; not is_end; ) {
      aoints_read_mat_value_block_(&ifile, &num, &is_end, vs, is, js, isyms);
      for(int k = 0; k < num; k++) {
	dcomplex v(vs[k].re, vs[k].im);
	int isym = isyms[k]-1;	
	int i = is[k]-1;
	int j = js[k]-1;
	MatrixXcd& T = (*tmat)(isym, isym);
	T(i, j) = v;
	if(j != i)
	  T(j, i) = v;
      }
    }
    for(bool is_end = false; not is_end; ) {
      aoints_read_mat_value_block_(&ifile, &num, &is_end, vs, is, js, isyms);
      for(int k = 0; k < num; k++) {
	dcomplex v(vs[k].re, vs[k].im);
	int isym = isyms[k]-1;	
	int i = is[k]-1;
	int j = js[k]-1;
	MatrixXcd& V = (*vmat)(isym, isym);
	V(i, j) = v;
	if(i != j)
	  V(j, i) = v;
      }
    }
    close_file_(&ifile);
    
  }
}


