#ifndef READ_AOINT_WRAPPER_H
#define READ_AOINT_WRAPPER_H

#include <string>
#include "../src_cpp/bmatset.hpp"

namespace cbasis {
  
  extern "C" {
  typedef struct {
    double re;
    double im;
  } for_complex;
  struct AoIntsHeader {
    char blabel[80];
    char ityp[8][4];
    char mtyp[10][4];
    double repfunc;
    int nst, ns, isfr;  
    short int nd[8], nso[8], ms[142], mnl[142], kstar[142];
    for_complex zscale;
  };
  
  

  void aoints_read_header_(int *, AoIntsHeader*);
  void aoints_read_mat_structure_(int *, int *);
  void aoints_read_mat_value_block_(int *, int *, bool *,
				    for_complex *, int *, int *, int *);

  void  open_file_binary_read_(int *, bool *, char *, long);
  void  close_file_(int *);
  } 

  void ReadAOINTS(char *filename, AoIntsHeader *header,
		  BMat *smat, BMat *tmat, BMat *vmat);  
}

#endif
