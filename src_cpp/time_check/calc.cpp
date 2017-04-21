#include <l2func.hpp>
#include <iostream>
#include <vector>
#include "timer.hpp"

using namespace std;
using namespace l2func;

template<class Prim>
void CreateSampleBasisSet(int num, vector<Prim>* us) {

  for(int i = 0; i < num; i++) 
    us->push_back( Prim(2, complex<double>(0.1 * (i + 1)), Normalized));
 }

template<class Prim>
void CalcMatrixByBasis(const vector<Prim>& us, Op<Prim> op) {

  complex<double> acc(0.0);
  typedef typename vector<Prim>::const_iterator IT;
  for(IT i = us.begin(), end = us.end(); i != end; ++i) {
    LinearComb<Prim> op_ui = op(*i);
    for(IT j = us.begin();  j != end; ++j){
      acc += CIP(*j, op_ui);
    }
  }

}

template<class Prim>
void CalcD1BasisTimes(const vector<Prim>& us) {

  int us_size(us.size());
  vector<LinearComb<Prim> > d_us(us_size);

  for(int j = 0; j < us_size; j++) 
    d_us[j].SetD1Normalized(us[j]);

}

template<class Prim>
void CalcVector(const vector<Prim>& us, const LinearComb<CSTO>& driv) {
  
  complex<double> acc(0.0);
  typedef typename vector<Prim>::const_iterator IT;
  for(IT it = us.begin(); it != us.end(); ++it) {
    acc += CIP(*it, driv);
  }
}

int main() {
  
  HLikeAtom<complex<double> > h_atom(1, 1.0, 0);
  Op<CSTO> op = h_atom.HMinusEnergy<CSTO>(0.5);
  Op<CGTO> op_g = h_atom.HMinusEnergy<CGTO>(0.5);
  Timer timer;

  LinearComb<CSTO> driv = h_atom.DipoleInitLength(1);

  vector<CSTO> us; CreateSampleBasisSet(5, &us);
  timer.Start("op_5STO");
  CalcMatrixByBasis<CSTO>(us, op);
  timer.End("op_5STO");
  timer.Start("driv_5STO");
  CalcVector(us, driv);
  timer.End("driv_5STO");

  vector<CSTO> us50; CreateSampleBasisSet(50, &us50);
  timer.Start("op_50STO");
  CalcMatrixByBasis<CSTO>(us50, op);
  timer.End("op_50STO");
  timer.Start("driv_50STO");
  CalcVector(us50, driv);
  timer.End("driv_50STO");
  

  vector<CGTO> g5; CreateSampleBasisSet(5, &g5);
  timer.Start("op_5GTO");
  CalcMatrixByBasis<CGTO>(g5, op_g);
  timer.End("op_5GTO");
  timer.Start("driv_5GTO");
  CalcVector(g5, driv);
  timer.End("driv_5GTO");

  vector<CGTO> g50; CreateSampleBasisSet(50, &g50);
  timer.Start("op_50GTO");
  CalcMatrixByBasis<CGTO>(g50, op_g);
  timer.End("op_50GTO");
  timer.Start("driv_50GTO");
  CalcVector(g50, driv);
  timer.End("driv_50GTO");

  timer.Start("d1_basis");
  CalcD1BasisTimes(g50);
  timer.End(  "d1_basis");
  timer.Display();
}



