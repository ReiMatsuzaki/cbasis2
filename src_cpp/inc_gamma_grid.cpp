//#include <stdio.h>
#include <iostream>
#include <fstream>
#include "molint.hpp"
using namespace std;
using namespace cbasis;

int main ()
{
  /*
  int m_max = 24;
  double hx = 0.1;
  double hy = 0.1;
  double x_max = 20.0;
  double y_max = 40.0;
  */
  int m_max = 3;
  double hx = 0.1;
  double hy = 0.1;
  double x_max = 0.4;
  double y_max = 0.4;

  ofstream fout;
  fout.open("inc_gamma_grid.dat", ios::out | ios::binary | ios::trunc);
  
  for(double x = 0.0; x < x_max; x += hx)
    for(double y = 0.0; y < y_max; y += hy) {
      dcomplex* fs = IncompleteGamma(m_max, dcomplex(x, y));
      for(int m = 0; m <= m_max; m++) {
	//cout << m << " " << x << " " << y << " " << real(fs[m]) << imag(fs[m]) << endl;
	//printf("%d %20.15f %20.15f %20.15f %20.15f \n", m, x, y, real(fs[m]), imag(fs[m]));
	printf("%20.15f %20.15f \n", real(fs[m]), imag(fs[m]));
	double rpart = real(fs[m]);
	double ipart = imag(fs[m]);
	fout.write((char*)&rpart, sizeof(double));
	fout.write((char*)&ipart, sizeof(double));
      }
    }

  fout.close();

  return 0;
}
