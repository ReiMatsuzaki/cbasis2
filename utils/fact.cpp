#include <string>
#include <stdexcept>
#include "fact.hpp"

namespace cbasis {
  
  int iabs(int a) {
    return a > 0 ? a : -a;
  }

  int ipow(int a, int n) {

    if(a == 1) {
      return 1;
    }

    if(a == -1) {
      if(n % 2 == 0)
	return 1;
      else
	return -1;
    }
    
    if(n < 0) {
      throw std::runtime_error("n must be non negative integer");
    }

    if(n == 0) {
      return 1;
    }

    int acc = 1;
    for (int i = 0; i < n; i++) {
      acc *= a;
    }
    
    return acc;
    
  }

  int Factorial(int n) {

    if(n < 0) {
      throw std::runtime_error("n must be zero or positive in Factorial.");
    }

    int acc = 1;
    for(int i = n; i > 0; i--)
      acc *= i;

    return acc;
  }
  double DFactorial(int n) {
    return Factorial(n) * 1.0;
  }
  int DoubleFactorial(int n) {

    if(n < 0) {
      throw std::invalid_argument("n must be zero or positive in DoubleFactorial.");
    }


    int acc = 1;
    for(int i = n; i > 0; i-=2) {
      acc *= i;
    }
    
    return acc;
  }
  double DDoubleFactorial(int n) {
    return DoubleFactorial(n) * 1.0;
  }
  int Combination(int n, int k) {
    return Factorial(n) / (Factorial(k) * Factorial(n-k));
  }
}
