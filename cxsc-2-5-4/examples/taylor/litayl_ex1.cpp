 /*---------------------------------------------------------------
   Sample program litayl_ex1.cpp for staggered Taylor arithmetic.
      Calculating guaranteed inclusions of the derivatives of the
      function  f = x / (1+x^2) in interval staggered arithmetic.
  ---------------------------------------------------------------*/

#include "litaylor.hpp"   // Header file for class l_itaylor
#include <l_interval.hpp> // Interval staggered arithmetic
#include <iostream>       // Input, output

using namespace std;
using namespace taylor;
using namespace cxsc;

int main()
{

    stagprec = 3;  // Provides a precision of about 3*16 
                   // decimal digits.
    int p = 80;    // Taylor expansion of order p

    l_itaylor x(p,succ(1.0));  // Constructor call
                               // Point of expansion: succ(1.0)
    l_itaylor f;
    f = x / (1+sqr(x));

    cout << SetDotPrecision(16*stagprec,16*stagprec) << Scientific;
    cout << "1.  derivative: " << get_j_derive(f,1)  << endl;
    cout << "80. derivative: " << get_j_derive(f,80) << endl;
    
    return 0;
}

