/*---------------------------------------------------------------
  Sample program litayl_ex2.cpp for staggered Taylor arithmetic.
  Calculating guaranteed inclusions of the derivatives of the
  function  f = atan(x^2/(1+x^2)) in interval staggered arithmetic.
-----------------------------------------------------------------*/

#include "litaylor.hpp"   // Header file for class l_itaylor
#include <l_interval.hpp> // Interval staggered arithmetic
#include <iostream>       // Input, output

using namespace std;
using namespace taylor;
using namespace cxsc;

l_itaylor f1(const l_itaylor& x)
{
    l_itaylor w;
    l_real S( abs( Sup(get_j_derive(x,0)) ) );
    if (S<1) {
	w = sqr(x);
	atan( w/(1+w) );
    }
    else 
	w = atan(1/(1+sqr(1/x)));  // to avoid overflow
    return w;
}

int main()
{
    stagprec = 3;  // Provides a precision of about 3*16=48 
                   // decimal digits.
    int p = 2;     // Taylor expansion of order p

    l_itaylor x(p,comp(0.5,134));   // Constructor call for  
                       // independentvariable x of order p=2;
    l_itaylor f,g,h;   // Default constructor call

    f = f1(x); // After assignment: f is Taylor object of order p.
    cout << SetDotPrecision(16*stagprec,16*stagprec) << Scientific
         << "Derivatives of function f1(x):" << endl;
    cout << "1. derivative: " << get_j_derive(f,1)  << endl;
    cout << "2. derivative: " << get_j_derive(f,2)  << endl;

    cout << "Derivatives of function g(x):" << endl;
    g = sqr(x); g = atan(g/(1+g)); // suitable for |(x)_0| < 1;
    cout << "1. derivative: " << get_j_derive(g,1)  << endl;
    cout << "2. derivative: " << get_j_derive(g,2)  << endl;

    cout << "Derivatives of function h(x):" << endl;
    h = atan(1/(1+sqr(1/x))); // suitable for |(x)_0| >= 1;
    cout << "1. derivative: " << get_j_derive(h,1)  << endl;
    cout << "2. derivative: " << get_j_derive(h,2)  << endl;
    
    return 0;
}

