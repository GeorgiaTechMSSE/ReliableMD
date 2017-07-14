/*---------------------------------------------------------
   Sample program citayl_ex1.cpp.  
   Calculating Taylor coefficients and derivatives up to 
   order p=10. Expansion point x.     
---------------------------------------------------------*/

#include "citaylor.hpp"  // Header file of class citaylor
#include <cinterval.hpp> // Complex interval arithmetic
#include <iostream>      // Input

using namespace cxsc;
using namespace std;
using namespace taylor;

int main()
{
    int p=10;    // order of expansion
 
    citaylor f;  // Default constructor for f(x) = e^(-x^2).
    cinterval c; // To include the point of expansion.
    cout << "Including the point of expansion z = x + i*y" << endl; 
    cout << "c = ([x,x],[y,y]) = ? "; cin >> c;
    citaylor x(p,c);      // Constructor call

    f = exp(-sqr(x));     // function f(x)

    cout << SetPrecision(15,15) << Scientific << endl;
    print_citaylor(f); // Output of Taylor coefficients

    civector derivative(0,p); // 11 components, indices: 0,1,...,10;
    for(int i=0;i<=p; i++) derivative[i] = get_j_derive(f,i); 
    cout << "f(x): Derivatives up to order 10 " << endl;
    for(int i=0; i<=p; i++) // Output derivatives
	cout  << i <<"-th. derivative: " << derivative[i] << endl;
    return 0;
} // main


