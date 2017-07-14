// Sample program itayl_ex1.cpp;

#include "itaylor.hpp" // Header file of class itaylor
#include <iostream>    // Input | output

using namespace cxsc;
using namespace std;
using namespace taylor;

int main()
{
 /*------------------------------------------------------------
   Inclusions of the Taylor coefficients and derivatives up to
   order p=5 for      P(x) = 2x^4 + x^3 + 4x^2 - 3x +2      at 
   points of expansion x0 included by the interval z.
  ------------------------------------------------------------*/

    int p = 5;    // Order of expansion
    interval z;   // Interval to include the point of expansion
    itaylor P;    // Default constructor No. 1.

    while(1) {
	cout << endl << "Inclusion of x0; [x0,x0] = ? ";
	cin >> z;
	itaylor x(p,z);  // Constructor No. 4, 

	P = (((real(2.0)*x + 1) // Polynomial of order 4.
	      *x + 4)
	      *x - real(3.0))
	      *x + real(2.0);

	cout << SetPrecision(16,16) << Scientific << endl;
	print_itaylor(P);  // Output of Taylor coefficients.

	ivector derivative(0,p); // interval vector with 
        // components from index 0 up to the order p.

	for(int i=0;i<=p; i++) 
	    derivative[i] = get_j_derive(P,i); // Derivatives
  
	cout << "Inclusions of the derivatives up to order 5 " 
             << endl;
	for(int i=0; i<=p; i++) // Output of the derivatives
	    cout  << i<<"th derivative:  " << derivative[i] 
		  << endl;
    }
    return 0;
} // main

/*Ausgabe 

Berechnung der ersten zehn Ableitungen der Funktion
f(x)= exp(x)*x+sin(x)+x an der Stelle x=1

Ausgabe itaylor der Ordnung 10 
i  0  component: [  4.559752,  4.559753]
i  1  component: [  6.976865,  6.976866]
i  2  component: [  3.656687,  3.656688]
i  3  component: [  1.722137,  1.722138]
i  4  component: [  0.601370,  0.601371]
i  5  component: [  0.140416,  0.140417]
i  6  component: [  0.025259,  0.025260]
i  7  component: [  0.004207,  0.004208]
i  8  component: [  0.000627,  0.000628]
i  9  component: [7.639748E-005,7.639749E-005]
i  10  component: [8.008054E-006,8.008055E-006]

Einschliessungen der Ableitungen bis Ordnung 10 
0-te  Ableitung:  [  4.559752,  4.559753]
1-te  Ableitung:  [  6.976865,  6.976866]
2-te  Ableitung:  [  7.313374,  7.313375]
3-te  Ableitung:  [ 10.332825, 10.332826]
4-te  Ableitung:  [ 14.432880, 14.432881]
5-te  Ableitung:  [ 16.849993, 16.849994]
6-te  Ableitung:  [ 18.186501, 18.186502]
7-te  Ableitung:  [ 21.205952, 21.205953]
8-te  Ableitung:  [ 25.306007, 25.306008]
9-te  Ableitung:  [ 27.723120, 27.723121]
10-te  Ableitung:  [ 29.059629, 29.059630]


*/
