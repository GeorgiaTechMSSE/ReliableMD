 /*---------------------------------------------------------
   Sample program itayl_ex2.cpp;
   Evaluation of the equivalent functions  
          g(x) = ln(1+e^x)     h(x) = x + ln(1+e^(-x));
   The selection of the appropriate expression is essential 
   for getting tight enclosures of the Taylor coefficients.
  ---------------------------------------------------------*/
#include "itaylor.hpp" // Header file of class itaylor
#include <iostream>    // Input | output

using namespace cxsc;
using namespace std;
using namespace taylor;

itaylor f1(const itaylor& u)
{
    itaylor w;
    interval z( get_j_coef(u,0) );
    if (Inf(z)>0) w = u + lnp1( exp(-u) );
    else w = lnp1( exp(u) );
    return w;
}

int main()
{
    int p = 4;    // Order of expansion
    itaylor f,g,h;  // Default constructor No. 1.
    itaylor x(p,interval(706)); 

    g = lnp1( exp(x) );
    h = x + lnp1( exp(-x) );
    f = f1(x);

    cout << SetPrecision(16,16) << Scientific << endl;
    cout << "Taylor coefficients of g(x) are included by:" << endl;
    print_itaylor(g);  // Output of Taylor coefficients.
    cout << "Taylor coefficients of h(x) are included by:" << endl;
    print_itaylor(h);
    cout << "Taylor coefficients of f1(x) are included by:" << endl;
    print_itaylor(f);

    return 0;
} // main


