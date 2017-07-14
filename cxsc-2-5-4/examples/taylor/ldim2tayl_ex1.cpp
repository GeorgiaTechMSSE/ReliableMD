/***************************************************************
*  Sample program ldim2tayl_ex1.cpp to calculate an enclosure  *
*  of the Taylor coefficient      f[2][1]     of the function  *
*  f(x,y) = sqrt( 1+x^2/y)  at (x,y)=(4,2) (Staggered Arithm.) *
***************************************************************/

    #include <iostream>
    #include "ldim2taylor.hpp"

    using namespace std;
    using namespace cxsc;
    using namespace taylor;

    ldim2taylor f(ldim2taylor_vector& x)
    {   //  f = sqrt( 1+x^2/y);
        ldim2taylor erg;
        erg = sqrt(l_real(1.0) + sqr(x[1])/x[2] ); // f(x,y)
        return erg;
    }

    int main()
    {
        int p = 3;     // Maximal order of Taylor expansion
	stagprec = 3;  // Desired precision of 3*16=48 digits

        l_ivector iv(2);  // 2 components with Lb=1 and Ub=2;
        iv[1] = l_interval(4);  // x-value,   x_0 = 4;
        iv[2] = l_interval(2);  // y-value,   y_0 = 2;

        ldim2taylor_vector tv;  // Default constructor call
        tv = init_var(p,iv);    // Initialization with vector iv
        ldim2taylor t;          // Default constructor call
        t = f(tv);              // function call
	cout << SetDotPrecision(16*stagprec,16*stagprec-3) 
             << Scientific;
        cout << "t[2][1] = " << t[2][1] << endl; 
    } 
