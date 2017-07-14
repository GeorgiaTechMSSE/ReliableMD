/***************************************************************
*  Sample program ldim2tayl_ex2.cpp to calculate an enclosure  *
*  of the Taylor coefficient      f[2][1]     of the function  *
*  f(x,y) = sqrt( 1+(x+y)^2)    at    (x,y) = (10^8,2.1*10^8)  *
***************************************************************/

    #include <iostream>
    #include "ldim2taylor.hpp"

    using namespace std;
    using namespace cxsc;
    using namespace taylor;

    ldim2taylor f(ldim2taylor_vector& x)
    {   //  f = sqrt( 1+(x+y)^2 );
        ldim2taylor erg;
	erg = sqrt( 1+sqr(x[1]+x[2]) );
        return erg;
    }

    ldim2taylor g(ldim2taylor_vector& x)
    {   //  g = sqrt( 1+(x+y)^2 );
        ldim2taylor erg;
	erg = sqrt1px2( x[1]+x[2] );
        return erg;
    }

    int main()
    {
        int p = 3;     // Maximal order of Taylor expansion
	stagprec = 3;  // Desired precision of 3*16=48 digits

	string string1 = "[1e8,1e8]";
	string string2 = "[2.1e8,2.1e8]";

	l_interval z1,z2;
	string1 >> z1;   string2 >> z2;  

        l_ivector iv(2);  // 2 components with Lb=1 and Ub=2;
        iv[1] = z1;  // x-value,   x_0 = 4;
        iv[2] = z2;  // y-value,   y_0 = 2;

        ldim2taylor_vector tv;  // Default constructor call
        tv = init_var(p,iv);    // Initialization with vector iv
        ldim2taylor t;          // Default constructor call
        t = f(tv);              // function call
	cout << SetDotPrecision(16*stagprec,16*stagprec-3) 
             << Scientific;
        cout << "Function f(x,y): t[2][1] = " << t[2][1] << endl; 
        t = g(tv);   
        cout << "Function g(x,y): t[2][1] = " << t[2][1] << endl; 
    } 
