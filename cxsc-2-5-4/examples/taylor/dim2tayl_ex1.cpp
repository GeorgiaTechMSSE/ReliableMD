/***********************************************************
*  Sample program dim2tayl_ex1.cpp to calculate enclosures *
*  of the Taylor coefficients of    f(x,y)=x^2-2xy^3    at *
*  (x,y)=(2,-3) up to the desired Taylor order p=3.        *
***********************************************************/

    #include <iostream>
    #include "dim2taylor.hpp"

    using namespace std;
    using namespace cxsc;
    using namespace taylor;

    dim2taylor f(dim2taylor_vector& x)
    {   //  f = x^2 - 2*x*y^3;
        dim2taylor erg;   // Default constructor call
        erg = sqr(x[1]) - 2*x[1]*x[2]*sqr(x[2]); // f(x,y)
        return erg;
    }

    int main()
    {
        int p = 3;      // Desired maximal order of Taylor expansion

        ivector iv(2);  // 2 interval components with Lb=1 and Ub=2;
        iv[1] = interval(2);   // x-value
        iv[2] = interval(-3);  // y-value

        dim2taylor_vector tv;  // Default constructor call
        tv = init_var(p,iv);   // Initialization with vector iv
        dim2taylor t;          // Default constructor call
        t = f(tv);             // function call

        cout << "t[1][2] = " << t[1][2] << endl; 
    } 
