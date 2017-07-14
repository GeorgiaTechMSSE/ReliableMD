/************************************************************
*  Sample program dim2tayl_ex2.cpp to calculate enclosures  *
*  of the Taylor coefficient     f[2][1]   of the function  *
*  f(x,y) = sqrt(1+(x+y)^2)    at (x,y) = (10^8, 2.1*10^8)  *
************************************************************/

    #include <iostream>
    #include "dim2taylor.hpp"

    using namespace std;
    using namespace cxsc;
    using namespace taylor;

    dim2taylor f(dim2taylor_vector& x)
    {   
        //  f(x,y) = sqrt( 1+sqr(x+y) );
        dim2taylor erg;    // Default constructor call
        erg = sqrt( 1 + sqr(x[1]+x[2]) ); // f(x,y)
        return erg;
    }

    dim2taylor g(dim2taylor_vector& x)
    {   
        //  g(x,y) = sqrt( 1+sqr(x+y) );
        dim2taylor erg;    // Default constructor call
        erg = sqrt1px2(x[1]+x[2]); // g(x,y)
        return erg;
    }

    int main()
    {
        int p = 3;      // Desired maximal order of Taylor expansion
        //char* string1 = "[1e8,1e8]";
        //char* string2 = "[2.1e8,2.1e8]";
        string string1 = "[1e8,1e8]"; // FIX
        string string2 = "[2.1e8,2.1e8]"; // FIX
        interval z1,z2;
        string1 >> z1;   string2 >> z2;
        ivector iv(2);  // 2 interval components with Lb=1 and Ub=2;
        iv[1] = z1;     // inclusion of x-coordinate
        iv[2] = z2;     // inclusion of y-coordinate

        dim2taylor_vector tv;
        tv = init_var(p,iv);   // Initialization with vector iv

        dim2taylor t;          // Default constructor call
        cout << SetPrecision(16,15) << Scientific;
        t = g(tv);             // function call
        cout << "Taylor coefficient t[2][1] with g(x,y):" << endl;
        cout << t[2][1] << endl;

        t = f(tv);             // function call
        cout << "Taylor coefficient t[2][1] with f(x,y):" << endl;
        cout << t[2][1] << endl;
    } 
