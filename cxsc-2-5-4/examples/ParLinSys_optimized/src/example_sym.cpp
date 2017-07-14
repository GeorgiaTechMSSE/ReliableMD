/*
 *  Example program for the usage of the symmetric linear system solver
 */


#include <iostream>
#include <string>
#include <symlinsys.hpp>
#include <ivector.hpp>      
#include <vector>

#include <fastlss.hpp>

using namespace std;
using namespace cxsc;

/*
 *  Example for solving the symmetric linear system (by Behnke)
 * 
 * A = ( 2    [1,2])  x  =  ([10,10.5])
 *     ([1,2]    2 )        ([10,10.5])
 * 
 */ 

int main() {
	int n, i, j, k, Err;

 	cout << SetPrecision(10, 16) << Scientific << endl; //OutputFormat

	cout << "Inner and outer bounds for the solution set hull of a symmetric linear system:" << endl << endl;

	n = 2;

	imatrix A(n,n);
	ivector b(n), x(n), y(n);

	A[1][1] = A[2][2] = 3.0;
	A[1][2] = A[2][1] = interval(1,2);
	
	b = interval(10,10.5);
	
	cout << "Starting symmetric system solver.." << endl << endl;
	struct parlinsysconfig cfg;
	cfg.refinement = true;
	cfg.msg = true;
        SymLinSolve(A, b, x, y, Err, cfg);
	cout << endl;
        if(Err != 0) 
	  cout << SymLinSolveErrMsg(Err) << endl;
	else {
          cout << "Verified outer enclosure:" << endl << x << endl;
	  cout << "Verified inner enclosure:" << endl << y << endl;
	}
	
	cout << endl << endl;

	cout << "Starting non symmetric linear system solver.." << endl << endl;
	struct lssconfig cfg2;
	cfg2.msg = true;
        ilss(A, b, x, y, Err, cfg2);
        cout << endl;
	
        if(Err != 0) 
	  cout << ParLinSolveErrMsg(Err) << endl;
	else {
          cout << "Verified outer enclosure:" << endl << x << endl;
	  cout << "Verified inner enclosure:" << endl << y << endl;
	}

	return 0;
}
/*
Inner and outer bounds for the solution set hull of a symmetric linear system:

Starting symmetric system solver..

Converting symmetric system to parametric system...
Using 8 thread(s)
Preparation: 0.0110469
Computation of approximate inverse: 0.0115681
Computing result for right hand side 1...
Defect iteration: 0.00123405
Computation of [C]: 0.000402927
Computation of [z]: 0.000174999
Verification: 4.50611e-05
Refinement Step: 1.00136e-05
Inner Estimation: 0.000594139

Verified outer enclosure:
[1.6480500143949043E+000,2.9075055411606523E+000]
[1.6480498014307762E+000,2.9075057541247800E+000]

Verified inner enclosure:
[2.0679339931443485E+000,2.4876215624112086E+000]
[2.0679339694816678E+000,2.4876215860738893E+000]



Starting non symmetric linear system solver..

Using 8 thread(s) for computations
Computing approximate inverse of A

Computing result for right hand side 1
Residual iteration...
 1. iteration step
Computing [z] = R*(b - A*x)
Verifying...
 Computing [C] = I-R*A
 1. verification step
 2. verification step
Computing inner enclosure
Verified solution found

Verified outer enclosure:
[8.3333333333333103E-001,3.7222222222222246E+000]
[8.3333333333333103E-001,3.7222222222222246E+000]

Verified inner enclosure:
[1.8333333333333348E+000,2.7222222222222206E+000]
[1.8333333333333348E+000,2.7222222222222206E+000]

*/
 
