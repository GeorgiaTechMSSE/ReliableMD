//============================================================================
//
// 	   SparseLinSys - Verified Solvers for Sparse Linear Systems
//
//	    Supplement to the C-XSC Toolbox for Verified Computing
//
// Author: Michael Zimmer
//
// This program is free software for non-commercial use.
//
// For details on theory see the papers:
//
// S. Rump: Verification methods for dense and sparse systems of equations,
// in: Topics in Validated Computations, ed. J. Herzberger, Elsevier Science
// B.V., 1994,	pp. 63-135.
//
// S. Rump: Validated Solution of Large Linear Systems. In: Albrecht, R., G. 
// Alefeld and H. Stetter (editors): Validated numerics: theory and 
// applications, part 9 of Computing Supplementum, pp. 191-212, Springer, 1993.
//
// S. Rump: Verification Methods: Rigorous results using floating-point arithmetic.
// Acta Numerica, 19:287-449, 2010.
//
// Kraemer, W., U. Kulisch and R. Lohner: Numerical Toolbox for Verified Computing II:
// Advanced Numerical Problems. Springer Verlag, 2009.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================

// Simple example for using the solvers

#include <iostream>
#include <sparselinsys.hpp>
#include <simatrix.hpp>
#include <timer.hpp>

using namespace std;
using namespace cxsc;

real RelErr(const ivector& x) {
    real ret(0.0);
    
    for(int i=1 ; i<=VecLen(x) ; i++) {
        real m = abs(mid(x[i]));
        real r = 0.5*diam(x[i]);
        if(m < 1.0) m=1.0;
        ret += r/m;
    }
    
    return ret/VecLen(x);
}


int main(int argc, char **argv) {

   int m,n;        

   //Dimension can be set via program argument
   if(argc == 2) {
     m = n = atoi(argv[1]);
   } else {
     m = n = 100000;
   }
   
   cout << "Simple example for using the verified sparse solvers" << endl;
   cout << "====================================================" << endl;

   //Generate a band matrix with upper,lower bandwidth 2 and bands 1,2,4,2,1
   simatrix A(m,n);
   imatrix B(m,5);
   B[Col(1)]=1;
   B[Col(2)]=-2;
   B[Col(3)]=4;
   B[Col(4)]=-2;
   B[Col(5)]=1;
   SetLb(B,COL,-2);
   A = simatrix(m,n,B);
   
   cout << "m=" << m << ",n=" << n << ",nnz=" << A.get_nnz() << endl; 

   cout << SetPrecision(10,16) << Scientific;

   ivector b(m);
   ivector xx(n);
   ivector yy(n);
   int err;

   //Right hand side is all ones
   b = 1.0;
  
   //Blow up A and b by factor 1e-12
   A = A + A * interval(-1e-12,1e-12);
   b = b + b * interval(-1e-12,1e-12);

   double start,end;

   //Use solver based on verified normwise bound
   cout << endl << "Solving system using solver based on verified normwise error bound" << endl;   
   slssnconfig cfgn;
   cfgn.msg = true;
   cfgn.spd = true; 
   start = GetTime();
   slssn(A,b,xx,err,cfgn);
   end = GetTime();
   if(err) {
     cout << "Error: " << SparseLinSolveErrMsg(err) << endl;
   } else {
     cout << "Verification successful" << endl;
     cout << "Time taken: " << end-start << endl;
     cout << "Relative error: " << RelErr(xx) << endl;
     cout << "x[1]=" << xx[1] << " , " << "x[n]=" << xx[n] << endl;
   }
   
   cout << endl << endl;
   
   //Use solver based on Krawczyk-operator
   cout << endl << "Solving system using solver based on the sparse Krawczyk operator and simple forward/backward substitution" << endl;
   slssconfig cfg;
   cfg.msg = true;           
   cfg.triMethod = FORWARD_BACKWARD;
   cfg.spd = true;
   start = GetTime();
   slss(A,b,xx,err,cfg);
   end = GetTime();
   if(err) {
     cout << "Error: " << SparseLinSolveErrMsg(err) << endl;
   } else {
     cout << "Verification successful" << endl;
     cout << "Time taken: " << end-start << endl;
     cout << "Relative error: " << RelErr(xx) << endl;
     cout << "x[1]=" << xx[1] << " , " << "x[n]=" << xx[n] << endl;
   }
   
   cout << endl << endl;
   
   //Use solver based on Krawczyk-operator
   cout << endl << "Solving system using solver based on the sparse Krawczyk operator and forward/backward substitution with parallel-epipeds" << endl;
   cfg.triMethod = BANDED;
   cfg.storeBase = true;
   start = GetTime();
   slss(A,b,xx,err,cfg);
   end = GetTime();
   if(err) {
     cout << "Error: " << SparseLinSolveErrMsg(err) << endl;
   } else {
     cout << "Verification successful" << endl;
     cout << "Time taken: " << end-start << endl;
     cout << "Relative error: " << RelErr(xx) << endl;
     cout << "x[1]=" << xx[1] << " , " << "x[n]=" << xx[n] << endl;
   }
   
   cout << endl << endl;

   return 0;
} 
