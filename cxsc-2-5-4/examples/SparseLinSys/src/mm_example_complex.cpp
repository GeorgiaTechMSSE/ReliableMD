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

// Example for using matrix market files

#include <iostream>
#include <sparselinsys.hpp>
#include <scimatrix.hpp>
#include <civector.hpp>
#include <fstream>
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

scmatrix conj(const scmatrix& A) {
  scmatrix B(A);
  
  vector<complex>& val = B.values();
  
  for(int i=0 ; i<val.size() ; i++)
    SetIm(val[i], -Im(val[i]));
  
  return B;
}

int main(int argc, char **argv) {

   int m,n;        
   char c;
   string filename;
   bool spd = false;
   
   cout << endl;
   cout << "Test program for sparse solvers using Matrix Market files" << endl;
   cout << "=========================================================" << endl;
   
   cout << "File name (matrix should be complex): "; cout.flush();
   cin >> filename;
   cout << "Reading file..." << endl;
   
   scmatrix A;
   ifstream file(filename.c_str());
   if(!file.good()) {
     cout << "File not found!" << endl;
     return 1;
   }
   cout << MatrixMarketInOut;
   file >> A;
   
   m = ColLen(A);
   n = RowLen(A);

   cout << "m=" << m << ",n=" << n << ",nnz=" << A.get_nnz() << endl << endl; 

   if(A == conj(transp(A))) {
     c = ' ';
     while(c!='j' && c!='n') {
       cout << "Assume matrix is positive definite (j/n)? "; cout.flush();
       cin >> c;
     }
     spd = (c == 'j');
     cout << endl;
   }

   cout << SetPrecision(10,16) << Scientific;

   civector bb(m);
   civector xx(n);
   civector yy(n);
   int err;

   bb = 1.0;
   
   scimatrix Ai = A + A*interval(-1e-12,1e-12);
   bb = bb + bb*interval(-1e-12,1e-12);
   
     
   double start,end;


   cout << "Starting normwise solver..." << endl;
   slssnconfig cfgn;
   cfgn.msg = true;
   cfgn.spd = spd;  
   start = GetTime();
   slssn(Ai,bb,xx,err,cfgn);
   end = GetTime();
   if(err) {
     cout << "Error: " << SparseLinSolveErrMsg(err) << endl;
   } else {
     cout << "Verification successful." << endl;
     cout << "Time taken: " << end-start << endl;
     cout << "Relative error: " << RelErr(Re(xx)) << ", " << RelErr(Im(xx)) << endl;
   }
   cout << endl << endl;
   
   cout << "Starting Krawczyk-solver with forward/backward substitution..." << endl;
   slssconfig cfg;
   cfg.msg = true;              
   cfg.triMethod = FORWARD_BACKWARD;   
   cfg.spd = spd;
   start = GetTime();
   slss(A,bb,xx,err,cfg);
   end = GetTime();
   if(err) {
     cout << "Error: " << SparseLinSolveErrMsg(err) << endl;
   } else {
     cout << "Verification successful." << endl;
     cout << "Time taken: " << end-start << endl;
     cout << "Relative error: " << RelErr(Re(xx)) << ", " << RelErr(Im(xx)) << endl;
   }
   cout << endl << endl; 

   return 0;
} 
