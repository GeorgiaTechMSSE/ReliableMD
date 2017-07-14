//============================================================================ 
// 
//		 Program/Module "Parametric Linear Solver" 
// 
//	    supplement to C++ Toolbox for Verified Computing 
// 
//		  Copyright (c) 2003   Evgenija Popova
// 
// This program/module is free software for non-commercial use. 
// For details on theory  see the papers 
// 
// S. Rump: Verification methods for dense and sparse systems of equations, 
// in: Topics in Validated Computations, ed. J. Herzberger, Elsevier Science 
// B.V., 1994,	pp. 63-135. 
//
// E. Popova, W. Kraemer: Parametric Fixed-point Iteration Implemented in 
// C-XSC,  Universitaet Wuppertal, Preprint BUW - WRSWT, 2003.
// 
// This program/module is distributed WITHOUT ANY WARRANTY. 
// 
//============================================================================ 
//---------------------------------------------------------------------------- 
// File: parlinsys (implementation) 
// 
// Purpose: Computing a verified inclusion for the parametric solution set 
// of a square parametric linear system of equations A(p)*x = b(p) with full 
// parametric matrix A(p) and parametric right-hand side b(p) which elements 
//             are affine-linear combinations of k given interval parameters. 
// Method: Transformation of A(p)*x = b(p) to fixed-point form and applying 
//  	    an interval residual iteration. 
// 
// Global functions: 
//    ParLinSolve()	 : to get a verified enclosure of the solution 
//    ParLinSolveErrMsg(): to get an error message text 
// 
//---------------------------------------------------------------------------- 
//---------------------------------------------------------------------------- 
// In the comments below the notations #*(...) and ##(...) are used to indi- 
// cate the evaluation of the expression specified in round brackets using 
// the exact scalar product. I.e. each component of the result is computed 
// with only one final rounding. The symbol #* holds for rounding to nearest 
// whereas ## holds for rounding to the smallest enclosing interval. An exact 
// scalar product may be implemented using dotprecision accumulators. 
//---------------------------------------------------------------------------- 

#include <string.h>         // String handling 
#include <imatrix.hpp>      // Interval matrix/vector arithmetic 
#include <mvi_util.hpp>     // Interval matrix/vector utility functions 
#include <matinv.hpp>       // Matrix inversion 
#include "parlinsys.hpp"    
 
using namespace cxsc; 
using namespace std; 
 
//---------------------------------------------------------------------------- 
// Error codes used in this module. 
//---------------------------------------------------------------------------- 
 
const int 
  NoError      = 0,   // No error occurred 
  NotSquare    = 1,   // System to be solved is not square 
  DimensionErr = 2,   // Dimensions of A and b are not compatible 
  InvFailed    = 3,   // System is probably singular 
  VerivFailed  = 4;   // Verification failed, system is probably 
                      // ill-conditioned 
 
//---------------------------------------------------------------------------- 
// Error messages depending on the error code. 
//---------------------------------------------------------------------------- 
 
char* ParLinSolveErrMsg (int Err) 
{ 
  static char Msg[80] = ""; 
 
  if (Err != NoError) { 
    char Hlp[60]; 
 
    switch (Err) { 
      case NotSquare: 
        strcpy(Hlp,"System to be solved is not square"); break; 
      case DimensionErr: 
        strcpy(Hlp,"Dimensions of A and b are not compatible"); break; 
      case InvFailed: 
        strcpy(Hlp,"System is probably singular"); break; 
     case VerivFailed: 
        strcpy(Hlp,"Verification failed, system is probably ill-conditioned"); 
        break; 
      default: 
        strcpy(Hlp,"Code not defined"); 
    } 
    sprintf(Msg,"Error: %s!",Hlp); 
  } 
  return(Msg); 
} // end  ParLinSolveErrMsg 
 
 
//---------------------------------------------------------------------------- 
// The vectors x and y are successive approximations for the solution of 
// a point linear system of equations computed by iterative refinement.  
// If a component of y is diminished by more than 'Factor', it is a good 
// candidate for a zero entry. Thus, it is set to zero. 
//---------------------------------------------------------------------------- 
 
static void CheckForZeros ( rvector& x, rvector& y ) 
{ 
  const real Factor = 1E+5; 
  int        i; 
 
  for (i = Lb(y); i <= Ub(y); i++) 
    if ( abs(y[i])*Factor < abs(x[i]) ) y[i] = 0.0; 
} 
 
//---------------------------------------------------------------------------- 
// The vectors x and y are successive iterates. The function returns true if 
// the relative error of all components x_i and y_i is <= 10^(-12), i.e. y_i 
// has about 12 correct decimals. If x_i or y_i vanishes, the relative error 
// is implicitly set to zero. 
//---------------------------------------------------------------------------- 
 
static int Accurate ( rvector& x, rvector& y ) 
{ 
  const real Delta = 1E-12;             // Relative error bound 
  int        i, n, ok; 
  real       abs_yi; 
 
  i = Lb(y); n = Ub(y); 
  do { 
    if (sign(x[i])*sign(y[i]) <= 0)     // Relative error set to 0 
      ok = 1; 
    else { 
     abs_yi = abs(y[i]);                // Relative error > Delta ? 
     ok = (abs(y[i] - x[i]) <= Delta * abs_yi ); 
    } 
    i++; 
  } while (ok && (i <= n)); 
  return ok; 
} // Accurate 
 
//---------------------------------------------------------------------------- 
// This function 'VerificationStep()' performs an iteration with  
// Einzelschrittverfahren: for k = 1,2,...  
// [y] = [x] = Blow([x],Epsilon), [x]_i = [z]_i + [C]_i*[x], i = 1,..n until 
// the new iterate [x] lies in [y] or until the maximum number of iterations  
// is exceeded. The factor 'Epsilon' for the inflation is specified by the
// user. The flag 'IsVerified' is set if an inclusion could be established. 
//---------------------------------------------------------------------------- 
 
static void VerificationStep(ivector& xx, ivector& zz, 
                              imatrix& C, real& Epsilon, int& IsVerified) 
{ 
  const int  MaxIter = 10;          // Maximal number of iteration steps 
 
  int     p, i; 
  ivector yy(Lb(xx), Ub(xx)); 
 
  xx = zz; p = 0;                   // Initialize:  [x] = [z] 
  do  
  {	   
    yy = xx = Blow(xx, Epsilon);          // Epsilon inflation 
    for(i = Lb(xx); i <= Ub(xx); i++)     // Einzelschrittferfahren  
      xx[i] = zz[i] + C[i] * xx;          // new iterate  
    IsVerified = in(xx, yy);              // Inclusion in the interior??? 
    p++; 
  }  
  while (!IsVerified && (p < MaxIter)); 

cout << "MaxIter = " << p << endl << endl;

}

// This function 'Dual' reverses the end-points of an interval.
// ------------------------------------------------------------
static void Dual(interval& xx)
{ real tmp;
    
  tmp = Inf(xx);
  UncheckedSetInf(xx, Sup(xx));
  UncheckedSetSup(xx, tmp);
} // End Dual
 
 
//---------------------------------------------------------------------------- 
// Purpose: The function 'ParLinSolve()' computes a verified enclosure for 
// the parametric solution set of a square parametric linear system of  
// equations A(p)*x=b(p), having a factorized by the parameters representation: 
//    (A0+[p1]A1+[p2]A2+ +[pk]Ak)x = (b0+[p1]b1+[p2]b2+ +[pk]bk).
// Parameters: 
//    In : 'Ap'         : coefficient matrix of the system in a factorized  
//                        representation Ap = (A0,A1,..,Ak) from R(n(k+1),n); 
//                        the matrix is passed as reference parameter to  
//                        save time for copying it; 
//         'bp'         : coefficients right-hand side of the system in a  
//                        factorized matrix representation 
//                        bp = (b0,b1,..bk) from R(n,k+1); matrix bp is passed  
//                        as reference parameter to save time for copying it; 
//         'ip'         : vector of the interval values for the parameters, 
//                        passed as reference parameter to save time for 
//                        copying it; 
//         'SharpC'     : flag signalling whether a sharp enclosure of the  
//                        contracting matrix C is to be computed; 
//         'Eps'        : constant for the epsilon inflation (e.g. 0.2);
//         'Inner'      : flag signalling whether inner estimation of the hull
//			  of the solution set is to be computed;
//    Out: 'xx'         : enclosure of the parametric solution set, resized to 
//                        standard index range with lower index bound 1; 
//	   'yy'		: inner estimation of the exact hull for the parametric
//			  solution set, resized to standard index range with 
//			  lower index bound 1;
//         'Err'        : error code; 
// Description: An approximate inverse 'R' of the mid-point matrix is computed  
//   by calling function 'MatInv()'. Then an approximate mid-point solution 'x' 
//   is computed applying a conventional real residual iteration. For the final 
//   verification, an interval residual iteration is performed. An enclosure of 
//   the unique solution is returned in the interval vector 'xx' and an inner  
//   approximation of the hull of the solution set is returned in the interval
//   vector 'yy'.
//---------------------------------------------------------------------------- 
 void ParLinSolve( rmatrix&  Ap,  
                   rmatrix&  bp, 
                   ivector&  ip,  
                       int&  SharpC, 
		      real&  Eps,
		       int&  Inner,
                   ivector&  xx, 
		   ivector&  yy,
                       int&  Err ) 
{ 
  const int MaxResCorr = 10;    // Maximal number of real residual corrections 
 
  int     p = VecLen(ip);            // Length of the parametric vector 'ip' 
				     // here it is of length p=p+1 for p parameters
  int     n = RowLen(Ap);            // Length of the rows of 'Ap' 
  int     m = ColLen(Ap);            // Length of the columns of 'Ap' 
  int     IsVerified, i, j, k; 
  int     ARowLow, AColLow,          // Lower index bounds of 'Ap', 'bp' 
          bRowLow, bColLow, ipLow;   // and ip 
   
  rmatrix R;                         // To store the inverse of 'mid(A([ip]))' 
   
   
  if (m != n*p)                      // Error: 'A(p)' is not square 
    { Err = NotSquare; return; } 
 
  if (n != ColLen(bp) &&  
         p != RowLen(bp))            // Error: Dimensions of 'A(p)' and  
    { Err = DimensionErr; return; }  //       'b(p)' are not compatible 
 
 
// Start algorithm 
//---------------- 
  rvector       x(n), y(n), d(n);  // Allocate dynamic vectors 
  imatrix       C(n, n);           // Allocate dynamic interval matrix C
  ivector       zz(n);	           // Allocate dynamic nonparametric vector 
  dotprecision  Accu;              // Allocate accumulators 
  idotprecision IAccu; 
 
  
// Normalize index range of 'Ap' and 'ip' to standard range 1..n 
// Save lower bounds of 'ip' and 'Ap' 
// to restore them before leaving 'ParLinSolve()'. 
//----------------------------------------------------------------- 
  ARowLow = Lb(Ap, ROW); SetLb(Ap, ROW, 1); 
  AColLow = Lb(Ap, COL); SetLb(Ap, COL, 1); 
  ipLow = Lb(ip); SetLb(ip, 1); 


 rmatrix    Am(n, n);               // Allocate dynamic mid-point matrix
 rvector    bm(n), pm(p);           // Allocate dynamic mid-point vectors
 
// Compute mid-point vectors/matrix: pm, Am, bm 
//---------------------------------------------------------------------- 
  pm = mid(ip); 
  bm = bp * pm;
  Am = Ap(1, n, 1, n); 
  for (i = 2; i <= p; i++)      
    Am += pm[i]* Ap(n*(i-1)+1, i*n, 1, n); 
//---------------------------------------------------------------------- 
 
  MatInv(Am, R, Err); 
  if (Err != NoError)                   // Error: Inversion failed 
    { Err = InvFailed; return; } 
  
// Normalize index range of 'R'  to standard index range 1..n 
//---------------------------------------------------------------------
 SetLb(R, ROW, 1); SetLb(R, COL, 1);
 
  x = R * bm; k = 0;   // Real residual iteration 
  do  
  {                
    y = x; 
    for (i = 1; i <= n; i++)  
    {                               // Compute: d = #*(bm-Am*y) 
      Accu = bm[i];                 //----------------------- 
      accumulate(Accu, -Am[i], y); 
      d[i] = rnd(Accu); 
    } 
    
    for (i = 1; i <= n; i++)        // Compute: x = #*(y+R*d) 
    {                               //----------------------  
      Accu = y[i];                
      accumulate(Accu, R[i], d); 
      x[i] = rnd(Accu); 
    } 
     
    CheckForZeros(y, x); 
    k++; 
  }  
  while (!Accurate(y, x) && (k < MaxResCorr)); 
 
//  Resize mid-point vectors and matrix to save memory
//  ----------------------------------------------------------------
    Resize(pm);  Resize(bm);  Resize(Am);
    Resize(y);


// Prepare verification step: 

imatrix T;	// Interval matrix for intermediate quantities
 
// compute enclosure [C]  of   I - R*A(p) for p in [ip]  
//--------------------------------------------------------------
 
  if (SharpC) 
    {                        // Sharp enclosure C = I - R*A(p) 
      Resize(T, 1, p);       // Allocate dynamic temporary matrix T
      for (i = 1; i <= n; i++) 
       for (j = 1; j <= n; j++)              
       {                                     
         Accu = (i == j) ? 1.0 : 0.0;         
         accumulate(Accu, -R[i], Ap(1, n, 1, n)[Col(j)]); 
         rnd(Accu, T[1][1]); 
         for (k = 2;k <= p; k++) 
          { 
           Accu = 0.0; 
	     accumulate(Accu, -R[i], Ap(n*(k-1)+1, n*k, 1, n)[Col(j)]); 
	     rnd(Accu, T[1][k]); 
          } 
         IAccu = 0.0; 
         accumulate(IAccu, T[Row(1)], ip); 
         C[i][j] = rnd(IAccu); 
       } 
    } 
  else 
    {                      // Rough enclosure C = I - R*A([ip])  
      Resize(T, n, n);     // Allocate dynamic temporary matrix T 
      T = Ap(1, n, 1, n);  //                  to contain A([ip])
      for (i = 2; i <= p; i++)
         T += ip[i]* Ap(n*(i-1)+1, i*n, 1, n); 

      for (i = 1; i <= n; i++) 
       for (j = 1; j <= n; j++)                
       {                                       
        IAccu = (i == j) ? 1.0 : 0.0;         
        accumulate(IAccu, -R[i], T[Col(j)]); 
        rnd(IAccu, C[i][j]); 
       } 
    };  


// Resize T to length (n,p).
// Normalize index range of 'bp' to standard range 1..n.
// Save lower bounds of 'bp'
// to restore them before leaving 'ParLinSolve()'.
//---------------------------------------------------------------------
   Resize(T, n, p);
   bRowLow = Lb(bp, ROW); SetLb(bp, ROW, 1);
   bColLow = Lb(bp, COL); SetLb(bp, COL, 1);   
 
// compute enclosure [z] of z = R*(b(p)-A(p)*x). 
//---------------------------------------------------------------------
  for (k = 1; k <= p; k++)  
  { j = n*(k-1); 
    for (i = 1; i <= n; i++)               // Compute: d = #*(b_k - A_k * x) 
    {                                      //------------------------------ 
      Accu = bp[Col(k)][i]; 
      m = j + i;                 
      accumulate(Accu, -Ap[m], x);  
      d[i] = rnd(Accu); 
    } 
 
    for (i = 1; i <= n; i++)                // Compute: [y] = ##(b_k - A_k*x - d) 
    {                                       //----------------------------------  
      Accu = bp[Col(k)][i];      
      m = j + i;                             
      accumulate(Accu, -Ap[m], x); 
      Accu -= d[i]; 
      rnd(Accu, zz[i]); 
    } 
 
                                //  Now b_k -A_k*x lies in d + [y].  
                                //  Thus R*(b_k-A_k*x) lies in [z]= ##(R*d + R*[y])     
    for (i = 1; i <= n; i++)  
    {   
      IAccu = 0.0; 
      accumulate(IAccu, R[i], d); 
      accumulate(IAccu, R[i], zz); 
      T[Col(k)][i] = rnd(IAccu); 
    } 
 
  }  // end for_k 
 
  for(i = 1; i <= n; i++) 
  { 
    IAccu = 0.0; 
    accumulate(IAccu, T[i], ip);   
    zz[i] = rnd(IAccu); 
  } 
 
 
// Resize 'xx' to length 'n'.
// ---------------------------------- 
   Resize(xx, n);  
  

// verified enclosure [x] for the absolute error of x. 
//-------------------------------------------------------------------------
   
//  VerificationStep(xx, zz, C, IsVerified);   // Attempt to compute [x] 

  VerificationStep(xx, zz, C, Eps, IsVerified); 
  if (!IsVerified) 
      Err = VerivFailed; 
  else 
  {
    if (Inner)        // Compute Inner Estimation of the outer solution
    {		      // --------------------------------------------------
      ivector cp(p);
      Resize(yy, n);  Resize(d);
    
      for (j = 1; j <= n; j++)
      {
         cp = ip;
         for (i=1; i<=p; i++)
         {
           if ( Inf(T[j][i]) >= 0 )  // T>=0
           { Dual(cp[i]);
             if ( Sup(ip[i]) <= 0 ) Dual(T[j][i]);  // ip <= 0
             else if ( Inf(ip[i]) < 0 ) UncheckedSetSup(T[j][i], Inf(T[j][i]));
           }
           else if ( Sup(T[j][i]) <= 0 ) // T<=0
           {
             if ( Sup(ip[i]) <= 0 ) Dual(T[j][i]);
             else if ( Inf(ip[i]) < 0 ) UncheckedSetInf(T[j][i], Sup(T[j][i]));
           }
              else   // 0 in T 
              { if ( Inf(ip[i]) < 0 && Sup(ip[i]) > 0 ) cp[i] = interval(0);
                else 
                  if ( Sup(ip[i]) <= 0 )
                    {  Dual(T[j][i]);
		       UncheckedSetInf(cp[i], Sup(ip[i]));
                    }
		  else UncheckedSetSup(cp[i], Inf(ip[i]));
              };
         };

         Accu = 0.0;
         accumulate(Accu, Inf(T[j]), Inf(cp));
         UncheckedSetInf(zz[j], rnd(Accu, RND_DOWN));   

         Accu = 0.0;
         accumulate(Accu, Sup(T[j]), Sup(cp));
         UncheckedSetSup(zz[j],  rnd(Accu, RND_UP));

         IAccu = x[j];
         IAccu += zz[j]; 
         accumulate(IAccu, C[j], xx);
         yy[j] = rnd(IAccu);
         Dual(yy[j]);

      };  // end for_j

    };    // end Inner Estimation of the solution enclosure. 
   xx = x + xx;             // The solution set lies in x+[x]
   };  
 
  // Restore lower index bounds of 'Ap', 'bp' and ip 
  //-------------------------------------------------   
  SetLb(bp, ROW, bRowLow); SetLb(bp, COL, bColLow); 
  SetLb(Ap, ROW, ARowLow); SetLb(Ap, COL, AColLow); 
  SetLb(ip, ipLow); 
 
} //ParLinSolve 
 
