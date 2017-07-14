//============================================================================
//
// 	      ParLinSys - A Solver for Parametric Linear Systems
//
//	    Supplement to the C-XSC Toolbox for Verified Computing
//
// Author: Michael Zimmer, based on a previous version by Evgenija Popova
//
// This program is free software for non-commercial use.
//
// For details on theory see the papers:
//
// S. Rump: Verification methods for dense and sparse systems of equations,
// in: Topics in Validated Computations, ed. J. Herzberger, Elsevier Science
// B.V., 1994,	pp. 63-135.
//
// E. Popova, W. Kraemer: Parametric Fixed-Point Iteration Implemented in
// C-XSC. Universitaet Wuppertal, Preprint BUW - WRSWT, 2003.
//
// M. Zimmer, W. Kraemer, E.Popova: Solvers for the verified solution of 
// parametric linear systems, Computing: Volume 94, Issue 2 (2012), Page 109-123,
// Springer.
//
// For versions using  the alternative algorithm:
// A. Neumaier, A. Pownuk: Linear systems wit hlarge uncertainties, with 
// applications to truss structures. Reliable Computing, 13(149-172), 2007.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================

/*
 * This file contains the implementation of the parametric linear solver 
 * based on the Krawczyk operator using methods by Rump and Popova
 */

#include <string>         
#include <simatrix.hpp>          
#include <matinv_aprx.hpp>       // Matrix inversion 
#include <parlinsys.hpp>
#include <timer.hpp>
#include <stdio.h>
#include <fenv.h>

#ifdef _OPENMP
#include <omp.h>
#endif
 
using namespace cxsc; 
using namespace std; 
 
namespace cxsc {
  
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
// string ParLinSolveErrMsg (int Err) 
// { 
//   string Msg;
//   if (Err != NoError) { 
//     switch (Err) { 
//       case NotSquare: 
//         Msg = "System to be solved is not square"; break; 
//       case DimensionErr: 
//         Msg = "Dimensions of A and b are not compatible"; break; 
//       case InvFailed: 
//         Msg = "System is probably singular"; break; 
//      case VerivFailed: 
//         Msg = "Verification failed, system is probably ill-conditioned"; 
//         break; 
//       default: 
//         Msg = "Code not defined"; 
//     } 
//   } 
//   return Msg; 
// } // end  ParLinSolveErrMsg 

 
// Tries to guess zero values in successive iterates
/* The vectors x and y are successive approximations for the solution of a
   linear system of equations computed by iterative refinement. If a compo-
   nent of y is diminished by more than 'Factor', it is a good candidate for
   a zero entry. Thus, it is set to zero.
   \param x An iterate
   \param y Iterate succeding x
*/
template<typename Tvec>
inline void CheckForZeros ( Tvec& x, Tvec& y )
{
  const real Factor = 1E+5;

  for (int i = Lb(y); i <= Ub(y); i++)
    if ( abs(y[i])*Factor < abs(x[i]) ) y[i] = 0.0;
}


// Checks if an iterate is sufficiently accurate
/*
   The vectors x and y are successive iterates. The function returns true if
   the relative error of all components x_i and y_i is <= 10^(-12), i.e. y_i
   has about 12 correct decimals. If x_i or y_i vanishes, the relative error
   is implicitly set to zero.
   \param x An iterate
   \param y Iterate succeding x
   \return true, if y hast about 12 correct decimals, else false
*/
template<typename Tvec>
inline int Accurate ( const Tvec& x, const Tvec& y )
{
  const real Delta = 1E-12;   // Relative error bound
  int        i, n, ok;

  i = Lb(y); n = Ub(y);
  do {
    // Relative error > Delta ?
    ok = (abs(y[i] - x[i]) <= Delta * abs(y[i]) );
    i++;
  } while (ok && (i <= n));

  return ok;
} // Accurate


//C-XSC Blow function for civector (epsilon inflation)
inline civector Blow(const civector &x, const real &eps) {
  civector y(Lb(x), Ub(x));
  
  for(int i=Lb(x) ; i<=Ub(x) ; i++)
    y[i] = Blow(x[i], eps);
  
  return y;
}
 
//---------------------------------------------------------------------------- 
// This function 'VerificationStep()' performs an iteration with  
// Einzelschrittverfahren: for k = 1,2,...  
// [y] = [x] = Blow([x],Epsilon), [x]_i = [z]_i + [C]_i*[x], i = 1,..n until 
// the new iterate [x] lies in [y] or until the maximum number of iterations  
// is exceeded. The factor 'Epsilon' for the inflation is specified by the
// user. The flag 'IsVerified' is set if an inclusion could be established. 
//---------------------------------------------------------------------------- 

template<class TVector, class TMatrix>
static void VerificationStep( TVector& xx, TVector& zz, TMatrix& C, 
			      int& IsVerified, int MaxIter, real Epsilon) 
{ 
  int     p, i; 
  TVector yy(Lb(xx), Ub(xx)); 
 
  xx = zz; p = 0;                   // Initialize:  [x] = [z] 
  do  
  {	 
    yy = xx = Blow(xx, Epsilon);          // Epsilon inflation 
    for(i = Lb(xx); i <= Ub(xx); i++) {     // Einzelschrittferfahren  
      xx[i] = zz[i] + C[i] * xx;       
    }
    IsVerified = in(xx, yy);              // Inclusion in the interior??? 
    p++; 
  }  
  while (!IsVerified && (p < MaxIter)); 

}


//----------------------------------------------------------------------------
//'maxDist()' computes the maximal distance between two 
// real/complex interval vectors.
// For interval x, y; dist(x,y)= Max(Abs(x.inf-y.inf), Abs(x.sup-y.sup))
//----------------------------------------------------------------------------
template<class TVec>
static real maxDist(const TVec& xx, const TVec& yy) {
  real dist = 0.0, tmp;
  int n = VecLen(xx);

  for(int i=1 ; i<=n ; i++) {
    tmp = abs(Inf(xx[i]) - Inf(yy[i]));
    if(tmp > dist) dist = tmp;
    tmp = abs(Sup(xx[i]) - Sup(yy[i]));
    if(tmp > dist) dist = tmp;
  }

  return dist;
}

//----------------------------------------------------------------------------
//'RefinementStep()' performs iterative refinement of a solution enclosure [x]: 
// p= 1,2,... [y]= [x], [x]= [z] + [C].[x], [x] &= [y], until 
// the distance between the two iterates [x] and [y] is > Epsilon 
// or until the maximum number of iterations is exceeded. 
// The factor 'Epsilon' for the distance between two iterates is user specified. 
// [z] and [C] are the same as in the Verification step.
//----------------------------------------------------------------------------
template<class TVec, class TMat> 
static void RefinementStep(TVec& xx, const TVec& zz, const TMat& C, int MaxIter, const real& Epsilon) { 
  int        p=0; 
  TVec yy(Lb(xx), Ub(xx)); 
 
  do {	   
    yy = xx;
    xx = zz + C*xx;
    xx &= yy;
    p++;
  } while (maxDist(xx, yy) > Epsilon && (p < MaxIter)); 

}

//Checks if Inf(c) <= Sup(c) for an interval c
static bool check(const cinterval& c) {
  return InfRe(c) <= SupRe(c) && InfIm(c) <= SupIm(c);
}

//Checks if Inf(c) <= Sup(c) for an interval c
static bool check(const interval& i) {
  return Inf(i) <= Sup(i);
}

//Reverse interval end points (Kaucher intervals)
void Dual(interval& i) {
  real tmp = Inf(i);
  UncheckedSetInf(i, Sup(i));
  UncheckedSetSup(i, tmp);
}
    
//Compute an inner estimation for a previously computed verified outer enclosure    
static void InnerEstimation(ivector& yy, ivector& zz, const imatrix& C, const imatrix& T, const ivector& xx, const ivector& ip, const rvector& x ) {
  int p = VecLen(ip);
  int n = VecLen(x);
  ivector cp(p);
  Resize(yy,n);
  dotprecision Accu;
  idotprecision IAccu;
 
  #ifdef _OPENMP
  #pragma omp parallel for private(cp,Accu,IAccu)
  #endif
  for (int j = 1; j <= n; j++) {

    cp = ip;

    for (int i=1; i<=p; i++)  {

      if ( Inf(T[j][i]) >= 0 )  { 

        // T>=0
        Dual(cp[i]);
        if ( Sup(ip[i]) <= 0 ) Dual(T[j][i]);  // ip <= 0
        else if ( Inf(ip[i]) < 0 ) UncheckedSetSup(T[j][i], Inf(T[j][i]));

      } else if ( Sup(T[j][i]) <= 0 )  {

        // T<=0
        if ( Sup(ip[i]) <= 0 ) Dual(T[j][i]);
        else if ( Inf(ip[i]) < 0 ) UncheckedSetInf(T[j][i], Sup(T[j][i]));

      } else { 

        // 0 in T 
        if ( Inf(ip[i]) < 0 && Sup(ip[i]) > 0 ) cp[i] = interval(0);
        else if ( Sup(ip[i]) <= 0 ) {  
          Dual(T[j][i]);
	  UncheckedSetInf(cp[i], Sup(ip[i]));
        } else {
          UncheckedSetSup(cp[i], Inf(ip[i]));
        }
      }
   }

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
   
   if(Inf(yy[j]) > Sup(yy[j])) yy[j] = interval(SignalingNaN);
  }
}


//Compute an inner estimation for a previously computed verified outer enclosure (complex version)
static void InnerEstimation(civector& yy, civector& zz, const cimatrix& C, const cimatrix& T, const civector& xx, const civector& ip, const cvector& x ) {
  int p = VecLen(ip);
  int n = VecLen(x);
  ivector cp(p),t(p);
  Resize(yy,n);
  dotprecision AccuInf, AccuSup;
  cidotprecision IAccu;

  #ifdef _OPENMP
  #pragma omp parallel for private(cp,t,AccuInf,AccuSup,IAccu)
  #endif
  for (int j = 1; j <= n; j++) {

    cp = Re(ip);
    t  = Re(T[j]);

    for (int i=1; i<=p; i++)  {

      if ( Inf(t[i]) >= 0 )  { 

        // T>=0
        Dual(cp[i]);
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);  // ip <= 0
        else if ( Inf(cp[i]) < 0 ) UncheckedSetSup(t[i], Inf(t[i]));

      } else if ( Sup(t[i]) <= 0 )  {

        // T<=0
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);
        else if ( Inf(cp[i]) < 0 ) UncheckedSetInf(t[i], Sup(t[i]));

      } else { 

        // 0 in T 
        if ( Inf(cp[i]) < 0 && Sup(cp[i]) > 0 ) cp[i] = interval(0);
        else if ( Sup(cp[i]) <= 0 ) {  
          Dual(t[i]);
	  UncheckedSetInf(cp[i], Sup(cp[i]));
        } else {
          UncheckedSetSup(cp[i], Inf(cp[i]));
        }
      }
   }

   AccuInf = 0.0;
   accumulate(AccuInf, Inf(t), Inf(cp));

   AccuSup = 0.0;
   accumulate(AccuSup, Sup(t), Sup(cp));

   cp = -Im(ip);
   t  = Im(T[j]);

    for (int i=1; i<=p; i++)  {

      if ( Inf(t[i]) >= 0 )  { 

        // T>=0
        Dual(cp[i]);
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);  // ip <= 0
        else if ( Inf(cp[i]) < 0 ) UncheckedSetSup(t[i], Inf(t[i]));

      } else if ( Sup(t[i]) <= 0 )  {

        // T<=0
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);
        else if ( Inf(cp[i]) < 0 ) UncheckedSetInf(t[i], Sup(t[i]));

      } else { 

        // 0 in T 
        if ( Inf(cp[i]) < 0 && Sup(cp[i]) > 0 ) cp[i] = interval(0);
        else if ( Sup(cp[i]) <= 0 ) {  
          Dual(t[i]);
	  UncheckedSetInf(cp[i], Sup(cp[i]));
        } else {
          UncheckedSetSup(cp[i], Inf(cp[i]));
        }
      }
   }
   
   accumulate(AccuInf, Inf(t), Inf(cp));
   UncheckedSetInf(Re(zz[j]), rnd(AccuInf, RND_DOWN));   

   accumulate(AccuSup, Sup(t), Sup(cp));
   UncheckedSetSup(Re(zz[j]), rnd(AccuSup, RND_UP));

   
    cp = Re(ip);
    t  = Im(T[j]);

    for (int i=1; i<=p; i++)  {

      if ( Inf(t[i]) >= 0 )  { 

        // T>=0
        Dual(cp[i]);
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);  // ip <= 0
        else if ( Inf(cp[i]) < 0 ) UncheckedSetSup(t[i], Inf(t[i]));

      } else if ( Sup(t[i]) <= 0 )  {

        // T<=0
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);
        else if ( Inf(cp[i]) < 0 ) UncheckedSetInf(t[i], Sup(t[i]));

      } else { 

        // 0 in T 
        if ( Inf(cp[i]) < 0 && Sup(cp[i]) > 0 ) cp[i] = interval(0);
        else if ( Sup(cp[i]) <= 0 ) {  
          Dual(t[i]);
	  UncheckedSetInf(cp[i], Sup(cp[i]));
        } else {
          UncheckedSetSup(cp[i], Inf(cp[i]));
        }
      }
   }

   AccuInf = 0.0;
   accumulate(AccuInf, Inf(t), Inf(cp));

   AccuSup = 0.0;
   accumulate(AccuSup, Sup(t), Sup(cp));

   cp = Im(ip);
   t  = Re(T[j]);

    for (int i=1; i<=p; i++)  {

      if ( Inf(t[i]) >= 0 )  { 

        // T>=0
        Dual(cp[i]);
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);  // ip <= 0
        else if ( Inf(cp[i]) < 0 ) UncheckedSetSup(t[i], Inf(t[i]));

      } else if ( Sup(t[i]) <= 0 )  {

        // T<=0
        if ( Sup(cp[i]) <= 0 ) Dual(t[i]);
        else if ( Inf(cp[i]) < 0 ) UncheckedSetInf(t[i], Sup(t[i]));

      } else { 

        // 0 in T 
        if ( Inf(cp[i]) < 0 && Sup(cp[i]) > 0 ) cp[i] = interval(0);
        else if ( Sup(cp[i]) <= 0 ) {  
          Dual(t[i]);
	  UncheckedSetInf(cp[i], Sup(cp[i]));
        } else {
          UncheckedSetSup(cp[i], Inf(cp[i]));
        }
      }
   }
   
   accumulate(AccuInf, Inf(t), Inf(cp));
   UncheckedSetInf(Im(zz[j]), rnd(AccuInf, RND_DOWN));   

   accumulate(AccuSup, Sup(t), Sup(cp));
   UncheckedSetSup(Im(zz[j]), rnd(AccuSup, RND_UP));
   
   
   IAccu = x[j];
   IAccu += zz[j]; 
   accumulate(IAccu, C[j], xx);
   yy[j] = rnd(IAccu);
   Dual(Re(yy[j]));
   Dual(Im(yy[j]));
   
   if(InfRe(yy[j]) > SupRe(yy[j])) SetRe(yy[j], interval(SignalingNaN));
   if(InfIm(yy[j]) > SupIm(yy[j])) SetIm(yy[j], interval(SignalingNaN));   
  }
}


//---------------------------------------------------------------------------- 
// Purpose: The function 'ParLinSolveMain()' computes a verified enclosure for 
// the parametric solution set of a square real or complex parametric linear system of  
// equations A(p)*x=b(p), having a factorized by the parameters representation: 
//    (A0+[p1]A1+[p2]A2+ +[pk]Ak)x = (b0+[p1]b1+[p2]b2+ +[pk]bk).
// Parameters: 
//    In : 'Ap'         : coefficient matrix of the system as a vector  
//                        of sparse or dense coefficient matrices (one for each parameter)
//         'bp'         : coefficients right-hand side of the system in a  
//                        factorized matrix representation 
//                        bp = (b0,b1,..bk) from R(n,k+1)
//         'ip'         : vector of the interval values for the parameters
//         'Inner'      : Also compute inner estimation of result
//         'cfg'        : Struct with configuration options
//    Out: 'xx'         : enclosure of the parametric solution set, resized to 
//                        standard index range with lower index bound 1; 
//         'yy'         : Inner estimation of result (only computed if Inner==true)
//         'Err'        : error code; 
// Description: An approximate inverse 'R' of the mid-point matrix is computed  
//   by calling function 'MatInv()'. Then an approximate mid-point solution 'x' 
//   is computed applying a conventional real residual iteration. For the final 
//   verification, an interval residual iteration is performed. An enclosure of 
//   the unique solution is returned in the interval vector 'xx'.
//----------------------------------------------------------------------------
template<typename TCoeff, typename TICoeff, typename Tb, typename Tp, typename Tx, 
         typename TInner, typename TR, typename TC, typename TVec, typename TIVec,
	 typename TDot, typename TIDot> 
 void ParLinSolveMain( vector<TCoeff>&  Ap, Tb&  bp, Tp&  ip,  
                       Tx&  xx, bool Inner, TInner& yy, int&  Err, 
		       struct parlinsysconfig cfg ) 
{ 
  double start=GetTime();
  double end;
 
  int     p = VecLen(ip);              // Length of the parametric vector 'ip' 
				       // here it is of length p=p+1 for p parameters
  int     n = RowLen(Ap[0]);           // Length of the rows of 'Ap' 
  int     m = ColLen(Ap[0])*Ap.size(); // Length of the columns of 'Ap' 
  int     rhs = RowLen(xx);            // Number of right hand sides
  int     IsVerified, i, j, k; 
//  vector<int> ARowLow, AColLow;          // Lower index bounds of 'Ap', 'bp' 
  int     bRowLow, bColLow, ipLow;   // and ip 
  int     oldprec = opdotprec;
  bool    C_computed = false;
  
  Err = NoError;
   
  TR R(n,n);                         // To store the inverse of 'mid(A([ip]))' 
      
  if (m != n*p)                      // Error: 'A(p)' is not square 
    { Err = NotSquare; return; } 

  if (n != ColLen(bp) || rhs*p != RowLen(bp))            // Error: Dimensions of 'A(p)' and  
    { Err = DimensionErr; return; }      // 'b(p)' are not compatible 
 
 
// Start algorithm 
//---------------- 
  TVec        x(n), y(n), d(n);  
  TC          C(n, n);           
  TIVec       zz(n);	         
  TDot        Accu;              
  TIDot       IAccu; 

  
#ifdef _OPENMP
  if(cfg.threads > 0) {
    omp_set_num_threads(cfg.threads);
  } else {
    cfg.threads = omp_get_max_threads();
    omp_set_num_threads(cfg.threads);
  }
    
  if(cfg.msg) cout << "Using " << cfg.threads << " thread(s)" << endl;
#endif

// Normalize index range of 'bp' and 'ip' to standard range 1..n 
// Save lower bounds of 'ip' and 'bp' 
// to restore them before leaving 'ParLinSolve()'. 
//----------------------------------------------------------------- 
/*  for(int i=0 ; i<p ; i++) {
    ARowLow.push_back(Lb(Ap[i], ROW));
    SetLb(Ap[i], ROW, 1); 
    AColLow.push_back(Lb(Ap[i], COL)); 
    SetLb(Ap[i], COL, 1); 
  }*/
  ipLow = Lb(ip); SetLb(ip, 1); 
  bRowLow = Lb(bp, ROW); SetLb(bp, ROW, 1);
  bColLow = Lb(bp, COL); SetLb(bp, COL, 1);   

  TR      Am(n, n);               // Allocate dynamic mid-point matrix
  TVec    bm(n), pm(p);           // Allocate dynamic mid-point vectors

  opdotprec = cfg.K;

// Compute mid-point vectors/matrix: pm, Am 
//---------------------------------------------------------------------- 
  pm = Inf(ip)+0.5*(Sup(ip)-Inf(ip)); 
  
#ifdef _OPENMP
  Am = 0.0;
  #pragma omp parallel private(i) default(shared)
  {
    TR MyAm(n,n);
    if(omp_get_thread_num() == 0) {
      if(Ap[0].isFull())
        MyAm = Ap[0].full();
      else
        MyAm = Ap[0].sparse();
    } else {
      MyAm = 0.0;
    }
    
    #pragma omp for
    for (i = 2; i <= p; i++)      
      MyAm += pm[i] * Ap[i-1]; 
    
    #pragma omp critical
    {
      Am += MyAm;
    }
  }
#else
  if(Ap[0].isFull())
    Am = Ap[0].full();
  else
    Am = Ap[0].sparse();
  for (i = 2; i <= p; i++)      
    Am += pm[i] * Ap[i-1]; 
#endif
    
    
  C = 0.0;
  
  end=GetTime();
  if(cfg.msg) cout << "Preparation: " << end-start << endl;

// Compute approximate inverse
//---------------------------------------------------------------------- 
  start = GetTime();
  MatInvAprx(Am, R, Err); 
  if (Err != NoError)                   // Error: Inversion failed 
    { Err = InvFailed; return; } 
  end = GetTime();
  if(cfg.msg) cout << "Computation of approximate inverse: " << end-start << endl;

// Normalize index range of 'R'  to standard index range 1..n 
//---------------------------------------------------------------------
  SetLb(R, ROW, 1); SetLb(R, COL, 1);


  for(int s=1 ; s<=rhs ; s++) {
  
  if(cfg.msg) cout << "Computing result for right hand side " << s << "..." << endl;

  Accu.set_k(cfg.K);
  IAccu.set_k(cfg.K);
  opdotprec = cfg.K;
  
// Compute mid-point vectors/matrix: bm 
//---------------------------------------------------------------------- 
  bm = bp(1,n,(s-1)*p+1,s*p) * pm;

  start = GetTime();
  #pragma omp parallel for firstprivate(Accu) private(i) default(shared)
  for (i = 1; i <= n; i++)  
  {        
    Accu = 0.0;
    accumulate(Accu, R[i], bm); 
    x[i] = rnd(Accu); 
  } 
  //x = R * bm; 
  k = 0;   // Real residual iteration 
  do  
  {                
    y = x; 
    #pragma omp parallel for firstprivate(Accu) private(i) default(shared)
    for (i = 1; i <= n; i++)  
    {                               // Compute: d = #*(bm-Am*y) 
      Accu = bm[i];                 //----------------------- 
      accumulate(Accu, -Am[i], y); 
      d[i] = rnd(Accu); 
    } 
    
    #pragma omp parallel for firstprivate(Accu) private(i) default(shared)
    for (i = 1; i <= n; i++)        // Compute: x = #*(y+R*d) 
    {                               //----------------------  
      Accu = y[i];                
      accumulate(Accu, R[i], d); 
      x[i] = rnd(Accu); 
    } 
     
    CheckForZeros(y, x); 
    k++; 
  }  
  while (!Accurate(y, x) && (k < cfg.maxIterResCorr)); 
 
  Resize(y);

  end = GetTime();
  if(cfg.msg) cout << "Defect iteration: " << end-start << endl;

  TC T(n,n);

  
  if(!C_computed) {
    // Prepare verification step: 
    start = GetTime();
    T = 0.0;

    //compute enclosure [C]  of   I - R*A(p) for p in [ip]  
    //--------------------------------------------------------------
    srmatrix I = Id(Ap[0]); 

    opdotprec = 1;

    if(cfg.SharpC) {
      #ifdef _OPENMP
        bool ompnested = omp_get_nested();
        omp_set_nested(false);
	
        #pragma omp parallel private(j,k) default(shared)
        {
          int size = n/omp_get_num_threads();

          #pragma omp for 
          for(j=1 ; j<=omp_get_num_threads() ; j++) {
            for(k=1 ; k<=p ; k++) {
              T((j-1)*size+1,j*size,1,n) -= ( R((j-1)*size+1, j*size, 1, n) * TICoeff(Ap[k-1]) ) * ip[k];
            } 
          }
	
          //Compute rest if matrix not evenly distributable among threads
          #pragma omp single
          {
            if(n % omp_get_num_threads() != 0) {
              for(k=1 ; k<=p ; k++) {
                T(size*omp_get_num_threads()+1,n,1,n) -= ( R(size*omp_get_num_threads()+1, n, 1, n) * TICoeff(Ap[k-1]) ) * ip[k];
              }
            }
          }

        }
        
        omp_set_nested(ompnested);
  
      #else 
 
        for(k=1 ; k<=p ; k++) {
          T -= (R * TICoeff(Ap[k-1])) * ip[k];
        }
   
      #endif

      C = I + T;

    } else {
      
      #ifdef _OPENMP
        T = 0.0;
       #pragma omp parallel private(i) default(shared)
       {
         TC MyT(n,n);
         if(omp_get_thread_num() == 0) {
           if(Ap[0].isFull())
             MyT = Ap[0].full(); 
         else 
             MyT = Ap[0].sparse();
         } else {
           MyT = 0.0;
         }
    
         #pragma omp for
         for (i = 2; i <= p; i++)
           MyT += ip[i]* Ap[i-1]; 

         #pragma omp critical
	 {
	   T += MyT;
	 }
       }

      #else
 
      if(Ap[0].isFull())
        T = Ap[0].full(); 
      else 
        T = Ap[0].sparse();
    
      for (i = 2; i <= p; i++)
         T += ip[i]* Ap[i-1]; 
      
      #endif

      #ifdef CXSC_USE_BLAS

        C = I - R*T;

      #else

        IAccu.set_k(1);
        #ifdef _OPENMP
        #pragma omp parallel for firstprivate(IAccu) private(i,j)  default(shared)
        #endif
        for (i = 1; i <= n; i++) 
         for (j = 1; j <= n; j++) {                                       
           IAccu = (i == j) ? 1.0 : 0.0;         
           accumulate(IAccu, -R[i], T[Col(j)]); 
           rnd(IAccu, C[i][j]); 
         } 
        IAccu.set_k(cfg.K);

      #endif
    }

    opdotprec = cfg.K;

    end=GetTime();
    if(cfg.msg) cout << "Computation of [C]: " << end-start << endl;

    C_computed = true;
  } //if(!C_computed)


// Resize T to length (n,p).
// Normalize index range of 'bp' to standard range 1..n.
// Save lower bounds of 'bp'
// to restore them before leaving 'ParLinSolve()'.
//---------------------------------------------------------------------
   start = GetTime();
   Resize(T, n, p);

// compute enclosure [z] of z = R*(b(p)-A(p)*x). 
//---------------------------------------------------------------------
  xx[Col(s)] = x;
  #ifdef _OPENMP
  #pragma omp parallel for private(k) default(shared)
  #endif
  for (k = 1; k <= p; k++)  {
    T[Col(k)] = bp[Col((s-1)*p+k)] - Ap[k-1] * xx[Col(s)];
  } 

  opdotprec = 1;
  T = R * T;
  zz = T * ip;

  end=GetTime();
  if(cfg.msg) cout << "Computation of [z]: " << end-start << endl;

// verified enclosure [x] for the absolute error of x. 
//-------------------------------------------------------------------------
  start = GetTime();
  TIVec xxtmp(xx[Col(s)]);
  VerificationStep(xxtmp, zz, C, IsVerified, cfg.maxIterVer, cfg.epsVer);   // Attempt to compute [x] 
  xx[Col(s)] = xxtmp;
  end=GetTime();
  if(cfg.msg) cout << "Verification: " << end-start << endl;

  if (!IsVerified) 
    Err = VerivFailed; 
  else {
    if(cfg.refinement) {
      start = GetTime();
      TIVec xxtmp(xx[Col(s)]);
      RefinementStep(xxtmp, zz, C, cfg.maxIterRef, cfg.epsRef);
      xx[Col(s)] = xxtmp;
      end=GetTime();
      if(cfg.msg) cout << "Refinement Step: " << end-start << endl;
    }

    opdotprec = cfg.K;
   
    if (Inner)        // Compute Inner Estimation of the outer solution
    {		      // --------------------------------------------------
      start = GetTime();
      TIVec yytmp(n);
      InnerEstimation(yytmp, zz, C, T, xx[Col(s)], ip, x);
      yy[Col(s)] = yytmp;
      opdotprec = cfg.K;
      end = GetTime();
      if(cfg.msg) cout << "Inner Estimation: " << end-start << endl;

    }    // end Inner Estimation of the solution enclosure. 

    xx[Col(s)] = x + xx[Col(s)];             // The solution set lies in x+[x]    
  }
  
  } //for... right hand sides
  
  // Restore lower index bounds of 'Ap', 'bp' and ip 
  //-------------------------------------------------   
  SetLb(bp, ROW, bRowLow); SetLb(bp, COL, bColLow); 
/*  for(int i=0 ; i<p ; i++) {
    SetLb(Ap[i], ROW, ARowLow[i]);
    SetLb(Ap[i], COL, AColLow[i]); 
  }*/
  SetLb(ip, ipLow); 
 
  opdotprec = oldprec;

} //ParLinSolve 



//The following functions are called by the user to start the solver and provide different combinations of data types for the system
//For alle functions the parameters have the following meaning:
//A  : vector of coefficient matrices
//bp : Right hand side(s)
//ip : Interval parameters
//xx : Solution (outer enclosure)
//yy : Inner enclosure
//err: An error code
//cfg: Struct with configuration options

void  ParLinSolve(std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::ivector& ip, cxsc::ivector& xx, cxsc::ivector& yy, int& err, struct parlinsysconfig cfg) {
  imatrix XX(VecLen(xx),1), YY(VecLen(yy),1);
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,XX,true,YY,err,cfg);
  xx = XX[Col(1)];
  yy = YY[Col(1)];
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::ivector& ip, cxsc::ivector& xx, int& err, struct parlinsysconfig cfg) {
  imatrix XX(VecLen(xx),1);
  imatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,XX,false,yy,err,cfg);
  xx = XX[Col(1)];
}

void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>& A, cxsc::scmatrix& bp, cxsc::civector& ip, cxsc::civector&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix yy;
  cimatrix XX(VecLen(xx),1);
  ParLinSolveMain<CoeffCMatrix,CoeffCIMatrix,scmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,XX,false,yy,err,cfg);
  xx = XX[Col(1)];
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::civector& ip, cxsc::civector&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix yy;
  cimatrix XX(VecLen(xx),1);
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,XX,false,yy,err,cfg);
  xx = XX[Col(1)];
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::civector& ip, cxsc::civector&xx, cxsc::civector& yy, int& err, struct parlinsysconfig cfg) {
  cimatrix XX(VecLen(xx),1), YY(VecLen(yy),1);
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,XX,true,YY,err,cfg);
  xx = XX[Col(1)];
  yy = YY[Col(1)];
}


void  ParLinSolve(std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::ivector& ip, cxsc::imatrix& xx, cxsc::imatrix& yy, int& err, struct parlinsysconfig cfg) {
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,xx,true,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::ivector& ip, cxsc::imatrix& xx, int& err, struct parlinsysconfig cfg) {
  imatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,xx,false,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>& A, cxsc::scmatrix& bp, cxsc::civector& ip, cxsc::cimatrix&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix yy;
  ParLinSolveMain<CoeffCMatrix,CoeffCIMatrix,scmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,xx,false,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::civector& ip, cxsc::cimatrix& xx, int& err, struct parlinsysconfig cfg) {
  cimatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,xx,false,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::srmatrix& bp, cxsc::civector& ip, cxsc::cimatrix&xx, cxsc::cimatrix& yy, int& err, struct parlinsysconfig cfg) {
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,srmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,xx,true,yy,err,cfg);
}



void  ParLinSolve(std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::ivector& ip, cxsc::ivector& xx, cxsc::ivector& yy, int& err, struct parlinsysconfig cfg) {
  imatrix XX(VecLen(xx),1), YY(VecLen(yy),1);  
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,XX,true,YY,err,cfg);
  xx = XX[Col(1)];
  yy = YY[Col(1)];
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::ivector& ip, cxsc::ivector& xx, int& err, struct parlinsysconfig cfg) {
  imatrix XX(VecLen(xx),1);  
  imatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,XX,false,yy,err,cfg);
  xx = XX[Col(1)];
}

void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>& A, cxsc::cmatrix& bp, cxsc::civector& ip, cxsc::civector&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix XX(VecLen(xx),1);
  cimatrix yy;
  ParLinSolveMain<CoeffCMatrix,CoeffCIMatrix,cmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,XX,false,yy,err,cfg);
  xx = XX[Col(1)]; 
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::civector& ip, cxsc::civector&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix XX(VecLen(xx),1);  
  cimatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,XX,false,yy,err,cfg);
  xx = XX[Col(1)]; 
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::civector& ip, cxsc::civector&xx, cxsc::civector& yy, int& err, struct parlinsysconfig cfg) {
  cimatrix XX(VecLen(xx),1), YY(VecLen(yy),1);  
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,XX,true,YY,err,cfg);
  xx = XX[Col(1)];
  yy = YY[Col(1)];
}


void  ParLinSolve(std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::ivector& ip, cxsc::imatrix& xx, cxsc::imatrix& yy, int& err, struct parlinsysconfig cfg) {
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,xx,true,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::ivector& ip, cxsc::imatrix& xx, int& err, struct parlinsysconfig cfg) {
  imatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,ivector,imatrix,imatrix,rmatrix,imatrix,rvector,ivector,dotprecision,idotprecision>
                 (A,bp,ip,xx,false,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>& A, cxsc::cmatrix& bp, cxsc::civector& ip, cxsc::cimatrix&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix yy;
  ParLinSolveMain<CoeffCMatrix,CoeffCIMatrix,cmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,xx,false,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::civector& ip, cxsc::cimatrix&xx, int& err, struct parlinsysconfig cfg) {
  cimatrix yy;
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,xx,false,yy,err,cfg);
}

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& A, cxsc::rmatrix& bp, cxsc::civector& ip, cxsc::cimatrix& xx, cxsc::cimatrix& yy, int& err, struct parlinsysconfig cfg) {
  ParLinSolveMain<CoeffMatrix,CoeffIMatrix,rmatrix,civector,cimatrix,cimatrix,cmatrix,cimatrix,cvector,civector,cdotprecision,cidotprecision>
                 (A,bp,ip,xx,true,yy,err,cfg);
}

} //namespace cxsc
