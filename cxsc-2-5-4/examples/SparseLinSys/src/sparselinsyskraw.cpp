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

/* Implementation of the sparse linear system solver slss based on the
 * Krawczyk-Operator. */

#include <scimatrix.hpp>   
#include <sparselinsys.hpp>
#include "cxsc_umfpack.hpp"
#include "cxsc_cholmod.hpp"
#include "trisolve.hpp"
#include "utility.hpp"
#include <timer.hpp>

#include <string>         
#include <fenv.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace cxsc;
using namespace std;

namespace cxsc {

//! Error codes used in this module.
static const int
  NoError      = 0,   // No error occurred
  NotSquare    = 1,   // System to be solved is not square
  DimensionErr = 2,   // Dimensions of A and b are not compatible
  DecompFailed = 3,   // System is probably singular
  VerivFailed  = 4,   // Verification failed, system is probably
                      // ill-conditioned
  NotPosDef    = 5;   // Matrix not positive definite

//! Error messages depending on the error code.
string SparseLinSolveErrMsg ( int Err )
{
  string msg;

  if (Err != NoError) {
    switch (Err) {
      case NotSquare:
        msg = "System to be solved is not square"; break;
      case DimensionErr:
        msg = "Dimensions of A and b are not compatible"; break;
      case DecompFailed:
        msg = "Matrix decomposition failed, system is probably badly conditioned or singular"; break;
      case VerivFailed:
        msg = "Verification failed, system is probably ill-conditioned";
        break;
      case NotPosDef:
	msg = "System matrix appears not to be symmetric/hermitian positive definite";
	break;
      default:
        msg = "Code not defined";
    }
  }
  return msg;
} // LinSolveErrMsg


//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
srmatrix* midpoint(const srmatrix& A) {
  return const_cast<srmatrix*>(&A);
}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
scmatrix* midpoint(const scmatrix& A) {
  return const_cast<scmatrix*>(&A);
}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
srmatrix* midpoint(const simatrix& A) {
  srmatrix *R = new srmatrix();
  *R = mid(A);
  return R;
}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
scmatrix* midpoint(const scimatrix& A) {
  scmatrix *R = new scmatrix();
  *R = mid(A);
  return R;
}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
rmatrix* midpoint(const rmatrix& A) {
  return const_cast<rmatrix*>(&A);
}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
cmatrix* midpoint(const cmatrix& A) {
  return const_cast<cmatrix*>(&A);

}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
rmatrix* midpoint(const imatrix& A) {
  rmatrix *R = new rmatrix();
  *R = mid(A);
  return R;
}

//! Return pointer to the midpoint of a matrix. For point matrices this is just a pointer to the matrix itself.
cmatrix* midpoint(const cimatrix& A) {
  cmatrix *R = new cmatrix();
  *R = mid(A);
  return R;
}


/*!
/* Compute aproximate LU-decomposition of banded matrix A with high precision without pivoting. 
/* Store factors in Lo,Up and store enclosure of defect LU - A in LU_A. All matrices are stored 
/* as dense matrices containing only the bands.
/* "A" and "LU_A" must be [1..n, -l..k]
/* Lo  must be [1..n, -1..0]
/* Up  must be [1..n,  0..k]
*/
template<typename TMat, typename TIMat, typename TDot>
void lu_decomp( TMat& A, TMat& Lo, TMat& Up, TIMat& LU_A, int& Err ) {
  int          i,j,m,n,k,l;
  TDot         dot;
  
  n = Ub(A,ROW);
  l = cxsc::abs( Lb(A,COL));
  k = Ub(A,COL);

  Lo   = 0;
  
  for(i = 1; i <= n; i++)   {
    Lo[i][0] = 1.0;

    for(j = i; j <= min( n,i + max(k,l) ); j++)  {
        if ( (j - i) <= k ) {
          dot = 0.0;
          for(m = max(1,i-l,j-k); m <= min(i-1,j); m++) 
            accumulate(dot,Lo[i][m-i],Up[m][j-m]);
          dot = A[i][j-i] - dot;
          Up[i][j-i] = rnd(dot);
          dot = Up[i][j-i] - dot;
          rnd(dot,LU_A[i][j-i]);
        }
    
	if ( (i != j) && (j-i <= l) ) {
          dot = 0.0;
          for(m = max(1,j-l,i-k); m <= min(j,i-1); m++) 
            accumulate(dot,Lo[j][m-j],Up[m][i-m]);
          dot = A[j][i-j] - dot;
	  if(Up[i][0] == 0.0) {
	    Err = 1;
	    return;
	  } else {
            Lo[j][i-j] = rnd(dot)/Up[i][0];
	  }
	  accumulate(dot,-Lo[j][i-j],Up[i][0]);
          rnd(-dot,LU_A[j][i-j]);
        }
    }
    
  }
  
  Err = 0;
} // end lu_decomp
    
    
/*!
/* Compute aproximate  cholesky decomposition of banded matrix A with high precision without pivoting. 
/* Store factor in Lo and store enclosure of defect LL' - A in LU_A. All matrices are stored 
/* as dense matrices containing only the bands.
/* "A" and "LU_A" must be [1..n, -l..k]
/* Lo  must be [1..n, -1..0]
*/
void chol( rmatrix& A, rmatrix& Lo, imatrix& LU_A, int& Err ) {
  int          i,j,m,n,k,l;
  dotprecision         dot;
        
  n = Ub(A,ROW);
  l = cxsc::abs( Lb(A,COL));
  k = Ub(A,COL);
        
  if(l != k) {
    //A not symmetric
    Err = 1;
    return;
  }
        
  Lo   = 0;
	
        
  for(i = 1; i <= n; i++)   {
    dot = A[i][0];
    
    for(j = -l; j < 0; j++)
      accumulate(dot,-Lo[i][j],Lo[i][j]);
      Lo[i][0] = rnd(dot);
      if(Lo[i][0] < 0.0) {
        Err = 2;
        return;
      }
      Lo[i][0] = sqrt(Lo[i][0]);
      rnd(dot,LU_A[i][0]);
      LU_A[i][0] = sqrt(LU_A[i][0]) - Lo[i][0];

      for(k = i+1 ; k<=min(n,i+l) ; k++) {
        dot = 0.0;
        for(m = -l; m < 0; m++) {
          if(k+m > 0  &&  m-(k-i)>=-l)
            accumulate(dot,Lo[k][m-(k-i)],Lo[i][m]);
	}
        dot = A[k][i-k] - dot;
        Lo[k][i-k] = rnd(dot)/Lo[i][0];                    
        accumulate(dot,-Lo[k][i-k],Lo[i][0]);
        rnd(-dot,LU_A[k][i-k]);
        LU_A[i][-(i-k)] = LU_A[k][i-k];
      }
            
  }        
        
  Err = 0;
} // end chol


/*!
/* Compute aproximate cholesky decomposition of complex banded matrix A with high precision without pivoting. 
/* Store factor in Lo and store enclosure of defect LL' - A in LU_A. All matrices are stored 
/* as dense matrices containing only the bands.
/* "A" and "LU_A" must be [1..n, -l..k]
/* Lo  must be [1..n, -1..0]
*/
void chol( cmatrix& A, cmatrix& Lo, cimatrix& LU_A, int& Err ) {
  int          i,j,m,n,k,l;
  cdotprecision         dot;
        
  n = Ub(A,ROW);
  l = cxsc::abs( Lb(A,COL));
  k = Ub(A,COL);
        
  if(l != k) {
    //A not symmetric
    Err = 1;
    return;
  }
        
  Lo   = 0;
        
  for(i = 1; i <= n; i++)   {
    dot = A[i][0];
    for(j = -l; j < 0; j++)
      accumulate(dot,-Lo[i][j],conj(Lo[i][j]));
    Lo[i][0] = rnd(dot);
    if(abs(Lo[i][0]) < 0.0) {
      Err = 2;
      return;
    }
    Lo[i][0] = sqrt(Lo[i][0]);
    rnd(dot,LU_A[i][0]);
    LU_A[i][0] = sqrt(LU_A[i][0]) - Lo[i][0];

    for(k = i+1 ; k<=min(n,i+l) ; k++) {
      dot = 0.0;
      for(m = -l; m < 0; m++) {
        if(k+m > 0  &&  m-(k-i)>=-l)
        accumulate(dot,Lo[k][m-(k-i)],conj(Lo[i][m]));
      }
      dot = A[k][i-k] - dot;
      Lo[k][i-k] = rnd(dot)/Lo[i][0];                    
      accumulate(dot,-Lo[k][i-k],Lo[i][0]);
      rnd(-dot,LU_A[k][i-k]);
      LU_A[i][-(i-k)] = LU_A[k][i-k];
    }
            
 }
        
 Err = 0;
} // end chol

//! Return the diagonal of a sparse matrix U as a dense vector
rvector diag(const srmatrix& U) {
  int n = RowLen(U);
  rvector d(n);
  const vector<int>& p = U.column_pointers();
  const vector<real>& val = U.values();
  
  for(int j=0 ; j<n ; j++) {
    d[j+1] = val[p[j+1]-1];      
  }
  
  return d;
}

/*! Compute enclosure C=[LU-A]. If parameter fast is true, a faster but less precise version needing
 *  only one matrix-matrix-product instead of two is used. Both versions call function defined in utility.cpp
 *  for the computation.
 */
void LUMinusA(const srmatrix& L, const srmatrix& U, const srmatrix& A, simatrix& C, bool fast) {
  if(fast) {

    int n = RowLen(L);
    srmatrix D(n,n);
    real denom = subd(1.0, mulu(n,Epsilon));
    real gamma = divu(mulu(n,Epsilon), denom);
    absLU(L,U,D);
    C = interval(-gamma,gamma) * D;

    //handle underflow
    rvector tmp(n);
    tmp = 1.0;
    real eta = power(2,-1074);
    fesetround(FE_UPWARD);
    tmp = n*tmp + diag(abs(U));
    tmp /= denom;
    fesetround(FE_TONEAREST);

    vector<int>& ind = C.row_indices();
    vector<interval>& val = C.values();

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0 ; i<ind.size() ; i++) {
      val[i] += tmp[ind[i]+1] * eta;
    }

  } else {

    srmatrix C_inf, C_sup;      
    LUMinusA(L,U,A,C_inf,C_sup);
    C = C_inf | C_sup;

  }
}

/*! Compute enclosure C=[LU-A]. If parameter fast is true, a faster but less precise version needing
 *  only one matrix-matrix-product instead of two is used. Both versions call function defined in utility.cpp
 *  for the computation.
 */
void LUMinusA(const scmatrix& L, const scmatrix& U, const scmatrix& A, scimatrix& C, bool fast) {
  if(fast) {

    int n = RowLen(L);
    srmatrix D(n,n);
    real denom = subd(1.0, mulu(2*n,Epsilon));
    real denom2 = subd(1.0, mulu(7*n,Epsilon));
    real gamma = divu(mulu(n,Epsilon), denom);
    absLU(L,U,D);
    C = interval(-gamma,gamma) * D;

    //handle underflow
    rvector tmp(n);
    real eta = power(2,-1074);
    fesetround(FE_UPWARD);
    tmp = ((2*n*eta) / denom);
    tmp += (diag(abs(U)) / denom2) * eta;
    fesetround(FE_TONEAREST);

    vector<int>& ind = C.row_indices();
    vector<cinterval>& val = C.values();

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif	
    for(unsigned int i=0 ; i<ind.size() ; i++) {
      val[i] += tmp[ind[i]+1];
      val[i] *= sqrt(2)/(1+Epsilon);
    }

    C *= complex(-1,1);

  } else {

    scmatrix C_inf, C_sup;
    LUMinusA(L,U,A,C_inf,C_sup);
    C = scimatrix(C_inf) | scimatrix(C_sup);

  }
}

/*! Compute enclosure C=[LU-A]. If parameter fast is true, a faster but less precise version needing
 *  only one matrix-matrix-product instead of two is used. Both versions call function defined in utility.cpp
 *  for the computation.
 */
void LUMinusA(const srmatrix& L, const srmatrix& U, const simatrix& A, simatrix& C, bool fast) {
  if(fast) {

    int n = RowLen(L);
    srmatrix D(n,n);
    real denom = subd(1.0, mulu(n,Epsilon));
    real gamma = divu(mulu(n,Epsilon), denom);
    absLU(L,U,D);
    C = interval(-gamma,gamma) * D;

    //handle underflow
    rvector tmp(n), eta(n);
    tmp = 1.0;
    eta = power(2,-1074);
    fesetround(FE_UPWARD);
    tmp = n*tmp + diag(abs(U));
    tmp /= denom;
    fesetround(FE_TONEAREST);

    vector<int>& ind = C.row_indices();
    vector<interval>& val = C.values();

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(unsigned int i=0 ; i<ind.size() ; i++) {
      val[i] += tmp[ind[i]+1] * eta[ind[i]+1];
    }

    C += 0.5 * (diam(A) * interval(-1,1));	
 
 } else {
    LUMinusA(L,U,A,C);
 }
}

/*! Compute enclosure C=[LU-A]. If parameter fast is true, a faster but less precise version needing
 *  only one matrix-matrix-product instead of two is used. Both versions call function defined in utility.cpp
 *  for the computation.
 */
void LUMinusA(const scmatrix& L, const scmatrix& U, const scimatrix& A, scimatrix& C, bool fast) {
  if(fast) {
    int n = RowLen(L);
    srmatrix D(n,n);
    real denom = subd(1.0, mulu(2*n,Epsilon));
    real denom2 = subd(1.0, mulu(7*n,Epsilon));
    real gamma = divu(mulu(n,Epsilon), denom);
    absLU(L,U,D);
    C = interval(-gamma,gamma) * D;
	
    //handle underflow
    rvector tmp(n);
    real eta = power(2,-1074);
    fesetround(FE_UPWARD);
    tmp = ((2*n*eta) / denom);
    tmp += (diag(abs(U)) / denom2) * eta;
    fesetround(FE_TONEAREST);
	
    vector<int>& ind = C.row_indices();
    vector<cinterval>& val = C.values();

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif	
    for(unsigned int i=0 ; i<ind.size() ; i++) {
      val[i] += tmp[ind[i]+1];
      val[i] *= sqrt(2)/(1+Epsilon);
    }

    C *= complex(-1,1);
    C += 0.5 * (diam(A) * interval(-1,1));
	
  } else {
    LUMinusA(L,U,A,C);
  }
}

//! Blow up the radius of the interval matrix C by the radius of A. Nothing to do for point matrices.
void addRadA(simatrix &C, const srmatrix& A) {
  return;
}

//! Blow up the radius of the interval matrix C by the radius of A.
void addRadA(simatrix &C, const simatrix& A) {
  fesetround(FE_UPWARD);
  srmatrix rad = 0.5 * diam(A);
  fesetround(FE_TONEAREST);
  C += rad * interval(-1,1);
  return;
}

//! Blow up the radius of the interval matrix C by the radius of A. Nothing to do for point matrices.
void addRadA(scimatrix &C, const scmatrix& A) {
  return;
}

//! Blow up the radius of the interval matrix C by the radius of A.
void addRadA(scimatrix &C, const scimatrix& A) {
  fesetround(FE_UPWARD);
  scmatrix rad = 0.5 * diam(A);
  fesetround(FE_TONEAREST);
  C += rad * interval(-1,1);
  return;
}

//! Compute final enclosure for fast verification step
void fastComputeXX(ivector& t1, ivector& t2, ivector& xx, bool& IsVerified) {
  int n = VecLen(t1);
  t1 = Sup(abs(t1));
  t2 = Sup(abs(t2));
  
  if( isLarger(t1,t2) ) {
	
    real max = 0.0, tmpm;
    for(int i=1 ; i<=n ; i++) {
      tmpm = Sup(t2[i]) / Inf((t1[i]-t2[i]));
      if(tmpm > max) max = tmpm;
    }

    xx = addu(1,max) * (t1 * interval(-1,1));
    IsVerified = true;
    
  } else {
    IsVerified = false;
  }
}

//! Compute final enclosure for fast verification step
void fastComputeXX(civector& t1, civector& t2, civector& xx, bool& IsVerified) {
  int n = VecLen(t1);
  ivector tt1, tt2;
  tt1 = Sup(abs(t1));
  tt2 = Sup(abs(t2));
  
  if( isLarger(tt1,tt2) ) {
	
    real max = 0.0, tmpm;
    for(int i=1 ; i<=n ; i++) {
      tmpm = Sup(tt2[i]) / Inf((tt1[i]-tt2[i]));
      if(tmpm > max) max = tmpm;
    }
    
    real k = addu(1,max);
    SetRe(xx, tt1 * interval(-k,k));
    SetIm(xx, tt1 * interval(-k,k));
    IsVerified = true;
    
  } else {
    IsVerified = false;
  }
}

/*!
// Main function of the solver. Solve the sparse system Ao*x=bo.
//
// Parameters:
//    In : 'Ao          : Sparse matrix of the system
//         'bo'         : Right-hand side(s) of the system
//         'cfg'        : Configuration for this solver run
//    Out: 'x'          : Enclosure of the unique solution
//         'Err'        : Error code
*/
template<typename TSysMat, typename TRhs, typename TSolution, typename TMidMat, typename TMidRhs, typename TC, 
         typename TRhsVec, typename TVec, typename TIVec, typename TElement, typename TDot, typename TIDot, typename TSIMat>
void slssmain ( const TSysMat  Ao, const TRhs  bo, TSolution&  x, int&  Err, slssconfig cfg )
{

  int     n = ColLen(Ao);                // Length of the columns of A
  int     m = RowLen(Ao);                // Length of the rows of A
  int     rhs = RowLen(bo);              // Number of right hand sides
  bool    IsVerified = true;             // Verification succesful?
  TMidMat L,U;                           // Stores decmoposition of A
  TC C(n,n);                             // C := [LU-A]
  intvector p(n), q(n);                  // Permutation vectors
  int oldprec = opdotprec;

  if (n != m)                            // Error: A is not square
    { Err = NotSquare; return; }

  if (n != ColLen(bo))                   // Error: Dimensions of A and b
    { Err = DimensionErr; return; }      //        are not compatible
    
  if(cfg.spd) {
    if(Ao != transpherm(Ao))
    { Err = NotPosDef; return; }
  }

  #ifdef _OPENMP
  if(cfg.threads > 0) {
     omp_set_num_threads(cfg.threads);
  } else {
     cfg.threads = omp_get_max_threads(); 
     omp_set_num_threads(cfg.threads);
  }
  if(cfg.msg) std::cout << "Using " << cfg.threads << " thread(s) for computations" << std::endl;    
  #endif

  double start, end;

  //Compute midpoint if necessary
  TMidMat& Am = *(midpoint(Ao));
  TMidRhs& bm = *(midpoint(bo));

  if(cfg.msg) std::cout << "Matrix decomposition..." << std::endl;
  
  if( cfg.forceHighPrecIterMatrix || (cfg.triMethod==BANDED && !cfg.forceSuiteSparse)) {
    if(cfg.msg) std::cout << "  Using high precision decomposition for banded systems..." << std::endl;
    start = GetTime();
    int Al = lower_bandwidth(Ao);
    int Au = upper_bandwidth(Ao);

    //Transform sparse banded matrices into dense matrices storing the bands
    TMidRhs AB(1,n,-Al,Au), LB(1,n,-Al,0), UB(1,n,0,Au);
    TSolution CB(1,n,-Al,Au);
    
    AB = 0.0;
    CB = 0.0;
    
    for(int j=1 ; j<=Al ; j++) {
      for(int i=-Al ; i<= Au ; i++) {
	if(j+i<1 || j+i>n)
	  AB[j][i] = 0.0;
	else {
	  AB[j][i] = Am(j,j+i);
	}
      }
    }
    
    for(int j=Al+1 ; j<=n-Au ; j++) {
      for(int i=-Al ; i<= Au ; i++) {
	AB[j][i] = Am(j,j+i);
      }
    }

    for(int j=n-Au+1 ; j<=n ; j++) {
      for(int i=-Al ; i<= Au ; i++) {
	if(j+i<1 || j+i>n)
	  AB[j][i] = 0.0;
	else
	  AB[j][i] = Am(j,j+i);
      }
    }

    //Decomposition
    if(cfg.spd) 
      chol(AB,LB,CB,Err);
    else
      lu_decomp<TMidRhs,TSolution,TDot>(AB, LB, UB, CB, Err);
    
    if(Err != 0) {
      Err = DecompFailed;
      opdotprec = oldprec;
      return;
    }

    //Retransform into true sparse matrices
    L = TMidMat(n,n,LB);
    if(cfg.spd)
      U = transpherm(L);
    else
      U = TMidMat(n,n,UB);
    C = TC(n,n,CB);

    //Add radius of Ao to computed C if Ao is interval matrix.
    addRadA(C,Ao);
    
    end = GetTime();
    if(cfg.msg) std::cout << "Time for matrix decomposition: " << end-start << std::endl;

  } else {
    
    //Compute decomposition with CHOLMOD or UMFPACK
    
    if(cfg.spd) {
      
      if(cfg.msg) std::cout << "  Using CHOLMOD for decomposition..." << std::endl;
      start = GetTime();
      chol(Am, L, p, Err);
      end = GetTime();
      U = transpherm(L);
      q = p;
      if(cfg.msg) std::cout << "Time for Cholesky-decomposition: " << end-start << std::endl;
      
    } else {
      
      if(cfg.msg) std::cout << "  Using UMFPACK for decomposition..." << std::endl;
      start = GetTime();
      lu_decomp(Am, L, U, p, q, Err);
      end = GetTime();
      if(cfg.msg) std::cout << "Time for LU-decomposition: " << end-start << std::endl;
      
    }

    if(Err != 0) {
      Err = DecompFailed;
      opdotprec = oldprec;
      return;
    }
    
  }
  

  //Perform permutation if necessary
  TSysMat A(Ao);
  TRhs    bmat(bo);

  if(!cfg.forceHighPrecIterMatrix && (cfg.triMethod!=BANDED || cfg.forceSuiteSparse)) {
    A    = A(p,q);
    bmat = bmat(p);
    Am   = Am(p,q);
    bm   = bm(p);
  }

  //Remaining allocations and preparations
  TVec        d(n);       
  TIVec       xxb(n), zz(n), tmp(n);
  TDot        Accu;       
  TIDot       IAccu;

  Resize(x,n,rhs);
  
  //Compute iteration matrix
  if(!cfg.forceHighPrecIterMatrix && (cfg.triMethod!=BANDED || cfg.forceSuiteSparse)) {
    if(cfg.msg) std::cout << "Computing iteration matrix [LU-A]..." << std::endl;
    start = GetTime();
    opdotprec = cfg.precIterMat;
    if(cfg.msg) std::cout << "  Using precision K=" << opdotprec << std::endl;
    
    if(opdotprec != 1) {
      C = L*TSIMat(U) - A;
    } else {
      LUMinusA(L,U,A,C,cfg.fastDecompErrorComputation);
    }
    end = GetTime();
    
    if(cfg.msg) std::cout << "Time for iteration matrix: " << end-start << std::endl;
  } 
  
  
  //Loop over all right hand sides
  for(int s=1 ; s<=rhs && IsVerified ; s++) {
    
    if(cfg.msg && rhs>1) std::cout << "Computing result for right hand side " << s << "..." << std::endl;
    TIVec xx(n);
    TRhsVec b(bmat[Col(s)]);

    //Approximate solution and residual
    if(cfg.msg) std::cout << "Computing approximate solution and residual..." << std::endl;
    start = GetTime();
    int prec = (cfg.apprxSolutionStagPrec > 0) ? cfg.apprxSolutionStagPrec : 1;
    TVec xapp(n);
    TMidRhs xappstag(n,prec);

    if(cfg.highPrecApprxSolution) {
      if(cfg.msg) std::cout << "  Computing high precision solution" << std::endl;    
      TVec h(n);
      d = bm[Col(s)];
    
      for(int p=1 ; p<prec ; p++) {
        lowtrisolve(L,d,h);
        uptrisolve(U,h,xapp);
        xappstag[Col(p)] = xapp;
        residual(d,bm[Col(s)],Am,xappstag,p);
      }
    
      lowtrisolve(L,d,h);
      uptrisolve(U,h,xapp);
      xappstag[Col(prec)] = xapp;
      residual(zz,b,A,xappstag,prec);
    
    } else {
    
      if(cfg.msg) std::cout << "  Computing floating point solution" << std::endl;    
      TVec h(bm[Col(s)]);
      lowtrisolve(L,h,d);
      uptrisolve(U,d,xapp);
      xx = xapp;
      opdotprec = cfg.precResidual;
      zz = b - A*xx;
          
    }
    end = GetTime();
    if(cfg.msg) std::cout << "Time for approximate solution and residual: " << end-start << std::endl;

    //Start verification
    if(cfg.msg) std::cout << "Starting verification..." << std::endl;
    xx = zz;
    opdotprec = 1;
    int steps = 0;
    real eps = 0.1;
    int lb = lower_bandwidth(L);
    int ub = upper_bandwidth(U);   

    TMidMat LT, UT;
    TSolution BasisL, BasisU;    
   
    if(cfg.triMethod == BANDED) {
      LT = transp(L);
      UT = transp(U);
    }    

    start = GetTime();
    if(isZeroVec(zz)) {  
      //Residual zero, exact solution (but not necessarily unique)
      IsVerified = true;
   
    } else if(cfg.fastVerification) {
   
      if(cfg.msg) std::cout << "  Using fast verification mode" << std::endl;
  
      TIVec t1(xx);
      TIVec t2(t1);
      
      if(cfg.msg) std::cout << "  Computing [t1] = L\\U\\[z]" << std::endl;      

      if(cfg.triMethod == FORWARD_BACKWARD)
        lowtrisolve(L, xx, tmp);  //forward substitution
      else if (cfg.triMethod == IMPLICIT_INVERSE)
        lowtrisolve_inv(L, xx, tmp); //implicit inverse
      else 
        lowtrisolve_band(L, LT, xx, tmp, lb, BasisL, cfg.storeBase, Err); //band      	    
      
      if(Err != 0) {
	IsVerified = false;
	continue;
      }
      
      if(cfg.triMethod == FORWARD_BACKWARD)
        uptrisolve(U, tmp, t1); //backward subst.
      else if (cfg.triMethod == IMPLICIT_INVERSE)
        uptrisolve_inv(U, tmp, t1); //implicit inverse
      else 
        uptrisolve_band(U, UT, tmp, t1, ub, BasisU, cfg.storeBase, Err); //band

      if(Err != 0) {
	IsVerified = false;
	continue;
      }	
	
      if(cfg.msg) std::cout << "  Computing [t2] = L\\U\\[LU-A]*[t1]" << std::endl;
      #ifdef _OPENMP
        if(cfg.threads == 1) {
	  xx = C * t1;
	} else {
          int size = n / (omp_get_max_threads() * 10);
          if(size <= 0) size = n;	  
          #pragma omp parallel for schedule(dynamic,1)
	  for(int i=0 ; i<n/size ; i++)
	    xx(i*size+1,(i+1)*size) = C(i*size+1,(i+1)*size,1,n) * t1;
	  int i = n/size;
          xx(i*size+1,n) = C(i*size+1,n,1,n) * t1;	
	}
      #else
          xx = C * t1;
      #endif
	      	
      if(cfg.triMethod == FORWARD_BACKWARD)
        lowtrisolve(L, xx, tmp);  //forward substitution
      else if (cfg.triMethod == IMPLICIT_INVERSE)
        lowtrisolve_inv(L, xx, tmp); //implicit inverse
      else  
        lowtrisolve_band(L, LT, xx, tmp, lb, BasisL, cfg.storeBase, Err); //band      	    

      if(Err != 0) {
	IsVerified = false;
	continue;
      }
	
      if(cfg.triMethod == FORWARD_BACKWARD)
        uptrisolve(U, tmp, t2); //backward subst.
      else if (cfg.triMethod == IMPLICIT_INVERSE)
        uptrisolve_inv(U, tmp, t2); //implicit inverse
      else 
        uptrisolve_band(U, UT, tmp, t2, ub, BasisU, cfg.storeBase, Err); //band

      if(Err != 0) {
	IsVerified = false;
	continue;
      }
	
      if(cfg.msg) std::cout << "  Trying to compute verified enclosure" << std::endl;
      fastComputeXX(t1,t2,xx,IsVerified);      
	
    } else {  //if(cfg.fastVerification)

      if(cfg.msg) std::cout << "  Using normal verification mode" << std::endl;

      do {
        if(cfg.msg) std::cout << "  Verification step " << steps << std::endl;
	steps++;	
        xxb = blow(xx,eps);
        #ifdef _OPENMP
	  if(cfg.threads == 1) {
	    xx = zz + C*xxb;
	  } else {
	    int size = n / (omp_get_max_threads() * 10);
	    if(size <= 0) size = n;
            #pragma omp parallel for schedule(dynamic,1)
	    for(int i=0 ; i<n/size ; i++)
	      xx(i*size+1,(i+1)*size) = zz(i*size+1,(i+1)*size) + C(i*size+1,(i+1)*size,1,n) * xxb;
	    int i = n/size;
	    if(n % size != 0)
              xx(i*size+1,n) = zz(i*size+1,n) + C(i*size+1,n,1,n) * xxb;	
	  }
        #else
          xx = zz + C*xxb;     
        #endif
	      
        if(cfg.msg) std::cout << "  Solving lower triangular system" << std::endl;
      
        if(cfg.triMethod == FORWARD_BACKWARD)
          lowtrisolve(L, xx, tmp);  //forward substitution
        else if (cfg.triMethod == IMPLICIT_INVERSE)
          lowtrisolve_inv(L, xx, tmp); //implicit inverse
        else 
	  lowtrisolve_band(L, LT, xx, tmp, lb, BasisL, cfg.storeBase, Err); //band      	    
	
        if(Err != 0) {
	  IsVerified = false;
	  continue;
        }

        if(cfg.msg) std::cout << "  Solving upper triangular system" << std::endl;
      
        if(cfg.triMethod == FORWARD_BACKWARD)
          uptrisolve(U, tmp, xx); //backward subst.
        else if (cfg.triMethod == IMPLICIT_INVERSE)
          uptrisolve_inv(U, tmp, xx); //implicit inverse
        else 
          uptrisolve_band(U, UT, tmp, xx, ub, BasisU, cfg.storeBase, Err); //band
	  
        if(Err != 0) {
	  IsVerified = false;
	  continue;
        }  
    
        IsVerified = in(xx,xxb);
      } while(!IsVerified && steps<=cfg.maxIterVer);
    
    } //if(cfg.FastVerification)
    
    end = GetTime();
    if(cfg.msg) std::cout << "Time for verification: " << end-start << std::endl;


    if(IsVerified) {
      if(cfg.msg) std::cout << "Verification successful " << std::endl;
    
      if(cfg.highPrecApprxSolution) {
        #ifdef _OPENMP
        #pragma omp parallel for private(IAccu)
        #endif
        for(int i=1 ; i<=n ; i++) {
          IAccu = xx[i];
	  for(int c=1 ; c<=prec ; c++)
            IAccu += xappstag[i][c];
          rnd(IAccu,xx[i]);
        }    
      } else {
        xx = xapp + xx;
      }
    
    } else {
    
      if(cfg.msg) std::cout << "Verification failed " << std::endl;
      if(cfg.highPrecApprxSolution) {
        xx = TVec(xappstag[Col(1)]);
      } else {
        xx = xapp;
      }
      Err = VerivFailed;
    
   }
 
   if(!cfg.forceHighPrecIterMatrix && (cfg.triMethod!=BANDED || cfg.forceSuiteSparse)) {
     xx = xx(perminv(q));
   }
 
   x[Col(s)] = xx;
  
   if(cfg.msg && rhs>1) std::cout << "Right hand side " << s << " finished..." << std::endl << std::endl;
   
 } //loop over all right hand sides


 opdotprec = oldprec; 
 
 if((TRhs*)&bm != &bo)  delete &bm;
 if((TSysMat*)&Am != &Ao)  delete &Am;
  
} //slss

/*! Central start function for the solver. Checks if system is square, transforms under- or overdetermined systems
 *  into corresponding square systems and starts the main solver routine. */
template<typename TSysMat, typename TRhs, typename TSolution, typename TMidMat, typename TMidRhs, typename TC, 
         typename TRhsVec, typename TVec, typename TIVec, typename TElement, typename TDot, typename TIDot, typename TSIMat>
inline void slssstart ( const TSysMat& A, const TRhs& b, TSolution& xx, int& Err, slssconfig cfg) {
  int m   = Ub(A,ROW) - Lb(A,ROW) + 1;
  int n   = Ub(A,COL) - Lb(A,COL) + 1;
  int rhs = Ub(b,COL) - Lb(b,COL) + 1;
  int dim = m+n;
  int olddotprec = opdotprec;
  
  if(m == n) {

    //Square
    slssmain<TSysMat, TRhs, TSolution, TMidMat, TMidRhs, TC, TRhsVec, TVec, TIVec, TElement, TDot, TIDot, TSIMat>(A, b, xx, Err, cfg);

  } else if (m > n) {

    //Overdetermined
    if(cfg.msg) std::cout << "Overdetermined system, generating equivalent square system" << std::endl; 

    TSysMat BIG_A(dim,dim);
    TRhs BIG_b(dim,rhs);
    TSolution  BIG_xx(dim,rhs);

    BIG_A( 1,m, 1,n ) = A;
    BIG_A( 1,m, n+1,n+m ) = -Id(srmatrix(m,m));
    BIG_A( m+1,m+n, n+1,n+m ) = transpherm(A);
    
    BIG_b( m+1,m+n,1,rhs ) = 0.0;    
    BIG_b( 1,m,1,rhs ) = b;
            
    slssmain<TSysMat, TRhs, TSolution, TMidMat, TMidRhs, TC, TRhsVec, TVec, TIVec, TElement, TDot, TIDot, TSIMat>( BIG_A, BIG_b, BIG_xx, Err, cfg); 
 
    xx = BIG_xx( 1,n,1,rhs );    

  } else if (m < n) {

    //Underdetermined
    if(cfg.msg) std::cout << "Underdetermined system, generating equivalent square system" << std::endl; 

    TSysMat BIG_A(dim,dim);
    TRhs BIG_b(dim,rhs);
    TSolution  BIG_xx(dim,rhs);

    BIG_A( 1,n, 1,m ) = transpherm(A);
    BIG_A( 1,n, m+1,m+n ) = -Id(srmatrix(n,n));
    BIG_A( n+1,n+m, m+1,m+n ) = A;
    
    BIG_b( 1,n,1,rhs ) = 0.0;    
    BIG_b( n+1,n+m,1,rhs ) = b;
    
    slssmain<TSysMat, TRhs, TSolution, TMidMat, TMidRhs, TC, TRhsVec, TVec, TIVec, TElement, TDot, TIDot, TSIMat>( BIG_A, BIG_b, BIG_xx, Err, cfg);    

    xx = BIG_xx( m+1,m+n,1,rhs );    

  }
  
  opdotprec = olddotprec;
}



//Below: Start functions of the solver for different data types

void  slss ( const srmatrix& A, const rmatrix& b, imatrix& x, int& err, slssconfig cfg ) {
  slssstart<srmatrix,rmatrix,imatrix,srmatrix,rmatrix,simatrix,rvector,rvector,ivector,real,dotprecision,idotprecision,simatrix>(A,b,x,err,cfg);
}

void  slss ( const srmatrix& A, const rvector& b, ivector& x, int& err, slssconfig cfg ) {
  rmatrix B(VecLen(b),1);
  B[Col(1)] = b;
  imatrix X;
  slssstart<srmatrix,rmatrix,imatrix,srmatrix,rmatrix,simatrix,rvector,rvector,ivector,real,dotprecision,idotprecision,simatrix>(A,B,X,err,cfg);
  x = X[Col(1)];
}

void  slss ( const srmatrix& A, const imatrix& b, imatrix& x, int& err, slssconfig cfg ) {
  slssstart<srmatrix,imatrix,imatrix,srmatrix,rmatrix,simatrix,ivector,rvector,ivector,real,dotprecision,idotprecision,simatrix>(A,b,x,err,cfg);
}

void  slss ( const srmatrix& A, const ivector& b, ivector& x, int& err, slssconfig cfg ) {
  imatrix B(VecLen(b),1);
  B[Col(1)] = b;
  imatrix X;
  slssstart<srmatrix,imatrix,imatrix,srmatrix,rmatrix,simatrix,ivector,rvector,ivector,real,dotprecision,idotprecision,simatrix>(A,B,X,err,cfg);
  x = X[Col(1)];
}



void  slss ( const simatrix& A, const imatrix& b, imatrix& x, int& err, slssconfig cfg ) {
  slssstart<simatrix,imatrix,imatrix,srmatrix,rmatrix,simatrix,ivector,rvector,ivector,real,dotprecision,idotprecision,simatrix>(A,b,x,err,cfg);
}

void  slss ( const simatrix& A, const ivector& b, ivector& x, int& err, slssconfig cfg ) {
  imatrix B(VecLen(b),1);
  B[Col(1)] = b;
  imatrix X;
  slssstart<simatrix,imatrix,imatrix,srmatrix,rmatrix,simatrix,ivector,rvector,ivector,real,dotprecision,idotprecision,simatrix>(A,B,X,err,cfg);
  x = X[Col(1)];
}


void  slss ( const scmatrix& A, const cmatrix& b, cimatrix& x, int& err, slssconfig cfg ) {
  slssstart<scmatrix,cmatrix,cimatrix,scmatrix,cmatrix,scimatrix,cvector,cvector,civector,complex,cdotprecision,cidotprecision,scimatrix>(A,b,x,err,cfg);
}

void  slss ( const scmatrix& A, const cvector& b, civector& x, int& err, slssconfig cfg ) {
  cmatrix B(VecLen(b),1);
  B[Col(1)] = b;
  cimatrix X;
  slssstart<scmatrix,cmatrix,cimatrix,scmatrix,cmatrix,scimatrix,cvector,cvector,civector,complex,cdotprecision,cidotprecision,scimatrix>(A,B,X,err,cfg);
  x = X[Col(1)];
}

void  slss ( const scmatrix& A, const cimatrix& b, cimatrix& x, int& err, slssconfig cfg ) {
  slssstart<scmatrix,cimatrix,cimatrix,scmatrix,cmatrix,scimatrix,civector,cvector,civector,complex,cdotprecision,cidotprecision,scimatrix>(A,b,x,err,cfg);
}

void  slss ( const scmatrix& A, const civector& b, civector& x, int& err, slssconfig cfg ) {
  cimatrix B(VecLen(b),1);
  B[Col(1)] = b;
  cimatrix X;
  slssstart<scmatrix,cimatrix,cimatrix,scmatrix,cmatrix,scimatrix,civector,cvector,civector,complex,cdotprecision,cidotprecision,scimatrix>(A,B,X,err,cfg);
  x = X[Col(1)];
}


void  slss ( const scimatrix& A, const cimatrix& b, cimatrix& x, int& err, slssconfig cfg ) {
  slssstart<scimatrix,cimatrix,cimatrix,scmatrix,cmatrix,scimatrix,civector,cvector,civector,complex,cdotprecision,cidotprecision,scimatrix>(A,b,x,err,cfg);
}

void  slss ( const scimatrix& A, const civector& b, civector& x, int& err, slssconfig cfg ) {
  cimatrix B(VecLen(b),1);
  B[Col(1)] = b;
  cimatrix X;
  slssstart<scimatrix,cimatrix,cimatrix,scmatrix,cmatrix,scimatrix,civector,cvector,civector,complex,cdotprecision,cidotprecision,scimatrix>(A,B,X,err,cfg);
  x = X[Col(1)];
}


} // namespace cxsc
