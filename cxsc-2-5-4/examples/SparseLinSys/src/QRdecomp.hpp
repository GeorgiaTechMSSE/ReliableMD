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

/* Header providing function for computing a QR decomposition of a matrix and
   the enclosure of the inverse of a nearly orthogonal matrix */

#ifndef _CXSC_QRDECOMP_HEADER_DEFINED
#define _CXSC_QRDECOMP_HEADER_DEFINED

#include <cimatrix.hpp>
#include <fenv.h>
#include "utility.hpp"

//Naming convention for LAPACK functions
#ifdef Unchanged
  #define dgeqrf_ dgeqrf
  #define dorgqr_ dorgqr
#elif defined(Add_)
  //Nothing to do
#elif defined(Uppercase)
  #define dgeqrf_ DGEQRF
  #define dorgqr_ DORGQR
#endif

//Declaration of LAPACK functions
extern "C" {
  void dgeqrf_(int *m, int *n, double* DA, int *lda, double* tau, double* work, int *lwork, int *info);  
  void dorgqr_(int *m, int *n , int *k, double* DA, int *lda, double* tau, double* work, int *lwork, int *info);
  void zgeqrf_(int *m, int *n, double* DA, int *lda, double* tau, double* work, int *lwork, int *info);  
  void zungqr_(int *m, int *n , int *k, double* DA, int *lda, double* tau, double* work, int *lwork, int *info);
}


namespace cxsc {
  

//! Compute enclosure of the inverse of an almost orthogonal matrix A
inline void QR_inv( const rmatrix& A, imatrix& A_inv, int& err ) {
  real         nrm; 
  interval     delta;
  imatrix      D(ColLen(A),RowLen(A));
  dotprecision Accu;     

  err = 0;
  if (Ub(A,ROW) - Lb(A,ROW) == 0) {
    
	A_inv = interval(1.0)/A[Lb(A,ROW)][Lb(A,COL)];
	
  } else if (Ub(A,ROW) - Lb(A,ROW) == 1) {
    
        Accu = 0.0;
        // delta = ##(a11*a22 - a12*a21)
        accumulate(Accu,A[Lb(A,ROW)][Lb(A,COL)],A[Ub(A,ROW)][Ub(A,COL)]);
        accumulate(Accu,-A[Lb(A,ROW)][Ub(A,COL)],A[Ub(A,ROW)][Lb(A,COL)]);
        rnd(Accu,delta);  

	if ( in(0.0, delta) ) 
            err = 1;
        else {
            A_inv[Lb(A_inv,ROW)][Lb(A_inv,COL)] =  
		      A[Ub(A,ROW)][Ub(A,COL)] / delta;
            A_inv[Lb(A_inv,ROW)][Ub(A_inv,COL)] = 
		      -A[Lb(A,ROW)][Ub(A,COL)] / delta;
            A_inv[Ub(A_inv,ROW)][Lb(A_inv,COL)] = 
			  -A[Ub(A,ROW)][Lb(A,COL)] / delta;
            A_inv[Ub(A_inv,ROW)][Ub(A_inv,COL)] =  
			  A[Lb(A,ROW)][Lb(A,COL)] / delta;
        }
        
  } else {
    
        rmatrix  At, TempMat; 
        At = transp(A);
	TempMat = Id(A) - At*A;

	nrm = Norm00(TempMat);  //nrm = Norm00 ( ##(Id(A) - At*A))
        if (nrm < 1.0) {
          nrm = divu(mulu(Norm00(At),nrm),subd(1.0,nrm)); 
          D   = interval(-nrm,nrm);
	  A_inv = At + D;
        } else {
          err = 1;
          A_inv = At;
        }     
        
  }    

} // end INV



//! Compute enclosure of the inverse of an almost orthogonal matrix A
inline void QR_inv( const cmatrix& A, cimatrix& A_inv, int& err ) {
  real         nrm; 
  cinterval     delta;
  cimatrix      D(ColLen(A),RowLen(A));
  cdotprecision Accu;     

  err = 0;
  if (Ub(A,ROW) - Lb(A,ROW) == 0) {
	A_inv = interval(1.0)/A[Lb(A,ROW)][Lb(A,COL)];
//   } else if (Ub(A,ROW) - Lb(A,ROW) == 1) {
//         Accu = 0.0;
//         // delta = ##(a11*a22 - a12*a21)
//         accumulate(Accu,A[Lb(A,ROW)][Lb(A,COL)],A[Ub(A,ROW)][Ub(A,COL)]);
//         accumulate(Accu,-A[Lb(A,ROW)][Ub(A,COL)],A[Ub(A,ROW)][Lb(A,COL)]);
//         rnd(Accu,delta);  
// 
// 	if ( in(0.0, Re(delta)) && in(0.0, Im(delta)) ) 
//             err = 1;
//         else {
//             A_inv[Lb(A_inv,ROW)][Lb(A_inv,COL)] =  
// 		      A[Ub(A,ROW)][Ub(A,COL)] / delta;
//             A_inv[Lb(A_inv,ROW)][Ub(A_inv,COL)] = 
// 		      -A[Lb(A,ROW)][Ub(A,COL)] / delta;
//             A_inv[Ub(A_inv,ROW)][Lb(A_inv,COL)] = 
// 			  -A[Ub(A,ROW)][Lb(A,COL)] / delta;
//             A_inv[Ub(A_inv,ROW)][Ub(A_inv,COL)] =  
// 			  A[Lb(A,ROW)][Lb(A,COL)] / delta;
//         }
  } else {
        cmatrix  At, TempMat; 
        At = transpherm(A);
	TempMat = Id(srmatrix(ColLen(A),RowLen(A))) - At*A;
	
	nrm = Norm00(TempMat);  //nrm = Norm00 ( ##(Id(A) - At*A))
        if (nrm < 1.0) {
          nrm = divu(mulu(Norm00(At),nrm),subd(1.0,nrm)); 
          D   = cinterval(interval(-nrm,nrm), interval(-nrm,nrm));
	  A_inv = At + D;
        } else {
          err = 1;
          A_inv = At;
        }     
  }    

} // end INV


//!Perform QR-decomposition of A by use of Householder matrices. Optionally use LAPACK for decomposition.
inline void QR( rmatrix& A, rmatrix& Q, int& err) {  
#ifdef CXSC_USE_LAPACK
  Q = transp(A);
  int m = ColLen(A);
  int n = RowLen(A);
  int lda = m;
  int info;
  int lwork = 256*n;
  double* tau = new double[n];
  double* work = new double[lwork];
  double* DA = Q.to_blas_array();
  
  dgeqrf_(&m, &n, DA, &lda, tau, work, &lwork, &info);
  
  if(info != 0) {
    err = 1;
    delete[] work;
    delete[] tau;
    return;
  } 
  
  dorgqr_(&m, &n , &n, DA, &lda, tau, work, &lwork, &info);
  
  if(info != 0) {
    err = 2;
    delete[] work;
    delete[] tau;
    return;
  }
  
  Q = transp(Q);
  
  delete[] work;
  delete[] tau;
  
  err = 0;
  
#else
  
  real         Accu;
  int          i, j, k, n;      
  real         r, s;              
  rvector      b(1,Ub(A,ROW)), c(1,Ub(A,ROW)), d(1,Ub(A,ROW)); 

  n = Ub(A,ROW);
  Q = Id(Q);
  
  b = 0.0;

  for ( k = 1; k <= (n-1); k++) {
    b(k,n) = A[Col(k)](k,n);
    if ( !isZeroVec(rvector(b(k,n))) ) {
      s = sqrt(b*b);
      if (A[k][k] >= 0.0) s = -s;
      b[k] = A[k][k] - s;
      if(b[k] == 0) { err=1; return; }
      s = 1 / (s * b[k]);
      
#ifdef _OPENMP
#pragma omp parallel for private(i,j,r,Accu) 
#endif
      for(i = k+1; i <= n; i++) {
        Accu = 0.0;
        // #*(for j=k to n sum(b[j]*A[i][j])) 
        for(j = k; j <= n; j++) Accu += b[j] * A[j][i];       
        r = s * Accu;
        for(j = k; j <= n; j++) A[j][i] += b[j] * r;
      }                           
              
#ifdef _OPENMP
#pragma omp parallel for private(i,j,r,Accu)
#endif
      for( i = 1; i <= n; i++) {
        Accu = 0.0;
        // #*(for j=k to n sum(b[j]*Q[i][j])) 
        for(j = k; j <= n; j++) Accu += b[j] * Q[i][j];            
        r = s * Accu;
        for(j = k; j <= n; j++) Q[i][j] += b[j] * r;              
      }
    }
    b[k] = 0;
  }
  
#endif
} // end QR


//!Perform QR-decomposition of A by use of Householder matrices. Optionally use LAPACK for decomposition.
inline void QR( const cmatrix& A, cmatrix& Q, int& err) {  
#ifdef CXSC_USE_LAPACK
  Q = transp(A);
  int m = ColLen(A);
  int n = RowLen(A);
  int lda = m;
  int info;
  int lwork = 256*n;
  double* tau = new double[2*n];
  double* work = new double[2*lwork];
  double* DA = Q.to_blas_array();
  
  zgeqrf_(&m, &n, DA, &lda, tau, work, &lwork, &info);
  
  if(info != 0) {
    err = 1;
    delete[] work;
    delete[] tau;
    return;
  }
  
  zungqr_(&m, &n , &n, DA, &lda, tau, work, &lwork, &info);
  
  if(info != 0) {
    err = 2;
    delete[] work;
    delete[] tau;
    return;
  }
  
  Q = transp(Q);
  
  delete[] work;
  delete[] tau;
  
  err = 0;

#else
  complex      Accu;
  int          i, j, k, n;      
  complex      r, s;              
  cvector      b(1,Ub(A,ROW)), c(1,Ub(A,ROW)), d(1,Ub(A,ROW)); 

  n = Ub(A,ROW);
  Q = Id(Q);
  
  b = 0.0;


  for ( k = 1; k <= (n-1); k++) {
    b(k,n) = A[Col(k)](k,n);
    if ( !isZeroVec(cvector(b(k,n))) ) {
      s = sqrt(b*b);
      if (Re(A[k][k]) >= 0.0) SetRe(s,-Re(s));
      if (Im(A[k][k]) >= 0.0) SetIm(s,-Im(s));
      b[k] = A[k][k] - s;
      if(b[k] == 0) { err=1; return; }
      s = 1 / (s * b[k]);

#ifdef _OPENMP
#pragma omp parallel for private(i,j,r,Accu)       
#endif
      for(i = k+1; i <= n; i++) {
        Accu = 0.0;
        // #*(for j=k to n sum(b[j]*A[i][j])) 
        for(j = k; j <= n; j++) Accu += b[j] * A[j][i];       
        r = s * Accu;
        for(j = k; j <= n; j++) A[j][i] += b[j] * r;
      }                           
              
#ifdef _OPENMP
#pragma omp parallel for private(i,j,r,Accu)               
#endif
      for( i = 1; i <= n; i++) {
        Accu = 0.0;
        // #*(for j=k to n sum(b[j]*Q[i][j])) 
        for(j = k; j <= n; j++) Accu += b[j] * Q[i][j];            
        r = s * Accu;
        for(j = k; j <= n; j++) Q[i][j] += b[j] * r;              
      }
    }
    b[k] = 0;
  }
  
  err = 0;
#endif
} // end QR


/*! Sort columns of A according to decreasing lenght 
/* and then do a QR-decomp. */
inline void QR_sorted( rmatrix& A, ivector& y, rmatrix& Q, int& err ) { 
  int     i, j, k, n;
  rvector hv(1,Ub(A,ROW)), laenge(1,Ub(A,ROW));
  real    hr;

  n = Ub(A,ROW);

  for(i = 1; i <= n; i++) {
    real coll = rvector(A[Col(i)])*rvector(A[Col(i)]);
      laenge[i] = sqr(diam(y[i]))*( coll );
  }
  

  for(j = 1; j <= (n-1); j++)
    for( k = j+1; k <= n; k++) 	{
      if (laenge[j] < laenge[k]) {
        hr = laenge[k]; laenge[k] = laenge[j]; laenge[j] = hr;
        hv = rvector(Col(A,k)); Col(A,k) = Col(A,j); Col(A,j) = hv;
      }
    }
  
  QR(A,Q,err);  
    
} // end QR


/*! Sort columns of A according to decreasing lenght 
/* and then do a QR-decomp. */
inline void QR_sorted( cmatrix& A, civector& y, cmatrix& Q, int& err ) { 
  int     i, j, k, n;
  cvector hv(1,Ub(A,ROW));
  rvector laenge(1,Ub(A,ROW));
  real    hr;

  n = Ub(A,ROW);
 
  for(i = 1; i <= n; i++) {
    laenge[i] = abs( (diam(y[i])*diam(y[i])) * ( cvector(A[Col(i)])*cvector(A[Col(i)]) ) );
  }
  
  for(j = 1; j <= (n-1); j++)
    for( k = j+1; k <= n; k++) 	{
      if (laenge[j] < laenge[k]) {
        hr = laenge[k]; laenge[k] = laenge[j]; laenge[j] = hr;
        hv = cvector(Col(A,k)); Col(A,k) = Col(A,j); Col(A,j) = hv;
      }
    }

  QR(A,Q,err);
} // end QR


} //namespace cxsc


#endif