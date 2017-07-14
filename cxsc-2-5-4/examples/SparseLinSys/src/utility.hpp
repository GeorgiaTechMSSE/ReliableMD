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

/* Header providing a number of utility functions used by the solvers */


#ifndef _CXSC_UTILITY_HEADER_INCLUDED
#define _CXSC_UTILITY_HEADER_INCLUDED

#include <scimatrix.hpp>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cxsc {

//! Return min a and b  
static inline int min(int a,int b) {
  return (a < b)? a : b;
} 

//! Return max a and b  
static inline int max(int a, int b) { 
  return (a > b) ? a : b; 
} 

//! Return min a and b and c
static inline int max(int a, int b, int c) { 
  return (a > b) ? max(a,c) : max(b,c); 
}

//! Return is r is NaN
static inline bool isNaN(const real& r) {
  return isnan(_double(r));
}

//! Return is r is NaN
static inline bool isNaN(const complex& r) {
  return ( isnan(_double(Re(r))) || isnan(_double(Im(r))) );
}

//! Return is r is NaN
static inline bool isNaN(const interval& r) {
  return ( isNaN(Inf(r)) || isNaN(Sup(r)) );
}

//! Return is r is NaN
static inline bool isNaN(const cinterval& r) {
  return ( isNaN(Inf(r)) || isNaN(Sup(r)) );
}

//! Return is i is empty
static inline bool isEmpty(const interval& i) {
  return Inf(i) > Sup(i);
}

//! Return is i is empty
static inline bool isEmpty(const cinterval& i) {
  return ( isEmpty(Re(i)) || isEmpty(Im(i)) );
}

//! Return is t1>t2
static inline bool isLarger(const ivector& t1, const ivector& t2) {
  return Inf(t1) > Sup(t2) ;
}

//! Return is re(t1)>re(t2) and im(t1)>im(t2)
static inline bool isLarger(const civector& t1, const civector& t2) {
  return ( InfRe(t1) > SupRe(t2)  &&  InfIm(t1) > SupIm(t2) );
}

//! Return if c lies in the interior of i
static inline bool in(const complex& c, const cinterval& i) {
  return ( in(Re(c),Re(i)) && in(Im(c),Im(i)) );
}

//! Check if r can be a diagonal element of a spd matrix
static inline bool checkPosDef(const real& r) {
  return r > 0.0;
}

//! Check if r can be a diagonal element of a hermitian positive definite matrix
static inline bool checkPosDef(const complex& r) {
  return (Re(r) > 0.0) && (Im(r) == 0.0);
}

//! Compute conjugate transposed of A
static inline srmatrix transpherm ( const srmatrix& A ) {                                            
  return transp(A);
}

//! Compute conjugate transposed of A
static inline simatrix transpherm ( const simatrix& A ) {                                            
  return transp(A);
}

//! Compute conjugate transposed of A
static inline scmatrix transpherm ( const scmatrix& A ) {                                   
  scmatrix res(transp(A));
  
  std::vector<complex>& val = res.values();
  
  for(unsigned int i=0 ; i<val.size() ; i++)
    SetIm(val[i], -Im(val[i]));
  
  return res;
}

//! Compute conjugate transposed of A
static inline scimatrix transpherm ( const scimatrix& A ) {                                   
  scimatrix res(transp(A));
  
  std::vector<cinterval>& val = res.values();
  
  for(unsigned int i=0 ; i<val.size() ; i++)
    SetIm(val[i], -Im(val[i]));
  
  return res;
}

//! Compute conjugate transposed of A
static inline cmatrix transpherm ( const cmatrix& A ) {                                   
  cmatrix res(transp(A));
  
  for(int i=Lb(A,ROW) ; i<=Ub(A,ROW) ; i++)
    for(int j=Lb(A,COL) ; j<=Ub(A,COL) ; j++)  
      SetIm(A[i][j], -Im(A[i][j]));
  
  return res;
}

//! Are all elements of v==0.0?
template<typename TVec>
static inline bool isZeroVec(const TVec& v) {
  bool zero = true;
  int i=Lb(v);

  while(zero && i<=Ub(v)) {
    if(v[i] != 0.0) 
      zero = false;
    i++;
  }

  return zero;
}

/*!
/* The vectors x and y are successive approximations for the solution of a
/* linear system of equations computed by iterative refinement. If a compo-
/* nent of y is diminished by more than 'Factor', it is a good candidate for
/* a zero entry. Thus, it is set to zero.
*/
template<typename TVec>
static inline void CheckForZeros ( const TVec& x, const TVec& y )
{
  const real Factor = 1E+5;
  int        i;

  for (i = Lb(y); i <= Ub(y); i++)
    if ( abs(y[i])*Factor < abs(x[i]) ) y[i] = 0.0;
}

//! Return the lower Bandwidth of A
template<typename TMat>
static inline int lower_bandwidth(const TMat& A) {
  int bw = 0, tmp = 0;
  int n = RowLen(A);

  const std::vector<int>& p   = A.column_pointers();
  const std::vector<int>& ind = A.row_indices();
  for(int j=0 ; j<n ; j++) {
    if(p[j]!=p[j+1]) tmp = ind[p[j+1]-1] - j;
    if(tmp > bw) bw = tmp;
  }

  return bw;
}

//! Return the upper Bandwidth of A
template<typename TMat>
static inline int upper_bandwidth(const TMat& A) {
  int bw = 0, tmp = 0;
  int n = RowLen(A);

  const std::vector<int>& p   = A.column_pointers();
  const std::vector<int>& ind = A.row_indices();

  for(int j=0 ; j<n ; j++) {
    if(p[j]!=p[j+1]) tmp = j - ind[p[j]];
    if(tmp > bw) bw = tmp;
  }

  return bw;
}

//! Epsilon inflation of x
civector blow(const civector& x, const real eps);
//! Epsilon inflation of x
ivector  blow(const ivector& x, const real eps);

//! Compute upper bound for the infinity norm of A
real Norm00(const srmatrix& A);
//! Compute upper bound for the infinity norm of A
real Norm00(const scmatrix& A);
//! Compute upper bound for the 1-norm of A
real Norm1(const srmatrix& A);
//! Compute upper bound for the 1-norm of A
real Norm1(const scmatrix& A);
//! Compute upper bound for the 2-norm of A
real Norm2(const rvector& a);
//! Compute upper bound for the 2-norm of A
real Norm2(const cvector& a);
//! Compute upper bound for the 2-norm of A
interval Norm2(const ivector& a);
//! Compute upper bound for the 2-norm of A
interval Norm2(const civector& a);

//! Compute upper bound for row-sum norm of a matrix A
inline real Norm00(const rmatrix &A)
{
  int i,j;
  real m,s;

  fesetround(FE_UPWARD);

  m = 0.0;
  for(i = Lb(A,ROW); i <= Ub(A,ROW); i++) 
  {
    s = 0.0;    
    for(j = Lb(A,COL); j <= Ub(A,COL);j++) {   
      s += abs(A[i][j]);
    }
    m = Max(m,s);
  }
  
  fesetround(FE_TONEAREST);
  
  return m;
} // end norm

//! Compute upper bound for row-sum norm of a matrix A
inline real Norm00(const cmatrix &A)
{
  int i,j;
  real m,s; 

  fesetround(FE_UPWARD);
  
  m = 0.0;
  for(i = Lb(A,ROW); i <= Ub(A,ROW); i++) 
  {
    s = 0.0;
    for(j = Lb(A,COL); j <= Ub(A,COL);j++) {   
      s += sqrt( Re(A[i][j])*Re(A[i][j]) + Im(A[i][j])*Im(A[i][j]) );  //abs(A[i][j]);
    }
    m = Max(m,s);
  }
  
  fesetround(FE_TONEAREST);
  
  return m;
} // end norm

//! Compute upper bound for row-sum norm of a matrix A
inline real Norm00(const imatrix &A)
{
  int i,j;
  real m,s;
  dotprecision Accu; 

  m = 0.0;
  for(i = Lb(A,ROW); i <= Ub(A,ROW); i++) 
  {
    Accu = 0.0; 
    for(j = Lb(A,COL); j <= Ub(A,COL); j++) 
      Accu += Sup( abs(A[i][j]) );
    s = rnd( Accu,RND_UP );
    m = Max(m,s);
  }
  return m;
} // end norm

//! Compute row-sum norm of A
inline real Norm00(const cimatrix &A)
{
  int i,j;
  real m,s;
  dotprecision Accu; 

  m = 0.0;
  for(i = Lb(A,ROW); i <= Ub(A,ROW); i++) 
  {
    Accu = 0.0; 
    for(j = Lb(A,COL); j <= Ub(A,COL); j++) 
      Accu += Sup( abs(A[i][j]) );
    s = rnd( Accu,RND_UP );
    m = Max(m,s);
  }
  return m;
} // end norm

//! Compute diagonal D1 and D2 for binormalization of A
void normbin(const srmatrix& A, srmatrix& D1, srmatrix& D2);
//! Compute diagonal D1 and D2 for binormalization of A
void normbin(const scmatrix& A, srmatrix& D1, srmatrix& D2);

//! High precision computation of def=b-A*xapp
void residual(rvector& def, const rvector& b, const srmatrix& A, const rmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
void residual(ivector& def, const rvector& b, const srmatrix& A, const rmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
void residual(ivector& def, const ivector& b, const srmatrix& A, const rmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
//! High precision computation of def=b-A*xapp
void residual(ivector& def, const ivector& b, const simatrix& A, const rmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
void residual(cvector& def, const cvector& b, const scmatrix& A, const cmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
void residual(civector& def, const cvector& b, const scmatrix& A, const cmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
void residual(civector& def, const civector& b, const scmatrix& A, const cmatrix& xapp, int cols, bool dotp = true);
//! High precision computation of def=b-A*xapp
void residual(civector& def, const civector& b, const scimatrix& A, const cmatrix& xapp, int cols, bool dotp = true);

//! Parallel computation of C=|L||U|
void absLU(const srmatrix& L, const srmatrix& U, srmatrix& C);
//! Parallel computation of C=|L||U|
void absLU(const scmatrix& L, const scmatrix& U, srmatrix& C);

//! Parallel computation of upper bound D=LU-A
void LUMinusA(srmatrix& D, const srmatrix& L, const srmatrix& U, const srmatrix& A);
//! Parallel computation of upper bound D=LU-A
void LUMinusA(scmatrix& D, const scmatrix& L, const scmatrix& U, const scmatrix& A);
//! Parallel computation of D=LU-A
void LUMinusA(const srmatrix& L, const srmatrix& U, const srmatrix& A, srmatrix& D_inf, srmatrix& D_sup);
//! Parallel computation of D=LU-A
void LUMinusA(const scmatrix& L, const scmatrix& U, const scmatrix& A, scmatrix& D_inf, scmatrix& D_sup);
//! Parallel computation of upper bound D=LU-A
void LUMinusA(const srmatrix& L, const srmatrix& U, const simatrix& A, simatrix& D);
//! Parallel computation of upper bound D=LU-A
void LUMinusA(const scmatrix& L, const scmatrix& U, const scimatrix& A, scimatrix& D);

//! Parallel computation of C=LU-A in higher precision
void LUMinusAHighPrec(const srmatrix& L, const srmatrix& U, const srmatrix& A, srmatrix& C);
//! Parallel computation of C=LU-A in higher precision
void LUMinusAHighPrec(const scmatrix& L, const scmatrix& U, const scmatrix& A, scmatrix& C);

//! Parallel computation of [C,D]=[A*B]
void spMatMulPar(const srmatrix& A, const srmatrix& B, srmatrix& C, srmatrix& D);
//! Parallel computation of [C,D]=[A*B]
void spMatMulPar(const scmatrix& A, const scmatrix& B, scmatrix& C, scmatrix& D);

} //namespace cxsc
#endif