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

/* Implementation of functions for the solution of triangular point systems
   or triangular systems with interval right hand side */

#ifndef _CXSC_TRISOLVE_HEADER_DEFINED
#define _CXSC_TRISOLVE_HEADER_DEFINED

#include <scmatrix.hpp>
#include <civector.hpp>
#include <sivector.hpp>

namespace cxsc {

//! Solve lower triangular system Lx=b with simple forward substitution
template<typename TRhs, typename Tx>
void lowtrisolve(const srmatrix& L, const TRhs& b, Tx& x);

//! Solve upper triangular system Ux=b with simple backward substitution
template<typename TRhs, typename Tx>
void uptrisolve(const srmatrix& U, const TRhs& b, Tx& x);

//! Solve lower triangular system Lx=b with simple forward substitution
template<typename TRhs, typename Tx>
void lowtrisolve(const scmatrix& L, const TRhs& b, Tx& x);

//! Solve upper triangular system Ux=b with simple backward substitution
template<typename TRhs, typename Tx>
void uptrisolve(const scmatrix& U, const TRhs& b, Tx& x);


//! Solve lower triangular system Lx=b by implicitly computing the inverse of L
void lowtrisolve_inv(const cxsc::srmatrix& L, const cxsc::ivector& b, cxsc::ivector& x);
//! Solve upper triangular system Ux=b by implicitly computing the inverse of U
void uptrisolve_inv(const cxsc::srmatrix& U, const cxsc::ivector& b, cxsc::ivector& x);
//! Solve lower triangular system Lx=b by implicitly computing the inverse of L
void lowtrisolve_inv(const cxsc::scmatrix& L, const cxsc::civector& b, cxsc::civector& x);
//! Solve upper triangular system Ux=b by implicitly computing the inverse of U
void uptrisolve_inv(const cxsc::scmatrix& U, const cxsc::civector& b, cxsc::civector& x);

//! Solve lower triangular system Lx=b by forward substitution with coordinate transformation using parallel epipeds
void lowtrisolve_band( cxsc::srmatrix& L,  cxsc::srmatrix& LT,  cxsc::ivector& b, cxsc::ivector& x, int l, cxsc::imatrix& Basis, bool baseStore, int& err);
//! Solve lower triangular system Lx=b by forward substitution with coordinate transformation using parallel epipeds
void lowtrisolve_band( cxsc::scmatrix& L,  cxsc::scmatrix& LT,  cxsc::civector& b, cxsc::civector& x, int l, cxsc::cimatrix& Basis,  bool baseStore, int& err);
//! Solve upper triangular system Ux=b by backward substitution with coordinate transformation using parallel epipeds
void uptrisolve_band( cxsc::srmatrix& U,  cxsc::srmatrix& UT,  cxsc::ivector& b, cxsc::ivector& x, int l, cxsc::imatrix& Basis, bool baseStore, int& err);
//! Solve upper triangular system Ux=b by backward substitution with coordinate transformation using parallel epipeds
void uptrisolve_band( cxsc::scmatrix& U,  cxsc::scmatrix& UT,  cxsc::civector& b, cxsc::civector& x, int l, cxsc::cimatrix& Basis,  bool baseStore, int& err);

}

#endif