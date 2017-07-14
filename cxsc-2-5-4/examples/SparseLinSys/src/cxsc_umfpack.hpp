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

/* Header providing interface functions to use UMFPACK with C-XSC data types */

#ifndef _CXSC_UMFPACK_HEADER
#define _CXSC_UMFPACK_HEADER

#include <scmatrix.hpp>
#include <intvector.hpp>

namespace cxsc { 

/*!  Compute LU-decomposition of A and fill-in reducing permutation p and q such that
 *  L*U=A(p,q). The error code err is 0 if no error occured and !=0 otherwise. */
void lu_decomp(const srmatrix& A, srmatrix& L, srmatrix& U, intvector& p, intvector& q, int& Err);

/*!  Compute complex LU-decomposition of A and fill-in reducing permutation p and q such that
 *  L*U=A(p,q). The error code err is 0 if no error occured and !=0 otherwise. */
void lu_decomp(const scmatrix& A, scmatrix& L, scmatrix& U, intvector& p, intvector& q, int& Err);

} //namespace cxsc

#endif