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

/* Header providing interface functions to use CHOLMOD with C-XSC data types */

#ifndef _CXSC_CHOLMOD_HEADER
#define _CXSC_CHOLMOD_HEADER

#include <scmatrix.hpp>
#include <intvector.hpp>

namespace cxsc {

/*! Compute Cholesky decomposition of A and fill-in reducing permutation p such that
 *  L*L'=A(p,p). The error code err is 0 if no error occured and !=0 otherwise. */  
void chol(srmatrix& A, srmatrix& L, intvector& p, int& err);
/*! Compute complex Cholesky decomposition of A and fill-in reducing permutation p such that
 *  L*L'=A(p,p). The error code err is 0 if no error occured and !=0 otherwise. */
void chol(scmatrix& A, scmatrix& L, intvector& p, int& err);


} //namespace cxsc

#endif