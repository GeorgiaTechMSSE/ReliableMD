//============================================================================
//
//		 Program/Module "Parametric Linear Solver"
//
//	       supplement to C++ Toolbox for Verified Computing
//
//		   Copyright (c) 2003   Evgenija Popova
//
// This program/module is free software for non-commercial use.
// For details on theory  see the papers:
//
// S. Rump: Verification methods for dense and sparse systems of equations,
// in: Topics in Validated Computations, ed. J. Herzberger, Elsevier Science
// B.V., 1994,	pp. 63-135.
//
// E. Popova, W. Kraemer: Parametric Fixed-Point Iteration Implemented in
// C-XSC. Universitaet Wuppertal, Preprint BUW - WRSWT, 2003.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================
//----------------------------------------------------------------------------
// File: parlinsys (header)
//
// Purpose: Computation of a verified inclusion for the parametric solution set
// of a square parametric linear system of equations A(p)*x = b(p) with full
// parametric matrix A(p) and parametric right-hand side b(p) which elements are
// affine-linear combinations of k given interval parameters.
//
// Global functions:
//
//    ParLinSolve()	 : to get a verified enclosure of the solution
//    ParLinSolveErrMsg(): to get an error message text
//
//----------------------------------------------------------------------------
#ifndef __PARLINSYS_HPP
#define __PARLINSYS_HPP

#include <rmatrix.hpp>	   // Real matrix/vector arithmetic
#include <ivector.hpp>	   // Interval vector arithmetic

using namespace cxsc;
using namespace std;

extern char* ParLinSolveErrMsg ( int );

extern void  ParLinSolve (rmatrix&, rmatrix&, ivector&,
				int&, real&, int&, ivector&, ivector&, int&);

#endif // end __PARLINSYS
