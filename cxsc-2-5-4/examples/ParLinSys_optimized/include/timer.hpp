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

#include "sys/time.h" 

//Helper function for time measurements
//Return current time value in seconds double format
inline double GetTime() {
   struct timeval _tp;

   gettimeofday(&_tp,0);
   
   return _tp.tv_sec + _tp.tv_usec / 1000000.0;
} 

