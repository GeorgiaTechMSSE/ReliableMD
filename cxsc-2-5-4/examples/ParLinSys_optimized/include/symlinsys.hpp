//============================================================================
//
//           Program/Module "Symmetric Linear Interval System Solver"
//
//	       supplement to the C-XSC Toolbox for Verified Computing
//
// Author: Michael Zimmer
//
// This program/module is free software for non-commercial use.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================

#ifndef __SYMLINSYS_HPP
#define __SYMLINSYS_HPP

#include "parlinsys.hpp"
#include <cimatrix.hpp>	   // Real matrix/vector arithmetic
#include <vector>

namespace cxsc {
  
//This function returns an error message corresponding to the error code returned by the solver
string SymLinSolveErrMsg ( int );

//Solver for symmetric linear interval systems
//The system Ax=b is transformed into an equivalent parametric linear system and solved with the parametric linear system solver
//(see header parlinsys.hpp)
//Parameters for all versions:
//A: system matrix
//b: Right hand side(s)
//x: Computed outer enclosure of solution
//y: (optional) computed inner enclosure of solution
//Err: error code (0 means no error). Meaning of error code can be queried by calling SymLinSolveErrMsg(Err)
//cfg: (optional) Configuration struct (the same as for the parametric solver, see parlinsys.hpp)
void  SymLinSolve (cxsc::imatrix& A, cxsc::ivector& b, cxsc::ivector& x, cxsc::ivector& y, int& Err, struct parlinsysconfig cfg = parlinsysconfig());
void  SymLinSolve (cxsc::imatrix&, cxsc::ivector&, cxsc::ivector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  SymLinSolve (cxsc::imatrix&, cxsc::ivector&, cxsc::imatrix&, cxsc::imatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  SymLinSolve (cxsc::imatrix&, cxsc::ivector&, cxsc::imatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());

void  SymLinSolve (cxsc::cimatrix&, cxsc::civector&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  SymLinSolve (cxsc::cimatrix&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  SymLinSolve (cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  SymLinSolve (cxsc::cimatrix&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());

} //namespace cxsc

#endif 
 
