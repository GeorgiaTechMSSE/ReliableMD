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
#ifndef __PARLINSYS_SPARSECOEFF_HPP
#define __PARLINSYS_SPARSECOEFF_HPP

#include "CoeffMatrix.hpp"
#include <scimatrix.hpp>	
#include <civector.hpp>	   
#include <vector>

namespace cxsc {
  
//This function returns an error message corresponding to the error code returned by the solver
std::string ParLinSolveErrMsg ( int );

//struct for optional configuration options when calling the following solver functions
//To use non default configuration, create an instance of the struct, change the desired
//values and pass it to the solver as the last argument
struct parlinsysconfig {
  int   K;              //Dot product precision for residual computation (K=0 maximum precision(slow), 
                        //K=1 pure floating point, K>=2 simulated K-fold double precision). Default: K=2
  bool  msg;            //Status message output? Default: false
  int   threads;        //Number of threads to use for OpenMP. If threads<=0, the max. available number of
                        //threads is used. Default: -1
  int   maxIterResCorr; //Maximum number of iterations during residual correction (not available for K=1)
                        //Default: 10
  int   maxIterVer;     //Maximum number of iterations during the verification step. Default: 5
                        //Increasing this number can help yield an outer enclusion for systems where
                        //no enclusion could be reached before
  bool  refinement;     //Perform an iterative refinement? This can improve the result a little with modest cost.
                        //Default: false
  int   maxIterRef;     //Maximum number of iterations during the refinement step. Default: 5
  real  epsVer;         //Epsilon for the verification step (for the epsilon inflation). Default: 0.1
  real  epsRef;         //Epsilon for the refinement step (stopping criterion). Default: 1e-5
  bool  SharpC;         //Sharp (but slower) enclosure of iteration matrix [C]? Sharp enclosure takes a lot more 
                        //time in general, but increases the class of solvable systems and leads to much better results
                        //Default: true
  //Default constructor, setting all values to the default values described above                        
  parlinsysconfig() : K(2),msg(false),threads(-1),maxIterResCorr(10), maxIterVer(5), 
                      refinement(false), maxIterRef(5), epsVer(0.1), epsRef(1e-5), SharpC(true) {}
};

//Versions using the standard algorithm.
//Parameters: 
//Apv: STL-vector of coefficient matrices for each parameter. CoeffMatrix is a wrapper class that can hold either a sparse or a full matrix
//bp:  Matrix of right hand sides. Each column is the right hand side for the according parameter. For more than one right hand side, the first p+1 columns are for right hand side 1, 
//     the next p+1 columns for right hand side 2 etc.
//ip:  Vector of parameter intervals, size p+1. The first entry corresponds to the parameter independent part of the system and should normally be 1.0
//x:   Computed outer enclosure of the system
//y:   (Optional) Computed inner enclosure of the system
//err: An error code (0 meaning "No error"). MEaning of the error code can be queried calling ParLinSolveErrMsg(Err)
//cfg: (Optional) Configuration struct (see above)
//See also the example programs for help on calling the solver
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>& Apv, cxsc::srmatrix& bp, cxsc::ivector& ip, cxsc::ivector& x, cxsc::ivector& y, int& err, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::ivector&, cxsc::ivector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>&, cxsc::scmatrix&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::civector&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::ivector&, cxsc::imatrix&, cxsc::imatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::ivector&, cxsc::imatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>&, cxsc::scmatrix&, cxsc::civector&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::civector&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::srmatrix&, cxsc::civector&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::ivector&, cxsc::ivector&, cxsc::ivector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::ivector&, cxsc::ivector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>&, cxsc::cmatrix&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::civector&, cxsc::civector&, cxsc::civector&, int&, struct parlinsysconfig cfg = parlinsysconfig());

void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::ivector&, cxsc::imatrix&, cxsc::imatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::ivector&, cxsc::imatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffCMatrix>&, cxsc::cmatrix&, cxsc::civector&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::civector&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());
void  ParLinSolve (std::vector<cxsc::CoeffMatrix>&, cxsc::rmatrix&, cxsc::civector&, cxsc::cimatrix&, cxsc::cimatrix&, int&, struct parlinsysconfig cfg = parlinsysconfig());


//============================================================================================================================================================================

//Version using alternative algorithm (see Neumaier/Pownuk paper mentioned above)

//struct for optional configuration options when calling the solver
struct parlinsysnpconfig {
  int   K;              //Precision for all operations (for possible values see above). K==1 should be sufficient for practical applications. Default: 1
  bool  msg;            //Status message output? Default: false
  int   threads;        //Number of threads for OpenMP. If <=0, max. available number of threads is used. Default: -1
  int   maxIter;        //maximum number of iteration steps for final iteration. Increasing this number could improve the result, especially for 
                        //low stopping criterion epsilons. Default: 20
  real  epsIter;        //Epsilon for stopping criterion for the final iteration. Iteration stops if difference between successive iterates < Epsilon. Default: 1e-6
  bool  verified;       //Compute verified result? Enabling this is much slower, but results are guaranteed enclosures (including rounding errors). Default: false
  
  //Default constructor using default values described above
  parlinsysnpconfig() : K(1),msg(false),threads(-1),maxIter(20), epsIter(1e-6), verified(false) {}
};

//Versions solving systems of the form (A^T * D * A)u = Fb
//Interval entries only in D and b, D should be a diagonal matrix
//Parameters A,D,u,F,b correspond to forumla above
//Err: error code (0 means no error). Meaning of error codes can be parse by calling ParLinSolveErrMsg(Err)
//cfg: (optional) Configuration struct, see above
void  ParLinSolve(const srmatrix& A, const simatrix& D, ivector& u, const srmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg=parlinsysnpconfig());
void  ParLinSolve(const rmatrix& A, const simatrix& D, ivector& u, const rmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg=parlinsysnpconfig());

//Versions solving systems of the form (K + B * D * A)u = a + Fb
//Interval entries only in D and b, D should be a diagonal matrix
//Parameters K,B,A,D,u,a,F,b correspond to forumla above
//Err: error code (0 means no error). Meaning of error codes can be parse by calling ParLinSolveErrMsg(Err)
//cfg: (optional) Configuration struct, see above
void  ParLinSolve(const srmatrix& K, const srmatrix& B, const simatrix& D, const srmatrix& A, ivector& u, 
                  const rvector& a, const srmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg=parlinsysnpconfig());
void  ParLinSolve(const rmatrix& K, const rmatrix& B, const simatrix& D, const rmatrix& A, ivector& u, 
                  const rvector& a, const rmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg=parlinsysnpconfig());

} //namespace cxsc
#endif // end __PARLINSYS
