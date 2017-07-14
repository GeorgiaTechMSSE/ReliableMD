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

/* General header file for the sparse solvers
 * There are two different versions:
 * slss - Solve system using algorithm based on Krawczyk-Operator
 *        Leads to componentwise good error bounds, range of solvable systems 
 *        is limited (primarily systems with M-matrices or diagonally dominant matrices)
 * 
 * slssn - Solve system using algorithm based on verified normwise error bound
 *         Applicable for many more systems than slss, especially useful for systems
 *         with s.p.d. matrices, but also useful for general systems if condition
 *         is moderate. Computes a normwise error bounds, enclosures might be
 *         qide for interval systems.
 */

#ifndef _CXSC_SPARSELINSYS_HPP
#define _CXSC_SPARSELINSYS_HPP

#include <srmatrix.hpp>  
#include <sivector.hpp>  

namespace cxsc {

/*! Enumeration for method of solving triangular systems with interval right hand side
 *  occuring in solver slss. The choice of method has a huge effect on run-time
 *  and applicability of the solver, see below.
 *  Methods:
 *  -FORWARD_BACKWARD: Solve triangular systems with forward/backward substitution.
 *                     Useful especially for diagonally dominant or M-matrices.
 *  -IMPLICIT_INVERSE: Solve triangular systems by implicitly computing the inverse.
 *                     Very slow in general, but should work for most systems.
 *  -BANDED:           Solve triangular system with forard/backward substitution using
 *                     a coordinate transformation. Useful for banded systems with small
 *                     bandwidth (<10).
*/
enum triSysMethod { FORWARD_BACKWARD, IMPLICIT_INVERSE, BANDED };

/*! Returns the error message corresponding to an error code returned by the solvers */
string SparseLinSolveErrMsg ( int );


/*! Configuration struct for solver slss. If you want to change the default configuration,
 * create an instance of this struct and overwrite the value you want to change. */
struct slssconfig {
  int   precResidual;              /*! Dot product precision for computation of residual */
  int   precIterMat;               /*! Dot product precision for computation of [C] */
  bool  msg;                       //! Status message output during solver run? The default is false
  int   threads;                   /*! Number of threads to use for OpenMP (only has an effect if OpenMP ist activated during compilation).
                                    *  Using a value <=0 means use all available threads (set for example by OMP_NUM_THREADS environment variable).
			            *  The default is -1.
			            */
  int   maxIterVer;                /*! Maximum number of iterations during the verification step.
                                    *  For systems with condition numbers approaching 10^15 for part one or 10^30 for part two
                                    *  higher than default settings might be necessary to find a verified result. The default is 5.
			            */
  real  epsVer;                    //! Epsilon for the verification step. This is used for the Epsilon inflation during verification. The default is 0.1.
  triSysMethod   triMethod;        //! Solving method for triangular systems (see above)
  bool  forceSuiteSparse;          //! Force usage of UMFPACK/CHOLMOD is triMethod==BANDED
  bool  forceHighPrecIterMatrix;   /*! Force usage of high precision computation of LU-A if triMethod!=BANDED. Should only be used for banded systems with
                                   /*  small bandwidth */
  bool  fastVerification;          //! Use faster verification method (in general leads to wider enclosures)
  bool  highPrecApprxSolution;     //! High precision approximate solution, leads to very accurate result especially for point systems at little cost
  int   apprxSolutionStagPrec;     //! Precision to use if highPrecApprxSolution==true
  bool  storeBase;                 //! If triMethod==BANDED, activates storage of occuring basis matrices. Much faster, but needs more memory.
  bool  fastDecompErrorComputation;//! Faster computation of iteration matrix [LU-A]. Much faster, but leads to wider enclosures and might hinder verification.
  bool  spd;                       //! Set to true if the system is symmetric/hermitian positive definite

  //! Default constructor setting all configuration variables to their default settings
  slssconfig() : precResidual(2),precIterMat(1),msg(false),threads(-1),maxIterVer(5), epsVer(0.1), triMethod(FORWARD_BACKWARD), 
                 forceSuiteSparse(false), fastVerification(false), highPrecApprxSolution(true), apprxSolutionStagPrec(3), storeBase(false),
		 fastDecompErrorComputation(false), spd(true), forceHighPrecIterMatrix(false) {}
}; 

/*
 * Start functions for solver based on the Krawczyk-Operator. Computes verified enclosure of the exact
 * solution of the sparse linear system Ax=b.
 * Parameters:
 * -A: System matrix
 * -b: Right hand side. Use a matrix type for multiple right hand sides.
 * -x: Interval vector/matrix containing verified enclosure of the true solution
 * -err: An error code (==0 means no error, !=0 means an error occured)
 * -cfg: (optional) Pass your own configuration if you do not want to use the standard configuration
 */
void  slss ( const srmatrix& A, const rmatrix& b, imatrix& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const srmatrix& A, const rvector& b, ivector& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const srmatrix& A, const imatrix& b, imatrix& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const srmatrix& A, const ivector& b, ivector& x, int& err, slssconfig cfg = slssconfig() );

void  slss ( const simatrix& A, const imatrix& b, imatrix& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const simatrix& A, const ivector& b, ivector& x, int& err, slssconfig cfg = slssconfig() );

void  slss ( const scmatrix& A, const cmatrix& b, cimatrix& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const scmatrix& A, const cvector& b, civector& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const scmatrix& A, const cimatrix& b, cimatrix& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const scmatrix& A, const civector& b, civector& x, int& err, slssconfig cfg = slssconfig() );

void  slss ( const scimatrix& A, const cimatrix& b, cimatrix& x, int& err, slssconfig cfg = slssconfig() );
void  slss ( const scimatrix& A, const civector& b, civector& x, int& err, slssconfig cfg = slssconfig() );


/*! Configuration struct for solver slssn. If you want to change the default configuration,
 * create an instance of this struct and overwrite the value you want to change. */
struct slssnconfig {
  bool  msg;                    //! Status message output during solver run? The default is false
  int   threads;                /*! Number of threads to use for OpenMP (only has an effect if OpenMP ist activated during compilation).
                                 *  Using a value <=0 means use all available threads (set for example by OMP_NUM_THREADS environment variable).
                                 *  The default is -1.
                                 */
  int   precDecompError;        /*! Dot product precision for verification step. Might help verification for systems wit higher condition, but is slower */
  int   maxIterPower;           /*! Maximum number of iterations for inverse power iteration */
  int   apprxSolutionStagPrec;  //! Precision for computation of approximate solution
  bool  highPrecApprxSolution;  //! Compute high precision approximate solution
  bool  scaleMatrix;            //! Apply scaling of system matrix?
  bool  spd;                    //! Set to true if the system is symmetric/hermitian positive definite
  

  //! Default constructor setting all configuration variables to their default settings
  slssnconfig() : msg(false),threads(-1),precDecompError(1),maxIterPower(4), highPrecApprxSolution(true),
                  apprxSolutionStagPrec(3), scaleMatrix(true), spd(false) {}
}; 

/*
 * Start functions for solver based on a verified normwise error bound. Computes verified enclosure of the exact
 * solution of the sparse linear system Ax=b.
 * Parameters:
 * -A: System matrix
 * -b: Right hand side. Use a matrix type for multiple right hand sides.
 * -x: Interval vector/matrix containing verified enclosure of the true solution
 * -err: An error code (==0 means no error, !=0 means an error occured)
 * -cfg: (optional) Pass your own configuration if you do not want to use the standard configuration
 */

void slssn(srmatrix& A, rvector& b, ivector& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(srmatrix& A, ivector& b, ivector& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(simatrix& A, ivector& b, ivector& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(srmatrix& A, rmatrix& b, imatrix& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(srmatrix& A, imatrix& b, imatrix& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(simatrix& A, imatrix& b, imatrix& x, int& err, slssnconfig cfg = slssnconfig());

void slssn(scmatrix& A, cvector& b, civector& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(scmatrix& A, civector& b, civector& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(scimatrix& A, civector& b, civector& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(scmatrix& A, cmatrix& b, cimatrix& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(scmatrix& A, cimatrix& b, cimatrix& x, int& err, slssnconfig cfg = slssnconfig());
void slssn(scimatrix& A, cimatrix& b, cimatrix& x, int& err, slssnconfig cfg = slssnconfig());


} //namespace cxsc

#endif




