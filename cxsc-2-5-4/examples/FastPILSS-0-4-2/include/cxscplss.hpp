/*
**  FastPLSS: A library of (parallel) verified linear (interval) system 
**  solvers using C-XSC (V 0.4.2)
**
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2012 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**  
**/

/*
**  Author: Michael Zimmer
**
**  This software is based on:
**    - Module LinSys of the C-XSC-Toolbox
**      Authors: Rolf Hammer, Matthias Hocks, Dietmar Ratz
**    - Self-verifying solver for a dense system of linear equations
**      Authors: Carlos Holbig, Walter Kraemer, Paulo Sergio Morandi Junior,
**               Bernardo Frederes Kramer Alcalde, 
**    - Parallel Interval Linear System Solver
**      Author: Markus Grimmer
**/
 

#ifndef _TEMPL_LSS_PAR_HPP
#define _TEMPL_LSS_PAR_HPP

#include <iostream> 
#include <fstream> 
#include <string>
#include <rmatrix.hpp>   
#include <imatrix.hpp>   
#include <cmatrix.hpp>   
#include <cimatrix.hpp>   
#include <rvector.hpp>
#include <ivector.hpp>
#include <cvector.hpp>
#include <civector.hpp>

#ifdef Unchanged
  #define sl_init_ sl_init
  #define blacs_gridinfo_ blacs_gridinfo
  #define numroc_ numroc
#elif defined(Add_)
  //Nothing to do
#elif defined(Uppercase)
  #define sl_init_ SL_INIT
  #define blacs_gridinfo_ BLACS_GRIDINFO
  #define numroc_ NUMROC
#endif

//Declaration of SCALAPACK routines
extern "C"  { 
   void sl_init_(int *ic, int *nr, int *nc);
   void blacs_gridinfo_(int *ic, int *nr, int *nc, int *mr, int *mc);       
   int numroc_ (int *n, int *nb, int *iproc, int *isrcproc, int *nprocs); 
} 

//Control constants
const int
  LSS_ONLY_PART_ONE = 0,
  LSS_ONLY_PART_TWO = 1,
  LSS_BOTH_PARTS    = 2;

//Function pointer type for distribution of matrix A
typedef void (*get_Element_of_rmatrix)(int,int,cxsc::real&);
typedef void (*get_Element_of_imatrix)(int,int,cxsc::interval&);
typedef void (*get_Element_of_cmatrix)(int,int,cxsc::complex&);
typedef void (*get_Element_of_cimatrix)(int,int,cxsc::cinterval&);

/*!
 * Configuration struct. The configuration of the solver can be changed by
 * creating an instance of this struct, changing the configuration options as needed
 * and passing it as an additional argument to the solver.
 */
struct plssconfig {
  int   K;              /*! Dot product precision for computation of the residual.
                            K=0: maximum precison using long accumulator (slow)
                            K=1: floating point computations, computes residual with splitting algorithm by Rump/Dekker
                            K>1: simulated K-fold double precision using DotK algorithm
                            Default is K=2. Selecting K=1 forces use of matrix mode (see below). */
  int   lssparts;       /*! Solver stages to use. The used algorithm hast two parts: one using a normal approximate inverse 
                            as preconditioner, one using an inverse R1+R2 of double length. The first part is much faster and works
                            for condition numbers of up to about 10^15. The second is slower, but can solve systems with higher 
			    condition numbers (up to about 10^30).
                            Possible values are:
                            LSS_ONLY_PART_ONE: Use only part one (automatically selected for K==1)
                            LSS_ONLY_PART_TWO: Use only part two
                            LSS_BOTH_PARTS: Try part one first, if it fails try part two (this is the default). */
  int   threads;        /*! Sets the number of threads used by OpenMP. Some BLAS libraries use this for thread-selection.
                            This setting has currently no further effect. A value <=0 sets the number of threads to all
                            available threads for OpenMP. The default is -1. */
  int   maxIterResCorr; /*! Maximum number of iteration steps during residual correction. This option has no effect for K==1. 
                            The default is 5.  */
  int   maxIterVer;     /*! Maximum number of iteration steps during the verification step. For systems with high condition numbers,
                            increasing the default value can sometimes help to find a verified solution. The default is 5. */
  bool  refinement;     /*! If set to true, a refinement step is added after finding a verified enclosure. This can sometimes improve
                            the solution at small cost. The default is false. */
  int   maxIterRef;     /*! Maximum number of iterations during the refinement step. The default is 5. */
  cxsc::real  epsVer;   /*! The Epsilon used for the epsilon inflation during the the verification step. The default is 0.1 */
  cxsc::real  epsRef;   /*! Sets the stopping criterion for the refinement step. If two successive iterates do not differ by more than epsRef,
                            the refinement iteration is stopped. The default is 1e-5. */
  int   nb;             /*! The blocksize used by ScaLAPACK. This setting can have a big effect on the performance of the solver. The optimal 
                            setting is system dependent. You can experiment with this setting to find a block size giving good performance.
                            The default value of 256 should give good results on modern machines. */
  bool  matrixMode;     /*! This mode can only be used for K==1. If activated for K!=1, K will automatically be set to one. In Matrix Mode, only 
                            part one of the algorithm can be used and no inner enclosure can be computed. This mode handles the right hand side as
                            a matrix and is much faster for many right hand sides. If you want to solve a system with many right hand sides 
                            (for example) when computing a verified inverse), it is suggested to use this mode. */

  /*! Default constructor setting all options to their default value. */
  plssconfig() : K(2),lssparts(LSS_BOTH_PARTS),threads(-1),maxIterResCorr(5), maxIterVer(5), 
                 refinement(false), maxIterRef(5), epsVer(0.1), epsRef(1e-5), nb(256), matrixMode(false) {}
};

//! Translates the error codes of the solver into corresponding error messages
std::string LinSolveErrMsg(int);

/** Entry Function to verified solver for linear real systems
 *  Computes the verified solution of the system AX=B
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void plss(cxsc::rmatrix& A, cxsc::rmatrix& B, cxsc::imatrix& X, int m, int n, 
          int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, Y will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y: Inner enclosure of solution
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void plss(cxsc::rmatrix& A, cxsc::rmatrix& B, cxsc::imatrix& X, cxsc::imatrix& Y, int m, int n, 
          int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real systems
 *  Computes the verified solution of the system AX=B
 *  A and B are defined by pointer to a function
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n : Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void plss(get_Element_of_rmatrix A, get_Element_of_rmatrix B, cxsc::imatrix& X, int m, int n, int rhs,
          int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are defined by pointer to a function
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X and Y will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y : inner enclosure
 *  m,n : Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void plss(get_Element_of_rmatrix A, get_Element_of_rmatrix B, cxsc::imatrix& X, cxsc::imatrix& Y, int m, int n, int rhs,
          int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real interval systems
 *  Computes the verified solution of the system AX=B
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pilss(cxsc::imatrix& A, cxsc::imatrix& B, cxsc::imatrix& X, int m, int n,
          int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real interval systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, Y will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y: Inner enclosure
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pilss(cxsc::imatrix& A, cxsc::imatrix& B, cxsc::imatrix& X, cxsc::imatrix& Y, int m, int n,
          int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real interval systems
 *  Computes the verified solution of the system AX=B
 *  A and B are defined by pointer to a function
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n: Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pilss(get_Element_of_imatrix A, get_Element_of_imatrix B, cxsc::imatrix& X, int m, int n, int rhs,
           int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear real interval systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are defined by pointer to a function
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X and Y will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n: Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pilss(get_Element_of_imatrix A, get_Element_of_imatrix B, cxsc::imatrix& X, cxsc::imatrix& Y, int m, int n, int rhs,
           int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex systems
 *  Computes the verified solution of the system AX=B
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pclss(cxsc::cmatrix& A, cxsc::cmatrix& B, cxsc::cimatrix& X, int m, int n,
           int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, Y will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y: Inner enclosure
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pclss(cxsc::cmatrix& A, cxsc::cmatrix& B, cxsc::cimatrix& X, cxsc::cimatrix& Y, int m, int n,
           int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex systems
 *  Computes the verified solution of the system AX=B
 *  A and B are defined by pointer to a function
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n : Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pclss(get_Element_of_cmatrix A, get_Element_of_cmatrix B, cxsc::cimatrix& X, int m, int n, int rhs,
           int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are defined by pointer to a function
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X and Y will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y: Inner enclosure
 *  m,n : Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pclss(get_Element_of_cmatrix A, get_Element_of_cmatrix B, cxsc::cimatrix& X, cxsc::cimatrix& Y, int m, int n, int rhs,
           int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);


/** Entry Function to verified solver for linear complex interval systems
 *  Computes the verified solution of the system AX=B
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pcilss(cxsc::cimatrix& A, cxsc::cimatrix& B, cxsc::cimatrix& X, int m, int n, 
            int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex interval systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are stored in process 0 at the beginning, will then be distributed and deleted
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, Y will be stored completely by process 0.
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y: Inner enclosure
 *  m,n : Dimension of system
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pcilss(cxsc::cimatrix& A, cxsc::cimatrix& B, cxsc::cimatrix& X, cxsc::cimatrix& Y, int m, int n, 
            int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex interval systems
 *  Computes the verified solution of the system AX=B
 *  A and B are defined by pointer to a function
 *  X will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  m,n: Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pcilss(get_Element_of_cimatrix A, get_Element_of_cimatrix B, cxsc::cimatrix& X, int m, int n, int rhs,
            int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

/** Entry Function to verified solver for linear complex interval systems
 *  Computes the verified solution of the system AX=B and an inner enclosure Y of the solution
 *  A and B are defined by pointer to a function
 *  X and Y will be cyclically distributed  by columns among the processes when not using matrix mode.
 *  When using matrix mode, X and Y will distributed in two dimensional block cyclic distribution according to
 *  selected block size (see ScaLAPACK reference for more information).
 *
 *  Parameters:
 *  A: System matrix
 *  B: Right hand side
 *  X: Solution
 *  Y: Inner enclosure
 *  m,n: Dimension of system
 *  rhs: Number of right hand sides
 *  procs: Number of processes used by MPI
 *  mypid: Own process ID
 *  errc: Error code (use function LinSolveErrMsg(int) to translate into error message)
 *  out: Output file stream for status messages per process
 *  cfg: Struct with configuration variables
 */
void pcilss(get_Element_of_cimatrix A, get_Element_of_cimatrix B, cxsc::cimatrix& X, cxsc::cimatrix& Y, int m, int n, int rhs,
            int procs, int mypid, int& errc, std::ofstream& out, struct plssconfig cfg);

#endif
 
