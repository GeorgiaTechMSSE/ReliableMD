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
**  
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
 

#include <mpi.h>
#include <iostream> 
#include <fstream> 
#include <iomanip>
#include <string>
#include <mvi_util.hpp>  

#ifdef _OPENMP
#include <omp.h>
#endif

#include <dot.hpp>
#include <idot.hpp>
#include <cdot.hpp>
#include <cidot.hpp>

#include "cxsc_mpicomm.hpp"
#include "matinv_aprx_par.hpp"    
#include "matmul_par.hpp"         
#include "cxsc_pblas.hpp"
#include "templplss.hpp"
#include <except.hpp>

#include "sys/time.h"

using namespace cxsc;
using namespace std;


/**
 * @name Constants
 */
//@{
const real zerotest = 1e6;
const real delta    = 1e-15;
const real eps1     = 1e-15;
const real limit    = 1e20; 
// Error codes
const int
  NoError      = 0,   // No error occurred
  NotSquare    = 1,   // System to be solved is not square
  DimensionErr = 2,   // Dimensions of A and b are not compatible
  InvFailed    = 3,   // System is probably singular
  VerivFailed  = 4,   // Verification failed, system is probably ill-conditioned
  DistFailed   = 5,   // Data Distribution failed
  CommErr      = 6;   // Communication Error
//@}


double time() {
   struct timeval _tp;

   gettimeofday(&_tp,0);
   
   return _tp.tv_sec + _tp.tv_usec / 1000000.0;
}

// ----------------------------------------------------------------------------

/*!
  This function translates the error codes used in the solver into
  a corresponding textmessage, which can be used as feedback for the
  user.
  \param Err The error code
  \return A textmessage corresponding to the error code
*/
string LinSolveErrMsg ( int Err )
{
  string ret;

  if (Err != NoError) {

    switch (Err) {
      case NotSquare:
        ret = "System to be solved is not square"; break;
      case DimensionErr:
        ret = "Dimensions of A and b are not compatible or Lb(A,1)!=Lb(A,2)"; break;
      case InvFailed:
        ret = "System is probably singular"; break;
      case VerivFailed:
        ret = "Verification failed, system is probably ill-conditioned";
        break;
      case DistFailed:
        ret = "Data Distribution failed";
        break;
      case CommErr:
        ret = "Communication Error";
        break;
      default:
        ret = "Error Code not defined";
    }
  }
  return ret;
} // LinSolveErrMsg

//----------------------------------------------------------------------------                      

//! Computes the transpose of a real matrix
/*
  \param A Matrix that is to be transposed
  \return The transposed of A
*/
inline rmatrix transpherm ( const rmatrix& A )
{                                            
  rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) 
      res[j][i] = A[i][j];
      
  return res;
}

//----------------------------------------------------------------------------                      

//! Computes the transpose of an interval matrix
/*
  \param A Matrix that is to be transposed
  \return The transposed of A
*/
inline imatrix transpherm ( const imatrix& A )
{                                            
  imatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) 
      res[j][i] = A[i][j];
      
  return res;
}

//----------------------------------------------------------------------------                      

//! Computes the hermitian of a complex matrix
/*
  \param A Matrix for which the hermitian is to be computed
  \return The hermitian of A
*/
inline cmatrix transpherm ( const cmatrix& A )
{                                            
  cmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) {
      res[j][i] = A[i][j];
      SetIm(res[j][i], -Im(res[j][i]));
    }
  
  return res;
}

//----------------------------------------------------------------------------                      

//! Computes the hermitian of a complex interval matrix
/*
  \param A Matrix for which the hermitian is to be computed
  \return The hermitian of A
*/
inline cimatrix transpherm ( const cimatrix& A )
{                                            
  cimatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  for (int i = Lb(A,1); i <= Ub(A,1); i++) 
    for (int j = Lb(A,2); j <= Ub(A,2); j++) {
      res[j][i] = A[i][j];
      SetIm(res[j][i], -Im(res[j][i]));
    }
  
  return res;
}

//----------------------------------------------------------------------------                      

//! Determines if all elements of an interval vector are 0
/*
  Determines if all elements of an interval vector are equal to 0
  \param zz Interval vector whose elements are to be checked
  \return true, if all elements of zz are zero
*/
inline bool isZero(ivector &zz) {
  return (Inf(zz) == 0 && Sup(zz) == 0);
}

//----------------------------------------------------------------------------                      

//! Determines if all elements of a complex interval vector are 0
/*
  Determines if all elements of a complex interval vector are equal to 0
  \param zz Interval vector whose elements are to be checked
  \return true, if all elements of zz are zero
*/
inline bool isZero(civector &zz) {
  return (InfRe(zz) == 0 && InfIm(zz) == 0 && SupRe(zz) == 0 && SupIm(zz) == 0);
}

//----------------------------------------------------------------------------                      

/**
 * RelativeError computes componentwise the maximum RelativeErrorative error of A w.r.t B. 
 * if A[i] and B[i] do not have the same sign or if B[i] = 0, then 
 * RelativeError. error = 0 for this component.
 * @param A The new value of an iteration 
 * @param B The old value of an iteration
 * @return RelativeErrorative error
 */
template<typename Tvec>
real RelativeError(const Tvec &A, const Tvec &B)
{
  real r,p=0.0;
  for(int i=Lb(A) ; i<=Ub(A) ; i++) {
    if( ( zerotest*abs(A[i]) < abs(B[i]) )  ||  sqr(abs(B[i])) == 0.0) 
      r = 0.0;
    else 
      r = abs( (A[i]-B[i])/B[i] );
    if (r>p) p = r;
  }
  return p;
} // end RelativeError

// ----------------------------------------------------------------------------

/**
 * Accurate = accuracy of A is far too bad, i.e. the diameter of A is extremely
 * big
 * 
 * note: 0 for false, 1 for true;
 *
 * @param A A vector
 * @return Too bad or not
 */
inline bool Accurate(const ivector &A) {
  bool bad = false;
  for(int i=Lb(A) ; i<=Ub(A) && !bad ; i++) {
    bad = Inf(A[i]) < -limit && Sup(A[i]) > limit;
  }
  return bad;
} // end Accurate

// ----------------------------------------------------------------------------

/**
 * Accurate = accuracy of A is far too bad
 * 
 * note: 0 for false, 1 for true;
 *
 * @param A A vector
 * @return Too bad or not
 */
inline bool Accurate(const civector &A)
{
  bool bad = false;
 
  for(int i=Lb(A) ; i<=Ub(A) && !bad ; i++)
  {
    bad = InfRe(A[i]) < -limit && SupRe(A[i]) > limit
          || InfIm(A[i]) < -limit && SupIm(A[i]) > limit;
  }
  return bad;
} // end Accurate

// ----------------------------------------------------------------------------

/**
 * x1 is the new, x0 the old value of an iteration. If a component of x1
 * has decreased by more than a factor of zerotest, then this component
 * is set to 0. The same is done if the sign of a component has changed.
 * 
 * @param  x0 old value of iteration
 * @param  x1 new value of iteration
 */
template<typename Tvec>
inline void GuessZeros(Tvec& x0, Tvec& x1) {
  for(int i=Lb(x1) ; i<=Ub(x1) ; i++)
    if( zerotest*abs(x1[i]) < abs(x0[i]) )
          x1[i] = 0.0;
} // end GuessZeros

// ----------------------------------------------------------------------------

/**
 * Computes the midpoint of the system matrix and the right hand side
 * Real point case: The matrices remain unchanged
 *
 * \param MyA System matrix
 * \param MyAm Mid of system matrix
 * \param b Right hand side
 * \param bm Mid of right hand side
 * \param numr Number of rows for this process
 * \param numc Number of columns for this process
 * \param n Dimension of the system
 */
inline void mid(rmatrix &MyA, rmatrix &MyAm, rmatrix &b, rmatrix &bm, 
                int numr, int numc, int n) {
  MyAm = MyA; bm = b;
}

// ----------------------------------------------------------------------------

/**
 * Computes the midpoint of the system matrix and the right hand side
 * Complex point case: The matrices remain unchanged
 *
 * \param MyA System matrix
 * \param MyAm Mid of system matrix
 * \param b Right hand side
 * \param bm Mid of right hand side
 * \param numr Number of rows for this process
 * \param numc Number of columns for this process
 * \param n Dimension of the system
 */
inline void mid(cmatrix &MyA, cmatrix &MyAm, cmatrix &b, cmatrix &bm, int numr, int numc, int n) {
  MyAm = MyA; bm = b;
}

// ----------------------------------------------------------------------------

/**
 * Computes the midpoint of the system matrix and the right hand side
 * Real inteval case: Midpoints are computed (in pure floating point)
 *
 * \param MyA System matrix
 * \param MyAm Mid of system matrix
 * \param b Right hand side
 * \param bm Mid of right hand side
 * \param numr Number of rows for this process
 * \param numc Number of columns for this process
 * \param n Dimension of the system
 */

inline void mid(imatrix &MyA, rmatrix &MyAm, imatrix &b, rmatrix &bm, int numr, int numc, int n) {
  MyAm = rmatrix(Lb(MyA,ROW), Ub(MyA,ROW), Lb(MyA,COL), Ub(MyA,COL));
  bm = rmatrix(Lb(b,ROW), Ub(b,ROW), Lb(b,COL), Ub(b,COL));
  for(int i=1 ; i<=numr ; i++) {
    for(int j=1 ; j<=numc ; j++) {
      MyAm[i][j] = (Inf(MyA[i][j]) + Sup(MyA[i][j])) / 2.0;
    }
  }

  for(int i=Lb(b,ROW) ; i<=Lb(b,COL) ; i++) 
    for(int j=Lb(b,COL) ; j<=Ub(b,COL) ; j++)
      bm[i][j] = (Inf(b[i][j]) + Sup(b[i][j])) / 2.0;
}

// ----------------------------------------------------------------------------

/**
 * Computes the midpoint of the system matrix and the right hand side
 * Complex inteval case: Midpoints are computed (in pure floating point)
 *
 * \param MyA System matrix
 * \param MyAm Mid of system matrix
 * \param b Right hand side
 * \param bm Mid of right hand side
 * \param numr Number of rows for this process
 * \param numc Number of columns for this process
 * \param n Dimension of the system
 */
inline void mid(cimatrix &MyA, cmatrix &MyAm, cimatrix &b, cmatrix &bm, int numr, int numc, int n) {
  for(int i=1 ; i<=numr ; i++) {
    for(int j=1 ; j<=numc ; j++) {
      MyAm[i][j] = (Inf(MyA[i][j]) + Sup(MyA[i][j])) / 2.0;
    }
  }

  for(int i=Lb(b,ROW) ; i<=Lb(b,COL) ; i++) 
    for(int j=Lb(b,COL) ; j<=Ub(b,COL) ; j++)
      bm[i][j] = (Inf(b[i][j]) + Sup(b[i][j])) / 2.0;
}

// ----------------------------------------------------------------------------

//----------------------------------------------------------------------------                      

//! Computes componentwise midpoint of a real matrix
/*
   Computes componentwise the midpoint of a matrix
   Real case: Just returns the matrix itself
   \param A A matrix
   \return Midpoint matrix
*/
inline rmatrix& midpoint(rmatrix &A, int numr, int numc) {
  return A;
}

//----------------------------------------------------------------------------                      

//! Computes componentwise midpoint of an interval matrix
/*
   Computes componentwise the midpoint of a matrix
   Interval case: Computes the midpoint and returns it in a
   new rmatrix
   \param A A matrix
   \return Midpoint matrix
*/
inline rmatrix& midpoint(imatrix &A, int numr, int numc) {
  rmatrix* Am = new rmatrix(Lb(A,ROW), Ub(A,ROW), Lb(A,COL), Ub(A,COL));

  for(int i=Lb(A,ROW) ; i<=Ub(A,ROW) ; i++) {
    for(int j=Lb(A,COL) ; j<=Ub(A,COL) ; j++) {
      (*Am)[i][j] = (Sup(A[i][j]) + Inf(A[i][j])) / 2.0;
    }
  }

  return *Am;
}

//----------------------------------------------------------------------------                      

//! Computes componentwise midpoint of a complex matrix
/*
   Computes componentwise the midpoint of a matrix
   Complex case: Just returns the matrix itself
   \param A A matrix
   \return Midpoint matrix
*/
inline cmatrix& midpoint(cmatrix &A, int numr, int numc) {
  return A;
}

//----------------------------------------------------------------------------                      

//! Computes componentwise midpoint of a complex interval matrix
/*
   Computes componentwise the midpoint of a matrix
   Complex interval case: Computes the midpoint and returns it in a
   new rmatrix
   \param A A matrix
   \return Midpoint matrix
*/
inline cmatrix& midpoint(cimatrix &A, int numr, int numc) {
  cmatrix* Am = new cmatrix(Lb(A,ROW), Ub(A,ROW), Lb(A,COL), Ub(A,COL));

  for(int i=Lb(A,ROW) ; i<=Ub(A,ROW) ; i++) {
    for(int j=Lb(A,COL) ; j<=Ub(A,COL) ; j++) {
      (*Am)[i][j] = complex( (SupRe(A[i][j]) + InfRe(A[i][j])) / 2.0,
                             (SupIm(A[i][j]) + InfIm(A[i][j])) / 2.0);
    }
  }

  return *Am;
}

/**
 * Distributes the Matrix A stored in Process 0 evenly to all Processes. 
 * The data is distributed according to the ScaLAPACK Standard in a
 * two-dimensional block cylic distribution scheme.
 *
 * The processes are organized in a grid of size nr x nc, which is initialized
 * in this function (if init==true), and each process receives blocks of size nb, which are 
 * stored cyclically in one matrix.
 *
 * \param[in]   A        The Matrix to be distributed (available only in 
 *                       Process 0), must be square
 * \param[out]  MyA      Matrix containing the distributed Data of each process
 * \param[in]   procs    Number of processes
 * \param[in]   mypid    Own process id
 * \param[out]  errc     Errorcode
 * \param[out]  commerrc  Errorcode for communication
 * \param[in]   out     Filestream for status messages of this process
 * \param[in]   nr       Number of rows in the process grid
 * \param[in]   nc       Number of columns in the process grid
 * \param[in]   nb       Blocksize of the cyclicaly distributed blocks
 * \param[out]  myr      Row coordinate of process in grid
 * \param[out]  myc      Column coordinate of process in grid
 * \param[out]  numr     Number of rows of own data
 * \param[out]  numc     Number of columns of own data
 * \param[out]  ic       Context for SCALAPACK
 * \param[in]   m        Rows of Matrix A
 * \param[in]   n        Columns of Matrix A
 * \param[in]   init     Must grid be initialised?
 *
 */
template<typename Tmat>
void distributeData(Tmat& A, Tmat& MyA, int procs, int mypid, int& errc, 
                    int& commerrc, ofstream& out, int nr, int nc, int nb,
                    int &myr, int &myc, int &numr, int &numc, int &ic, int &m, int &n, bool init=true) {
  MPI_Status status;

  //ID of process in upper left corner of grid  
  int irsrc = 0, icsrc = 0;

  //Init grid  
  if(init) {
    sl_init_(&ic, &nr, &nc);
    blacs_gridinfo_(&ic, &nr, &nc, &myr, &myc);    
  }

  //Process is not part of the grid
  if(mypid >= nr*nc) {
     return;
  }  
  
  //Determine size of MyA
  numr = numroc_(&m, &nb, &myr, &irsrc, &nr);
  numc = numroc_(&n, &nb, &myc, &icsrc, &nc);  

  //Data Distribution
  if(mypid == 0) {
    int lb1=Lb(A,1), lb2=Lb(A,2);
    for(int i=nr-1 ; i>=0 ; i--) {
      for(int j=nc-1 ; j>=0 ; j--) {
         //Size of MyA of current process
         int row = numroc_(&m, &nb, &i, &irsrc, &nr);
         int col = numroc_(&n, &nb, &j, &icsrc, &nc);  
                 
         //check if process gets any data        
         if(row == 0  ||  col == 0) continue;
         
         //Prepare MyA for current process
         MyA = Tmat(row,col);         
         int currow=1, curcol=1;
         for(int k=0 ; k<m ; k++) {
            for(int l=0 ; l<n ; l++) {
               if((k/nb)%nr == i  &&  (l/nb)%nc == j) {
                  MyA[currow][curcol] = A[lb1+k][lb2+l];
                  curcol++;
                  if(curcol > col) {
                     curcol = 1;
                     currow++;
                  }
               }
            }
         }         
         
         //Send data if root is not current process         
         if(i!=0 || j!=0) {
           //senden 
           commerrc=MPI_Send(MyA,1,row,1,col,i*nc+j, 0, MPI_COMM_WORLD);           
         } 
      }
    }
  } else if(mypid < nr*nc) {
     //Receive own MyA if not empty
     if(numr != 0  &&  numc != 0) {
       MyA = Tmat(numr, numc);
       commerrc=MPI_Recv(MyA, 1, numr, 1, numc, 0, 0, MPI_COMM_WORLD, &status);
     }
  }
  
  if(commerrc != 0)
    errc = DistFailed;

}



/**
 * Collects the Matrix A, distributed in a two-dimensional block cyclic distribution scheme, in
 * process 0.
 *
 *
 * \param[in]   A        The Matrix to be distributed (available only in 
 *                       Process 0), must be square
 * \param[out]  MyA      Matrix containing the distributed Data of each process
 * \param[in]   procs    Number of processes
 * \param[in]   mypid    Own process id
 * \param[out]  errc     Errorcode
 * \param[out]  commerrc  Errorcode for communication
 * \param[in]   out     Filestream for status messages of this process
 * \param[in]   nr       Number of rows in the process grid
 * \param[in]   nc       Number of columns in the process grid
 * \param[in]   nb       Blocksize of the cyclicaly distributed blocks
 * \param[out]  myr      Row coordinate of process in grid
 * \param[out]  myc      Column coordinate of process in grid
 * \param[out]  numr     Number of rows of own data
 * \param[out]  numc     Number of columns of own data
 * \param[out]  ic       Context for SCALAPACK
 * \param[in]   m        Rows of Matrix A
 * \param[in]   n        Columns of Matrix A
 *
 */
template<typename Tmat>
void collectData(Tmat& A, Tmat& MyA, int procs, int mypid, int& errc, 
                    int& commerrc, ofstream& out, int nr, int nc, int nb,
                    int &myr, int &myc, int &numr, int &numc, int &ic, int &m, int &n) {
  MPI_Status status;
  out << "Starting Data Collection" << endl;

  //ID of process in upper left corner of grid  
  int irsrc = 0, icsrc = 0;

  //Process is not part of the grid
  if(mypid >= nr*nc) {
     out << "Data Collection finished" << endl;
     return;
  }  
  
  //Determine size of MyA
  numr = numroc_(&m, &nb, &myr, &irsrc, &nr);
  numc = numroc_(&n, &nb, &myc, &icsrc, &nc);  

  //Data collection
  if(mypid == 0) {
    A = Tmat(m,n);
    for(int i=0 ; i<nr ; i++) {
      for(int j=0 ; j<nc ; j++) {
         //Size of MyA of current process
         int row = numroc_(&m, &nb, &i, &irsrc, &nr);
         int col = numroc_(&n, &nb, &j, &icsrc, &nc);  

         if(row==0 || col==0)
	   continue;
	 
	 //Receive data if root is not current process         
         if(i!=0 || j!=0) {
           commerrc=MPI_Recv(MyA, i*nc+j, 0, MPI_COMM_WORLD, &status);
         } 
          
         int currow=1, curcol=1;
         for(int k=0 ; k<m ; k++) {
            for(int l=0 ; l<n ; l++) {
               if((k/nb)%nr == i  &&  (l/nb)%nc == j) {
                  A[k+1][l+1] = MyA[currow][curcol];
                  curcol++;
                  if(curcol > col) {
                     curcol = 1;
                     currow++;
                  }
               }
            }
         }         
	 
      }
    }

	
  } else {
    //Send data
    if(numr != 0  &&  numc != 0 && RowLen(MyA)>0 && ColLen(MyA)>0) {
      commerrc=MPI_Send(MyA,0, 0, MPI_COMM_WORLD);
    }    
  }
  
    
  if(commerrc != 0)
    errc = DistFailed;

  out << "Data Collection complete" << endl;
}

// ----------------------------------------------------------------------------

/*!
 * Determines the number of rows and columns in the process grid used
 * as logical topography. Tries to compute both values such that they are
 * about equal. However, all processes are used always, so that if the number
 * of processes p is prime, a px1 grid will be used.
 *
 * \param[in] Number of processes
 * \param[out] Smaller dimension of grid
 * \param[out] Larger dimension of grid
 */
void determine_grid(int procs, int &p1, int &p2) {  
  if(procs == 1  ||  procs == 2  ||  procs == 3) {
    p1 = 1; p2 = procs;
    return;
  }
  
  p1 = (int)floor(sqrt((double)procs));
  p2 = p1;
  int p1_tmp, p2_tmp = p2;
  int maxp = 0;
  
  for(p1_tmp = p1 ; p1_tmp>=1 ; p1_tmp--) {
    while (p1_tmp*(p2_tmp+1)<=procs) 
      p2_tmp++;
    
    if(p1_tmp * p2_tmp > maxp) {
       p1 = p1_tmp;
       p2 = p2_tmp;
       maxp = p1 * p2;
    }
    
    if(p1 * p2 == procs)
      break;
  }     
  
}

// ----------------------------------------------------------------------------

//Blow-Function for class civector
civector Blow(const civector &x, const real &eps) {
  civector y(Lb(x), Ub(x));
  
  for(int i=Lb(x) ; i<=Ub(x) ; i++)
    y[i] = Blow(x[i], eps);
  
  return y;
}


//----------------------------------------------------------------------------                      

//! Determines if an interval is empty
/*
  Determines if an interval x is empty, i.e. Inf(x) > Sup(x)
  \param zz Interval to check
  \return true if Inf(x) > Sup(x)
*/
inline bool Empty(const interval &i) {
  return Inf(i) > Sup(i);
}

//----------------------------------------------------------------------------                      

//! Determines if a complex interval is empty
/*
  Determines if a complex interval x is empty, i.e. Inf > Sup for the real or imaginary part
  \param zz Interval to check
  \return true if interval is empty
*/
inline bool Empty(const cinterval &i) {
  return InfRe(i) > SupRe(i) || InfIm(i) > SupIm(i);
}

//----------------------------------------------------------------------------           

//! Compute the max distance between two interval vectors
/*'maxDist()' computes the maximal distance between two 
   real/complex interval vectors.
   For interval x, y; dist(x,y)= Max(Abs(x.inf-y.inf), Abs(x.sup-y.sup))
   \param xx First interval
   \param yy Second interval
   \return Maximal distance between xx and yy
*/
template<class TVec1, class TVec2>
static real maxDist(const TVec1& xx, const TVec2& yy) {
  real dist = 0.0, tmp;
  int n = VecLen(xx);

  for(int i=1 ; i<=n ; i++) {
    tmp = abs(Inf(xx[i]) - Inf(yy[i]));
    if(tmp > dist) dist = tmp;
    tmp = abs(Sup(xx[i]) - Sup(yy[i]));
    if(tmp > dist) dist = tmp;
  }

  return dist;
}

//! Compute the max distance between two interval matrices
/*'maxDist()' computes the maximal distance between two 
   real/complex interval matrices.
   For interval x, y; dist(x,y)= Max(Abs(x.inf-y.inf), Abs(x.sup-y.sup))
   \param xx First interval
   \param yy Second interval
   \return Maximal distance between xx and yy
*/
template<class TMat1, class TMat2>
static real maxDistMat(const TMat1& xx, const TMat2& yy) {
  real dist = 0.0, tmp;
  int n = RowLen(xx);

  for(int i=1 ; i<=n ; i++) {
    tmp = maxDist(xx[Col(i)], yy[Col(i)]);
    if(tmp > dist) dist = tmp;
  }

  return dist;
}

// ----------------------------------------------------------------------------

/*!
 * Broadcasts own part of a parallel computed vector in the row of the
 * process in the process grid. All parts are the accumulated. 
 * After this function, every process in this row has the parts
 * of the vector, that correspond to the rows of MyA.
 *
 * \param[in] myr Own row coordinate in process grid
 * \param[in] myc Own column coordinate in process grid
 * \param[in] numr Number of rows in MyA
 * \param[in] nc Number of columns in process grid
 * \param[in] row_comms MPI-Communicators for own row of process grid
 * \param[in/out] dotvec Accumulator for end result
 * \paramm[in/out] tmpvec Temp accumulator for communication
 * \param[in/out] mydotvec Accumulator conatining own part of computation
 * 
 */
template<typename TDot>
void bcast_iteration_parts(int myr, int myc, int numr, int nc,
                                  MPI_Comm *row_comms, 
                                  vector<TDot> &dotvec, 
                                  vector<TDot> &mydotvec, 
                                  vector<TDot> &tmpdotvec) {
  //Send own parts of dot product of each row as vector of accumulators
  //broadcast only in own row of processes!
  for(int l=0 ; l<nc ; l++) {
    if(myc == l)
      MPI_Bcast(mydotvec, l, row_comms[myr]);
    else
      MPI_Bcast(tmpdotvec, l, row_comms[myr]);
    for(int i=0 ; i<numr ; i++) {
      if(l == 0) dotvec[i] = 0.0;
      if(l == myc)
        dotvec[i] += mydotvec[i];
      else
        dotvec[i] += tmpdotvec[i];
    }
  }
}

// ----------------------------------------------------------------------------

/*!
 * Exchanges own part of a parallel computed vector in the own process column.
 * After this function every process has the complete vector.
 *
 * \param[in] myr Own row coordinate in process grid
 * \param[in] myc Own column coordinate in process grid
 * \param[in] numr Number of rows in MyA
 * \param[in] nr Number of rows in process grid
 * \param[in] nc Number of columns in process grid
 * \param[in] nb Blocksize for data distribution
 * \param[in] n Dimension of global matrix A
 * \param[in] irsrc row coordinate of upper left process
 * \param[in] col_comms MPI-Communicators for own column of process grid
 * \param[in/out] x Vector to be communicated (at entry only own part, after 
 *                  this function completly)
 */
template<typename Tvec>
void bcast_iteration_vector(int myr, int myc, int numr, int nr, int nc, 
                                   int nb, int n, int irsrc, 
                                   MPI_Comm *col_comms, Tvec &x) {
  //Exchange own part of x1 with other processes in the same grid column
  //(same will be done throughout the programm for parallel computed vectors)   
  int i,j,l;  
  for(i=0 ; i<nc ; i++) {   
    if(myc != i) continue;
    Tvec tmpvec;
    for(j=0 ; j<nr ; j++) {               
      tmpvec = Tvec(numroc_(&n, &nb, &j, &irsrc, &nr));
      int start = 1;

      if(myr == j) {
        for(l=j*nb+1 ; l<=n ; l+=nr*nb) {
          int end = l+nb-1;
          if(end>n) end=n;  
          if(end == n)
            tmpvec(start,numr) = x(l,end);
          else
            tmpvec(start,start+nb-1) = x(l,end);          
          start += nb;
        }
      }
      
      MPI_Bcast(tmpvec, j, col_comms[myc]);
      
      start = 1;
      for(l=j*nb+1 ; l<=n ; l+=nr*nb) {
        int end = l+nb-1;
        if(end>n) end=n;  
        if(end == n)
          x(l,end) = tmpvec(start,Ub(tmpvec));
        else
          x(l,end) = tmpvec(start,start+nb-1);          
        start += nb;
      }

    }
  }
    
}    

// ----------------------------------------------------------------------------

//! Compute the residual using a splitting algorithm by Rump/Dekker
template<typename TB>
inline TB FastResidual(TB MyB, rmatrix& MyA, rmatrix& MyX, int r, int s, int t, int nb, int ic, int numr, ofstream& out) {
  static real factor = 68719476737;
  bool empty = false;
  rmatrix MyT, MyT2;
  if(RowLen(MyA)>0 && ColLen(MyA)>0) {
    MyT = factor * MyA;
    MyT2 = MyT - MyA;
    MyT = MyT - MyT2;
    MyT2 = MyA - MyT;
  } else {
    empty = true;
  }
  
  rmatrix MyY, MyY2;
  if(RowLen(MyX)>0 && ColLen(MyX)>0) {
    MyY = factor * MyX;
    MyY2 = MyY - MyX;
    MyY = MyY - MyY2;
    MyY2 = MyX - MyY;
  } else {
    empty = true;
  }
  
  TB MyTmp(MyX), MyTmp2(MyT2), MyRes(MyX);
  matmul(MyTmp2,MyX,MyTmp,r,s,t,nb,ic,numr,out);    
  MyTmp2 = TB(MyT);
  matmul(MyTmp2,MyY2,MyRes,r,s,t,nb,ic,numr,out);  
  if(!empty)
    MyTmp += MyRes;
  matmul(MyTmp2,MyY,MyRes,r,s,t,nb,ic,numr,out);  
  
  if(!empty) 
    return (MyB - MyRes) - MyTmp;
  else
    return TB();
}

//! Compute the residual using a splitting algorithm by Rump/Dekker
template<typename TB>
inline TB FastResidual(TB MyB, cmatrix& MyA, cmatrix& MyX, int r, int s, int t, int nb, int ic, int numr, ofstream& out) {
  static real factor = 68719476737;
  bool empty = false;
  cmatrix MyT, MyT2;
  if(RowLen(MyA)>0 && ColLen(MyA)>0) {  
    MyT = factor * MyA;
    MyT2 = MyT - MyA;
    MyT = MyT - MyT2;
    MyT2 = MyA - MyT;
  } else {
    empty = true;
  }
  
  cmatrix MyY, MyY2, MyMX;
  if(RowLen(MyX)>0 && ColLen(MyX)>0) {
    MyMX = -MyX;
    MyY = factor * MyMX;
    MyY2 = MyY - MyMX;
    MyY = MyY - MyY2;
    MyY2 = MyMX - MyY;
  } else {
    empty = true;
  }
  
  TB MyTmp(MyX), MyTmp2(MyT2), MyRes(MyX);
  matmul(MyTmp2,MyMX,MyTmp,r,s,t,nb,ic,numr,out);    
  MyTmp2 = TB(MyT);
  matmul(MyTmp2,MyY2,MyRes,r,s,t,nb,ic,numr,out);  
  if(!empty)
    MyTmp += MyRes;
  matmul(MyTmp2,MyY,MyRes,r,s,t,nb,ic,numr,out);  
  
  if(!empty) 
    return (MyRes + MyB) + MyTmp;
  else
    return TB();
}

//! Compute the residual
template<typename TB>
inline TB FastResidual(TB MyB, imatrix& MyA, rmatrix& MyX, int r, int s, int t, int nb, int ic, int numr, ofstream& out) {
  imatrix MyTmp(ColLen(MyB),RowLen(MyB));
  matmul(MyA,MyX,MyTmp,r,s,t,nb,ic,numr,out); 
  return MyB - MyTmp;
}

//! Compute the residual
template<typename TB>
inline TB FastResidual(TB MyB, cimatrix& MyA, cmatrix& MyX, int r, int s, int t, int nb, int ic, int numr, ofstream& out) {
  cimatrix MyTmp(ColLen(MyB),RowLen(MyB));
  matmul(MyA,MyX,MyTmp,r,s,t,nb,ic,numr,out);      
  return MyB - MyTmp;
}

  

// ----------------------------------------------------------------------------


//   typename TSysMat,      //Type of the system matrix
//   typename TRhs,         //Type of the right hand side
//   typename TSolution,    //Type of the solution vector
//   typename TInverse,     //Type of the approximate inverse
//   typename TC,           //Type of [C]=I-RA
//   typename TMidMat,      //Type of midpoint matrix
//   typename TMidVec,      //Type of midpoint vector
//   typename TDot,         //Type for point dot products
//   typename TIDot,        //Type for interval dot products
//   typename TFunc,        //Type of function delivering matrix elements
//   typename Tb,		 //Type of one column of right hand side
//   typename TVerivVec	    //Type of solution vectors
//   typename TInterval     //Interval type (interval or cinterval)

/**
//! Computes a verified solution of a square linear system of equations A*x=b
/*! 
 *   An approximate inverse 'R' of 'A' is computed by using the parallel 
 *   ScaLAPACK routine. Then an approximate solution 'x' is computed
 *   applying a residual iteration. For the final verification, 
 *   an interval residual iteration is performed. An enclosure of the
 *   unique solution is returned in  'Y'. All iterations
 *   will also be done in parallel.
 *   
 *   If this first step fails, the solver can try to compute a verified
 *   solution by using an approximate Inverse of double length. This second
 *   step takes considerably longer than the first step, because all computations
 *   must me performed using high precision scalar products. It is configurable
 *   if both steps should be used, or only one of two.
 *   
 *   This solver uses PBLAS/ScaLAPACK routines for the computation of the 
 *   approximate inverse and for the matrix-matrix multiplication in step 1. 
 *   All additional matrix-matrix multiplications and scalar products are 
 *   computed using the DotK-Algorithm, which allows computation in K-fold
 *   double precision.
 *
 *
 * @param[in]  MyA        The matrix of the system
 * @param[in]  b          Right hand side(s)
 * @param[out] Y          Solution enclosure
 * @param[in]  procs      Number of parallel processes
 * @param[in]  mypid      Process ID
 * @param[in]  n          Matrix dimensions
 * @param[in]  lb1        Start of row index
 * @param[in]  lb2        Start of column index
 * @param[in]  nr         Number of rows in process grid
 * @param[in]  nc         Number of columns in proces grid
 * @param[in]  myr        Own row coordinate in process grid
 * @param[in]  myc        Own column in process grid
 * @param[in]  numr       Rows of own matrix
 * @param[in]  numc       Columns of row matrix
 * @param[in]  ic         Blacs context
 * @param[in]  irsrc      Row of first process in grid
 * @param[in]  icsrc      Column of first process in grid
 * @param[out] errc       Error code
 * @param[out] commerrc   Error code
 * @param[in]  out        Output file stream
 * @param[in]  cfg        Struct containing configuration information for solver
 * @param[in]  inner      Compute inner solution?
 * @param[out] ienc       Computed inner solution
 */
template <typename TSysMat, typename TRhs, typename TSolution, typename TInverse, typename TC,
          typename TMidMat, typename TMidVec, typename TDot, typename TIDot, typename TFunc, typename Tb, 
          typename TVerivVec, typename TInterval>
void LSSMain(TSysMat &MyA, TRhs &Myb, TSolution &Y, int procs, int mypid, int n, int rhs,
               int lb1, int lb2, int nr, int nc, int myr, int myc, int numr, int numc, 
	       int ic, int irsrc, int icsrc, int &errc, int &commerrc, ofstream &out, 
	       struct plssconfig cfg, bool inner, TSolution& ienc) {

  //Variable Declaration  
  TInverse      MyR;
  TC            MyC;
  Tb		b;
  TMidVec       bm;
  TMidVec       x0(n), x1(n), y0(n), x1temp(n); 
  TVerivVec     Y1(n), YA(n), Z(n),  Y1temp(n);
  int           i,j,k,l;
  bool          ready;
  real          p,eps1;
  TDot          dot;
  TIDot         idot;
  int           currow, curcol=1;
  bool		C_computed = false;
  int           nb = cfg.nb;

  dot.set_k(cfg.K);
  idot.set_k(cfg.K);
  opdotprec = cfg.K;

  double timer = time();

  MPI_Status status;
  int root = 0;   
  ready = false;
  errc = NoError;
  
  #ifdef _OPENMP
    if(cfg.threads > 0)
      omp_set_num_threads(cfg.threads);
  #endif

  //Create MPI communicators for each row in process grid
  MPI_Group MPI_GROUP_WORLD;
  MPI_Group row_groups[nr];
  MPI_Comm row_comms[nr];
  MPI_Group col_groups[nc];
  MPI_Comm col_comms[nc];
  int process_ranks_row[nr][nc];
  int process_ranks_col[nc][nr];  
  
  for(i=0 ; i<nr ; i++)
    for(j=0 ; j<nc ; j++) {
      process_ranks_row[i][j] = i*nc + j;
      process_ranks_col[j][i] = j + i*nc;
    }
  
  MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  
  for(i=0 ; i<nr ; i++) {
    MPI_Group_incl(MPI_GROUP_WORLD, nc, process_ranks_row[i], &(row_groups[i]));
    MPI_Comm_create(MPI_COMM_WORLD, row_groups[i], &(row_comms[i]));
  }

  for(i=0 ; i<nc ; i++) {
    MPI_Group_incl(MPI_GROUP_WORLD, nr, process_ranks_col[i], &(col_groups[i]));
    MPI_Comm_create(MPI_COMM_WORLD, col_groups[i], &(col_comms[i]));
  }
  
  
//  MyAm = mid(MyA);
//  bm = mid(b);
  //MyAm = TMidMat(numr,numc);
  //Mybm = TMidMat(n,procs/rhs+1);
  TMidMat&  MyAm = midpoint(MyA, numr, numc);
  TMidMat&  Mybm = midpoint(Myb, numr, numc);


  //mid(MyA, MyAm, Myb, Mybm, numr, numc, n);

  bm   = TMidVec(n);
  b    = Tb(n);


  //Compute approximate Inverse of A
  timer = time();
  MatInv( MyR, MyAm, procs, mypid, errc, commerrc, out, nr, nc, nb, myr, myc,
          numr, numc, ic, n);
  if(errc != 0) {
    errc = InvFailed;
    if((TMidMat*)&MyA != &MyAm) delete &MyAm;
    if((TMidMat*)&Myb != &Mybm) delete &Mybm;
    return;    
  }
  out << "Time for inversion:" << time() - timer << endl;
 
  if(cfg.lssparts != LSS_ONLY_PART_TWO) {
    out << endl << "Starting stage one of the solver..." << endl;
    //Start of part 1               
    vector<TDot> mydotvec(numr);
    vector<TIDot> myidotvec(numr);
    vector<TDot> tmpdotvec(numr);
    vector<TIDot> tmpidotvec(numr);
    vector<TDot> dotvec(numr);
    vector<TIDot> idotvec(numr);
    
    if(cfg.matrixMode) {
      out << "Using matrix mode" << endl;      
      TInverse  MyX(ColLen(Myb),RowLen(Myb));
      TSolution MyZ(ColLen(Myb),RowLen(Myb));
      TSolution MyTmp(ColLen(Myb),RowLen(Myb));      

      //X=R*bm;
      double tmp = time();	
      out << "Computing approximate solution x" << endl;
      matmul(MyR,Mybm,MyX,n,n,rhs,cfg.nb,ic,numr,out);
      
      out << "Defect iteration" << endl;
      TInverse MyTmpP2 = FastResidual(Mybm, MyAm, MyX, n, n, rhs, cfg.nb, ic, numr, out);
      TInverse MyTmpP(MyTmpP2);
      matmul(MyR,MyTmpP,MyTmpP2,n,n,rhs,cfg.nb,ic,numr,out);
      if(RowLen(MyX)>0 && ColLen(MyX)>0)
        MyX += MyTmpP2;
      out << "Time for computation of x: " << time()-tmp<<endl;

      tmp = time();
      //[Z]=[R*(B-A*X)];
      out << "Computing [z]" << endl;
      MyTmp = FastResidual(TSolution(Myb), MyA, MyX, n, n, rhs, cfg.nb, ic, numr, out);
      matmul(MyR,MyTmp,MyZ,n,n,rhs,cfg.nb,ic,numr,out);
      out << "Time for computation of [z]: " << time()-tmp<<endl;

      //[C]=[I-R*A]
      if(numr>0 && numc>0)
        Resize(MyC, numr, numc);
      out << "Computing [C]" << endl;
      double timerCMM = time();
      IminusAB(MyR, MyA, MyC, n, n, n, nb, ic, numr, myr, myc, nr, nc, out);    
      out << "Time for computation of [C]: " << time()-timerCMM << endl;
      
      //Verification
      tmp = time();
      out << "Verification" << endl;      
      Y = MyX;
      k = 0;
      do {
	MyTmp = Y;
        out << " Step " << k+1 << endl;      
	for(int i=1 ; i<=ColLen(MyTmp) ; i++)
	  for(int j=1 ; j<=RowLen(MyTmp) ; j++) 
	    MyTmp[i][j] = Blow(MyTmp[i][j],cfg.epsVer);
	  
        matmul(MyC,MyTmp,Y,n,n,rhs,cfg.nb,ic,numr,out);	

	ready = true;
	if(RowLen(Y)>0 && ColLen(Y)>0) {
	  Y = MyZ + Y;

	  for(int i=1 ; i<=ColLen(MyTmp) ; i++)
	    for(int j=1 ; j<=RowLen(MyTmp) ; j++)
              ready = ready && in(Y[i][j],MyTmp[i][j]);
	}

	int myready;
	for(int p=0 ; p<procs ; p++) {
	  myready = ready;
	  MPI_Bcast(&myready, 1, MPI_INT, p, MPI_COMM_WORLD);
	  if(!myready) {
	    ready = false;
	    break;
	  }
	}
	
	k++;
      } while(!ready && k<cfg.maxIterVer);
      out << "Time for verification: " << time()-tmp << endl;

      k = 0;
      if(cfg.refinement) {
	bool ref_criterion = true;
	out << "Refinement step" << endl;
	do {
          out << " Step " << k+1 << endl;   
	
	  MyTmp = Y;
          matmul(MyC,MyTmp,Y,n,n,rhs,cfg.nb,ic,numr,out);	
	  if(RowLen(Y)>0 && ColLen(Y)>0) {
	    Y = MyZ + Y;
	    Y &= MyTmp;
	    ref_criterion = maxDistMat(Y,MyTmp) > cfg.epsRef;
	  }
	  
	  int mycriterion;
	  for(int p=0 ; p<procs ; p++) {
	    mycriterion = ref_criterion;
	    MPI_Bcast(&mycriterion, 1, MPI_INT, p, MPI_COMM_WORLD);
	    if(!mycriterion) {
	      ref_criterion = false;
	      break;
	    }
	  }	  
	
	  k++;
        } while(ref_criterion && k<cfg.maxIterRef);
      }
      
      if(ready) {
	if(RowLen(Y)>0 && ColLen(Y)>0)
	  Y += MyX;
      } else {
	errc = VerivFailed;
      }
    }

    for(int s=1 ; s<=rhs && errc!=VerivFailed && !cfg.matrixMode ; s++) {
       
        out << endl << "Computing result for right hand side " << s << endl;

        //distribute current right hand side
        if((s-1)%procs == mypid) {
          b  = Myb[Col((s-1)/procs+1)];
          bm = Mybm[Col((s-1)/procs+1)];
        }
        MPI_Bcast(b, (s-1)%procs, MPI_COMM_WORLD);
        MPI_Bcast(bm, (s-1)%procs, MPI_COMM_WORLD);
  
	if (errc == 0) { 
	  out << "Starting defect iteration..." << endl;
	  //------------------------------------------------------
	  // floating point defect iteration: result is x1
	  //------------------------------------------------------
	  timer = time();
	  k = 0;
	
	  //x1 = R1*b;
	  currow = 0;

	  for(i=1 ; i<=n ; i++) {
	    //Check if row is part of own blocks
	    if(((i-1)/nb)%nr == myr) {
	      currow++;
	      dot = 0.0;
	
	      if(nc == 1) {
		//Complete dot product can be computed with own Data
		accumulate_approx(dot, MyR[Row(currow)], bm);
		x1[i] = rnd(dot);  
	      } else {
		//Can only compute parts of the dot product, save them in accumulators
		for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
		  int end = j+nb-1;
		  if(end>n) end=n;
		  accumulate_approx(dot, MyR[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), bm(j,end));
		}
		
		mydotvec[currow-1] = dot;
	
	      }
	    }
	  }
	
	  if(nc != 1) {
	    //Send own parts of dot product of each row as vector of accumulators
	    //broadcast only in own row of processes!
	    bcast_iteration_parts(myr, myc, numr, nc, row_comms, dotvec, mydotvec, 
				tmpdotvec);
	
	    //compute final result of dot products for own rows
	    currow = 0;
	    for(i=1 ; i<=n ; i++) {
	      if(((i-1)/nb)%nr == myr) {
		x1[i] = rnd(dotvec[currow]);
		currow++;
	      }
	    }    
	  }
			
	  //Exchange own part of x1 with other processes in the same grid column
	  //(same will be done throughout the programm for parallel computed vectors)     
	  bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, x1);
	
	  // iterate x = x + R*(b-Ax)	    
	  do {
	    k = k + 1;
	    x0 = x1;
	
	    out << " " << k << ". iteration step" << endl;
	    // x1 := #*(b - A*x0)     
	    currow = 0;
	    for(i=1 ; i<=n ; i++) {
	      //Check if row is part of own blocks
	      if(((i-1)/nb)%nr == myr) {
	        currow++;
		if(nc == 1) {
		  dot = bm[i];
		  accumulate_approx(dot,-MyAm[Row(currow)], x0);
		  x1[i] = rnd(dot);
		} else {
		  dot = 0.0;
		  for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			int end = j+nb-1;
			if(end>n) end=n;
			accumulate_approx(dot, -MyAm[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), x0(j,end));
		  }
	
		  mydotvec[currow-1] = dot;
	
		}
	      }
	    }

	    if(nc != 1) {
	      //Send own parts of dot product of each row as vector of accumulators
	      //broadcast only in own row of processes!
	      bcast_iteration_parts(myr, myc, numr, nc, row_comms, dotvec, mydotvec, 
					tmpdotvec);
	
	      //compute final result of dot products for own rows
	      currow = 0;
	      for(i=1 ; i<=n ; i++) {
		if(((i-1)/nb)%nr == myr) {
		  dotvec[currow] += bm[i];                            
		  x1[i] = rnd(dotvec[currow]);
		  currow++;
		}
	      }    
	    }

	    //Exchange own part of x1 with other processes in the same grid column
	    //(same will be done throughout the programm for parallel computed vectors)     
	    bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, x1);
		
	    // x1 := #*(x0 + R1*x1)
	    currow = 0;
	    for(i=1 ; i<=n ; i++) {
	      //Check if row is part of own blocks
	      if(((i-1)/nb)%nr == myr) {
		currow++;
	        if(nc == 1) {
		  dot = x0[i];
		  accumulate_approx(dot,MyR[Row(currow)], x1);
		  x1temp[i] = rnd(dot);
		} else {
		  dot = 0.0;
		  for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			int end = j+nb-1;
			if(end>n) end=n;
			accumulate_approx(dot, MyR[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), x1(j,end));
		  }
	
		  mydotvec[currow-1] = dot;
	
		}
	      }
	    }
	
	    if(nc != 1) {
	      //Send own parts of dot product of each row as vector of accumulators
	      //broadcast only in own row of processes!
	      bcast_iteration_parts(myr, myc, numr, nc, row_comms, dotvec, mydotvec, 
					tmpdotvec);
	
	      //compute final reault of dot products for own rows
	      currow = 0;
	      for(i=1 ; i<=n ; i++) {
		if(((i-1)/nb)%nr == myr) {
		  dotvec[currow] += x0[i];                            
		  x1temp[i] = rnd(dotvec[currow]);
		  currow++;
		}
	      }    
	    }
			
	    //Exchange own part of x1 with other processes in the same grid column
	    //(same will be done throughout the programm for parallel computed vectors)     
	    bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, 
				x1temp);
		
	    x1 = x1temp;
	    p = RelativeError(x1,x0);
	    GuessZeros(x1,x0);
	} while (k<=cfg.maxIterResCorr && p>=delta);

	
	out << "Time for defect iteration: " << time() - timer << endl;
	
	timer = time();  
	out << "Computing [z]" << endl; 
	// y0 := #*(b-A*x1)
	currow = 0;
	for(i=1 ; i<=n ; i++) {
	  //Check if row is part of own blocks
	  if(((i-1)/nb)%nr == myr) {
	    currow++;
	    if(nc == 1) {
	      idot = b[i];
	      accumulate(idot,-MyA[Row(currow)], x1);
	      y0[i] = mid(rnd(idot));
	    } else {
	      idot = 0.0;
	      for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
	        int end = j+nb-1;
	        if(end>n) end=n;
	        accumulate(idot, -MyA[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), x1(j,end));
	     }
	     myidotvec[currow-1] = idot;
	    }
	  }
	}
	
	if(nc != 1) {
	  //Send own parts of dot product of each row as vector of accumulators
	  //broadcast only in own row of processes!
	  bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, 
				tmpidotvec);
	
	  //compute final reault of dot products for own rows
	  currow = 0;
	  for(i=1 ; i<=n ; i++) {
	    if(((i-1)/nb)%nr == myr) {
	      idotvec[currow] += b[i];                            
	      y0[i] = mid(rnd(idotvec[currow]));
	      currow++;
	    }
	  }    
	}
		
	//Exchange own part of x1 with other processes in the same grid column
	//(same will be done throughout the programm for parallel computed vectors)   
	bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, y0);
	
	// Y1 := ##(b-A*x1-y0)
	currow = 0;
	for(i=1 ; i<=n ; i++) {
	  //Check if row is part of own blocks
	  if(((i-1)/nb)%nr == myr) {
	    currow++;
	    if(nc == 1) {
	      idot = b[i];
	      accumulate(idot, -MyA[Row(currow)], x1);
	      idot += -y0[i];
	      Y1[i] = rnd(idot);
	     } else {
	      idot = 0.0;
	      for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
	        int end = j+nb-1;
		if(end>n) end=n;
		accumulate(idot, -MyA[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), x1(j,end));
	      }
	
	      myidotvec[currow-1] = idot;
	    }
	  }
	}
	
	if(nc != 1) {
	  //Send own parts of dot product of each row as vector of accumulators
	  //broadcast only in own row of processes!
	  bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, 
				tmpidotvec);
	
	  //compute final reault of dot products for own rows
	  currow = 0;
	  for(i=1 ; i<=n ; i++) {
	    if(((i-1)/nb)%nr == myr) {
	      idotvec[currow] += b[i];                            
	      idotvec[currow] += -y0[i];                                     
	      Y1[i] = rnd(idotvec[currow]);
	      currow++;
	    }
	  }    
	}
		
	//Exchange own part of x1 with other processes in the same grid column
	//(same will be done throughout the programm for parallel computed vectors)     
	bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, Y1); 
	
	// Y1 := ##( R1*y0 + R1*Y1 );
	currow = 0;
	for(i=1 ; i<=n ; i++) {
	  //Check if row is part of own blocks
	    if(((i-1)/nb)%nr == myr) {
	      currow++;
	      idot = 0.0;
	      if(nc == 1) {
		accumulate(idot,MyR[Row(currow)], y0);
		accumulate(idot,MyR[Row(currow)], Y1);
		Y1temp[i] = rnd(idot);
	      } else {
		for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
		  int end = j+nb-1;
		  if(end>n) end=n;
		  accumulate(idot, MyR[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), y0(j,end));
		  accumulate(idot, MyR[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), Y1(j,end));         
		}
	
		myidotvec[currow-1] = idot;
	      }
	    }
	}
	
	if(nc != 1) {
	  //Send own parts of dot product of each row as vector of accumulators
	  //broadcast only in own row of processes!
	  bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, 
				tmpidotvec);
	
	  //compute final reault of dot products for own rows
	  currow = 0;
	  for(i=1 ; i<=n ; i++) {
	    if(((i-1)/nb)%nr == myr) {
	      Y1temp[i] = rnd(idotvec[currow]);
	      currow++;
	    }
	  } 
	}
		
	//Exchange own part of x1 with other processes in the same grid column
	//(same will be done throughout the programm for parallel computed vectors)     
	bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, Y1temp);     
	
	Y1 = Y1temp;
	Z = Y1;
	
	out << "Time for computation of [z]: " << time() - timer << endl;
	
	out << "Verifying..." << endl;      
	timer = time();    
	if (isZero(Z)) {  // exact solution! (however, not necessarily unique!)
                if((s-1)%procs == mypid) {
                   Y[Col(curcol)] = x1;
                   curcol++;
                }
		errc = NoError;
		ready = true;
	} else {
		// interval iteration until inclusion is obtained (or max. iteration count)
		k = 0;
		eps1 = cfg.epsVer;
	
		if(!C_computed) {
                  //MyC = TC(numr, numc);
                  Resize(MyC, numr, numc);
		  out << " Computing [C]" << endl;
		  double timerC = time();
		  IminusAB(MyR, MyA, MyC, n, n, n, nb, ic, numr, myr, myc, nr, nc, out);    
		  out << " Time for computation of [C]: " << time()-timerC << endl;
		  C_computed = true;
		}
	
		do {         
		  k++;
		  out << " " << k << ". verification step" << endl;
		  YA = Blow(Y1,eps1);

		  //Y1 = Z + C*YA
		  currow = 0;
		  for(i=1 ; i<=n ; i++) {
		    //Check if row is part of own blocks
		    if(((i-1)/nb)%nr == myr) {
		      currow++;
		      if(nc == 1) {
			idot = Z[i];
			accumulate(idot, MyC[Row(currow)], YA);
			Y1[i] = rnd(idot);
		      } else {
			idot = 0.0;
			for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			  int end = j+nb-1;
			  if(end>n) end=n;
			  accumulate(idot, MyC[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), YA(j,end));
			}
	
			myidotvec[currow-1] = idot;
		      }
		    }
		  }


		  if(nc != 1) {
		    //Send own parts of dot product of each row as vector of accumulators
		    //broadcast only in own row of processes!
		    bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, tmpidotvec);
		    //compute final reault of dot products for own rows
		    currow = 0;
		    for(i=1 ; i<=n ; i++) {
		      if(((i-1)/nb)%nr == myr) {
			idotvec[currow] += Z[i];
			Y1[i] = rnd(idotvec[currow]);
			currow++;
		      }
		    }    
		  }


		  //Exchange own part of x1 with other processes in the same grid column
		  //(same will be done throughout the programm for parallel computed vectors)     
		  bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, Y1);      
		  ready = in(Y1,YA);


		} while ( (!ready && k<cfg.maxIterVer && !Accurate(Y1)) );


		// result
		if (ready) {
                    
                  TVerivVec D, zzi, zze;

                  if(cfg.refinement) {
		    out << " Starting iterative refinement" << endl;
		    TVerivVec yy(n), xx(n);
		    int p = 0;
		    real Epsilon = cfg.epsRef;

		    do {
		      yy = Y1;
		      
		      currow = 0;
		      for(i=1 ; i<=n ; i++) {
		        //Check if row is part of own blocks
		        if(((i-1)/nb)%nr == myr) {
		          currow++;
		          if(nc == 1) {
		            idot = Z[i];
			    accumulate(idot, MyC[Row(currow)], Y1);
			    rnd(idot, xx[i]);
		          } else {
			    idot = Z[i];
			    for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			      int end = j+nb-1;
			      if(end>n) end=n;
			      accumulate(idot, MyC[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), Y1(j,end));
			    }
	
			    myidotvec[currow-1] = idot;
		          }
		        } 
		      }


		      if(nc != 1) {
		        //Send own parts of dot product of each row as vector of accumulators
		        //broadcast only in own row of processes!
		        bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, tmpidotvec);
		        //compute final reault of dot products for own rows
		        currow = 0;
		        for(i=1 ; i<=n ; i++) {
		          if(((i-1)/nb)%nr == myr) {
			    rnd(idotvec[currow], xx[i]);
			    currow++;
		          }
		        }   
		      }
		    
		      bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, xx); 		
		      
		      Y1 = xx & yy;
		      p++;
		      
		    } while (maxDist(Y1,yy) > Epsilon  &&  (p < cfg.maxIterRef));
		  }


                  if(inner) {
		    out << " Computing inner enclosure" << endl;
		    D = TVerivVec(n);
		    zzi = TVerivVec(n);
		    zze = TVerivVec(n);
  		    idot.set_dotprec(0);

		    currow = 0;
		    for(i=1 ; i<=n ; i++) {
		      //Check if row is part of own blocks
		      if(((i-1)/nb)%nr == myr) {
		        currow++;
		        if(nc == 1) {
		          idot = 0.0;
			  accumulate(idot, MyC[Row(currow)], Y1);
			  rnd(idot, D[i]);
		        } else {
			  idot = 0.0;
			  for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			    int end = j+nb-1;
			    if(end>n) end=n;
			    accumulate(idot, MyC[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), Y1(j,end));
			  }
	
			  myidotvec[currow-1] = idot;
		        }
		      } 
		    }


		    if(nc != 1) {
		      //Send own parts of dot product of each row as vector of accumulators
		      //broadcast only in own row of processes!
		      bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, tmpidotvec);
		      //compute final reault of dot products for own rows
		      currow = 0;
		      for(i=1 ; i<=n ; i++) {
		        if(((i-1)/nb)%nr == myr) {
			  rnd(idotvec[currow], D[i]);
			  currow++;
		        }
		      }    
		      
		    }
		    
		    bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, D); 		    
		    

		    currow = 0;
		    for(i=1 ; i<=n ; i++) {
		      //Check if row is part of own blocks
		      if(((i-1)/nb)%nr == myr) {
		        currow++;
		        if(nc == 1) {
		          idot = b[i];
			  accumulate(idot, -MyA[Row(currow)], x1);
			  UncheckedSetInf(zzi[i], rnd(Inf(idot),RND_UP));
			  UncheckedSetSup(zzi[i], rnd(Sup(idot),RND_DOWN));
		        } else {
			  idot = b[i];
			  for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			    int end = j+nb-1;
			    if(end>n) end=n;
			    accumulate(idot, -MyA[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), x1(j,end));
			  }
	
			  myidotvec[currow-1] = idot;
		        }
		      } 
		    }


		    if(nc != 1) {
		      //Send own parts of dot product of each row as vector of accumulators
		      //broadcast only in own row of processes!
		      bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, tmpidotvec);
		      //compute final reault of dot products for own rows
		      currow = 0;
		      for(i=1 ; i<=n ; i++) {
		        if(((i-1)/nb)%nr == myr) {
			  UncheckedSetInf(zzi[i], rnd(Inf(idotvec[currow]),RND_UP));
			  UncheckedSetSup(zzi[i], rnd(Sup(idotvec[currow]),RND_DOWN));
			  currow++;
		        }
		      }    
		    }

    		    bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, zzi); 		    

		    currow = 0;
		    for(i=1 ; i<=n ; i++) {
		      //Check if row is part of own blocks
		      if(((i-1)/nb)%nr == myr) {
		        currow++;
		        if(nc == 1) {
		          idot = 0.0;
			  accumulate(idot, MyR[Row(currow)], zzi);
			  UncheckedSetInf(zze[i], rnd(Inf(idot),RND_UP));
			  UncheckedSetSup(zze[i], rnd(Sup(idot),RND_DOWN));
		        } else {
			  idot = 0.0;
			  for(j=myc*nb+1 ; j<=n ; j+=nc*nb) {
			    int end = j+nb-1;
			    if(end>n) end=n;
			    accumulate(idot, MyR[Row(currow)](j/(nb*nc)*nb+1,j/(nb*nc)*nb+(end-j)+1), zzi(j,end));
			  }
	
			  myidotvec[currow-1] = idot;
		        }
		      } 
		    }


		    if(nc != 1) {
		      //Send own parts of dot product of each row as vector of accumulators
		      //broadcast only in own row of processes!
		      bcast_iteration_parts(myr, myc, numr, nc, row_comms, idotvec, myidotvec, tmpidotvec);
		      //compute final reault of dot products for own rows
		      currow = 0;
		      for(i=1 ; i<=n ; i++) {
		        if(((i-1)/nb)%nr == myr) {
			  UncheckedSetInf(zze[i], rnd(Inf(idotvec[currow]),RND_UP));
			  UncheckedSetSup(zze[i], rnd(Sup(idotvec[currow]),RND_DOWN));
			  currow++;
		        }
		      }    
		    }
		    idot.set_dotprec(cfg.K);	

    		    bcast_iteration_vector(myr, myc, numr, nr, nc, nb, n, irsrc, col_comms, zze); 
                  } //if(inner)
                  

                  if((s-1)%procs == mypid) {

                     if(inner) {
		       for(i=1 ; i<=n ; i++) {
	                 UncheckedSetInf(ienc[i][curcol], addu(addu(x1[i],Inf(zze[i])), Sup(D[i])));
		         UncheckedSetSup(ienc[i][curcol], addd(addd(x1[i],Sup(zze[i])), Inf(D[i])));
			 if(Empty(ienc[i][curcol]))
			   ienc[i][curcol] = TInterval(SignalingNaN);
		       }
                     }

		     Y[Col(curcol)] = x1 + Y1;
                     curcol++;
                   }

		   errc = NoError;

		} else {

		   errc = VerivFailed;

		}
	}
	
	out << "Time for verification: " << time() - timer << endl;
      } // end USE_SINGLE_R 

   }

   if(!ready  &&  (errc != VerivFailed || cfg.lssparts==LSS_ONLY_PART_ONE)) {
      out << "LSS finished in error." << endl;
      if((TMidMat*)&MyA != &MyAm) delete &MyAm;
      if((TMidMat*)&Myb != &Mybm) delete &Mybm;
      return;
   } else if(ready) {     
      out << "LSS finished." << endl;
      if((TMidMat*)&MyA != &MyAm) delete &MyAm;
      if((TMidMat*)&Myb != &Mybm) delete &Mybm;
      return;
   }
 
 } //end lss_cfg!=LSS_ONLY_PART_TWO


//------------------------------------------------------------------------------
// Part II
//------------------------------------------------------------------------------
// if no success: try again with approximate inverse R = R1+R2 of double length
// USE_DOUBLE_R to try again with approximate inverse R = R1+R2 of double length 
//------------------------------------------------------------------------------ 
 if(cfg.lssparts == LSS_BOTH_PARTS)
   out << endl << "Stage 1 not successful. Starting stage 2..." << endl;
 else
   out << "Starting stage 2..." << endl;
 
 errc = NoError;
 C_computed = false; 
 curcol = 1;

 if(cfg.lssparts != LSS_ONLY_PART_TWO) {
   out << "Reorganizing data for stage 2" << endl;
   timer = time();
   //Rearrange A and R for second part
   int x_global, y_global;
   int nb_new = n / procs;
   if(n%procs != 0) nb_new++;
   TSysMat A, Atmp;
   TInverse R, Rtmp;

   if(numr != 0  &&  numc != 0) {
     A = MyA;
     R = MyR;       
   }

   if(n-mypid*nb_new >= nb_new) {
      MyA = TSysMat(nb_new, n);
      MyR = TInverse(nb_new, n);         
      MyC = TC(nb_new, n);
   } else {
      int size = n-mypid*nb_new;
      if(size > 0) {
        MyA = TSysMat(size, n);
        MyR = TInverse(size, n);         
        MyC = TC(size, n);
      } else {
        MyA = TSysMat(0,0);
        MyR = TInverse(0,0);
      }
   }           

   for(int r=0 ; r<nr ; r++) {
     for(int c=0 ; c<nc ; c++) {      
       //Communicate block of current process
       int rows = numroc_(&n, &nb, &r, &irsrc, &nr);
       int cols = numroc_(&n, &nb, &c, &icsrc, &nc);   

       if(rows == 0  || cols == 0) continue;          

       Atmp = TSysMat(rows,cols);
       Rtmp = TInverse(rows,cols);
       if(myr == r  && myc == c) {
         MPI_Bcast(A, r*nc+c, MPI_COMM_WORLD); 
         MPI_Bcast(R, r*nc+c, MPI_COMM_WORLD); 
       } else {
         MPI_Bcast(Atmp, r*nc+c, MPI_COMM_WORLD); 
         MPI_Bcast(Rtmp, r*nc+c, MPI_COMM_WORLD);          
       }

       //Rearrange own MyA and MyR, copy data from received block
       int nbr = (nr==1) ? n : nb;
       int nbc = (nc==1) ? n : nb;          
 
       for(i=1 ; i<=rows ; i++) {
         for(j=1 ; j<=cols ; j++) {
           x_global = ((i-1)/nbr)*nr*nbr+r*nbr+((i-1)%nbr)+1;
           y_global = ((j-1)/nbc)*nc*nbc+c*nbc+((j-1)%nbc)+1;
           if(x_global > mypid*nb_new  &&  x_global <= (mypid+1)*nb_new) {
             if(myr == r && myc == c) {
               MyA[(x_global-1)%nb_new+1][y_global] = A[i][j];
               MyR[(x_global-1)%nb_new+1][y_global] = R[i][j];               
             } else {
               MyA[(x_global-1)%nb_new+1][y_global] = Atmp[i][j];
               MyR[(x_global-1)%nb_new+1][y_global] = Rtmp[i][j];               
             }
           }
         }
       }
     
     }
   }
   
   //delete tmp matrices
   A = TSysMat();
   Atmp = TSysMat();
   R = TInverse();
   Rtmp = TInverse();
   
   //Reset data distribution attributes according to rearrangement
   nb = nb_new;
   nr = procs;
   nc = 1;
   sl_init_(&ic, &nr, &nc);
   blacs_gridinfo_(&ic, &nr, &nc, &myr, &myc);    
   numr = numroc_(&n, &nb, &myr, &irsrc, &nr);
   numc = numroc_(&n, &nb, &myc, &icsrc, &nc);   
   out << "Time for Reorganization: " << time() - timer << endl;    
 }

 //Set indices of MyA and MyR in relation to A
 int lm = lb1 + mypid * nb;
 int um = lm + numr - 1;
 int ln = lb2;
 int un = ln + n - 1;

 if(cfg.lssparts == LSS_ONLY_PART_TWO)
   MyC = TC(lm,um,ln,un);


 if(numr == 0 || numc == 0) {
   //Process has no data
   lm=-1; um=-2;
 } else {  
   //set index ranges
   SetLb(MyA, ROW, lm); SetUb(MyA, ROW, um);
   SetLb(MyA, COL, ln); SetUb(MyA, COL, un);
   //SetLb(MyAm, ROW, lm); SetUb(MyAm, ROW, um);
   //SetLb(MyAm, COL, ln); SetUb(MyAm, COL, un);
   SetLb(MyR, ROW, lm); SetUb(MyR, ROW, um);
   SetLb(MyR, COL, ln); SetUb(MyR, COL, un);
   SetLb(MyC, ROW, lm); SetUb(MyC, ROW, um);
   SetLb(MyC, COL, ln); SetUb(MyC, COL, un);

   if((TMidMat*)&MyA != &MyAm) {
     delete &MyAm;
   }
   if((TMidMat*)&Myb != &Mybm) {
     delete &Mybm;
   }

//   MyAm = rmatrix(lm,um,ln,un);
 }

 TMidMat&  MyAm2 = midpoint(MyA, numr, numc);
 TMidMat&  Mybm2 = midpoint(Myb, numr, numc);
 
 //Determine row index ranges for every process 
 //( MyA of process p: Lb(MyA,1)==ind1[p] and Ub(MyA,1)==ind2[p] )
 int ind1[procs];
 int ind2[procs];
 for(i=0 ; i<procs ; i++) {          
   if(mypid == i) {
     ind1[i] = lm;
     ind2[i] = um;
   } else {
     ind1[i] = lb1 + i * nb;
     ind2[i] = ind1[i] + nb - 1;
     if(ind1[i]>n) {
       ind1[i]=-1;
       ind2[i]=-2;
     } else if(ind2[i]>n) {
       ind2[i] = n;             
     }
   }
 }

 
 //----------------------------------------------------------------------
 // Parallel: R1A=R1*A 
 //----------------------------------------------------------------------
 out << "Computing R*A" << endl;
 timer = time();
 TInverse MyRA(lm, um, ln, un);
 MatMul(MyR, MyAm2, MyRA, mypid, procs, ind1, ind2, errc, commerrc, out, cfg.K, cfg.threads);
 out << "Time for computation of R*A: " << time() - timer << endl;
 

 //----------------------------------------------------------------------------
 // Parallel computation of approximate inverse
 //----------------------------------------------------------------------------
 out << "Computing approximate inverse of R*A" << endl;
 timer = time();
 TInverse MyR2;
 MatInv( MyR2, MyRA, procs, mypid, errc, commerrc, out, nr, nc, nb, myr, myc,
         numr, numc, ic, n);
 SetLb(MyR2, ROW, lm); SetUb(MyR2, ROW, um);
 SetLb(MyR2, COL, ln); SetUb(MyR2, COL, un);         
//delete MyRA
 MyRA = TInverse(); 
 out << "Time for inversion of R*A: " << time() - timer << endl;
      
 if(errc!=0) {
   out << "LSS finished in error" << endl;
   errc = InvFailed;
   return;            
 }

 //----------------------------------------------------------------------
 // Parallel: "D := #*(R2*R1)"
 // Parallel: "R2 := #* (R2*R1 - D);"
 //----------------------------------------------------------------------
 out << "Computing R2=#*(R2*R1 - #*(R2*R1))" << endl;
 timer = time();
 TInverse MyR1temp(lm, um, ln, un);
 TInverse MyR2temp(lm, um, ln, un); 
 ABminusRndAB(MyR2,MyR,MyR1temp,MyR2temp,mypid,procs,ind1,ind2,errc,commerrc,out,cfg.K,cfg.threads);

 MyR  = MyR1temp;
 MyR2 = MyR2temp;

 //Delete temp matrices
 MyR1temp = TInverse();
 MyR2temp = TInverse();     

 out << "Time for computation of R2: " << time() - timer << endl;

 if (errc!=0) {
   out << "LSS finished in error" << endl;
   errc = VerivFailed;
   return;
 }


 for(int s=1 ; s<=rhs && errc!=VerivFailed ; s++) {

        out << endl << "Computing result for right hand side " << s << endl;

        //distribute current right hand side
        if((s-1)%procs == mypid) {
          b  = Myb[Col((s-1)/procs+1)];
          bm = Mybm2[Col((s-1)/procs+1)];
        }
        MPI_Bcast(b, (s-1)%procs, MPI_COMM_WORLD);
        MPI_Bcast(bm, (s-1)%procs, MPI_COMM_WORLD);

	//------------------------------------------------------
	// floating point defect iteration: result is x1+x0
	//------------------------------------------------------
	out << "Starting defect iteration..." << endl;
	timer = time();
	
	
	// x1 := #*( R1*b + R2*b );
	for (i=lm ; i<=um ; i++) {
	  dot = 0.0;
	  accumulate_approx(dot,MyR[Row(i)],bm);
	  accumulate_approx(dot,MyR2[Row(i)],bm);
	  x1[i] = rnd(dot);
	}
	
	for(i=0 ; i<procs ; i++) {   
	  if(ind1[i] != -1) {                   
	    MPI_Bcast(x1, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	  }
	}  
	
	// x0 = #*( R1*b + R2*b - x1):
	for (i=lm ; i<=um ; i++) {
	  dot = -x1[i];
	  accumulate_approx(dot,MyR[Row(i)],bm);
	  accumulate_approx(dot,MyR2[Row(i)],bm);
	  x0[i] = rnd(dot);
	}
	
	for(i=0 ; i<procs ; i++) {   
	  if(ind1[i] != -1) {                   
	    MPI_Bcast(x0, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	  }
	}  
	
	// iteration x = x + (R1+R2)*(b-Ax), x = x1 + x0)
	k = 0;
	do {
	  k++;
	  out << " " << k << ". iteration step" << endl;
	  // y0 = #*(b - A*x1 - A*x0)
	  for (i=lm ; i<=um ; i++) {
	    dot = bm[i];
	    accumulate_approx(dot,-MyAm2[Row(i)],x1);
	    accumulate_approx(dot,-MyAm2[Row(i)],x0);
	    y0[i] = rnd(dot);
	  }
	
	  for(i=0 ; i<procs ; i++) {   
	    if(ind1[i] != -1) {                   
	      MPI_Bcast(y0, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	    }
	  }   
	
	  // y0 := #*(x0 + R1*y0 + R2*y0);
	  for (i=lm ; i<=um ; i++) {
	    dot = x0[i];
	    accumulate_approx(dot,MyR[Row(i)],y0);
	    accumulate_approx(dot,MyR2[Row(i)],y0);
	    x1temp[i] = rnd(dot);
	  }
	
	  for(i=0 ; i<procs ; i++) {   
	    if(ind1[i] != -1) {                   
	      MPI_Bcast(x1temp, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	    }
	  }   
	
	  y0 = x1temp;
	  p = RelativeError (x1+y0,x1+x0);
	  y0 = x1 + y0;
	
	  // x0 := #*(x1 + x0 - y0)
	  for (i=Lb(x1);i<=Ub(x1);i++) {
	    dot = x1[i];
	    dot += x0[i];
	    dot += -y0[i];
	    x0[i] = rnd(dot);
	  }
	  x1 = y0;

	} while ( k<=cfg.maxIterResCorr && p>=delta );
	
	out << "Time for iteration: " << time() - timer << endl;
	
	// compute enclosure y0+Y1 of the residuum b-A*x1 of the approximation x1 and
	// initialize Y1:= Z:= (R1+R2)*(b-A*x1), C:= I-(R1+R2)*A
	out << "Computing [z]" << endl;
	timer = time();
	// y0 := #*(b-A*x1) 
	for (i=lm ; i<=um ; i++) {
	  dot = bm[i];
	  accumulate_approx(dot,-MyAm2[Row(i)],x1);
	  y0[i] = rnd(dot);
	}
	
	for(i=0 ; i<procs ; i++) {   
	  if(ind1[i] != -1) {                   
	    MPI_Bcast(y0, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	  }
	}   
	
	// Y1 := ##(b-A*x1-y0)
	for (i=lm ; i<=um ; i++) {
	  idot = b[i];
	  accumulate(idot,-MyA[Row(i)],x1);
	  idot += -y0[i];
	  Y1[i] = rnd(idot);
	}
	
	for(i=0 ; i<procs ; i++) {   
	  if(ind1[i] != -1) {                   
	    MPI_Bcast(Y1, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	  }
	}   
	
	// Y1 := ##(R1*y0 + R2*y0 + R1*Y1 + R2*Y1 );
	for (i=lm ; i<=um ; i++) {
	  idot = 0.0;
	  accumulate(idot,MyR[Row(i)],y0);
	  accumulate(idot,MyR2[Row(i)],y0);
	  accumulate(idot,MyR[Row(i)],Y1);
	  accumulate(idot,MyR2[Row(i)],Y1);
	  Y1temp[i] = rnd(idot);
	}
	
	for(i=0 ; i<procs ; i++) {   
	  if(ind1[i] != -1) {                   
	    MPI_Bcast(Y1temp, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	  }
	}   
	
	Y1 = Y1temp;
	Z = Y1;
	
	out << "Time for computation of [z]: " << time() - timer << endl;
	
	if (errc!=0) {
	  out << "LSS finished in error" << endl;
	  errc = VerivFailed;
	  return;
	}
	
	out << "Verifying..." << endl;
	timer = time();
	if (isZero(Z)) {  // exact solution! (however, not necessarily unique!)
          if((s-1)%procs == mypid) {
	    Y[Col(curcol)] = x1;
            curcol++;
          }
	  errc = NoError;
	  ready = true;
	} else {
          if(!C_computed) {
  	    //------------------------------------------------------
	    // Compute C:= ## (ID(A) - R1*A - R2*A);
	    //------------------------------------------------------
	    out << " Computing [C]" << endl;
	    double timerC = time();
	    IminusA1A2B(MyR,MyR2,MyA,MyC,mypid,procs,ind1,ind2,errc,commerrc,out,cfg.K,cfg.threads);
	    out << " Time for computation of [C]: " << time() - timerC << endl;
            C_computed = true;
          }

	  // interval iteration until inclusion is obtained (or max. iteration count)
	  k = 0;
	  eps1 = cfg.epsVer;

	  do {
	    k++;
	    out << " " << k << ". verification step" << endl;
	
	    YA = Blow( Y1, eps1);
	
	    //Y1 = Z + C*YA
	    for(i=lm ; i<=um ; i++) {   
		idot = Z[i];
		accumulate(idot,MyC[Row(i)], YA);
		Y1[i] = rnd(idot);
	    }
	
	    for(i=0 ; i<procs ; i++) { 
		if(ind1[i] != -1)                          
		MPI_Bcast(Y1, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	    }           
	
	    ready = in(Y1,YA);
	  } while ( !ready && k<cfg.maxIterVer && !Accurate(Y1) );
	
	  // output of the result
	  if (ready) {
             TVerivVec D, zzi, zze;
	     
             if(cfg.refinement) {
	      out << " Starting iterative refinement" << endl;
	      TVerivVec yy(n), xx(n);
	      int p = 0;
	      real Epsilon = cfg.epsRef;

	      do {
		yy = Y1;
		      
	        for(i=lm ; i<=um ; i++) {
	          idot = Z[i];
		  accumulate(idot, MyC[Row(i)], Y1);
		  rnd(idot, xx[i]);
		}

	        for(i=0 ; i<procs ; i++) {   
	          if(ind1[i] != -1) {                   
	            MPI_Bcast(xx, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	          }
	        }  

		
		Y1 = xx & yy;
		p++;
		      
	      } while (maxDist(Y1,yy) > Epsilon  &&  (p < cfg.maxIterRef));
             }
            
             if(inner) {
		out << " Computing inner enclosure" << endl;
		D   = TVerivVec(n);
                zzi = TVerivVec(n);
		zze = TVerivVec(n);
  		idot.set_dotprec(0);

	        for(i=lm ; i<=um ; i++) {
	          idot = 0.0;
		  accumulate(idot, MyC[Row(i)], Y1);
		  rnd(idot, D[i]);
		}

	        for(i=0 ; i<procs ; i++) {   
	          if(ind1[i] != -1) {                   
	            MPI_Bcast(D, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	          }
	        }  

	        for(i=lm ; i<=um ; i++) {
	          idot = b[i];
		  accumulate(idot, -MyA[Row(i)], x1);
		  UncheckedSetInf(zzi[i], rnd(Inf(idot),RND_UP));
		  UncheckedSetSup(zzi[i], rnd(Sup(idot),RND_DOWN));	  
		}

	        for(i=0 ; i<procs ; i++) {   
	          if(ind1[i] != -1) {                   
	            MPI_Bcast(zzi, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	          }
	        }  

	        for(i=lm ; i<=um ; i++) {
	          idot = 0.0;
		  accumulate(idot, MyR[Row(i)], zzi);
		  accumulate(idot, MyR2[Row(i)], zzi);		  
		  UncheckedSetInf(zze[i], rnd(Inf(idot),RND_UP));
		  UncheckedSetSup(zze[i], rnd(Sup(idot),RND_DOWN));	  
		}

	        for(i=0 ; i<procs ; i++) {   
	          if(ind1[i] != -1) {                   
	            MPI_Bcast(zze, ind1[i], ind2[i], i, MPI_COMM_WORLD);
	          }
	        }  

		idot.set_dotprec(cfg.K);	

            } //if(inner)

            if((s-1)%procs == mypid) {
              if(inner) {
		  for(i=1 ; i<=n ; i++) {
		    UncheckedSetInf(ienc[i][curcol], addu(addu(x1[i],Inf(zze[i])), Sup(D[i])));
		    UncheckedSetSup(ienc[i][curcol], addd(addd(x1[i],Sup(zze[i])), Inf(D[i])));
		    if(Empty(ienc[i][curcol]))
	              ienc[i][curcol] = TInterval(SignalingNaN);
		  }
              }

	      Y[Col(curcol)] = x1 + Y1;
              curcol++;
            }
	    errc = NoError;
	  } else errc = VerivFailed;
     }
 }
 out << "Time for verification: " << time() - timer << endl;

 if((TMidMat*)&MyA != &MyAm2) {
   delete &MyAm2;
 }
 if((TMidMat*)&Myb != &Mybm2) {
   delete &Mybm2;
 }

 if(ready)
   out << "LSS finished. " << endl;
 else
   out << "LSS finished in error." << endl;  
  
}  // end LSS


// ----------------------------------------------------------------------------


/*!
 * Starting function for verified parallel linear system solver. This is the entry
 * function for the case that A is not available on process 0. In this case, a 
 * function pointer to a function returning the value at a position (i,j) of A must 
 * be supplied. With the help of this function, every process fills his own part of the
 * matrix himself.
 *
 * Two dimensional block cyclic distribution is used according to the
 * blocksize paramter nb.
 *
 * The result Y is an enclosure of the solution
 *
 *
 * @param[in]  Get_A      Pointer to a function f(i,j,r), that computes r=A(i,j)
 * @param[in]  b          Right hand side
 * @param[out] Y          Solution enclosure
 * @param[in]  m,n        Dimensions of Matrix A
 * @param[in]  procs      Number of parallel processes
 * @param[in]  mypid      Process ID
 * @param[out] errc       Error code
 * @param[out] commerrc   Error code
 * @param[in]  out        Output file stream
 * @param[in]  cfg        Struct containing configuration information for solver
 * @param[in]  inner      Compute inner solution?
 * @param[out] ienc       Computed inner solution
 */
template <typename TSysMat, typename TRhs, typename TSolution, typename TInverse, typename TC,
          typename TMidMat, typename TMidVec, typename TDot, typename TIDot, typename TFunc, typename Tb, 
          typename TVerivVec, typename TInterval>
void SolverStart(TFunc Get_A, TFunc Get_b, TSolution& Y, int m, int n, int rhs, 
                   int procs, int mypid, int& errc, int& commerrc,
                   ofstream& out, struct plssconfig cfg, bool inner, TSolution& ienc) {
  int dim;
  
  if(m!=n) 
    dim = m+n;
  else
    dim = n;
            
  //Indices of A           
  int lb1=1, ub1=m, lb2=1, ub2=n; 
  
  int breakrow = -1;
  
  if(cfg.K != 1 && cfg.matrixMode) {
    out << "Warning: Activating Matrix Mode automatically sets K=1" << endl;
    cfg.K = 1;  
  }
  
  if(cfg.K == 1) {
    if(!cfg.matrixMode) {
      out << "Remark: Using K==1 automatically activates Matrix Mode" << endl;
      cfg.matrixMode = true;
    }
    
    if(cfg.lssparts != LSS_ONLY_PART_ONE) {
      out << "Warning: Only part one of solver can be used with K==1, adjusting configuration" << endl;
      cfg.lssparts = LSS_ONLY_PART_ONE;
    }
    
    if(inner) {
      out << "Warning: Inner enclosure can not be computed with K==1" << endl;
      inner = false;
    }
  }

  out << endl;

  double timer = time();

  out << "Starting parallel solver..." << endl;
  out << "Dimension of system: " << m << "x" << n << endl;
  out << "Precision K: " << cfg.K << endl;
  out << "Number of processes: " << procs << endl;
  out << "ID of this process: " << mypid << endl << endl;

  int root = 0;   
  errc = NoError;
  
  //Values for data distribution (two dimensional block cyclic with
  //block dimension nb)
  int nr, nc;
  int myr, myc, numr, numc, numcrhs, ic, irsrc=0, icsrc=0;
  
  if(cfg.nb<=0) cfg.nb = 256;

  if(cfg.lssparts != LSS_ONLY_PART_TWO) {
    determine_grid(procs,nc,nr);
  } else {
    nc = 1; nr = procs;
    cfg.nb = dim / procs;
    if(dim%procs != 0) cfg.nb++;
  }
  
  TSysMat MyA;
  TRhs    Myb;
  
  //Init grid
  sl_init_(&ic, &nr, &nc);
  blacs_gridinfo_(&ic, &nr, &nc, &myr, &myc);     

  out << "Starting data distribution" << endl;
  
  //Distribute A and b
  if(n > m) {
    out << "Underdetermined system, creating equivalent square system" << endl;    
    //Determine size of MyA
    numr = numroc_(&dim, &(cfg.nb), &myr, &irsrc, &nr);
    numc = numroc_(&dim, &(cfg.nb), &myc, &icsrc, &nc);  
    numcrhs = numroc_(&rhs, &(cfg.nb), &myc, &icsrc, &nc);   
    
    MyA = TSysMat(numr,numc);
    
    //Fill own MyA
    int currow=1, curcol=1;
    for(int k=0 ; k<dim ; k++) {
      for(int l=0 ; l<dim ; l++) {
         if((k/cfg.nb)%nr == myr  &&  (l/cfg.nb)%nc == myc) {
            if(k>=n && l>=m) {
              Get_A(k-m+1,l-n+1,MyA[currow][curcol]);
            } else if(k<n && l>=m) {
              if(k == l-m)
                MyA[currow][curcol] = -1.0;
              else 
                MyA[currow][curcol] = 0.0;
            } else if(k>=n && l<m) {
              MyA[currow][curcol] = 0.0;
            } else if(k<n && l<m) {
              Get_A(l+1,k+1,MyA[currow][curcol]);              
            }
            curcol++;
            if(curcol > numc) {
               curcol = 1;
               currow++;
            }
         }
      }
    }        


    int mybsize = rhs/procs;
    if(rhs%procs != 0  &&  mypid < rhs%procs)
      mybsize++; 
    TRhs BIG_b(1,dim,1,mybsize);    
    
    if(cfg.matrixMode) {

      BIG_b = TRhs(numr,numcrhs);
      Y = TSolution(numr,numcrhs);      
      currow=1; curcol=1;
      for(int k=0 ; k<m+n ; k++) {
        for(int l=0 ; l<rhs ; l++) {
           if((k/cfg.nb)%nr == myr  &&  (l/cfg.nb)%nc == myc) {
	      if(k>=m) {
		if(breakrow == -1) breakrow = currow;
	      }
	      
              if(k<n) {
		BIG_b[currow][curcol] = 0.0;
	      } else {
                Get_b(k-n+1,l+1,BIG_b[currow][curcol]);
	      }
              curcol++;
              if(curcol > numcrhs) {
                curcol = 1;
                currow++;
              }
            }
        }
      }         
      Myb = BIG_b;
      //delete old Data to save memory    
      BIG_b = TRhs();

      
    } else {

      int mybsize = rhs/procs;
      if(rhs%procs != 0  &&  mypid < rhs%procs)
        mybsize++; 
      BIG_b( 1, n, 1, mybsize ) = 0.0;
      curcol = 1;
      for(int j=1 ; j<=rhs ; j++) {
         if((j-1)%procs == mypid) {
            for(int i=1 ; i<=m ; i++)
              Get_b(i,j,BIG_b[n+i][curcol]);
            curcol++;
         }
      }
      Myb = BIG_b;
      //delete old Data to save memory    
      BIG_b = TRhs();

      Y = TSolution(1,dim,1,RowLen(Y));
      if(inner)
        ienc = TSolution(1,dim,1,RowLen(Y)); 
    }
           
  } else if(m > n) {
    //over determined case
    out << "Overdetermined system, creating equivalent square system" << endl;    
    //Determine size of MyA
    numr = numroc_(&dim, &(cfg.nb), &myr, &irsrc, &nr);
    numc = numroc_(&dim, &(cfg.nb), &myc, &icsrc, &nc);  
    numcrhs = numroc_(&rhs, &(cfg.nb), &myc, &icsrc, &nc);   
    
    MyA = TSysMat(numr,numc);
    
    //Fill own MyA
    int currow=1, curcol=1;
    for(int k=0 ; k<dim ; k++) {
      for(int l=0 ; l<dim ; l++) {
         if((k/cfg.nb)%nr == myr  &&  (l/cfg.nb)%nc == myc) {
            if(k<m && l<n) {
              Get_A(k+1,l+1,MyA[currow][curcol]);
            } else if(k<m && l>=n) {
              if(k == l-n)
                MyA[currow][curcol] = -1.0;
              else
                MyA[currow][curcol] = 0.0;
            } else if(k>=m && l<n) {
              MyA[currow][curcol] = 0.0;
            } else if(k>=m && l>=n) {
              Get_A(l-n+1,k-m+1,MyA[currow][curcol]);              
            }
            curcol++;
            if(curcol > numc) {
               curcol = 1;
               currow++;
            }
         }
      }
    }        

    int mybsize = rhs/procs;
    if(rhs%procs != 0  &&  mypid < rhs%procs)
      mybsize++; 
    TRhs BIG_b(1,dim,1,mybsize);
    
    if(cfg.matrixMode) {

      BIG_b = TRhs(numr,numcrhs);
      Y = TSolution(numr,numcrhs);      
      currow=1; curcol=1;
      for(int k=0 ; k<m+n ; k++) {
        for(int l=0 ; l<rhs ; l++) {
           if((k/cfg.nb)%nr == myr  &&  (l/cfg.nb)%nc == myc) {
	      if(k>=n) {
     		if(breakrow == -1) breakrow = currow - 1;
	      }

	      if(k>=m) {
		BIG_b[currow][curcol] = 0.0;
	      } else {
                Get_b(k+1,l+1,BIG_b[currow][curcol]);
	      }
	      
              curcol++;
              if(curcol > numcrhs) {
                curcol = 1;
                currow++;
              }
            }
        }
      }         
      
      Myb = BIG_b;
      BIG_b = TRhs();

      
    } else {

      curcol = 1;
      for(int j=1 ; j<=rhs ; j++)
         for(int i=1 ; i<=m ; i++)
           if((j-1)%procs == mypid) {
             Get_b(i,j,BIG_b[i][curcol]);
             curcol++;
           }

      int mybsize = rhs/procs;
      if(rhs%procs != 0  &&  mypid < rhs%procs)
        mybsize++;            
      BIG_b( m+1,m+n,1,mybsize ) = 0.0;
      Myb = BIG_b;
      //delete old Data to save memory    
      BIG_b = TRhs();

      Y = TSolution(1,dim,1,RowLen(Y));
      if(inner)
        ienc = TSolution(1,dim,1,RowLen(Y));
    }
    
  } else {
    //Square system
  
    //Determine size of MyA
    numr = numroc_(&n, &(cfg.nb), &myr, &irsrc, &nr);
    numc = numroc_(&n, &(cfg.nb), &myc, &icsrc, &nc);  
    numcrhs = numroc_(&rhs, &(cfg.nb), &myc, &icsrc, &nc);   

    MyA = TSysMat(numr,numc);

    //Fill own MyA
    int currow=1, curcol=1;
    for(int k=0 ; k<n ; k++) {
      for(int l=0 ; l<n ; l++) {
         if((k/cfg.nb)%nr == myr  &&  (l/cfg.nb)%nc == myc) {
            Get_A(k+1,l+1,MyA[currow][curcol]);
            curcol++;
            if(curcol > numc) {
               curcol = 1;
               currow++;
            }
          }
      }
    }         

    //Fill own Myb
    if(cfg.matrixMode) {
      
      Myb = TRhs(numr,numcrhs);
      Y = TSolution(numr,numcrhs);

      currow=1; curcol=1;
      for(int k=0 ; k<n ; k++) {
        for(int l=0 ; l<rhs ; l++) {
         if((k/cfg.nb)%nr == myr  &&  (l/cfg.nb)%nc == myc) {
            Get_b(k+1,l+1,Myb[currow][curcol]);
            curcol++;
            if(curcol > numcrhs) {
               curcol = 1;
               currow++;
            }
         }
        }
      }      
      
    } else {

      int mybsize = rhs/procs;
      if(rhs%procs != 0  &&  mypid < rhs%procs)
        mybsize++;      
      Myb = rmatrix(n,mybsize);
      curcol = 1;
      for(int j=1 ; j<=rhs ; j++) {
         if((j-1)%procs == mypid) {
           for(int i=1 ; i<=n ; i++)
              Get_b(i,j,Myb[i][curcol]);
           curcol++;
         }
      }
      
    }
         
  }

  out << "Data distribution finished" << endl;    
  

  out << "Time for data distribution and initialization: " << time() - timer << endl;  
  
  LSSMain<TSysMat,TRhs,TSolution,TInverse,TC,TMidMat,TMidVec,TDot,TIDot,TFunc,Tb,TVerivVec,TInterval>
         (MyA, Myb, Y, procs, mypid, dim, rhs, lb1, lb2, nr, nc, myr, myc, numr, numc,
          ic, irsrc, icsrc, errc, commerrc, out, cfg, inner, ienc);

  if(m > n) {
    if(cfg.matrixMode) {
      if(breakrow != -1) 
        Y = Y(1,breakrow,1,RowLen(Y));
    } else {
      Y = Y( 1,n,1,RowLen(Y) );         
      if(inner)
        ienc = ienc(1,n,1,RowLen(ienc));
    }
  } else if(m < n) {
    if(cfg.matrixMode) {
      if(breakrow != -1)
        Y = Y(breakrow,ColLen(Y),1,RowLen(Y));
    } else {
      Y = Y( m+1,m+n,1,RowLen(Y) );     
      if(inner)
        ienc = ienc(m+1,m+n,1,RowLen(ienc));
    }
  }                   
}


/**
 * Start function for parallel linear system solver. This is the entry function for the case
 * that A and B are avaiable in process 0 and must be distributed.
 *
 * Two dimensional block cyclic distribution is used according to the
 * blocksize paramter nb.
 *
 * The result Y is an enclosure of the solution
 *
 *
 * @param[in]  A          The matrix of the system
 * @param[in]  b          Right hand side
 * @param[out] Y          Solution enclosure
 * @param[in]  m,n        Dimensions of Matrix A
 * @param[in]  procs      Number of parallel processes
 * @param[in]  mypid      Process ID
 * @param[out] errc       Error code
 * @param[out] commerrc   Error code
 * @param[in]  out        Output file stream
 * @param[in]  cfg        Struct containing configuration information for solver
 * @param[in]  inner      Compute inner solution?
 * @param[out] ienc       Computed inner solution
 */
template <typename TSysMat, typename TRhs, typename TSolution, typename TInverse, typename TC,
          typename TMidMat, typename TMidVec, typename TDot, typename TIDot, typename TFunc, typename Tb, 
          typename TVerivVec, typename TInterval>
void SolverStart(TSysMat& A, TRhs& b, TSolution& Y, int m, int n, 
                   int procs, int mypid, int& errc, int& commerrc,
                   ofstream& out, struct plssconfig cfg, bool inner, TSolution& ienc) {
  int dim;
  if(m!=n)
    dim = m+n;
  else
    dim = n; 

  int rhs = Ub(b,COL)-Lb(b,COL)+1;
  
  if(cfg.K != 1 && cfg.matrixMode) {
    out << "Warning: Activating Matrix Mode automatically sets K=1" << endl;
    cfg.K = 1;  
  }
  
  if(cfg.K == 1) {
    if(!cfg.matrixMode) {
      out << "Remark: Using K==1 automatically activates Matrix Mode" << endl;
      cfg.matrixMode = true;
    }
    
    if(cfg.lssparts != LSS_ONLY_PART_ONE) {
      out << "Warning: Only part one of solver can be used with K==1, adjusting configuration" << endl;
      cfg.lssparts = LSS_ONLY_PART_ONE;
    }
    
    if(inner) {
      out << "Warning: Inner enclosure can not be computed with K==1" << endl;
      inner = false;
    }
  }

  if(mypid == 0) {  
           
    if(n > m) {
        //under determined case
        out << "Underdetermined system, creating equivalent square system" << endl;
        TSysMat BIG_A(1,dim,1,dim);
        TRhs BIG_b(1,dim,1,rhs);
        TSolution BIG_y(1,dim,1,RowLen(Y));
        TSolution BIG_ienc;
        if(inner) BIG_ienc = TSolution(1,dim,1,RowLen(Y));
        int     i, j;
      
        BIG_A( 1,n, 1,m ) = transpherm(A);
          
        //BIG_A( 1,n, m+1,m+n ) = -Id(m);
        for( i = 1; i <= n; i++)
          for( j = m+1; j <= n+m; j++) 
            (j == i+2)? BIG_A[i][j]=-1 : BIG_A[i][j] = 0;   
      
        BIG_A( n+1,n+m, 1,m ) = 0.0;
        BIG_A( n+1,n+m, m+1,m+n ) = A;
        BIG_b( 1,n,1,rhs ) = 0.0;
        BIG_b( n+1,n+m,1,rhs ) = b;
        
        //delete old Data to save memory
        A = BIG_A;
        BIG_A = TSysMat();
        b = BIG_b;
        BIG_b = TRhs();
        Y = BIG_y;
        BIG_y = TSolution();
        if(inner) {
          ienc = BIG_ienc;
          BIG_ienc = TSolution();
        }
    } else if(m > n) {
        //over determined case
        out << "Overdetermined system, creating equivalent square system" << endl;      
        TSysMat BIG_A(1,dim,1,dim);
        TRhs BIG_b(1,dim,1,rhs);
        TSolution BIG_y(1,dim,1,RowLen(Y));
        TSolution BIG_ienc;
        if(inner) BIG_ienc = TSolution(1,dim,1,RowLen(Y));
        int     i, j;
        
        BIG_A( 1,m, 1,n ) = A;
        
        // BIG_A( 1,m, n+1,n+m ) = -Id(m);
        for( i = 1; i <= m; i++)  
          for( j = n+1; j <= n+m; j++) 
            (j == i+2)? BIG_A[i][j]=-1 : BIG_A[i][j] = 0;  
      
        BIG_A( m+1,m+n, 1,n ) = 0.0;
        BIG_A( m+1,m+n, n+1,n+m ) = transpherm(A);
        BIG_b( 1,m,1,rhs ) = b;
        BIG_b( m+1,m+n,1,rhs ) = 0.0;

        //delete old Data to save memory
        A = BIG_A;
        BIG_A = TSysMat();
        b = BIG_b;
        BIG_b = TRhs();
        Y = BIG_y;
        BIG_y = TSolution();  
        if(inner) {
          ienc = BIG_ienc;
          BIG_ienc = TSolution();
        }
    }  
    
  } //end mypid==0       
  
  //Indices of A           
  int lb1, ub1, lb2, ub2; 
  
  //Broadcast Matrix Dimensions
  if (mypid==0) {
    lb1 = Lb(A,ROW); ub1 = Ub(A,ROW);
    lb2 = Lb(A,COL); ub2 = Ub(A,COL);
  }
  commerrc=MPI_Bcast(&lb1,1,MPI_INT,0,MPI_COMM_WORLD);
  commerrc=MPI_Bcast(&ub1,1,MPI_INT,0,MPI_COMM_WORLD);
  commerrc=MPI_Bcast(&lb2,1,MPI_INT,0,MPI_COMM_WORLD);
  commerrc=MPI_Bcast(&ub2,1,MPI_INT,0,MPI_COMM_WORLD);

  double timer = time();
  MPI_Status status;

  if (lb1 != lb2) {
    errc = DimensionErr;
    return;
  }  

  out << endl;
  out << "Starting parallel solver..." << endl;
  out << "Dimension of system: " << m << "x" << n << endl;
  out << "Precision K: " << cfg.K << endl;
  out << "Number of processes: " << procs << endl;
  out << "ID of this process: " << mypid << endl << endl;

  int root = 0;   
  errc = NoError;
  
  //Values for data distribution (two dimensional block cyclic with
  //block dimension nb)
  int nr, nc;
  int myr, myc, numr, numc, ic, irsrc=0, icsrc=0;
  
  if(cfg.nb<=0) cfg.nb = 256;

  if(cfg.lssparts != LSS_ONLY_PART_TWO) {
    determine_grid(procs,nc,nr);
  } else {
    nc = 1; nr = procs;
    cfg.nb = dim / procs;
    if(dim%procs != 0) cfg.nb++;
  }
  
  TSysMat MyA;
  TSolution MyY;
  int numcrhs;

  //Distribute A and b 
  out << "Starting data distribution" << endl;
  
  int mybsize = rhs/procs;
  if(rhs%procs != 0  &&  mypid < rhs%procs)
    mybsize++;
  TRhs Myb(dim, mybsize);
  if(cfg.matrixMode) {
    distributeData(b, Myb, procs, mypid, errc, commerrc, out, nr, nc, cfg.nb,
                   myr, myc, numr, numcrhs, ic, dim, rhs);
    MyY = TSolution(numr,numcrhs);
  } else {
    int currcol = 1;  
    for(int s=1 ; s<=rhs ; s++) {
      if(mypid == 0) {
        if((s-1)%procs == 0) {
          Myb[Col(currcol)] = b[Col(s)];
          currcol++;
        } else {
          commerrc=MPI_Send(b,Lb(b,ROW),Ub(b,ROW),s,s,(s-1)%procs, 0, MPI_COMM_WORLD);
        }
      } else {
        if((s-1)%procs == mypid) {
          commerrc=MPI_Recv(Myb, Lb(Myb,ROW), Ub(Myb,ROW), currcol, currcol, 0, 0, MPI_COMM_WORLD, &status);
          currcol++;
        }
      }
    }
  }
   
  distributeData(A, MyA, procs, mypid, errc, commerrc, out, nr, nc, cfg.nb,
                 myr, myc, numr, numc, ic, dim, dim, !cfg.matrixMode);  
  
  //delete A and b to save memory
  A = TSysMat();
  b = TRhs();

  out << "Data distribution complete" << endl;  

  LSSMain<TSysMat,TRhs,TSolution,TInverse,TC,TMidMat,TMidVec,TDot,TIDot,TFunc,Tb,TVerivVec,TInterval>
         (MyA, Myb, Y, procs, mypid, dim, rhs, lb1, lb2, nr, nc, myr, myc, numr, numc,
          ic, irsrc, icsrc, errc, commerrc, out, cfg, inner, ienc);

  if(cfg.matrixMode) {
    collectData(MyY, Y, procs, mypid, errc, commerrc, out, nr, nc, cfg.nb,
                myr, myc, numr, numcrhs, ic, dim, rhs);  
    if(mypid == 0) {
      if(m == n) Y = MyY;
      if(m > n) {
        Y = MyY( 1,n,1,RowLen(Y) );         
        if(inner)
          ienc = ienc( 1,n,1,RowLen(Y) );
      } else if(m < n) {
        Y = MyY( m+1,m+n,1,RowLen(Y) );
        if(inner)
          ienc = ienc( m+1,m+n,1,RowLen(Y) );     
      }
    } else {
      Y = TSolution();
    }
      
  } else {
    
    if(m > n) {
      Y = Y( 1,n,1,RowLen(Y) );         
      if(inner)
        ienc = ienc( 1,n,1,RowLen(Y) );
    } else if(m < n) {
      Y = Y( m+1,m+n,1,RowLen(Y) );
      if(inner)
        ienc = ienc( m+1,m+n,1,RowLen(Y) );     
    }
    
  }

  
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real point system solver
 */
void plss(rmatrix& A, rmatrix& b, imatrix& X, int m, int n, 
          int procs, int mypid, int& errc, ofstream& out,
          struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  imatrix ienc;
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_rmatrix, rvector, ivector, interval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real point system solver
 */
void plss(rmatrix& A, rmatrix& b, imatrix& X, imatrix& Y, int m, int n, 
          int procs, int mypid, int& errc, ofstream& out,
          struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_rmatrix, rvector, ivector, interval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real point system solver
 */
void plss(get_Element_of_rmatrix A, get_Element_of_rmatrix b, imatrix& X, int m, int n, int rhs,
          int procs, int mypid, int& errc, ofstream& out,
          struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  imatrix ienc;
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_rmatrix, rvector, ivector, interval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real point system solver
 */
void plss(get_Element_of_rmatrix A, get_Element_of_rmatrix b, imatrix& X, imatrix& Y, int m, int n, int rhs,
          int procs, int mypid, int& errc, ofstream& out,
          struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<rmatrix, rmatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_rmatrix, rvector, ivector, interval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real interval system solver
 */
void pilss(imatrix& A, imatrix& b, imatrix& X, int m, int n, 
           int procs, int mypid, int& errc, ofstream& out,
           struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  imatrix ienc;
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_imatrix, ivector, ivector, interval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real interval system solver
 */
void pilss(imatrix& A, imatrix& b, imatrix& X, imatrix& Y, int m, int n, 
           int procs, int mypid, int& errc, ofstream& out,
           struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_imatrix, ivector, ivector, interval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real interval system solver
 */
void pilss(get_Element_of_imatrix A, get_Element_of_imatrix b, imatrix& X, int m, int n, int rhs,
           int procs, int mypid, int& errc, ofstream& out, 
           struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  imatrix ienc;
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_imatrix, ivector, ivector, interval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for real interval system solver
 */
void pilss(get_Element_of_imatrix A, get_Element_of_imatrix b, imatrix& X, imatrix& Y, int m, int n, int rhs,
           int procs, int mypid, int& errc, ofstream& out, 
           struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<imatrix, imatrix, imatrix, rmatrix, imatrix, rmatrix, rvector, dotprecision, idotprecision, get_Element_of_imatrix, ivector, ivector, interval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex point system solver
 */
void pclss(cmatrix& A, cmatrix& b, cimatrix& X, int m, int n, 
           int procs, int mypid, int& errc, ofstream& out, 
           struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  cimatrix ienc;
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cmatrix, cvector, civector, cinterval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex point system solver
 */
void pclss(cmatrix& A, cmatrix& b, cimatrix& X, cimatrix& Y, int m, int n, 
           int procs, int mypid, int& errc, ofstream& out, 
           struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cmatrix, cvector, civector, cinterval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex point system solver
 */
void pclss(get_Element_of_cmatrix A, get_Element_of_cmatrix b, cimatrix& X, int m, int n, int rhs,
           int procs, int mypid, int& errc, ofstream& out, struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  cimatrix ienc;
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cmatrix, cvector, civector, cinterval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex point system solver
 */
void pclss(get_Element_of_cmatrix A, get_Element_of_cmatrix b, cimatrix& X, cimatrix& Y, int m, int n, int rhs,
           int procs, int mypid, int& errc, ofstream& out, struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<cmatrix, cmatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cmatrix, cvector, civector, cinterval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex interval system solver
 */
void pcilss(cimatrix& A, cimatrix& b, cimatrix& X, int m, int n, 
            int procs, int mypid, int& errc, ofstream& out, struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  cimatrix ienc;
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cimatrix, civector, civector, cinterval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex interval system solver
 */
void pcilss(cimatrix& A, cimatrix& b, cimatrix& X, cimatrix& Y, int m, int n, 
            int procs, int mypid, int& errc, ofstream& out, struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cimatrix, civector, civector, cinterval>(A, b, X, m, n, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex interval system solver
 */
void pcilss(get_Element_of_cimatrix A, get_Element_of_cimatrix b, cimatrix& X, int m, int n, int rhs,
            int procs, int mypid, int& errc, ofstream& out, struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  cimatrix ienc;
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cimatrix, civector, civector, cinterval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, false, ienc);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}

// ----------------------------------------------------------------------------

/**
 * Entry function for complex interval system solver
 */
void pcilss(get_Element_of_cimatrix A, get_Element_of_cimatrix b, cimatrix& X, cimatrix& Y, int m, int n, int rhs,
            int procs, int mypid, int& errc, ofstream& out, struct plssconfig cfg) {
  if(m!=n) 
    init_CXSC_MPI();
  else
    init_CXSC_MPI();
  int commerrc = 0;
  SolverStart<cimatrix, cimatrix, cimatrix, cmatrix, cimatrix, cmatrix, cvector, cdotprecision, cidotprecision, get_Element_of_cimatrix, civector, civector, cinterval>(A, b, X, m, n, rhs, procs, mypid, errc, commerrc, out, cfg, true, Y);
  if(commerrc != 0) errc = CommErr;
  finish_CXSC_MPI();
}
