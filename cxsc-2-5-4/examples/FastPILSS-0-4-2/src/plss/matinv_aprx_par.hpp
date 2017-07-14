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

#ifndef __MATINV_APRX_PAR_HPP
#define __MATINV_APRX_PAR_HPP

#include <rmatrix.hpp>   

using namespace cxsc;
using namespace std;

#ifdef Unchanged
  #define pdgetri_ pdgetri
  #define pdgetrf_ pdgetrf
  #define pzgetri_ pzgetri
  #define pzgetrf_ pzgetrf
  #define descinit_ descinit
#elif defined(Add_)
  //Nothing to do
#elif defined(Uppercase)
  #define pdgetri_ PDGETRI
  #define pdgetrf_ PDGETRF
  #define pzgetri_ PZGETRI
  #define pzgetrf_ PZGETRF
  #define descinit_ DESCINIT
#endif

//Declaration of SCALAPACK routines
extern "C"  { 
   void descinit_ (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info); 
   void pdgetrf_( int *m, int *n, double *A, int *ia, int *ja, int *desca, int *ipiv, int *info );   
   void pdgetri_( int *n, double *A, int *ia, int *ja, int *descA, int *ipiv, double *work, int *lwork, int *iwork, int *liwork, int *info );   
   void pzgetrf_( int *m, int *n, double *A, int *ia, int *ja, int *desca, int *ipiv, int *info );   
   void pzgetri_( int *n, double *A, int *ia, int *ja, int *descA, int *ipiv, double *work, int *lwork, int *iwork, int *liwork, int *info );   
} 

//! Parallel inversion of a real matrix
void MatInv( rmatrix& MyR, rmatrix& MyA, int procs, int mypid, int& errc, 
             int& commerrc, ofstream& ausg, int nr, int nc, int nb,
             int myr, int myc, int numr, int numc, int ic, int n) ;

//! Parallel inversion of a complex matrix
void MatInv( cmatrix& MyR, cmatrix& MyA, int procs, int mypid, int& errc, 
             int& commerrc, ofstream& ausg, int nr, int nc, int nb,
             int myr, int myc, int numr, int numc, int ic, int n) ;
        
#endif
