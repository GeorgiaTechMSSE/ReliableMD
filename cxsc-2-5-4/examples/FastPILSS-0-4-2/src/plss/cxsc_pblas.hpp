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
 
#ifndef CXSC_PBLAS_HEADER_INCLUDED
#define CXSC_PBLAS_HEADER_INCLUDED

#include <rmatrix.hpp>
#include <imatrix.hpp>
#include <cmatrix.hpp>
#include <cimatrix.hpp>

using namespace cxsc;
using namespace std;

#ifdef Unchanged
  #define pdgemm_ pdgemm
  #define pzgemm_ pzgemm
  #define descinit_ descinit
#elif defined(Add_)
  //Nothing to do
#elif defined(Uppercase)
  #define dgemm_ DGEMM
  #define zgemm_ ZGEMM
  #define descinit_ DESCINIT
#endif


//Declaration of PBLAS-routines
extern "C" {
   void pdgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, 
                double *A, int *IA, int *JA, int *DESCA, 
                double *B, int *IB, int *JB, int *DESCB, 
 		        double *BETA, double *C, int *IC, int *JC, int *DESCC);
   void pzgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, double *ALPHA, 
                double *A, int *IA, int *JA, int *DESCA, 
                double *B, int *IB, int *JB, int *DESCB, 
 		        double *BETA, double *C, int *IC, int *JC, int *DESCC);
   void descinit_ (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info); 
 		        
}

//! Multiplication of two real matrices
void matmul(rmatrix &A, rmatrix &B, rmatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Multiplication of two real matrices, computes enclosure of result
void matmul(rmatrix &A, rmatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product real matrix times interval matrix
void matmul(rmatrix &A, imatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product real interval matrix times point matrix
void matmul(imatrix &A, rmatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product real interval matrix times interval matrix
void matmul(imatrix &A, imatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);

//! Computes enclosure of a product complex matrix times complex matrix
void matmul(cmatrix &A, cmatrix &B, cmatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product complex matrix times complex matrix
void matmul(cmatrix &A, cmatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product complex matrix times cinterval matrix
void matmul(cmatrix &A, cimatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product cinterval matrix times complex matrix
void matmul(cimatrix &A, cmatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);
//! Computes enclosure of a product cinterval matrix times cinterval matrix
void matmul(cimatrix &A, cimatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg);

//! Computes enclosure of (I-A*B), for A and B of type rmatrix
void IminusAB(rmatrix &A, rmatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg);
//! Computes enclosure of (I-A*B), for A of type rmatrix and B of type imatrix
void IminusAB(rmatrix &A, imatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg);
//! Computes enclosure of (I-A*B), for A and B of type cmatrix
void IminusAB(cmatrix &A, cmatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg);
//! Computes enclosure of (I-A*B), for A of type cmatrix and B of type cimatrix
void IminusAB(cmatrix &A, cimatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg);



#endif
