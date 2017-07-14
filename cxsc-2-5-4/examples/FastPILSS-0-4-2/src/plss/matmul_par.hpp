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
 
#ifndef __MATMUL_PAR_HPP
#define __MATMUL_PAR_HPP

#include <fstream> 
#include <vector>      
#include <rmatrix.hpp> 
#include <imatrix.hpp> 

using namespace cxsc;
using namespace std;

//! Parallel matrix multiplication using higher precision scalar product
void MatMul(rmatrix& A, rmatrix& B, rmatrix& R, int mypid, int procs, int *ind1, 
            int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np); 

//! Parallel matrix multiplication using higher precision scalar product
void MatMul(cmatrix& A, cmatrix& B, cmatrix& R, int mypid, int procs, int *ind1, 
            int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np); 

//! Parallel computation of #*(A*B-#*(A*B)) higher precision scalar product
void ABminusRndAB(rmatrix& A, rmatrix& B, rmatrix& C1, rmatrix &C2, int mypid, int procs, int *ind1,
                  int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np);

//! Parallel computation of #*(A*B-#*(A*B)) using higher precision scalar product
void ABminusRndAB(cmatrix& A, cmatrix& B, cmatrix& C1, cmatrix &C2, int mypid, int procs, int *ind1,
                  int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np);

//! Parallel computation of I-A1*B-A2*B using higher precision scalar product                  
void IminusA1A2B(rmatrix& A1, rmatrix& A2, rmatrix& B, imatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np);                                

//! Parallel computation of I-A1*B-A2*B using higher precision scalar product                  
void IminusA1A2B(cmatrix& A1, cmatrix& A2, cmatrix& B, cimatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np);                                

//! Parallel computation of I-A1*B-A2*B using higher precision scalar product                  
void IminusA1A2B(rmatrix& A1, rmatrix& A2, imatrix& B, imatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np);                                

//! Parallel computation of I-A1*B-A2*B using higher precision scalar product                  
void IminusA1A2B(cmatrix& A1, cmatrix& A2, cimatrix& B, cimatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np);                                
            

#endif
