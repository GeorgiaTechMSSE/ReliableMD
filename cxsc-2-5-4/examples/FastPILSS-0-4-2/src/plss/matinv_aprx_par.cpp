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
 */

#include <mpi.h>

#include <iostream>
#include <fstream>

#include <rmatrix.hpp>
#include <cmatrix.hpp>
#include "cxsc_mpicomm.hpp"   // MPI Communication for C-XSC
#include "matinv_aprx_par.hpp"      // Local header


using namespace cxsc;
using namespace std;



// ----------------------------------------------------------------------------

/**
 * Parallel computation of an approximate inverse of a real matrix. Expects
 * the Data to be distributed according to ScaLAPACK standards.
 * 
 * @param[out] MyR        Approximate inverse
 * @param[out] MyA        Input Matrix
 * @param[in]  procs      Number of parallel processes
 * @param[in]  mypid      Process ID
 * @param[out] errc       Error code
 * @param[out] commerrc   Error code
 * @param[in]  ausg       Output file stream
 * @param[in]  nr         Number of rows of processor grid
 * @param[in]  nc         Number of columns of processor grid
 * @param[in]  nb         Blocksize of distribution
 * @param[in]  myr        Own row in process grid
 * @param[in]  myc        Own column in process grid
 * @param[in]  numr       Row size of own data
 * @param[in]  numc       Column size of own data
 * @param[in]  ic         SCALAPACK context
 * @param[in]  n          Dimension of global matrix
 */
void MatInv( rmatrix& MyR, rmatrix& MyA, int procs, int mypid, int& errc, 
             int& commerrc, ofstream& ausg, int nr, int nc, int nb,
             int myr, int myc, int numr, int numc, int ic, int n) 

{
  MPI_Status status;
  ausg << "Starting matrix inversion" << endl;
     
  int irsrc = 0, icsrc = 0;

  //Process not part of process grid
  if(mypid >= nr*nc) {
     ausg << "Inversion finished" << endl;
     return;
  }  

  //Copy into double array for SCALAPACK (transposed!)
  double *DA = new double[numr*numc];  
  for(int i=0 ; i<numc ; i++)
    for(int j=0 ; j<numr ; j++)
       DA[i*numr+j] = _double(MyA[Lb(MyA,1)+j][Lb(MyA,2)+i]);

  //Initialize descriptor for SCALAPACK
  int descA[9];
  int lld = numr;
  int info;
  if(lld == 0) lld = 1;
  descinit_(descA, &n, &n, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 

  //LU-Decomp.
  int ia=1, ja=1;
  int *ipiv = new int[numr + nb];

  pdgetrf_( &n, &n, DA, &ia, &ja, descA, ipiv, &info );  

  if(info != 0) {
    errc = 1;
    ausg << "Error during Inversion!" << endl;
    return;
  }
  
  //Inversion
  int lwork = n*nb;
  double *work = new double[lwork];
  int liwork; 
  if((nr==nc) || (nb>=sqrt((double)n*nc)))
    liwork = 4*(n+nb);
  else
    liwork = 4*(n+(int)(sqrt((double)n*nc))+1);
  int *iwork = new int[liwork];
  pdgetri_( &n, DA, &ia, &ja, descA, ipiv, work, &lwork, iwork, &liwork, &info);  

  if(info != 0) {
    errc = 1;
    ausg << "Error during Inversion!" << endl;    
    return;
  }

  delete[] ipiv;
  delete[] work;
  delete[] iwork;

  //Copy result into MyR  
  Resize(MyR,numr,numc);
  for(int i=0 ; i<numc ; i++)
    for(int j=0 ; j<numr ; j++)
       MyR[j+1][i+1] = DA[i*numr+j];
  delete[] DA;

      
  ausg << "Inversion finished" << endl;

}


/**
 * Parallel computation of an approximate inverse of a complex matrix. Expects
 * the Data to be distributed according to ScaLAPACK standards.
 * 
 * @param[out] MyR        Approximate inverse
 * @param[out] MyA        Input Matrix
 * @param[in]  procs      Number of parallel processes
 * @param[in]  mypid      Process ID
 * @param[out] errc       Error code
 * @param[out] commerrc   Error code
 * @param[in]  ausg       Output file stream
 * @param[in]  nr         Number of rows of processor grid
 * @param[in]  nc         Number of columns of processor grid
 * @param[in]  nb         Blocksize of distribution
 * @param[in]  myr        Own row in process grid
 * @param[in]  myc        Own column in process grid
 * @param[in]  numr       Row size of own data
 * @param[in]  numc       Column size of own data
 * @param[in]  ic         SCALAPACK context
 * @param[in]  n          Dimension of global matrix
 */
void MatInv( cmatrix& MyR, cmatrix& MyA, int procs, int mypid, int& errc, 
             int& commerrc, ofstream& ausg, int nr, int nc, int nb,
             int myr, int myc, int numr, int numc, int ic, int n) 

{
  //----------------------------------------------------------------------------
  // MPI preparation
  //----------------------------------------------------------------------------
  MPI_Status status;
  ausg << "Starting matrix inversion" << endl;
     
  int irsrc = 0, icsrc = 0;

  if(mypid >= nr*nc) {
     ausg << "Inversion finished" << endl;
     return;
  }  
  
  //Copy into double array for SCALAPACK (transposed!)
  double *DA = new double[2*numr*numc];  
  for(int i=0 ; i<numc ; i++)
    for(int j=0 ; j<numr ; j++) {
       int i1 = 2*(i*numr+j);
       int i2 = i1+1;
       DA[i1] = _double(Re(MyA[Lb(MyA,1)+j][Lb(MyA,2)+i]));
       DA[i2] = _double(Im(MyA[Lb(MyA,1)+j][Lb(MyA,2)+i]));       
    }

  //Initialize SCALAPACK descriptor
  int descA[9];
  int lld = numr;
  int info;
  if(lld == 0) lld = 1;

  descinit_(descA, &n, &n, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
  //LU-Decomp.
  int ia=1, ja=1;
  int *ipiv = new int[2*(numr + nb)];
  pzgetrf_( &n, &n, DA, &ia, &ja, descA, ipiv, &info );  

  if(info != 0) {
    errc = 1;
    ausg << "Error during Inversion!" << endl;
    return;
  }
  
  //Inversion
  int lwork = 2*n*nb;
  double *work = new double[lwork];
  int liwork;
  if((nr==nc) || (nb>=sqrt((double)n*nc)))
    liwork = 4*(2*(n+nb));
  else
    liwork = 4*(n+(int)(sqrt((double)n*nc))+1);
  int *iwork = new int[liwork];
  pzgetri_( &n, DA, &ia, &ja, descA, ipiv, work, &lwork, iwork, &liwork, &info);  
  
  if(info != 0) {
    errc = 1;
    ausg << "Error during Inversion!" << endl;    
    return;
  }

  delete[] ipiv;
  delete[] work;
  delete[] iwork;
   
  //Copy result into MyR
  Resize(MyR,numr,numc);
  for(int i=0 ; i<numc ; i++)
    for(int j=0 ; j<numr ; j++) {
       int i1 = 2*(i*numr+j);
       int i2 = i1+1;
       SetRe(MyR[j+1][i+1], DA[i1]);
       SetIm(MyR[j+1][i+1], DA[i2]);       
    }

  delete[] DA;
  

  
      
  ausg << "Inversion finished" << endl;

}


