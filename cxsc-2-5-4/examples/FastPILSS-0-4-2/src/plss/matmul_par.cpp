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
#include <rmatrix.hpp>  
#include <imatrix.hpp>  
#include <cmatrix.hpp>  
#include <cimatrix.hpp>  
#include <cmath>
#include <iostream>     
#include "cxsc_mpicomm.hpp"
#include "matmul_par.hpp"        
#include <dot.hpp>
#include <idot.hpp>
#include <cdot.hpp>
#include <cidot.hpp>


using namespace cxsc;
using namespace std;



// ----------------------------------------------------------------------------

/**
 * Parallel Matrix multiplication rmatrix*rmatrix using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A         Input matrix
 * @param[in]   B         Input matrix
 * @param[out]  R         Result
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */
void MatMul(rmatrix& A, rmatrix& B, rmatrix& R, int mypid, int procs, int *ind1,
            int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting matrix multiplication..." << endl; 
  
  dotprecision dot(0.0);
  dot.set_k(K);
  int n = Ub(A,2) - Lb(A,2) + 1;
  rmatrix tmp;  
    
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;
     
     tmp = rmatrix(1,n,ind1[p],ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }
     
     //Compute own result of current block
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
           dot = 0.0;
           accumulate_approx(dot, A[Row(i)], tmp[Col(j)]);
           R[i][j] = rnd(dot);
        }
     }
  }

  ausg << "Matrix multiplication ok. " << endl;
}

/**
 * Parallel Matrix multiplication cmatrix*cmatrix using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A         Input matrix
 * @param[in]   B         Input matrix
 * @param[out]  R         Result
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */
void MatMul(cmatrix& A, cmatrix& B, cmatrix& R, int mypid, int procs, int *ind1,
            int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting matrix multiplication..." << endl;
  
  cdotprecision dot;
  dot.set_k(K);
  int n = Ub(A,2) - Lb(A,2) + 1;
  cmatrix tmp;  
    
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;
     
     tmp = cmatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }
     
     //Compute own result of current block     
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
           dot = 0.0;
           accumulate_approx(dot, A[Row(i)], tmp[Col(j)]);
           R[i][j] = rnd(dot);
        }
     }
  }

  ausg << "Matrix multiplication ok. " << endl;
}

/**
 * Parallel calculation of #*(A*B-#*(A*B)) using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A         Input matrix
 * @param[in]   B         Input matrix
 * @param[out]  C1        Result of #*(A*B)
 * @param[out]  C2        Result of #*(A*B-#*(A*B))
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */
void ABminusRndAB(rmatrix& A, rmatrix& B, rmatrix& C1, rmatrix &C2, int mypid, int procs, int *ind1,
                  int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting ABMinusRndAB..." << endl;

  
  dotprecision dot;
  dot.set_k(K);
  int n = Ub(A,2) - Lb(A,2) + 1;
  rmatrix tmp;  
  
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;

     tmp = rmatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;             
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }

     //Compute own result of current block          
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
          dot = 0.0;
          accumulate_approx(dot, A[Row(i)],tmp[Col(j)]);
          C1[i][j] = rnd(dot);
          dot -= C1[i][j];
          C2[i][j] = rnd(dot);
        }
     }
  }

  ausg << "ABMinusRndAB ok. " << endl;

}

/**
 * Parallel calculation of #*(A*B-#*(A*B)) using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A         Input matrix
 * @param[in]   B         Input matrix
 * @param[out]  C1        Result of #*(A*B)
 * @param[out]  C2        Result of #*(A*B-#*(A*B))
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */
void ABminusRndAB(cmatrix& A, cmatrix& B, cmatrix& C1, cmatrix &C2, int mypid, int procs, int *ind1,
                  int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting ABMinusRndAB..." << endl;
  
  cdotprecision dot;
  dot.set_k(K);
  int n = Ub(A,2) - Lb(A,2) + 1;
  cmatrix tmp;  
  
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;

     tmp = cmatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;             
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }

     //Compute own result of current block          
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
          dot = 0.0;
          accumulate_approx(dot, A[Row(i)],tmp[Col(j)]);
          C1[i][j] = rnd(dot);
          dot += -C1[i][j];
          C2[i][j] = rnd(dot);
        }
     }
  }

  ausg << "ABMinusRndAB ok. " << endl;

}

/**
 * Parallel calculation of I-A1*B-A2*B using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A1        Input matrix
 * @param[in]   A2        Input matrix
 * @param[out]  B         Input matrix
 * @param[out]  R         Result
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */

void IminusA1A2B(rmatrix& A1, rmatrix& A2, rmatrix& B, imatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting IminusA1A2B..." << endl;
  
  idotprecision idot;
  idot.set_k(K);
  int n = Ub(A1,2) - Lb(A1,2) + 1;
  rmatrix tmp;  

  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;

     tmp = rmatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;             
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }

     //Compute own result of current block          
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
          idot = 0.0;
          accumulate(idot,A1[Row(i)],tmp[Col(j)]);
          accumulate(idot,A2[Row(i)],tmp[Col(j)]);          
          R[i][j] = -rnd(idot);
          if(i==j) R[i][j] += 1.0;
        }
     }
  }

  ausg << "IminusA1A2B ok. " << endl;

}

/**
 * Parallel calculation of I-A1*B-A2*B using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A1        Input matrix
 * @param[in]   A2        Input matrix
 * @param[out]  B         Input matrix
 * @param[out]  R         Result
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */

void IminusA1A2B(cmatrix& A1, cmatrix& A2, cmatrix& B, cimatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting IminusA1A2B..." << endl;
  
  cidotprecision idot;
  idot.set_k(K);
  int n = Ub(A1,2) - Lb(A1,2) + 1;
  cmatrix tmp;  
  
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;

     tmp = cmatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;             
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }

     //Compute own result of current block          
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
          idot = 0.0;
          accumulate(idot,A1[Row(i)],tmp[Col(j)]);
          accumulate(idot,A2[Row(i)],tmp[Col(j)]);          
          R[i][j] = -rnd(idot);
          if(i==j) R[i][j] += 1.0;
        }
     }
  }

  ausg << "IminusA1A2B ok. " << endl;

}

/**
 * Parallel calculation of I-A1*B-A2*B using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A1        Input matrix
 * @param[in]   A2        Input matrix
 * @param[out]  B         Input matrix
 * @param[out]  R         Result
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */
void IminusA1A2B(rmatrix& A1, rmatrix& A2, imatrix& B, imatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting IminusA1A2B..." << endl;
  
  idotprecision idot;
  idot.set_k(K);
  int n = Ub(A1,2) - Lb(A1,2) + 1;
  imatrix tmp;  
  
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;

     tmp = imatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;             
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }

     //Compute own result of current block          
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
          idot = 0.0;
          accumulate(idot,A1[Row(i)],tmp[Col(j)]);
          accumulate(idot,A2[Row(i)],tmp[Col(j)]);          
          R[i][j] = -rnd(idot);
          if(i==j) R[i][j] += 1.0;
        }
     }
  }

  ausg << "IminusA1A2B ok. " << endl;

}

/**
 * Parallel calculation of I-A1*B-A2*B using higher precision scalar
 * product. Each processor owns a block (ind1[p],ind2[p],1,n) of the input
 * matrices an computes the corresponding block of the result R
 * 
 * @param[in]   A1        Input matrix
 * @param[in]   A2        Input matrix
 * @param[out]  B         Input matrix
 * @param[out]  R         Result
 * @param[in]   mypid     Process ID
 * @param[in]   procs     Number of processes
 * @param[in]   ind1      Start indices for each processor
 * @param[in]   ind2      End indices for each processor
 * @param[out]  errc      Error code
 * @param[out]  commerrc  Error code
 * @param[in]   ausg      Output file stream
 * @param[in]   K         Precision for DotK
 * @param[in]   np        Number of threads to be used by DotK
 */

void IminusA1A2B(cmatrix& A1, cmatrix& A2, cimatrix& B, cimatrix &R, int mypid, int procs, int *ind1,
                 int *ind2, int& errc, int& commerrc, ofstream& ausg, int K, int np) 
{

  MPI_Status status; 
  ausg << "Starting IminusA1A2B..." << endl;
  
  cidotprecision idot;
  idot.set_k(K);
  int n = Ub(A1,2) - Lb(A1,2) + 1;
  cimatrix tmp;  
  
  for(int p=0 ; p<procs ; p++) {
     if(ind1[p] == -1) continue;

     tmp = cimatrix(1,n,ind1[p], ind2[p]);

     //Broadcast own part of current vertical block 
     for(int i=0 ; i<procs ; i++) {
       if(ind1[i] == -1) continue;             
       if(mypid == i) {
         tmp(ind1[i],ind2[i],ind1[p],ind2[p]) = B(ind1[i],ind2[i],ind1[p],ind2[p]);
       }
       MPI_Bcast(tmp, ind1[i],ind2[i],ind1[p],ind2[p], i, MPI_COMM_WORLD);
     }

     //Compute own result of current block          
     for(int i=ind1[mypid] ; i<=ind2[mypid] ; i++) {
        for(int j=Lb(tmp,2) ; j<=Ub(tmp,2) ; j++) {
          idot = 0.0;
          accumulate(idot,A1[Row(i)],tmp[Col(j)]);
          accumulate(idot,A2[Row(i)],tmp[Col(j)]);          
          R[i][j] = -rnd(idot);
          if(i==j) R[i][j] += 1.0;
        }
     }
  }

  ausg << "IminusA1A2B ok. " << endl;

}
