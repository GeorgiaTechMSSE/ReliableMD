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
 
#include "cxsc_pblas.hpp"
#include <fenv.h>
#include <iomanip>
#include <fstream> 

#define ROUND_NEAR  FE_TONEAREST
#define ROUND_UP    FE_UPWARD
#define ROUND_DOWN  FE_DOWNWARD

using namespace cxsc;
using namespace std;

//Upper bound of sqrt(2)*0.5
static const real cr = mulu( Sup(sqrt(interval(2.0,2.0))) , 0.5 );

/*void ausgabe(double* A, int m , int n, ofstream& ausg) {
    for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         ausg << A[(i-1)*n+(j-1)] << " ";
      }
      ausg << endl;
    }
    ausg << endl;
} */

//! Sets the rounding mode of the processor
/*
  To be able to compute rigid enclosures using BLAS routines, it is necessary
  to control the rounding mode of the processor. This function can set the
  rounding mode to "round up", "round down" or "round to nearest".
  
  \param rnd Sets the rounding mode (-1 for round down, 0 for round to nearest,
             1 for round up)
*/
void setRound( int rnd ) {

  switch (rnd) {
    /* round up */
    case  1 :  fesetround(ROUND_UP);   break;
    /* round to nearest */
    case  0 :  fesetround(ROUND_NEAR); break;
    /* round down */
    case -1 :  fesetround(ROUND_DOWN); break;
    /* round to nearest */
    default :  fesetround(ROUND_NEAR); break;
  }

} 


/*!
  The product C=A*B for A and B of type rmatrix is computed using BLAS. The
  result ist NOT accurate and will be very bad for ill conditioned matrices.
  
  \param A first operator
  \param B second operator
  \param C The result C=A*B, with rounding errors
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
*/

void matmul(rmatrix &A, rmatrix &B, rmatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
 // ausg  << "matmul rr" << endl;
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     

   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;


   //Copy A and B into double arrays   
   double *DA = new double[m*n];
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*m+(i-1)] = _double(A[lbA_i+i-1][lbA_j+j-1]);
      }
   }        
   A = rmatrix();

   double *DB = new double[o*p];
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         DB[(j-1)*o+(i-1)] = _double(B[lbB_i+i-1][lbB_j+j-1]);         
      }
   }
   B = rmatrix();
  
  int descA[9], descB[9], descC[9];
  int lld = numr;
  int info;
  int irsrc = 0, icsrc = 0;
  if(lld == 0) lld = 1;

  descinit_(descA, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descB, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
  int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
  char *TRANS=(char*)"No Transpose";
  double alpha = 1.0, beta = 0.0;

   //Multiply using PBLAS 
  double *DC = new double[o*p];
  pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DB, &ib, &jb, 
          descB, &beta, DC, &icc, &jc, descC);

  A = rmatrix(m,n);
  for(int i=1 ; i<=m ; i++) {
     for(int j=1 ; j<=n ; j++) {
       A[lbA_i+i-1][lbA_j+j-1] = DA[(j-1)*m+(i-1)];
     }
  }        
  delete[] DA;

  B = rmatrix(o,p);
  for(int i=1 ; i<=o ; i++) {
     for(int j=1 ; j<=p ; j++) {
        B[lbB_i+i-1][lbB_j+j-1] = DB[(j-1)*o+(i-1)];
     }
  }
  delete[] DB;

  //Copy result into C

  for(int i=1 ; i<=o ; i++) {
     for(int j=1 ; j<=p ; j++) {
        C[lbC_i+i-1][lbC_j+j-1] = DC[(j-1)*o+(i-1)];
     }
  }
  delete[] DC;

 // ausg  << "matmul rr ende" << endl;
}


/*!
  The product C=A*B for A and B of type rmatrix is computed using BLAS. C is
  an interval matrix enclosing the correct result. However, the diameter of
  C can get very large depending of the condition of A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(rmatrix &A, rmatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
  //ausg << "matmul rr=i" << endl;  
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     

   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;

   
   C = imatrix();

   //Copy A and B into double array    
   double *DA = new double[m*n];
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*m+(i-1)] = _double(A[lbA_i+i-1][lbA_j+j-1]);
      }
   }
   A = rmatrix();

   double *DB = new double[o*p];
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         DB[(j-1)*o+(i-1)] = _double(B[lbB_i+i-1][lbB_j+j-1]);         
      }
   }
   B= rmatrix();

  int descA[9], descB[9], descC[9];
  int lld = numr;
  int info;
  int irsrc = 0, icsrc = 0;
  if(lld == 0) lld = 1;

  descinit_(descA, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descB, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
  int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
  char *TRANS=(char*)"No Transpose";
  double alpha = 1.0, beta = 0.0;

   //Compute Infimum
  double *DC = new double[o*p];
  setRound(-1);
  pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DB, &ib, &jb, 
          descB, &beta, DC, &icc, &jc, descC);

  
  //Copy infimum into C
  C = imatrix(o,p);
  for(int i=1 ; i<=p ; i++) {
     for(int j=1 ; j<=o ; j++) {
        C[lbC_i+j-1][lbC_j+i-1] = DC[(i-1)*o+(j-1)];
     }
  }        

  //Compute supremum
  setRound(1);
  pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DB, &ia, &ja, 
          descB, &beta, DC, &ia, &ja, descC);

  //copy supremum into C
  setRound(0);
  for(int i=1 ; i<=p ; i++) {
     for(int j=1 ; j<=o ; j++) {
        SetSup(C[lbC_i+j-1][lbC_j+i-1], DC[(i-1)*o+(j-1)]);
     }
  }        
  delete[] DC;

  //Copy A and B into double array    
  A = rmatrix(m,n);
  for(int i=1 ; i<=m ; i++) {
     for(int j=1 ; j<=n ; j++) {
        A[lbA_i+i-1][lbA_j+j-1] = DA[(j-1)*m+(i-1)];
     }
  }
  delete[] DA;

  B = rmatrix(o,p);
  for(int i=1 ; i<=o ; i++) {
     for(int j=1 ; j<=p ; j++) {
        B[lbB_i+i-1][lbB_j+j-1] = DB[(j-1)*o+(i-1)];
     }
  }
  delete[] DB;
  
 // ausg << "matmul rr=i ende" << endl;  
}


/*!
  The product C=A*B for A of type rmatrix B of type imatrix is computed using 
  BLAS. C is an interval matrix enclosing the correct result. However, the 
  diameter of C can get very large depending of the condition of A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(rmatrix &A, imatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
 // ausg << "matmul ri" << endl;
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;
   
   C = imatrix();

   int descA[9], descBmid[9], descBrad[9], descC1[9], descC2[9];
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descA, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBmid, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBrad, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC1, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC2, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;

   //copy A into double array
   double *DA = new double[m*n];
   setRound(1);
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*m+(i-1)] = _double(A[lbA_i+i-1][lbA_j+j-1]);
      }
   }         

   //compute mid(B) and rad(B) and copy into double arrays
   double *DBmid = new double[o*p];
   double *DBrad = new double[o*p];   
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         DBmid[(j-1)*m+(i-1)] = _double(Inf(B[lbB_i+i-1][lbB_j+j-1]) + 0.5*(Sup(B[lbB_i+i-1][lbB_j+j-1]) - Inf(B[lbB_i+i-1][lbB_j+j-1])));         
         DBrad[(j-1)*m+(i-1)] = _double(DBmid[(j-1)*o+(i-1)] - Inf(B[lbB_i+i-1][lbB_j+j-1]));    
      }
   }          

   //compute lower bound for mid of result   
   double *DC1   = new double[o*p];
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DBmid, &ib, &jb, 
           descBmid, &beta, DC1, &icc, &jc, descC1);

   //compute upper bound for mid of result
   double *DC2   = new double[o*p];   
   setRound(1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DBmid, &ib, &jb, 
           descBmid, &beta, DC2, &icc, &jc, descC2);

   delete[] DBmid;

   
   rmatrix Cmid(o,p);   
   
   //compute mid of result
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         Cmid[j][i] = DC1[(i-1)*o+(j-1)] + 0.5 * (DC2[(i-1)*o+(j-1)] - DC1[(i-1)*o+(j-1)]); 
      }
   }       


   //copy abs(A) into double array
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*m+(i-1)] = _double(abs(A[lbA_i+i-1][lbA_j+j-1]));
      }
   }   


   //compute upper bound for rad of result
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DBrad, &ib, &jb, 
           descBrad, &beta, DC2, &icc, &jc, descC2);

   delete[] DA;
   delete[] DBrad;  
   
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         DC2[(i-1)*o+(j-1)] += _double(Cmid[j][i]) - DC1[(i-1)*o+(j-1)];
      }
   }

   setRound(-1);

   //compute infimum of result and copy into C
   C = imatrix(o,p);
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         C[lbC_i+j-1][lbC_j+i-1] = Cmid[j][i] - DC2[(i-1)*o+(j-1)]; 
      }
   }        
   
   setRound(1);

   //compute supremum of result and copy into C
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         SetSup(C[lbC_i+j-1][lbC_j+i-1], Cmid[j][i] + DC2[(i-1)*o+(j-1)]); 
      }
   }        
   
   setRound(0);
   
 
   delete[] DC1;          
   delete[] DC2;             

  //   ausg << "matmul ri ende" << endl;
}



/*!
  The product C=A*B for A of type imatrix, B of type imatrix is computed using 
  PBLAS. C is an interval matrix enclosing the correct result. However, the 
  diameter of C can get very large depending of the condition of A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(imatrix &A, imatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
//  ausg << "matmul ii" << endl;
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;
   
   C = imatrix();

   int descAmid[9], descAabs[9], descArad[9], descBmid[9], descBabs[9], descBrad[9], descC1[9], descC2[9];
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descAmid, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAabs, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descArad, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBmid, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBabs, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBrad, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC1, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC2, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;
   
   //copy A into double array
   double *DAmid = new double[m*n];
   double *DArad = new double[m*n];   
   double *DAabs = new double[m*n];     
   
   setRound(1);
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DAmid[(j-1)*m+(i-1)] = _double(Inf(A[lbA_i+i-1][lbA_j+j-1]) + 0.5*(Sup(A[lbA_i+i-1][lbA_j+j-1]) - Inf(A[lbA_i+i-1][lbA_j+j-1])));         
         DArad[(j-1)*m+(i-1)] = _double(DAmid[(j-1)*m+(i-1)] - Inf(A[lbA_i+i-1][lbA_j+j-1]));    
	 DAabs[(j-1)*m+(i-1)] = abs(DAmid[(j-1)*m+(i-1)]);
      }
   }         

   //compute mid(B) and rad(B) and copy into double arrays
   double *DBmid = new double[o*p];
   double *DBrad = new double[o*p];   
   double *DBabs = new double[o*p];
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         DBmid[(j-1)*o+(i-1)] = _double(Inf(B[lbB_i+i-1][lbB_j+j-1]) + 0.5*(Sup(B[lbB_i+i-1][lbB_j+j-1]) - Inf(B[lbB_i+i-1][lbB_j+j-1])));         
         DBrad[(j-1)*o+(i-1)] = _double(DBmid[(j-1)*o+(i-1)] - Inf(B[lbB_i+i-1][lbB_j+j-1]));    
	 DBabs[(j-1)*o+(i-1)] = abs(DBmid[(j-1)*o+(i-1)]);
      }
   }          

   //compute lower bound for mid of result   
   double *DC1   = new double[o*p];
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAmid, &ia, &ja, descAmid, DBmid, &ib, &jb, 
           descBmid, &beta, DC1, &icc, &jc, descC1);

   //compute upper bound for mid of result
   double *DC2   = new double[o*p];   
   setRound(1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAmid, &ia, &ja, descAmid, DBmid, &ib, &jb, 
           descBmid, &beta, DC2, &icc, &jc, descC2);

   delete[] DBmid;
   delete[] DAmid;
   
   rmatrix Cmid(o,p);   
   rmatrix Crad(o,p);
   
   //compute mid of result
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         Cmid[j][i] = DC1[(i-1)*o+(j-1)] + 0.5 * (DC2[(i-1)*o+(j-1)] - DC1[(i-1)*o+(j-1)]); 
	 Crad[j][i] = Cmid[j][i] - DC1[(i-1)*o+(j-1)];
      }
   }       

   //copy abs(B) into double array
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         DBabs[(i-1)*o+(j-1)] += DBrad[(i-1)*o+(j-1)];
      }
   }


   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DArad, &ia, &ja, descArad, DBabs, &ib, &jb, 
           descBabs, &beta, DC1, &icc, &jc, descC1);

   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAabs, &ia, &ja, descAabs, DBrad, &ib, &jb, 
           descBrad, &beta, DC2, &icc, &jc, descC2);

   delete[] DAabs;
   delete[] DArad;
   delete[] DBabs;   
   delete[] DBrad;  
   
   C = imatrix(o,p);

   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         Crad[j][i] += DC1[(i-1)*o+(j-1)] + DC2[(i-1)*o+(j-1)];
	 UncheckedSetSup(C[j][i], Cmid[j][i] + Crad[j][i]);
      }
   }
   
   setRound(-1);

   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
	 UncheckedSetInf(C[j][i], Cmid[j][i] - Crad[j][i]);
      }
   }
   
  
   setRound(0);
   
 
   delete[] DC1;          
   delete[] DC2;             
//ausg << "matmul ii ende" << endl; 
}

/*!
  The product C=A*B for A of type rmatrix B of type imatrix is computed using 
  BLAS. C is an interval matrix enclosing the correct result. However, the 
  diameter of C can get very large depending of the condition of A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(imatrix &A, rmatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
//ausg << "matmul ir" << endl;  
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;
     
   C = imatrix();

   int descB[9], descAmid[9], descArad[9], descC1[9], descC2[9];
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descB, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAmid, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descArad, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC1, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC2, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;

   //copy B into double array
   double *DB = new double[o*p];
   setRound(1);
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         DB[(j-1)*o+(i-1)] = _double(B[lbB_i+i-1][lbB_j+j-1]);
      }
   }         

   //compute mid(A) and rad(A) and copy into double arrays
   double *DAmid = new double[m*n];
   double *DArad = new double[m*n];   
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DAmid[(j-1)*m+(i-1)] = _double(Inf(A[lbA_i+i-1][lbA_j+j-1]) + 0.5*(Sup(A[lbA_i+i-1][lbA_j+j-1]) - Inf(A[lbA_i+i-1][lbA_j+j-1])));         
         DArad[(j-1)*m+(i-1)] = _double(DAmid[(j-1)*m+(i-1)] - Inf(A[lbA_i+i-1][lbA_j+j-1]));    
      }
   }      
   
   //compute lower bound for mid of result   
   double *DC1   = new double[o*p];
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAmid, &ia, &ja, descAmid, DB, &ib, &jb, 
           descB, &beta, DC1, &icc, &jc, descC1);

   //compute upper bound for mid of result
   double *DC2   = new double[o*p];   
   setRound(1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAmid, &ia, &ja, descAmid, DB, &ib, &jb, 
           descB, &beta, DC2, &icc, &jc, descC2);

   delete[] DAmid;

   rmatrix Cmid(o,p);   
   
   //compute mid of result
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         Cmid[j][i] = DC1[(i-1)*o+(j-1)] + 0.5 * (DC2[(i-1)*o+(j-1)] - DC1[(i-1)*o+(j-1)]); 
      }
   }       

   //copy abs(B) into double array
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         DB[(i-1)*o+(j-1)] = _double(abs(B[lbB_i+j-1][lbB_j+i-1])); //jetzt ist DA=abs(A)
      }
   }

   //compute upper bound for rad of result
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DArad, &ia, &ja, descArad, DB, &ib, &jb, 
           descB, &beta, DC2, &icc, &jc, descC2);

   delete[] DB;
   delete[] DArad;  
   
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         DC2[(i-1)*o+(j-1)] += _double(Cmid[j][i]) - DC1[(i-1)*o+(j-1)];
      }
   }

   setRound(-1);

   //compute infimum of result and copy into C
   C = imatrix(o,p);
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         C[lbC_i+j-1][lbC_j+i-1] = Cmid[j][i] - DC2[(i-1)*o+(j-1)]; 
      }
   }        
   
   setRound(1);

   //compute supremum of result and copy into C
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         SetSup(C[lbC_i+j-1][lbC_j+i-1], Cmid[j][i] + DC2[(i-1)*o+(j-1)]); 
      }
   }        
   
   setRound(0);
 
   delete[] DC1;          
   delete[] DC2;          
   
//   ausg << "matmul ir ende" << endl;
 
}

//void matmul(imatrix &A, imatrix &B, imatrix &C);

/*!
  The product C=A*B for A and B of type cmatrix is computed using 
  BLAS. C is a complex interval matrix enclosing the correct result. 
  However, the diameter of C can get very large depending of the condition of 
  A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(cmatrix &A, cmatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
//  ausg << "matmul cc=ci" << endl;
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;
   
   C = cimatrix();

   
   //Copy A and B into corresponding double arrays   
   double *DAR  = new double[m*n]; //A.re
   double *DAI  = new double[m*n]; //A.Im
   double *DAIm = new double[m*n]; //-A.Im
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         DAR[ind] = _double(Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] =  -DAI[ind];         
      }
   }
   A = cmatrix();

   double *DBR = new double[p*o];  //B.re
   double *DBI = new double[p*o];  //B.Im
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((j-1)*p+(i-1));
         DBR[ind] = _double(Re(B[lbB_i+i-1][lbB_j+j-1]));         
         DBI[ind] = _double(Im(B[lbB_i+i-1][lbB_j+j-1]));                  
      }
   }
   B = cmatrix();

   int descAR[9], descAI[9], descAIm[9], descBR[9], descBI[9], descCR[9], descCI[9];;
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descAR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAI, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAIm, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBR, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descCR, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descCI, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;
   
   //Compute lower bound of result     
   double *DCR = new double[o*p];  //C.re
   double *DCI = new double[o*p];  //C.Im
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DCR, &icc, &jc, descCR);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DCR, &icc, &jc, descCR);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DCI, &icc, &jc, descCI);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DCI, &icc, &jc, descCI);   

   //Copy lower bound into C
   C = cimatrix(o,p);
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         C[lbC_i+j-1][lbC_j+i-1] = complex(DCR[ind],DCI[ind]);
      }
   }        

   //Compute upper bound of result 
   setRound(1);
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DCR, &icc, &jc, descCR);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DCR, &icc, &jc, descCR);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DCI, &icc, &jc, descCI);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DCI, &icc, &jc, descCI);   

   delete[] DAIm;

   //Copy upper bound into C
   setRound(0);
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         SetSup(C[lbC_i+j-1][lbC_j+i-1], complex(DCR[ind],DCI[ind]));
      }
   }
   delete[] DCR;    
   delete[] DCI;    

   A = cmatrix(lbA_i,ubA_i,lbA_j,ubA_j);
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         SetRe(A[lbA_i+i-1][lbA_j+j-1],DAR[ind]);
         SetIm(A[lbA_i+i-1][lbA_j+j-1],DAI[ind]);
      }
   }
   delete[] DAR;
   delete[] DAI;

   B = cmatrix(lbB_i,ubB_i,lbB_j,ubB_j);
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         int ind = ((j-1)*o+(i-1));
         SetRe(B[lbB_i+i-1][lbB_j+j-1],DBR[ind]);         
         SetIm(B[lbB_i+i-1][lbB_j+j-1],DBI[ind]);                  
      }
   }
   delete[] DBR;
   delete[] DBI;   
   
//  ausg << "matmul cc=ci ende" << endl;   
}


/*!
  The product C=A*B for A and B of type cmatrix is computed using 
  BLAS. C is a complex result.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(cmatrix &A, cmatrix &B, cmatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
//  ausg << "matmul cc" << endl;  
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;
   
   C = cmatrix();

   
   //Copy A and B into corresponding double arrays   
   double *DAR  = new double[m*n]; //A.re
   double *DAI  = new double[m*n]; //A.Im
   double *DAIm = new double[m*n]; //-A.Im
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         DAR[ind] = _double(Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] =  -DAI[ind];         
      }
   }
   A = cmatrix();

   double *DBR = new double[o*p];  //B.re
   double *DBI = new double[o*p];  //B.Im
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         int ind = ((j-1)*o+(i-1));
         DBR[ind] = _double(Re(B[lbB_i+i-1][lbB_j+j-1]));         
         DBI[ind] = _double(Im(B[lbB_i+i-1][lbB_j+j-1]));                  
      }
   }
   B = cmatrix();

   int descAR[9], descAI[9], descAIm[9], descBR[9], descBI[9], descCR[9], descCI[9];;
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descAR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAI, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAIm, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBR, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descCR, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descCI, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;
   
   double *DCR = new double[o*p];  //C.re
   double *DCI = new double[o*p];  //C.Im

   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DCR, &icc, &jc, descCR);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DCR, &icc, &jc, descCR);   
   
   delete[] DAIm;
   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DCI, &icc, &jc, descCI);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DCI, &icc, &jc, descCI);   

   //Copy result into C
   C = cmatrix(o,p);
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         C[lbC_i+j-1][lbC_j+i-1] = complex(DCR[ind],DCI[ind]);
      }
   }      
   delete[] DCR;
   delete[] DCI;


   A = cmatrix(m,n);
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         SetRe(A[lbA_i+i-1][lbA_j+j-1],DAR[ind]);
         SetIm(A[lbA_i+i-1][lbA_j+j-1],DAI[ind]);
      }
   }
   delete[] DAR;
   delete[] DAI;

   B = cmatrix(o,p);
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         int ind = ((j-1)*o+(i-1));
         SetRe(B[lbB_i+i-1][lbB_j+j-1],DBR[ind]);         
         SetIm(B[lbB_i+i-1][lbB_j+j-1],DBI[ind]);                  
      }
   }
   delete[] DBR;
   delete[] DBI;   
 // ausg << "matmul cc ende" << endl;  
   
}

/*!
  The product C=A*B for A of type cmatrix and B of type cimatrix is computed 
  using BLAS. C is a complex interval matrix enclosing the correct result. 
  However, the diameter of C can get very large depending of the condition of 
  A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(cmatrix &A, cimatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
// ausg << "matmul cci" << endl;  
  
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   int o = ubB_i - lbB_i + 1;
   int p = ubB_j - lbB_j + 1;

   C = cimatrix();
   
   //copy Re(A), Im(A) and -Im(A) into double arrays   
   double *DAR = new double[m*n];   
   double *DAI = new double[m*n];   
   double *DAIm = new double[m*n]; 
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         DAR[ind] = _double(Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] = -DAI[ind];
      }
   }        

   //copy Re(mid(B)) and Im(mid(B)) into double arrays
   double *DBR = new double[o*p];   
   double *DBI = new double[o*p];  
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         int ind = ((j-1)*o+(i-1));
         complex c = Inf(B[lbB_i+i-1][lbB_j+j-1]) + 0.5*(Sup(B[lbB_i+i-1][lbB_j+j-1]) - Inf(B[lbB_i+i-1][lbB_j+j-1]));
         DBR[ind] = _double(Re(c));         
         DBI[ind] = _double(Im(c));                  
      }
   }        

  int descAR[9], descAI[9], descAIm[9], descBR[9], descBI[9], 
      descC1R[9], descC1I[9];
  int lld = numr;
  int info;
  int irsrc = 0, icsrc = 0;
  if(lld == 0) lld = 1;

  descinit_(descAR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descAI, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descAIm, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descBR, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descBI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC1R, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC1I, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
  int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
  char *TRANS=(char*)"No Transpose";
  double alpha = 1.0, beta = 0.0;

  double *DC1R = new double[o*p];  
  double *DC1I = new double[o*p];   

   //Compute lower bound of mid of result   
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DC1I, &icc, &jc, descC1I);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DC1I, &icc, &jc, descC1I);  

   //mid of result   
   cmatrix Cmid(o,p);

   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         Cmid[j][i] = complex(DC1R[ind],DC1I[ind]);
      }
   }    
                     
   //Compute upper bound of mid of result
   setRound(1);
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DC1I, &icc, &jc, descC1I);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DC1I, &icc, &jc, descC1I);  

   delete[] DAI;
   delete[] DAIm;
   delete[] DBI; 

   //rad of result   
   cmatrix Crad(o,p);
   
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         complex c = Cmid[j][i];
         //setting mid(C) to mid between lower and upper bound of mid of result         
         Cmid[j][i] += 0.5 * (complex(DC1R[ind],DC1I[ind]) - c);
         //First part of rad(C): difference between lower bound of midpoint and 
         //actually used midpoint         
         Crad[j][i] = Cmid[j][i] - c;
      }
   }       

   delete[] DC1I;    
   
   //Compute the radius of a disc with midpoint mid(B) that encloses B, which
   //is a rectangle in the complex plane. To save computing time, the radius
   //of the disc is computed as max(sqrt(2)*0.5*diam(Re(B)), sqrt(2)*0.5*diam(Im(B)))
   //which overestimates the real radius of the closest enclosing disc if
   //B is not a sqaure in the complex plane   
   for(int i=1 ; i<=o ; i++) {
      for(int j=1 ; j<=p ; j++) {
         int ind = ((j-1)*o+(i-1));
         real diam_re = diam(Re(B[lbB_i+i-1][lbB_j+j-1]));
         real diam_im = diam(Im(B[lbB_i+i-1][lbB_j+j-1]));
         real rad = 0;

         //(roundin mode still set to round up!)
         rad = sqrt((0.5*diam_re*0.5*diam_re + 0.5*diam_im*0.5*diam_im));
                                   
         //DBR nor rad(B)
         DBR[ind] = _double(rad);         
      }
   }    
   
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));	
         //DAR now abs(A)
         DAR[ind] = sqrt(_double(Re(A[lbA_i+i-1][lbA_j+j-1])*Re(A[lbA_i+i-1][lbA_j+j-1]) +
                         Im(A[lbA_i+i-1][lbA_j+j-1])*Im(A[lbA_i+i-1][lbA_j+j-1])));             
      }
   }    
   

   //Compute an upper bound for the radius of the result
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DC1R, &icc, &jc, descC1R);   

   delete[] DAR;
   delete[] DBR;
 
   C = cimatrix(o,p);  
   for(int i=1 ; i<=p ; i++) {
      for(int j=1 ; j<=o ; j++) {
         int ind = ((i-1)*o+(j-1));
         //With this radius C reliably encloses the correct result         
         Crad[j][i] += DC1R[ind];
         //Convert midpoint and radius into infimum-supremum-representation         
         C[lbC_i+j-1][lbC_j+i-1] = cinterval(Cmid[j][i] - Crad[j][i], 
                                             Cmid[j][i] + Crad[j][i]);
      }
   } 
   delete[] DC1R;    

   setRound(0);   
//  ausg << "matmul cci ende" << endl;  

}


/*!
  The product C=A*B for A of type cmatrix and B of type cimatrix is computed 
  using BLAS. C is a complex interval matrix enclosing the correct result. 
  However, the diameter of C can get very large depending of the condition of 
  A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param ausg Output stream
  
*/

void matmul(cimatrix &A, cmatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) {
//  ausg << "matmul cic" << endl;  

  imatrix tmp(ColLen(A),RowLen(A)), tmp1;
  rmatrix tmp2;
  tmp1 = Re(A); tmp2 = Re(B);
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  C = tmp;
  tmp1 = Im(A); tmp2 = Im(B);
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  C -= tmp;
  tmp1 = Re(A); tmp2 = Im(B);  
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  SetIm(C,tmp);
  tmp1 = Im(A); tmp2 = Re(B);  
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  SetIm(C, Im(C)+tmp);
  
//    int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
//    int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
//    int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
//    int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
//    int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
//    int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
//    
//    int m = ubA_i - lbA_i + 1;
//    int n = ubA_j - lbA_j + 1;
//    int o = ubB_i - lbB_i + 1;
//    int p = ubB_j - lbB_j + 1;
// 
//    C = cimatrix();
//    
//    //copy Re(A), Im(A) and -Im(A) into double arrays   
//    double *DBR = new double[o*p];   
//    double *DBI = new double[o*p];   
//    double *DBIm = new double[o*p]; 
//    for(int i=1 ; i<=o ; i++) {
//       for(int j=1 ; j<=p ; j++) {
//          int ind = ((j-1)*o+(i-1));
//          DBR[ind] = _double(Re(B[lbB_i+i-1][lbB_j+j-1]));
//          DBI[ind] = _double(Im(B[lbB_i+i-1][lbB_j+j-1]));
//          DBIm[ind] = -DBI[ind];
//       }
//    }        
// 
//    //copy Re(mid(B)) and Im(mid(B)) into double arrays
//    double *DAR = new double[m*n];   
//    double *DAI = new double[m*n];  
//    for(int i=1 ; i<=m ; i++) {
//       for(int j=1 ; j<=n ; j++) {
//          int ind = ((j-1)*m+(i-1));
//          complex c = Inf(A[lbA_i+i-1][lbA_j+j-1]) + 0.5*(Sup(A[lbA_i+i-1][lbA_j+j-1]) - Inf(A[lbA_i+i-1][lbA_j+j-1]));
//          DAR[ind] = _double(Re(c));         
//          DAI[ind] = _double(Im(c));                  
//       }
//    }        
// 
//   int descBR[9], descBI[9], descBIm[9], descAR[9], descAI[9], 
//       descC1R[9], descC1I[9];
//   int lld = numr;
//   int info;
//   int irsrc = 0, icsrc = 0;
//   if(lld == 0) lld = 1;
// 
//   descinit_(descBR, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   descinit_(descBI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   descinit_(descBIm, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   descinit_(descAR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   descinit_(descAI, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   descinit_(descC1R, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   descinit_(descC1I, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
//   
//   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
//   char *TRANS=(char*)"No Transpose";
//   double alpha = 1.0, beta = 0.0;
// 
//   double *DC1R = new double[o*p];  
//   double *DC1I = new double[o*p];   
// 
//    //Compute lower bound of mid of result   
//    setRound(-1);
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
//            descBR, &beta, DC1R, &icc, &jc, descC1R);   
//    beta = 1.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBIm, &ib, &jb, 
//            descBIm, &beta, DC1R, &icc, &jc, descC1R);   
//    beta = 0.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
//            descBI, &beta, DC1I, &icc, &jc, descC1I);   
//    beta = 1.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
//            descBR, &beta, DC1I, &icc, &jc, descC1I);  
// 
//    //mid of result   
//    cmatrix Cmid(o,p);
// 
//    for(int i=1 ; i<=p ; i++) {
//       for(int j=1 ; j<=o ; j++) {
//          int ind = ((i-1)*o+(j-1));
//          Cmid[j][i] = complex(DC1R[ind],DC1I[ind]);
//       }
//    }    
//                      
//    //Compute upper bound of mid of result
//    setRound(1);
//    beta = 0.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
//            descBR, &beta, DC1R, &icc, &jc, descC1R);   
//    beta = 1.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBIm, &ib, &jb, 
//            descBIm, &beta, DC1R, &icc, &jc, descC1R);   
//    beta = 0.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
//            descBI, &beta, DC1I, &icc, &jc, descC1I);   
//    beta = 1.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
//            descBR, &beta, DC1I, &icc, &jc, descC1I);  
// 
//    delete[] DAI;
//    delete[] DBIm;
//    delete[] DBI; 
// 
//    //rad of result   
//    cmatrix Crad(o,p);
//    
//    for(int i=1 ; i<=p ; i++) {
//       for(int j=1 ; j<=o ; j++) {
//          int ind = ((i-1)*o+(j-1));
//          complex c = Cmid[j][i];
//          //setting mid(C) to mid between lower and upper bound of mid of result         
//          Cmid[j][i] += 0.5 * (complex(DC1R[ind],DC1I[ind]) - c);
//          //First part of rad(C): difference between lower bound of midpoint and 
//          //actually used midpoint         
//          Crad[j][i] = Cmid[j][i] - c;
//       }
//    }       
// 
//    delete[] DC1I;    
// 
//    //Compute the radius of a disc with midpoint mid(B) that encloses B, which
//    //is a rectangle in the complex plane. To save computing time, the radius
//    //of the disc is computed as max(sqrt(2)*0.5*diam(Re(B)), sqrt(2)*0.5*diam(Im(B)))
//    //which overestimates the real radius of the closest enclosing disc if
//    //B is not a sqaure in the complex plane   
//    for(int i=1 ; i<=m ; i++) {
//       for(int j=1 ; j<=n ; j++) {
//          int ind = ((j-1)*m+(i-1));
//          real diam_re = diam(Re(A[lbA_i+i-1][lbA_j+j-1]));
//          real diam_im = diam(Im(A[lbA_i+i-1][lbA_j+j-1]));
//          real rad = 0;
// 
//          //(roundin mode still set to round up!)
//          rad = sqrt((0.5*diam_re*0.5*diam_re + 0.5*diam_im*0.5*diam_im)); 
//                                    
//          //DAR now rad(A)
//          DAR[ind] = _double(rad);         
//       }
//    }    
//    
//    for(int i=1 ; i<=o ; i++) {
//       for(int j=1 ; j<=p ; j++) {
//          int ind = ((j-1)*o+(i-1));	
//          //DBR now abs(B)
//          DBR[ind] = sqrt(_double(Re(B[lbB_i+i-1][lbB_j+j-1])*Re(B[lbB_i+i-1][lbB_j+j-1]) +
//                          Im(B[lbB_i+i-1][lbB_j+j-1])*Im(B[lbB_i+i-1][lbB_j+j-1])));             
//       }
//    }    
//    
//    //Compute an upper bound for the radius of the result
//    beta = 0.0;
//    pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
//            descBR, &beta, DC1R, &icc, &jc, descC1R);   
// 
//    delete[] DAR;
//    delete[] DBR;
//  
//    C = cimatrix(o,p);  
// 
//    for(int i=1 ; i<=o ; i++) {
//       for(int j=1 ; j<=p ; j++) {
//          int ind = ((i-1)*p+(j-1));
//          //With this radius C reliably encloses the correct result         
//          Crad[j][i] += DC1R[ind];
//          //Convert midpoint and radius into infimum-supremum-representation         
//          UncheckedSetInf(C[lbC_i+j-1][lbC_j+i-1], Cmid[j][i] - Crad[j][i]); 
//          UncheckedSetSup(C[lbC_i+j-1][lbC_j+i-1], Cmid[j][i] + Crad[j][i]);
//       }
//    } 
//    delete[] DC1R;    
// 
//    setRound(0);   

 //  ausg << "matmul cic ende" << endl;  

}


void matmul(cimatrix &A, cimatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, ofstream& ausg) { 
 // ausg << "matmul cici" << endl;  
  
  imatrix tmp(ColLen(A),RowLen(A)), tmp1, tmp2;
  tmp1 = Re(A); tmp2 = Re(B);
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  C = tmp;
  tmp1 = Im(A); tmp2 = Im(B);
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  C -= tmp;
  tmp1 = Re(A); tmp2 = Im(B);  
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  SetIm(C,tmp);
  tmp1 = Im(A); tmp2 = Re(B);  
  matmul(tmp1,tmp2,tmp,r,s,t,nb,ic,numr,ausg);
  SetIm(C, Im(C)+tmp);
  
 // ausg << "matmul cici ende" << endl;    
}


//void matmul(cimatrix &A, cimatrix &B, cimatrix &C);

/*!
  (I-A*B) for A and B of type rmatrix is computed using BLAS. C is
  an interval matrix enclosing the correct result. However, the diameter of
  C can get very large depending of the condition of A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param myr Own row coordinate in process grid
  \param myc Own column coordinate in process grid
  \param nr Number of rows of process grid
  \param nc Number of columns of process grid
  \param ausg Output stream  
*/

void IminusAB(rmatrix &A, rmatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg) {
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     

   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   
   //Delete C
   C = imatrix();

   double *DA = new double[m*n];
   //Copy A and B into double array    
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*m+(i-1)] = _double(-A[lbA_i+i-1][lbA_j+j-1]);
      }
   }
   //delete original A
   A = rmatrix();

   double *DB = new double[m*n];
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DB[(j-1)*m+(i-1)] = _double(B[lbB_i+i-1][lbB_j+j-1]);         
      }
   }
   //delete original B
   B = rmatrix();

  int descA[9], descB[9], descC[9];
  int lld = numr;
  int info;
  int irsrc = 0, icsrc = 0;
  if(lld == 0) lld = 1;

  descinit_(descA, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descB, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
  int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
  char *TRANS=(char*)"No Transpose";
  double alpha = 1.0, beta = 0.0;

   //Compute Infimum
  double *DC = new double[m*n];
  setRound(-1);
  pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DB, &ib, &jb, 
          descB, &beta, DC, &icc, &jc, descC);

  
  //Copy infimum into C
  //C = imatrix(lbC_i,ubC_i,lbC_j,ubC_j);
  if(m>0 && n>0)
    Resize(C,lbC_i,ubC_i,lbC_j,ubC_j);
  int nbr = (nr==1) ? r : nb;
  int nbc = (nc==1) ? r : nb;
  for(int i=1 ; i<=n ; i++) {
     for(int j=1 ; j<=m ; j++) {
        if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
          C[lbC_i+j-1][lbC_j+i-1] = DC[(i-1)*m+(j-1)] + 1.0;
        else
          C[lbC_i+j-1][lbC_j+i-1] = DC[(i-1)*m+(j-1)];
     }
  }        

  //Compute supremum
  setRound(1);
  pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DB, &ia, &ja, 
          descB, &beta, DC, &ia, &ja, descC);

  //copy supremum into C
  for(int i=1 ; i<=n ; i++) {
     for(int j=1 ; j<=m ; j++) {
        if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
          SetSup(C[lbC_i+j-1][lbC_j+i-1], DC[(i-1)*m+(j-1)] + 1.0);
        else
          SetSup(C[lbC_i+j-1][lbC_j+i-1], DC[(i-1)*m+(j-1)]);
     }
  }        
  setRound(0);

  delete[] DC;     

  //Copy back A
  //A = rmatrix(lbA_i,ubA_i,lbA_j,ubA_j);
  if(m>0 && n>0)
    Resize(A,lbA_i,ubA_i,lbA_j,ubA_j);
  for(int i=1 ; i<=m ; i++) {
     for(int j=1 ; j<=n ; j++) {
        A[lbA_i+i-1][lbA_j+j-1] = -DA[(j-1)*m+(i-1)];
     }
  }
  delete[] DA;

  //B = rmatrix(lbB_i,ubB_i,lbB_j,ubB_j);
  if(m>0 && n>0)
    Resize(B,lbB_i,ubB_i,lbB_j,ubB_j);
  for(int i=1 ; i<=m ; i++) {
     for(int j=1 ; j<=n ; j++) {
        B[lbB_i+i-1][lbB_j+j-1] = DB[(j-1)*m+(i-1)];         
     }
  }
  delete[] DB;
}


/*!
  (I-A*B) for A of type rmatrix and B of type imatrix is computed using 
  BLAS. C is an interval matrix enclosing the correct result. However, the 
  diameter of C can get very large depending of the condition of A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param myr Own row coordinate in process grid
  \param myc Own column coordinate in process grid
  \param nr Number of rows of process grid
  \param nc Number of columns of process grid
  \param ausg Output stream  
  
*/

void IminusAB(rmatrix &A, imatrix &B, imatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg) {
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;

   C = imatrix(); 

   int descA[9], descBmid[9], descBrad[9], descC1[9], descC2[9];
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descA, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBmid, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBrad, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC1, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descC2, &r, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;
   
   //copy A into double array
   setRound(1);
   double *DA    = new double[m*n];
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DA[(j-1)*m+(i-1)] = _double(-A[lbA_i+i-1][lbA_j+j-1]);
      }
   }         

   //compute mid(B) and rad(B) and copy into double arrays
   double *DBmid = new double[m*n];
   double *DBrad = new double[m*n];   
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         DBmid[(j-1)*m+(i-1)] = _double(Inf(B[lbB_i+i-1][lbB_j+j-1]) + 0.5*(Sup(B[lbB_i+i-1][lbB_j+j-1]) - Inf(B[lbB_i+i-1][lbB_j+j-1])));         
         DBrad[(j-1)*m+(i-1)] = _double(DBmid[(j-1)*m+(i-1)] - Inf(B[lbB_i+i-1][lbB_j+j-1]));    
      }
   }          

   //compute lower bound for mid of result   
   double *DC1   = new double[m*n];
   double *DC2   = new double[m*n];   
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DBmid, &ib, &jb, 
           descBmid, &beta, DC1, &icc, &jc, descC1);

   //compute upper bound for mid of result
   setRound(1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DBmid, &ib, &jb, 
           descBmid, &beta, DC2, &icc, &jc, descC2);

   delete[] DBmid;
   
   rmatrix Cmid(m,n); 
   
   //compute mid of result
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         Cmid[j][i] = DC1[(i-1)*m+(j-1)] + 0.5 * (DC2[(i-1)*m+(j-1)] - DC1[(i-1)*m+(j-1)]); 
      }
   }       

   //copy abs(A) into double array
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         DA[(i-1)*m+(j-1)] = _double(abs(A[lbA_i+j-1][lbA_j+i-1])); //jetzt ist DA=abs(A)
      }
   }


   //compute upper bound for rad of result
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DA, &ia, &ja, descA, DBrad, &ib, &jb, 
           descBrad, &beta, DC2, &icc, &jc, descC2);

   delete[] DA;
   delete[] DBrad;  
   
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         DC2[(i-1)*m+(j-1)] += _double(Cmid[j][i]) - DC1[(i-1)*m+(j-1)];
      }
   }
   delete[] DC1;
   
   setRound(-1);

   C = imatrix(lbC_i,ubC_i,lbC_j,ubC_j);

   //compute infimum of result and copy into C
  int nbr = (nr==1) ? r : nb;
  int nbc = (nc==1) ? r : nb;   
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
        if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
           C[lbC_i+j-1][lbC_j+i-1] = Cmid[j][i] - DC2[(i-1)*m+(j-1)] + 1.0; 
         else              
           C[lbC_i+j-1][lbC_j+i-1] = Cmid[j][i] - DC2[(i-1)*m+(j-1)]; 
      }
   }        
   
   setRound(1);

   //compute supremum of result and copy into C
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
        if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
           SetSup(C[lbC_i+j-1][lbC_j+i-1], Cmid[j][i] + DC2[(i-1)*m+(j-1)] + 1.0); 
         else
           SetSup(C[lbC_i+j-1][lbC_j+i-1], Cmid[j][i] + DC2[(i-1)*m+(j-1)]); 
      }
   }        
   
   setRound(0);
    
   delete[] DC2;             
 
}


/*!
  (I-A*B) for A and B of type cmatrix is computed using 
  BLAS. C is a complex interval matrix enclosing the correct result. 
  However, the diameter of C can get very large depending of the condition of 
  A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param myr Own row coordinate in process grid
  \param myc Own column coordinate in process grid
  \param nr Number of rows of process grid
  \param nc Number of columns of process grid
  \param ausg Output stream  
  
*/

void IminusAB(cmatrix &A, cmatrix &B, cimatrix &C, int r, int s, int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, ofstream& ausg) {
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;
   
   C = cimatrix();
   
   //Copy A and B into corresponding double arrays   
   double *DAR  = new double[m*n]; //A.re
   double *DAI  = new double[m*n]; //A.Im
   double *DAIm = new double[m*n]; //-A.Im
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         DAR[ind] = _double(-Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(-Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] =  -DAI[ind];         
      }
   }
   A = cmatrix();

   double *DBR = new double[m*n];  //B.re
   double *DBI = new double[m*n];  //B.Im
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         DBR[ind] = _double(Re(B[lbB_i+i-1][lbB_j+j-1]));         
         DBI[ind] = _double(Im(B[lbB_i+i-1][lbB_j+j-1]));                  
      }
   }
   B = cmatrix();

   int descAR[9], descAI[9], descAIm[9], descBR[9], descBI[9], descCR[9], descCI[9];;
   int lld = numr;
   int info;
   int irsrc = 0, icsrc = 0;
   if(lld == 0) lld = 1;

   descinit_(descAR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descAIm, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descBI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descCR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
   descinit_(descCI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
   int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
   char *TRANS=(char*)"No Transpose";
   double alpha = 1.0, beta = 0.0;

   double *DCR = new double[m*n];  //C.re
   double *DCI = new double[m*n];  //C.Im

   //Compute lower bound of result     
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DCR, &icc, &jc, descCR);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DCR, &icc, &jc, descCR);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DCI, &icc, &jc, descCI);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DCI, &icc, &jc, descCI);   

   //Copy lower bound into C  
  C = cimatrix(lbC_i,ubC_i,lbC_j,ubC_j);
  int nbr = (nr==1) ? r : nb;
  int nbc = (nc==1) ? r : nb;        
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
        if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
           C[lbC_i+j-1][lbC_j+i-1] = complex(DCR[ind] + 1.0, DCI[ind]);
         else
           C[lbC_i+j-1][lbC_j+i-1] = complex(DCR[ind],DCI[ind]);
      }
   }        

   //Compute upper bound of result 
   setRound(1);
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DCR, &icc, &jc, descCR);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DCR, &icc, &jc, descCR);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DCI, &icc, &jc, descCI);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DCI, &icc, &jc, descCI);   

   delete[] DAIm;

   //Copy upper bound into C
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
        if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
           SetSup(C[lbC_i+j-1][lbC_j+i-1], complex(DCR[ind] + 1.0, DCI[ind]));         
         else
           SetSup(C[lbC_i+j-1][lbC_j+i-1], complex(DCR[ind],DCI[ind]));
      }
   }       

   delete[] DCR;    
   delete[] DCI;    

   A = cmatrix(lbA_i,ubA_i,lbA_j,ubA_j);
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         SetRe(A[lbA_i+i-1][lbA_j+j-1],-DAR[ind]);
         SetIm(A[lbA_i+i-1][lbA_j+j-1],-DAI[ind]);
      }
   }
   delete[] DAR;
   delete[] DAI;

   B = cmatrix(lbC_i,ubC_i,lbC_j,ubC_j);
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         SetRe(B[lbB_i+i-1][lbB_j+j-1],DBR[ind]);         
         SetIm(B[lbB_i+i-1][lbB_j+j-1],DBI[ind]);                  
      }
   }
   delete[] DBR;
   delete[] DBI;   


   setRound(0);
      
}


/*!
  (I-A*B) for A of type cmatrix and B of type cimatrix is computed 
  using BLAS. C is a complex interval matrix enclosing the correct result. 
  However, the diameter of C can get very large depending of the condition of 
  A and B.
  
  \param A first operator
  \param B second operator
  \param C An enclosure of A*B
  \param r Number of rows of A
  \param s Number of columns of A / rows of B
  \param t Number of columns of B
  \param nb Blocksize for ScaLAPACK
  \param ic Blacs context
  \param numr Number ofr rows of own matrix
  \param myr Own row coordinate in process grid
  \param myc Own column coordinate in process grid
  \param nr Number of rows of process grid
  \param nc Number of columns of process grid
  \param ausg Output stream  
  
*/

void IminusAB(cmatrix &A, cimatrix &B, cimatrix &C, int r, int s, 
              int t, int nb, int ic, int numr, int myr, int myc, int nr, int nc, 
              ofstream& ausg) {
   int lbA_i = Lb(A,1); int lbA_j = Lb(A,2);
   int lbB_i = Lb(B,1); int lbB_j = Lb(B,2); 
   int ubA_i = Ub(A,1); int ubA_j = Ub(A,2);
   int ubB_i = Ub(B,1); int ubB_j = Ub(B,2);   
   int lbC_i = Lb(C,1); int lbC_j = Lb(C,2); 
   int ubC_i = Ub(C,1); int ubC_j = Ub(C,2);     
   
   int m = ubA_i - lbA_i + 1;
   int n = ubA_j - lbA_j + 1;

   C = cimatrix();
   
   //copy Re(A), Im(A) and -Im(A) into double arrays   
   double *DAR = new double[m*n];   
   double *DAI = new double[m*n];   
   double *DAIm = new double[m*n]; 
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         DAR[ind] = _double(-Re(A[lbA_i+i-1][lbA_j+j-1]));
         DAI[ind] = _double(-Im(A[lbA_i+i-1][lbA_j+j-1]));
         DAIm[ind] = -DAI[ind];
      }
   }        

   //copy Re(mid(B)) and Im(mid(B)) into double arrays
   double *DBR = new double[m*n];   
   double *DBI = new double[m*n];  
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         complex c = Inf(B[lbB_i+i-1][lbB_j+j-1]) + 0.5*(Sup(B[lbB_i+i-1][lbB_j+j-1]) - Inf(B[lbB_i+i-1][lbB_j+j-1]));
         DBR[ind] = _double(Re(c));         
         DBI[ind] = _double(Im(c));                   
      }
   }        

  int descAR[9], descAI[9], descAIm[9], descBR[9], descBI[9], 
      descC1R[9], descC1I[9];
  int lld = numr;
  int info;
  int irsrc = 0, icsrc = 0;
  if(lld == 0) lld = 1;

  descinit_(descAR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descAI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descAIm, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descBR, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descBI, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC1R, &r, &s, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  descinit_(descC1I, &s, &t, &nb, &nb, &irsrc, &icsrc, &ic, &lld, &info); 
  
  int ia=1, ja=1, ib=1, jb=1, icc=1, jc=1;
  char *TRANS=(char*)"No Transpose";
  double alpha = 1.0, beta = 0.0;

  double *DC1R = new double[m*n];  
  double *DC1I = new double[m*n];   

   //Compute lower bound of mid of result   
   setRound(-1);
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DC1I, &icc, &jc, descC1I);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DC1I, &icc, &jc, descC1I);  

   //mid of result   
   cmatrix Cmid(m,n);

   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         Cmid[j][i] = complex(DC1R[ind],DC1I[ind]);
      }
   }    
                     
   //Compute upper bound of mid of result
   setRound(1);
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAIm, &ia, &ja, descAIm, DBI, &ib, &jb, 
           descBI, &beta, DC1R, &icc, &jc, descC1R);   
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBI, &ib, &jb, 
           descBI, &beta, DC1I, &icc, &jc, descC1I);   
   beta = 1.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAI, &ia, &ja, descAI, DBR, &ib, &jb, 
           descBR, &beta, DC1I, &icc, &jc, descC1I);  

   delete[] DAI;
   delete[] DAIm;
   delete[] DBI; 

   //rad of result   
   cmatrix Crad(m,n);

  int nbr = (nr==1) ? r : nb;
  int nbc = (nc==1) ? r : nb;      
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         complex c = Cmid[j][i];
         //setting mid(C) to mid between lower and upper bound of mid of result         
         Cmid[j][i] += 0.5 * (complex(DC1R[ind],DC1I[ind]) - c);
         //First part of rad(C): difference between lower bound of midpoint and 
         //actually used midpoint         
         Crad[j][i] = abs(Cmid[j][i] - c);
         if(((j-1)/nbr)*nr*nbr+myr*nbr+((j-1)%nbr)+1  ==  ((i-1)/nbc)*nc*nbc+myc*nbc+((i-1)%nbc)+1) 
           Cmid[j][i] += 1.0;         
      }
   }       

   delete[] DC1I;    
   
   //Compute the radius of a disc with midpoint mid(B) that encloses B, which
   //is a rectangle in the complex plane. To save computing time, the radius
   //of the disc is computed as max(sqrt(2)*0.5*diam(Re(B)), sqrt(2)*0.5*diam(Im(B)))
   //which overestimates the real radius of the closest enclosing disc if
   //B is not a sqaure in the complex plane   
   for(int i=1 ; i<=m ; i++) {
      for(int j=1 ; j<=n ; j++) {
         int ind = ((j-1)*m+(i-1));
         real diam_re = diam(Re(B[lbB_i+i-1][lbB_j+j-1]));
         real diam_im = diam(Im(B[lbB_i+i-1][lbB_j+j-1]));
         real rad = 0;

         //(roundin mode still set to round up!)
         if(diam_re > diam_im)
           rad = cr * diam_re;
         else  
           rad = cr * diam_im;  
                             
         //DAR now abs(A)
         DAR[ind] = sqrt(_double(Re(A[lbA_i+i-1][lbA_j+j-1])*Re(A[lbA_i+i-1][lbA_j+j-1]) +
                         Im(A[lbA_i+i-1][lbA_j+j-1])*Im(A[lbA_i+i-1][lbA_j+j-1])));             
         //DBR nor rad(B)
         DBR[ind] = _double(rad);         
      }
   }    

   //Compute an upper bound for the radius of the result
   beta = 0.0;
   pdgemm_(TRANS, TRANS, &r, &t, &s, &alpha, DAR, &ia, &ja, descAR, DBR, &ib, &jb, 
           descBR, &beta, DC1R, &icc, &jc, descC1R);   

   delete[] DAR;
   delete[] DBR;

   C = cimatrix(lbC_i,ubC_i,lbC_j,ubC_j); 
   
   for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=m ; j++) {
         int ind = ((i-1)*m+(j-1));
         //With this radius C reliably encloses the correct result         
         Crad[j][i] += DC1R[ind];
         //Convert midpoint and radius into infimum-supremum-representation         
         C[lbC_i+j-1][lbC_j+i-1] = cinterval(Cmid[j][i] - Crad[j][i], 
                                             Cmid[j][i] + Crad[j][i]);                                             
      }
   } 

   setRound(0);   

   delete[] DC1R;         
}

