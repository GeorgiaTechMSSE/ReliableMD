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

/* Implementation of interface functions to use UMFPACK with C-XSC data types */

#include "cxsc_umfpack.hpp"
#include <srmatrix.hpp>
#include <intvector.hpp>
#include <umfpack.h>
#include <iostream>

namespace cxsc { 

/*  Compute LU-decomposition of A and fill-in reducing permutation p and q such that
 *  L*U=A(p,q). The error code err is 0 if no error occured and !=0 otherwise. */
void lu_decomp(const srmatrix& A, srmatrix& L, srmatrix& U, intvector& p, intvector& q, int& Err) {
  int n = RowLen(A);    

  void *Symbolic, *Numeric;

  const std::vector<int>&  pv   = A.column_pointers();
  const std::vector<int>&  indv = A.row_indices();
  const std::vector<real>& xv   = A.values();

  int*    Ap = const_cast<int*>(&pv[0]);
  int*    Ai = const_cast<int*>(&indv[0]);
  double* Ax = const_cast<double*>(reinterpret_cast<const double*>(&xv[0]));

  double info[UMFPACK_INFO];
  double control[UMFPACK_CONTROL];

  //Set default control parameters
  umfpack_di_defaults(control);
  //Deactivate scaling
  control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
  
  //Symbolic analysis
  umfpack_di_symbolic(n,n,Ap,Ai,Ax,&Symbolic,control,info);
  
  if(info[UMFPACK_STATUS] != UMFPACK_OK) {
    Err = 1;
    return;
  }

  //Numeric decomposition
  umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,control,info);

  if(info[UMFPACK_STATUS] != UMFPACK_OK) {
    Err = 2;
    return;
  }

  //Delete data of symbolic analysis
  umfpack_di_free_symbolic(&Symbolic);

  //Get nnz(L) and nnz(U)
  int lnz = static_cast<int>(info[UMFPACK_LNZ]);
  int unz = static_cast<int>(info[UMFPACK_UNZ]);

  //Read out data for L,U,p,q,R...
  p = intvector(n);
  q = intvector(n);

  int* Lp    = new int[n+1];
  int* Li    = new int[lnz];
  double* Lx = new double[lnz];
  int* Up    = new int[n+1];
  int* Ui    = new int[unz];
  double* Ux = new double[unz];
  int* pi    = static_cast<int*>(&p[1]);
  int* qi    = static_cast<int*>(&q[1]);

  umfpack_di_get_numeric(Lp,Li,Lx,Up,Ui,Ux,pi,qi,NULL,NULL,NULL,Numeric); 

  //Copy L
  L = srmatrix(n,n,lnz,Lp,Li,reinterpret_cast<real*>(Lx),compressed_row);

  delete[] Lp;
  delete[] Li;
  delete[] Lx;

  //Copy U
  U = srmatrix(n,n,unz,Up,Ui,reinterpret_cast<real*>(Ux),compressed_column);
  
  delete[] Up;
  delete[] Ui;
  delete[] Ux;

  //Delete numeric data
  umfpack_di_free_numeric(&Numeric);
  
  Err = 0;
}



/*  Compute complex LU-decomposition of A and fill-in reducing permutation p and q such that
 *  L*U=A(p,q). The error code err is 0 if no error occured and !=0 otherwise. */
void lu_decomp(const scmatrix& A, scmatrix& L, scmatrix& U, intvector& p, intvector& q, int& Err) {
  int n = RowLen(A);  

  void *Symbolic, *Numeric;

  const std::vector<int>&  pv   = A.column_pointers();
  const std::vector<int>&  indv = A.row_indices();
  const std::vector<complex>& xv   = A.values();

  int*    Ap = const_cast<int*>(&pv[0]);
  int*    Ai = const_cast<int*>(&indv[0]);
  double* Ax = const_cast<double*>(reinterpret_cast<const double*>(&xv[0]));

  double info[UMFPACK_INFO];
  double control[UMFPACK_CONTROL];

  //Set default control parameters
  umfpack_zi_defaults(control);
  //Deactivate scaling
  control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE;
  
  //Symbolic analysis
  umfpack_zi_symbolic(n,n,Ap,Ai,Ax,NULL,&Symbolic,control,info);
  
  if(info[UMFPACK_STATUS] != UMFPACK_OK) {
    Err = 1;
    return;
  }

  //Numeric decomposition
  umfpack_zi_numeric(Ap,Ai,Ax,NULL,Symbolic,&Numeric,control,info);

  if(info[UMFPACK_STATUS] != UMFPACK_OK) {      
    Err = 2;
    return;
  }

  //Delete data of symbolic analysis
  umfpack_zi_free_symbolic(&Symbolic);

  //Get nnz(L) and nnz(U)
  int lnz = static_cast<int>(info[UMFPACK_LNZ]);
  int unz = static_cast<int>(info[UMFPACK_UNZ]);

  //Read out data for L,U,p,q,R...
  p = intvector(n);
  q = intvector(n);

  int* Lp    = new int[n+1];
  int* Li    = new int[lnz];
  double* Lx = new double[2*lnz];
  int* Up    = new int[n+1];
  int* Ui    = new int[unz];
  double* Ux = new double[2*unz];
  int* pi    = static_cast<int*>(&p[1]);
  int* qi    = static_cast<int*>(&q[1]);
  int do_recip(0);

  umfpack_zi_get_numeric(Lp,Li,Lx,NULL,Up,Ui,Ux,NULL,pi,qi,NULL,NULL,&do_recip,NULL,Numeric); 

  //Copy L
  L = scmatrix(n,n,lnz,Lp,Li,reinterpret_cast<complex*>(Lx),compressed_row);

  delete[] Lp;
  delete[] Li;
  delete[] Lx;

  //Copy U
  U = scmatrix(n,n,unz,Up,Ui,reinterpret_cast<complex*>(Ux),compressed_column);
  
  delete[] Up;
  delete[] Ui;
  delete[] Ux;

  //Delete numeric data
  umfpack_zi_free_numeric(&Numeric);
  
  Err = 0; 
  
}

} //namespace cxsc
