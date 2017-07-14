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

/* Implementation of interface functions to use CHOLMOD with C-XSC data types */

#include "cxsc_cholmod.hpp"
#include <scmatrix.hpp>
#include <intvector.hpp>
#include <cholmod.h>
#include <cholmod_core.h>

namespace cxsc {

/*  Compute Cholesky decomposition of A and fill-in reducing permutation p such that
 *  L*L'=A(p,p). The error code err is 0 if no error occured and !=0 otherwise. */
void chol(srmatrix& A, srmatrix& L, intvector& p, int& err) {
  ::cholmod_common cm;

  std::vector<int>&  Ap  = A.column_pointers();
  std::vector<int>&  Ai  = A.row_indices();
  std::vector<real>& Ax  = A.values();

  cholmod_start(&cm);
  cm.print = 0;
  cm.prefer_upper = false;

  ::cholmod_sparse *A_chol = cholmod_allocate_sparse(ColLen(A), RowLen(A), RowLen(A)+0.5*(A.get_nnz()-RowLen(A)),
 						     true, true, -1, CHOLMOD_REAL, &cm);


  int nnz = 0;
  for(int j=0 ; j<RowLen(A) ; j++) {
    ((int*)(A_chol->p))[j] = nnz;
    for(int i=Ap[j] ; i<Ap[j+1] ; i++) {  
      if(Ai[i] >= j) {
        ((int*)(A_chol->i))[nnz] = Ai[i];
        ((double*)(A_chol->x))[nnz] = _double(Ax[i]);
	nnz++;
      }
    }
  }
  ((int*)(A_chol->p))[RowLen(A)] = nnz;
 
  ::cholmod_factor* factor = (::cholmod_factor*) cholmod_analyze(A_chol,&cm);

  if(cm.status != 0) {    
    err = 1;
    cholmod_free_sparse(&A_chol,&cm);
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);
    return;
  }

  if(!cholmod_factorize(A_chol,factor,&cm)) {
    err = 2;
    cholmod_free_sparse(&A_chol,&cm);  
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);
    return;
  }

  if(cm.status == CHOLMOD_NOT_POSDEF) {
    err = 3;
    cholmod_free_sparse(&A_chol,&cm);
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);    
    return;
  }

  if(!(factor->is_ll)) {   
    cholmod_change_factor(CHOLMOD_REAL, 1, 0, 1, 1, factor, &cm);
  }

  if(cm.status == CHOLMOD_NOT_POSDEF) {
    err = 4;    
    cholmod_free_sparse(&A_chol,&cm);
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);    
    return;
  }


  Resize(p, A_chol->nrow);
  int *perm_int = (int*)(factor->Perm);
  for(int i=1 ; i<=A_chol->nrow ; i++)
    p[i] = perm_int[i-1];

  
  cholmod_free_sparse(&A_chol,&cm);
  
  ::cholmod_sparse* L_chol = (::cholmod_sparse*) cholmod_factor_to_sparse(factor,&cm);
  

  L = srmatrix(L_chol->nrow, L_chol->ncol, ((int*)(L_chol->p))[L_chol->ncol], (int*)L_chol->p, 
               (int*)L_chol->i, reinterpret_cast<real*>(L_chol->x), compressed_column);

  cholmod_free_sparse(&L_chol,&cm);
  cholmod_free_factor(&factor,&cm);
  cholmod_finish(&cm);

  err = 0;
}



/*  Compute complex Cholesky decomposition of A and fill-in reducing permutation p such that
 *  L*L'=A(p,p). The error code err is 0 if no error occured and !=0 otherwise. */
void chol(scmatrix& A, scmatrix& L, intvector& p, int& err) {
  ::cholmod_common cm;

  std::vector<int>&  Ap  = A.column_pointers();
  std::vector<int>&  Ai  = A.row_indices();
  std::vector<complex>& Ax  = A.values();

  cholmod_start(&cm);
  cm.print = 0;

  ::cholmod_sparse *A_chol = cholmod_allocate_sparse(ColLen(A), RowLen(A), RowLen(A)+0.5*(A.get_nnz()-RowLen(A)),
						     true, true, -1, CHOLMOD_COMPLEX, &cm);

  int nnz = 0;
  for(int j=0 ; j<RowLen(A) ; j++) {
    ((int*)(A_chol->p))[j] = nnz;
    for(int i=Ap[j] ; i<Ap[j+1] ; i++) {  
      if(Ai[i] >= j) {
        ((int*)(A_chol->i))[nnz] = Ai[i];
        ((double*)(A_chol->x))[2*nnz] = _double(Re(Ax[i]));
        ((double*)(A_chol->x))[2*nnz+1] = _double(Im(Ax[i]));
	nnz++;
      }
    }
  }
  ((int*)(A_chol->p))[RowLen(A)] = nnz;

 
  ::cholmod_factor* factor = (::cholmod_factor*) cholmod_analyze(A_chol,&cm);

  if(cm.status != 0) {    
    err = 1;
    cholmod_free_sparse(&A_chol,&cm);
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);
    return;
  }

  if(!cholmod_factorize(A_chol,factor,&cm)) {
    err = 2;
    cholmod_free_sparse(&A_chol,&cm);  
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);
    return;
  }

  if(cm.status == CHOLMOD_NOT_POSDEF) {
    err = 3;
    cholmod_free_sparse(&A_chol,&cm);
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);    
    return;
  }

  if(!(factor->is_ll)) {   
    cholmod_change_factor(CHOLMOD_COMPLEX, 1, 0, 1, 1, factor, &cm);
  }

  if(cm.status == CHOLMOD_NOT_POSDEF) {
    err = 4;    
    cholmod_free_sparse(&A_chol,&cm);
    cholmod_free_factor(&factor,&cm);
    cholmod_finish(&cm);    
    return;
  }

  Resize(p, A_chol->nrow);
  int *perm_int = (int*)(factor->Perm);
  for(int i=1 ; i<=A_chol->nrow ; i++)
    p[i] = perm_int[i-1];

  
  cholmod_free_sparse(&A_chol,&cm);
  
  ::cholmod_sparse* L_chol = (::cholmod_sparse*) cholmod_factor_to_sparse(factor,&cm);
  

  L = scmatrix(L_chol->nrow, L_chol->ncol, ((int*)(L_chol->p))[L_chol->ncol], (int*)L_chol->p, 
               (int*)L_chol->i, reinterpret_cast<complex*>(L_chol->x), compressed_column);


  cholmod_free_sparse(&L_chol,&cm);
  cholmod_free_factor(&factor,&cm);
  cholmod_finish(&cm);

  err = 0;
}


} //namespace cxsc

