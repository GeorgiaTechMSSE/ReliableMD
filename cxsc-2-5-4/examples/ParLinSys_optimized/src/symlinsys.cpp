//============================================================================
//
//           Program/Module "Symmetric Linear Interval System Solver"
//
//	       supplement to the C-XSC Toolbox for Verified Computing
//
// Author: Michael Zimmer
//
// This program/module is free software for non-commercial use.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================

/* 
 * Containt the implementation of the symmetric linear interval system solver.
 * Symmetric systems are transformed into equivalent parametric systems and solved
 * with the parametric linear systems solver. During the transformation, one parameter
 * for each pair a_ij=a_ji of non-diagonal entries of A is introduced.
 * Other symmetries, for example if whole diagonals are symmetric, should be handled manually.
 */

#include <symlinsys.hpp>
#include <timer.hpp>

using namespace std;
using namespace cxsc;

namespace cxsc {

 //Return the error message for an error code
string SymLinSolveErrMsg ( int err ) {
  return ParLinSolveErrMsg(err);
}

//Tranform symmetric linear interval system into corresponding parametric system
//and compute an inner and outer enclosure using the parametric solver
template<class TMat, class TRhs, class TVec, class TSMat, class TCoeff>
void  SymLinSolveInnerOuter (TMat& A, TRhs& b, TRhs& x, TRhs& y, int& Err, struct parlinsysconfig cfg) {
  int k=0;
  int n = RowLen(A);
  
  if(cfg.msg) cout << "Converting symmetric system to parametric system..." << endl;

  for(int i=1 ; i<=n; i++) k+=i;

  TVec p(k+n+1);
  p = 0.0;
  p[1] = 1.0;
  int pind = 2;

  vector<TCoeff> coeff;
  coeff.push_back(TSMat(n,n));
  
  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=i ; j++) {
      TSMat C(n,n);
      if(i==j) {
	if(A[Lb(A,1)+i-1][Lb(A,2)+j-1] != 0) {
          C[i][j] = 1.0;
          p[pind++] = A[Lb(A,1)+i-1][Lb(A,2)+j-1];
          coeff.push_back(TCoeff(C));  
	}
      } else {
	if(A[Lb(A,1)+i-1][Lb(A,2)+j-1] != 0) {
          C[i][j] = 1.0;
          C[j][i] = 1.0;
          p[pind++] = A[Lb(A,1)+i-1][Lb(A,2)+j-1];
          coeff.push_back(TCoeff(C));  
	}
      }
    }
  }

  TMat B(b);  
  k = pind - 1;  
  TSMat rhs(n,k + RowLen(B)*n);
  for(int j=1 ; j<=RowLen(B) ; j++) {
    for(int i=1 ; i<=n ; i++) {
      rhs[i][k+i] = 1.0;
      p[pind++] = B[i][j];
      coeff.push_back(TSMat(n,n));
    }
  }
  
  Resize(p,pind-1);

  ParLinSolve(coeff, rhs, p, x, y, Err, cfg);
}


//Tranform symmetric linear interval system into corresponding parametric system
//and compute an outer enclosure using the parametric solver
template<class TMat, class TRhs, class TVec, class TSMat, class TCoeff>
void  SymLinSolveOuter (TMat& A, TRhs& b, TRhs& x, int& Err, struct parlinsysconfig cfg) {
  int k=0;
  int n = RowLen(A);

  if(cfg.msg) cout << "Converting symmetric system to parametric system..." << endl;

  for(int i=1 ; i<=n; i++) k+=i;

  TVec p(k+n+1);
  p = 0.0;
  p[1] = 1.0;
  int pind = 2;

  vector<TCoeff> coeff;
  coeff.push_back(TSMat(n,n));
  
  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=i ; j++) {
      TSMat C(n,n);
      if(i==j) {
	if(A[Lb(A,1)+i-1][Lb(A,2)+j-1] != 0) {
          C[i][j] = 1.0;
          p[pind++] = A[Lb(A,1)+i-1][Lb(A,2)+j-1];
          coeff.push_back(TCoeff(C));  
	}
      } else {
	if(A[Lb(A,1)+i-1][Lb(A,2)+j-1] != 0) {
          C[i][j] = 1.0;
          C[j][i] = 1.0;
          p[pind++] = A[Lb(A,1)+i-1][Lb(A,2)+j-1];
          coeff.push_back(TCoeff(C));  
	}
      }
    }
  }

  TMat B(b);  
  k = pind - 1;  
  TSMat rhs(n,k + RowLen(B)*n);
  for(int j=1 ; j<=RowLen(B) ; j++) {
    for(int i=1 ; i<=n ; i++) {
      rhs[i][k+i] = 1.0;
      p[pind++] = B[i][j];
      coeff.push_back(TSMat(n,n));
    }
  }
  
  Resize(p,pind-1);
  
  ParLinSolve(coeff, rhs, p, x, Err, cfg);
}




//Entry function for the solver called by the user

void SymLinSolve(imatrix& A, ivector& b, ivector& x, ivector& y, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveInnerOuter<imatrix,ivector,ivector,srmatrix,CoeffMatrix>(A,b,x,y,Err,cfg);
}

void SymLinSolve(imatrix& A, imatrix& b, imatrix& x, imatrix& y, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveInnerOuter<imatrix,imatrix,ivector,srmatrix,CoeffMatrix>(A,b,x,y,Err,cfg);
}

void SymLinSolve(cimatrix& A, civector& b, civector& x, civector& y, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveInnerOuter<cimatrix,civector,civector,srmatrix,CoeffMatrix>(A,b,x,y,Err,cfg);
}

void SymLinSolve(cimatrix& A, cimatrix& b, cimatrix& x, cimatrix& y, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveInnerOuter<cimatrix,cimatrix,civector,srmatrix,CoeffMatrix>(A,b,x,y,Err,cfg);
}


void SymLinSolve(imatrix& A, ivector& b, ivector& x, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveOuter<imatrix,ivector,ivector,srmatrix,CoeffMatrix>(A,b,x,Err,cfg);
}

void SymLinSolve(imatrix& A, imatrix& b, imatrix& x, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveOuter<imatrix,imatrix,ivector,srmatrix,CoeffMatrix>(A,b,x,Err,cfg);
}

void SymLinSolve(cimatrix& A, civector& b, civector& x, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveOuter<cimatrix,civector,civector,srmatrix,CoeffMatrix>(A,b,x,Err,cfg);
}

void SymLinSolve(cimatrix& A, cimatrix& b, cimatrix& x, int& Err, struct parlinsysconfig cfg) {
  SymLinSolveOuter<cimatrix,cimatrix,civector,srmatrix,CoeffMatrix>(A,b,x,Err,cfg);
}

} //namespace cxsc