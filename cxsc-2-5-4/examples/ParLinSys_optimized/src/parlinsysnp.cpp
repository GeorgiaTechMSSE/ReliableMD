//============================================================================
//
// 	      ParLinSys - A Solver for Parametric Linear Systems
//
//	    Supplement to the C-XSC Toolbox for Verified Computing
//
// Author: Michael Zimmer, based on a previous version by Evgenija Popova
//
// This program is free software for non-commercial use.
//
// For details on theory see the papers:
//
// S. Rump: Verification methods for dense and sparse systems of equations,
// in: Topics in Validated Computations, ed. J. Herzberger, Elsevier Science
// B.V., 1994,	pp. 63-135.
//
// E. Popova, W. Kraemer: Parametric Fixed-Point Iteration Implemented in
// C-XSC. Universitaet Wuppertal, Preprint BUW - WRSWT, 2003.
//
// M. Zimmer, W. Kraemer, E.Popova: Solvers for the verified solution of 
// parametric linear systems, Computing: Volume 94, Issue 2 (2012), Page 109-123,
// Springer.
//
// For versions using  the alternative algorithm:
// A. Neumaier, A. Pownuk: Linear systems wit hlarge uncertainties, with 
// applications to truss structures. Reliable Computing, 13(149-172), 2007.
//
// This program/module is distributed WITHOUT ANY WARRANTY.
//
//============================================================================

/*
 * This file contains the implementation of the parametric linear solver 
 * based on an alternative algorithm by Neumaier/Pownuk (see paper above)
 * 
 * For systems of the form (K+B*D*A)u=a+Fb with interval entries in D and b
 * only and D diagonal, this algorithm is much faster and gives better enclosures
 * than the algorithm based on the Krawczyk operator.
 * A specialized version for systems of the form (A^T*D*A)u=Fb is also included. Such
 * systems arise for example when using the FEM-method for truss structures.
 */


#include <parlinsys.hpp>
#include <simatrix.hpp>
#include <matinv_aprx.hpp>
#include <string>
#include <fastlss.hpp>

using namespace std;
using namespace cxsc;

namespace cxsc {

//---------------------------------------------------------------------------- 
// Error codes used in this module. 
//----------------------------------------------------------------------------  
 const int 
//   NoError      = 0,   // No error occurred 
//   NotSquare    = 1,   // System to be solved is not square 
//   DimensionErr = 2,   // Dimensions of A and b are not compatible 
//   InvFailed    = 3,   // System is probably singular 
//   VerivFailed  = 4;   // Verification failed, system is probably 
//                       // ill-conditioned 
     EncFailed    = 5;   // Could not compute initial enclosure (other error codes are re-used from fastlss module)
//---------------------------------------------------------------------------- 
// Error messages depending on the error code. 
//---------------------------------------------------------------------------- 
string ParLinSolveErrMsg (int Err) 
{ 
  static string msg;
  if (Err != NoError) { 
    switch (Err) { 
      case NotSquare: 
        msg = "System to be solved is not square"; break; 
      case DimensionErr: 
        msg = "Dimensions of matrices/right hand side are not compatible"; break; 
      case InvFailed: 
        msg = "System is probably singular"; break; 
      case VerivFailed:
        msg = "Verification failed, (parametric) system matrix is not strongly regular"; break; 
      case EncFailed: 
        msg = "Failed to compute initial enclosure"; 
        break; 
      default: 
        msg = "Code not defined"; 
    } 
    msg = "Error: " + msg;
  } 
  return(msg); 
} // end  ParLinSolveErrMsg 

//Compute 2-norm of vector x
real norm(const rvector& x) {
  real s = 0.0;
  
  for(int i=1 ; i<=VecLen(x) ; i++)
    s += x[i]*x[i];
  
  return sqrt(s);
}

//Checks stopping criterion. Returns true if the difference of the componentwise diameter of two successive iterates
//is smaller than Epsilon
bool criterion(const ivector& d, const ivector& dold, real Epsilon = 1e-5) {
  real sumd(0), sumdold(0);
  
  for(int i=1 ; i<=VecLen(d) ; i++) {
    sumd += diam(d[i]);
    sumdold += diam(dold[i]);
  }
  
  return sumd < (1-Epsilon) * sumdold;
}

//! Return the lower bandwidth of A
template<typename TMat>
static inline int lower_bandwidth(const TMat& A) {
  int bw = 0, tmp = 0;
  int n = RowLen(A);

  const std::vector<int>& p   = A.column_pointers();
  const std::vector<int>& ind = A.row_indices();
  for(int j=0 ; j<n ; j++) {
    if(p[j]!=p[j+1]) tmp = ind[p[j+1]-1] - j;
    if(tmp > bw) bw = tmp;
  }

  return bw;
}

//! Return the upper bandwidth of A
template<typename TMat>
static inline int upper_bandwidth(const TMat& A) {
  int bw = 0, tmp = 0;
  int n = RowLen(A);

  const std::vector<int>& p   = A.column_pointers();
  const std::vector<int>& ind = A.row_indices();

  for(int j=0 ; j<n ; j++) {
    if(p[j]!=p[j+1]) tmp = j - ind[p[j]];
    if(tmp > bw) bw = tmp;
  }

  return bw;
}

//Tries to compute an initial enclosure of the solution
template <typename TMat, typename TVec>
bool initialEnclosure(rvector &w, real& alpha, const simatrix &D, const TMat& ACB, const srmatrix &D0mD, const TVec& ACFb, simatrix& DD) {
    
    int m = VecLen(w);
    rvector ww(m), www(m);     
    
    if(diam(D).get_nnz() == RowLen(D)*ColLen(D)) {
      w = 0.0;
      return true;
    }
    
    w   = 1.0;
    DD = mid(D) - D;
        
    rmatrix M = Sup(abs(DD)) * ACB;
    int k = 0;
    ww = w;
    do {
        w = ww / norm(ww);
        ww = M*w;
        k++;
    } while(!(ww < w) && k<10);
    w = ww / norm(ww);
    
    ww  = w - D0mD * ACB * w;
    www = Inf(simatrix(D0mD) * Sup(abs(ACFb)));    

    for(int i=Lb(ww) ; i<=Ub(ww) ; i++) {
      if(ww[i] <= 0)
	return false;
    }
                      
    alpha = www[1] / ww[1];
    real tmp;
    for(int i=2 ; i<=m ; i++) {
         tmp = www[i] / ww[i];
         if(tmp > alpha) alpha = tmp;
    }
        
    return true;
}    

//Solves a system of the form (K+B*D*A)z=a+Fb
//non-verified version
template<typename TK, typename TB, typename TD, typename TA, typename TF>
void ParLinSolveNPMain(const TK& K, const TB& B, const TD& D, const TA& A, ivector& u, 
		       const rvector& a, const TF& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  int n = RowLen(K);
  int m = RowLen(D);
  int k;
  opdotprec = cfg.K;
  
  Err = NoError;
  
  //Dimension check
  if(ColLen(K) != n || ColLen(B) != n || RowLen(B) != ColLen(D) || RowLen(D) != ColLen(D) || 
     RowLen(D) != ColLen(A) || RowLen(A) != n || VecLen(a) != n || ColLen(F) != n || RowLen(F) != VecLen(b) ) {
    Err = DimensionErr;
    return ;
  }
  
  #ifdef _OPENMP
    if(cfg.threads > 0) {
       omp_set_num_threads(cfg.threads);
    } else {
       cfg.threads = omp_get_max_threads(); 
       omp_set_num_threads(cfg.threads);
    }
    if(cfg.msg) cout << "Using " << cfg.threads << " thread(s) for computations" << endl;    
  #endif  

  rmatrix Am(n,n);
  rmatrix C(n,n);

  if(cfg.msg) cout << "Computing midpoint matrix" << endl;
  Am = K + B * mid(D) * A;

  if(cfg.msg) cout << "Computing inverse matrix" << endl;
  MatInvAprx(Am, C, Err); 
  if (Err != NoError)                   // Error: Inversion failed 
    { Err = InvFailed; return; } 
  
  //Compute necessary expressions
  if(cfg.msg) cout << "Computing recurring expressions" << endl;
  rmatrix AC = A * C;
  rmatrix ACB  = AC * B;
  rvector ACa  = AC * a;
  ivector ACFb = (AC * F) * b;
  srmatrix D0mD = Sup(abs(mid(D) - D));
  AC = rmatrix();
  simatrix DD;
    
  
  //Compute initial enclosure
  if(cfg.msg) cout << "Computing initial enclosure" << endl;

  rvector w(m);  
  ivector d(m), dold(m), v(m);
  real alpha;
  
    
  if(!initialEnclosure(w, alpha, D, ACB, D0mD, ACFb, DD)) {
    Err = EncFailed;
    return;
   }
  
  d = (alpha * w) * interval(-1,1);
  v = ACa +ACFb + ACB*d;
  
  if(cfg.msg) cout << "Starting iteration" << endl;
  k = 0;
  do {
    if(cfg.msg) cout << " Step " << k+1 << endl;
    dold = d;
    v = ( ACa + ACFb + ACB * d ) & v;
    d = (DD * v) & d;
    k++;
  } while(k<cfg.maxIter && criterion(d,dold,cfg.epsIter));
  
  if(cfg.msg) cout << "Computing final enclosure" << endl;
  u = C*a + (C*F)*b + (C*B)*d;
  
  return;
}


//Solves a system of the form (K+B*D*A)z=a+Fb
//verified version
template<typename TK, typename TB, typename TD, typename TA, typename TF>
void ParLinSolveNPMainVer(const TK& K, const TB& B, const TD& D, const TA& A, ivector& u, 
		          const rvector& a, const TF& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  int n = RowLen(K);
  int m = RowLen(D);
  int k;

  opdotprec = cfg.K;
  
  Err = NoError;  

  #ifdef _OPENMP
    if(cfg.threads > 0) {
       omp_set_num_threads(cfg.threads);
    } else {
       cfg.threads = omp_get_max_threads(); 
       omp_set_num_threads(cfg.threads);
    }
    if(cfg.msg) cout << "Using " << cfg.threads << " thread(s) for computations" << endl;    
  #endif  

  //Dimension check
  if(ColLen(K) != n || ColLen(B) != n || RowLen(B) != ColLen(D) || RowLen(D) != ColLen(D) || 
     RowLen(D) != ColLen(A) || RowLen(A) != n || VecLen(a) != n || ColLen(F) != n || RowLen(F) != VecLen(b) ) {
    Err = DimensionErr;
    return ;
  }
    
  rmatrix Am(n,n);
  imatrix C(n,n);

  if(cfg.msg) cout << "Computing midpoint matrix" << endl;
  Am = K + B * mid(D) * A;

  if(cfg.msg) cout << "Computing inverse matrix" << endl;
  rmatrix I = Id(Am);
  lssconfig cfglss;
  cfglss.K = 1;
  cfglss.matrixMode = true;  
  lss(Am,I,C,Err,cfglss);
  if (Err != NoError)                   // Error: Inversion failed 
    { Err = InvFailed; return; } 
  I = rmatrix();  
  
  //Compute necessary expressions
  if(cfg.msg) cout << "Computing recurring expressions" << endl;
  imatrix AC = A * C;
  imatrix ACB  = AC * B;
  ivector ACa  = AC * a;
  ivector ACFb = (AC * F) * b;
  srmatrix D0mD = Sup(abs(mid(D) - D));
  simatrix DD;
    
  AC = imatrix();
  
  //Compute initial enclosure
  if(cfg.msg) cout << "Computing initial enclosure" << endl;

  rvector w(m);
  real alpha;

  if(!initialEnclosure(w, alpha, D, Sup(abs(ACB)), D0mD, ACFb, DD)) {
    Err = EncFailed;
    return;
  }
  
  ivector d(m), dold(m), v(m);
  d = (alpha * w) * interval(-1,1);
  v = ACa +ACFb + ACB*d;

  if(cfg.msg) cout << "Starting iteration" << endl;
  k = 0;
  do {
    if(cfg.msg) cout << " Step " << k+1 << endl;
    dold = d;
    v = ( ACa + ACFb + ACB * d ) & v;
    d = (DD * v) & d;
    k++;
  } while(k<cfg.maxIter && criterion(d,dold,cfg.epsIter));
  
  if(cfg.msg) cout << "Computing final enclosure" << endl;
  u = C*a + (C*F)*b + (C*B)*d;
  
  return;
}

/*
 *  Solve system (A^T * D * A)* u = F*b with intervals only in D and b.
 *  D must be positive definite and should be a diagonal matrix
 *  verified version
 */
template<typename TA, typename TD, typename TF>
void ParLinSolveNPMainVer(const TA& A, const TD& D, ivector& u, 
		          const TF& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  
  int n = RowLen(A);
  int m = RowLen(D);
  int k;
    
  opdotprec = cfg.K;
  
  Err = NoError;  

  #ifdef _OPENMP
    if(cfg.threads > 0) {
       omp_set_num_threads(cfg.threads);
    } else {
       cfg.threads = omp_get_max_threads(); 
       omp_set_num_threads(cfg.threads);
    }
    if(cfg.msg) cout << "Using " << cfg.threads << " thread(s) for computations" << endl;    
  #endif  

  //Dimension check
  if(RowLen(D) != ColLen(D) || RowLen(D) != ColLen(A) || RowLen(A) != n || ColLen(F) != n || RowLen(F) != VecLen(b) ) {
    Err = DimensionErr;
    return ;
  }
    
  rmatrix Am(n,n);
  imatrix C(n,n);

  if(cfg.msg) cout << "Computing midpoint matrix" << endl;
  TA B = transp(A);
  Am = B * mid(D) * A;

  if(cfg.msg) cout << "Computing inverse matrix" << endl;
  rmatrix I = Id(Am);
  lssconfig cfglss;
  cfglss.K = 1;
  cfglss.matrixMode = true;     
  lss(Am, I, C, Err, cfglss); 
  if (Err != NoError)                   // Error: Inversion failed 
    { Err = InvFailed; return; } 
  
  //Compute necessary expressions
  if(cfg.msg) cout << "Computing recurring expressions" << endl;
  imatrix AC = A * C;
  imatrix ACB  = AC * B;
  imatrix ACF = AC * F;
  ivector ACFb = ACF * b;
  srmatrix D0mD = Sup(abs(mid(D) - D));
  simatrix DD;
    
  
  //Compute initial enclosure
  if(cfg.msg) cout << "Computing initial enclosure" << endl;

  rvector w(m);
  real alpha;
    
  ivector d(m), dold(m), v(m);
  ivector c(m), z(m);
    
  if(!initialEnclosure(w, alpha, D, Sup(abs(ACB)), D0mD, ACFb, DD)) {
      if(cfg.msg) cout << "  Method 1 failed, trying method 2" << endl;
      
      if(lower_bandwidth(D) != 0 || upper_bandwidth(D) != 0) {
	Err = EncFailed;
	return;
      } else {
	const vector<interval>& val = D.values();
	for(int i=0 ; i<val.size() ; i++) {
	  if(Inf(val[i]) <= 0) {
	    Err = EncFailed;
	    if(cfg.msg) cout << "  Method 2 failed" << endl;
	    return;
	  }
	}
      }

      c = (mid(D)*ACF)*b;
      
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int i=1 ; i<=m ; i++) {
          real sum = 0;
          for(int k=1 ; k<=m ; k++)
              sum += (Sup(abs(c(k)))*Sup(abs(c(k)))) / Inf(D(k,k)); 
          z[i] = 0.5 * (c[i] + interval(-1,1) * sqrt( Sup(abs(c[i]))*Sup(abs(c[i]))*(1-Sup(D(i,i))/Inf(D(i,i))) + Sup(D(i,i))*sum ) );  
          v[i] = z[i] / D(i,i);
          d[i] = (mid(D(i,i))/D(i,i) - 1) * z[i];
      }

  } else {
    
    d = (alpha * w) * interval(-1,1);
    v = ACFb + ACB*d;
  }
    

  if(cfg.msg) cout << "Starting iteration" << endl;
  k = 0;
  do {
    if(cfg.msg) cout << " Step " << k+1 << endl;
    dold = d;
    v = ( ACFb + ACB * d ) & v;
    d = (DD * v) & d;
    k++;
  } while(k<cfg.maxIter && criterion(d,dold,cfg.epsIter));
  
  if(cfg.msg) cout << "Computing final enclosure" << endl;
  u = (C*F)*b + (C*B)*d;
  
  return;
}


/*
 *  Solve system (A^T * D * A)* u = F*b with intervals only in D and b.
 *  D must be positive definite and should be a diagonal matrix
 *  non-verified version
 */
template<typename TA, typename TD, typename TF>
void ParLinSolveNPMain(const TA& A, const TD& D, ivector& u, 
		       const TF& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {

  int n = RowLen(A);
  int m = RowLen(D);
  int k;  

  opdotprec = cfg.K;
  
  Err = NoError;  

  #ifdef _OPENMP
    if(cfg.threads > 0) {
       omp_set_num_threads(cfg.threads);
    } else {
       cfg.threads = omp_get_max_threads(); 
       omp_set_num_threads(cfg.threads);
    }
    if(cfg.msg) cout << "Using " << cfg.threads << " thread(s) for computations" << endl;    
  #endif  
  
  //Dimension check
  if(RowLen(D) != ColLen(D) || RowLen(D) != ColLen(A) || RowLen(A) != n || ColLen(F) != n || RowLen(F) != VecLen(b) ) {
    Err = DimensionErr;
    return ;
  }

  rmatrix Am(n,n);
  rmatrix C(n,n);

  if(cfg.msg) cout << "Computing midpoint matrix" << endl;
  TA B = transp(A);
  Am = B * mid(D) * A;

  if(cfg.msg) cout << "Computing inverse matrix" << endl;
  MatInvAprx(Am, C, Err); 
  if (Err != NoError)                   // Error: Inversion failed 
    { Err = InvFailed; return; } 
  
  //Compute necessary expressions
  if(cfg.msg) cout << "Computing recurring expressions" << endl;
  rmatrix AC = A * C;
  rmatrix ACB  = AC * B;
  rmatrix ACF = AC * F;
  ivector ACFb = ACF * b;
  srmatrix D0mD = Sup(abs(mid(D) - D));
  simatrix DD;
    
  
  //Compute initial enclosure
  if(cfg.msg) cout << "Computing initial enclosure" << endl;

    rvector w(m);
    real alpha;
    
    ivector d(m), dold(m), v(m);
    ivector c(m), z(m);
    
    if(!initialEnclosure(w, alpha, D, ACB, D0mD, ACFb, DD)) {
        if(cfg.msg) cout << "  Method 1 failed,trying method 2" << endl;
      
        if(lower_bandwidth(D) != 0 || upper_bandwidth(D) != 0) {
	  Err = EncFailed;
	  return;
        } else {
	  const vector<interval>& val = D.values();
	  for(int i=0 ; i<val.size() ; i++) {
	    if(Inf(val[i]) <= 0) {
	      Err = EncFailed;
              if(cfg.msg) cout << "  Method 2 failed" << endl;	      
	      return;
	    }
	 }
        }
      
        c = (mid(D)*ACF)*b;
        
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i=1 ; i<=m ; i++) {
            real sum = 0;
            for(int k=1 ; k<=m ; k++)
                sum += (Sup(abs(c(k)))*Sup(abs(c(k)))) / Inf(D(k,k)); 
            z[i] = 0.5 * (c[i] + interval(-1,1) * sqrt( Sup(abs(c[i]))*Sup(abs(c[i]))*(1-Sup(D(i,i))/Inf(D(i,i))) + Sup(D(i,i))*sum ) );  
            v[i] = z[i] / D(i,i);
            d[i] = (mid(D(i,i))/D(i,i) - 1) * z[i];
        }
        
    } else {
        
        d = (alpha * w) * interval(-1,1);
        v = ACFb + ACB*d;
    }


  if(cfg.msg) cout << "Starting iteration" << endl; 
  k = 0;
  do {
    if(cfg.msg) cout << " Step " << k+1 << endl;
    dold = d;
    v = ( ACFb + ACB * d ) & v;
    d = (DD * v) & d;
    k++;
  } while(k<cfg.maxIter && criterion(d,dold,cfg.epsIter));
  
  if(cfg.msg) cout << "Computing final enclosure" << endl;
  u = (C*F)*b + (C*B)*d;
  
  return;
}



//Entry functions for the solver called by the user

void ParLinSolve(const srmatrix& A, const simatrix& D, ivector& u, const srmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  if(cfg.verified) {
    ParLinSolveNPMainVer(A,D,u,F,b,Err,cfg);
  } else {
    ParLinSolveNPMain(A,D,u,F,b,Err,cfg);
  }
}

void ParLinSolve(const rmatrix& A, const simatrix& D, ivector& u, const rmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  if(cfg.verified) {
    ParLinSolveNPMainVer(A,D,u,F,b,Err,cfg);
  } else {
    ParLinSolveNPMain(A,D,u,F,b,Err,cfg);
  }
}

void ParLinSolve(const srmatrix& K, const srmatrix& B, const simatrix& D, const srmatrix& A, ivector& u, 
		   const rvector& a, const srmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  if(cfg.verified) {
    ParLinSolveNPMainVer(K,B,D,A,u,a,F,b,Err,cfg);
  } else {
    ParLinSolveNPMain(K,B,D,A,u,a,F,b,Err,cfg);
  }		     
}

void ParLinSolve(const rmatrix& K, const rmatrix& B, const simatrix& D, const rmatrix& A, ivector& u, 
		   const rvector& a, const rmatrix& F, const ivector& b, int& Err, parlinsysnpconfig cfg) {
  if(cfg.verified) {
    ParLinSolveNPMainVer(K,B,D,A,u,a,F,b,Err,cfg);
  } else {
    ParLinSolveNPMain(K,B,D,A,u,a,F,b,Err,cfg);
  }		     
}


} //namespace cxsc
