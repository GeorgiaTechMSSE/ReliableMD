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

/* Implementation of utility functions used by the solvers */

#include "utility.hpp"
#include <scimatrix.hpp>
#include <vector>
#include <fenv.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cxsc {
  
//! Epsilon inflation for complex intervals
civector blow(const civector& x, const real eps) {
  civector y(x);
  SetRe( y, Blow(Re(x),eps) );
  SetIm( y, Blow(Im(x),eps) );
  return y;
}

//! Epsilon inflation for real intervals
ivector blow(const ivector& x, const real eps) {
  return Blow(x,eps);
}

//! Computes upper bounds for the infinity norm of A
real Norm00(const srmatrix& A) {
  real maxi = 0.0;
  int n = RowLen(A);
  rvector tmp(n);

  const std::vector<int>& p    = A.column_pointers();
  const std::vector<int>& ind  = A.row_indices();
  const std::vector<real>& val = A.values();

  tmp = 0.0;

  fesetround(FE_UPWARD);

  for(int j=0 ; j<n ; j++) {
    for(int k=p[j] ; k<p[j+1] ; k++) {
      tmp[ind[k]+1] += abs(val[k]);
    }
  }

  for(int i=1 ; i<=n ; i++) {
    if(tmp[i] > maxi) maxi = tmp[i];
  }

  fesetround(FE_TONEAREST);

  return maxi;
}

//! Computes upper bounds for the infinity norm of A
real Norm00(const scmatrix& A) {
  real maxi = 0.0;
  int n = RowLen(A);
  rvector tmp(n);

  const std::vector<int>& p    = A.column_pointers();
  const std::vector<int>& ind  = A.row_indices();
  const std::vector<complex>& val = A.values();

  tmp = 0.0;

  fesetround(FE_UPWARD);

  for(int j=0 ; j<n ; j++) {
    for(int k=p[j] ; k<p[j+1] ; k++) {
      tmp[ind[k]+1] += abs(val[k]);
    }
  }

  for(int i=1 ; i<=n ; i++) {
    if(tmp[i] > maxi) maxi = tmp[i];
  }

  fesetround(FE_TONEAREST);

  return maxi;
}

//! Computes upper bounds for the 1-norm of A
real Norm1(const srmatrix& A) {
  real maxj=0.0, sumi;
  int n = RowLen(A);

  const std::vector<int>& p    = A.column_pointers();
  const std::vector<real>& val = A.values();

  fesetround(FE_UPWARD);

  for(int j=0 ; j<n ; j++) {

    sumi = 0.0;
    for(int k=p[j] ; k<p[j+1] ; k++) {
      sumi += abs(val[k]);
    }

    if(sumi > maxj) maxj = sumi;

  }
 
  fesetround(FE_TONEAREST);

  return maxj;
}

//! Computes upper bounds for the 1-norm of A
real Norm1(const scmatrix& A) {
  real maxj=0.0, sumi;
  int n = RowLen(A);

  const std::vector<int>& p    = A.column_pointers();
  const std::vector<complex>& val = A.values();

  fesetround(FE_UPWARD);

  for(int j=0 ; j<n ; j++) {

    sumi = 0.0;
    for(int k=p[j] ; k<p[j+1] ; k++) {
      sumi += abs(val[k]);
    }

    if(sumi > maxj) maxj = sumi;

  }
 
  fesetround(FE_TONEAREST);

  return maxj;
}

//! Computes upper bounds for the 2-norm of a
real Norm2(const rvector& a) {
  int prec = opdotprec;
  opdotprec = 1;

  fesetround(FE_UPWARD);
  real max = sqrt(a*a);
  fesetround(FE_TONEAREST);

  opdotprec = prec;

  return max;
}

//! Computes upper bounds for the 2-norm of a
real Norm2(const cvector& a) {
  int prec = opdotprec;
  opdotprec = 1;

  fesetround(FE_UPWARD);
  real max = sqrt(Re(a)*Re(a) + Im(a)*Im(a));
  fesetround(FE_TONEAREST);

  opdotprec = prec;

  return max;
}

//! Computes upper bounds for the 2-norm of v
real Norm2(const srvector& v) {
  real norm = 0.0;

  const std::vector<real>& val = v.values();

  fesetround(FE_UPWARD);

  for(unsigned int i=0 ; i<val.size() ; i++) {
    norm += sqr(val[i]);
  }
  
  real ret = sqrt(norm);

  fesetround(FE_TONEAREST);

  return ret;
}

//! Computes upper bounds for the 2-norm of v
real Norm2(const scvector& v) {
  real norm = 0.0;

  const std::vector<complex>& val = v.values();

  fesetround(FE_UPWARD);

  for(unsigned int i=0 ; i<val.size() ; i++) {
    norm += sqr(Re(val[i])) + sqr(Im(val[i]));
  }
  
  real ret = sqrt(norm);

  fesetround(FE_TONEAREST);

  return ret;
}

//! Computes upper bounds for the 2-norm of a
interval Norm2(const ivector& a) {
  interval ret(0.0);
  for(int i=Lb(a) ; i<=Ub(a) ; i++)
    ret += sqr(a[i]);
  return sqrt(ret);
}

//! Computes upper bounds for the 2-norm of a
interval Norm2(const civector& a) {
  interval ret(0.0);
  for(int i=Lb(a) ; i<=Ub(a) ; i++)
    ret += sqr(Re(a[i])) + sqr(Im(a[i]));
  return sqrt(ret);
}

inline real sum(const rvector& x) {
  real s = 0.0;
  for(int i=1 ; i<=VecLen(x) ; i++)
    s += x[i];
  return s;
}

inline real sum(const srvector& x) {
  const std::vector<real>& val = x.values();
  real s = 0.0;
  for(int i=0 ; i<val.size() ; i++)
    s += val[i];
  return s;
}

/*! Binormalization of A: Compue diagonal matrices D1 and D2 such that D1*A*D2
 *  is cloeser to being doubly stochastic */
void normbin(const srmatrix& A, srmatrix& D1, srmatrix& D2) {
  int n = RowLen(A);
  srmatrix T = transp(A);

  D1 = Id(A);
  D2 = Id(A);

  std::vector<real>& D1_val = D1.values();
  std::vector<real>& D2_val = D2.values();

  for(int i=0 ; i<n ; i++) {
    D1_val[i] = 1.0 / sum(abs(srvector(T[Col(i+1)])));
  }

  T = D1 * A;

  for(int i=0 ; i<n ; i++) {
    D2_val[i] = 1.0 / sum(abs(srvector(T[Col(i+1)])));
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++) {
    D1_val[i-1] = ::pow(2, static_cast<int>(_double(q_log2(D1_val[i-1]))) );
    D2_val[i-1] = ::pow(2, static_cast<int>(_double(q_log2(D2_val[i-1]))) );
  }
}

/*! Binormalization of A: Compue diagonal matrices D1 and D2 such that D1*A*D2
 *  is cloeser to being doubly stochastic */
void normbin(const scmatrix& A, srmatrix& D1, srmatrix& D2) {
  int n = RowLen(A);
  scmatrix T = transp(A);

  D1 = Id(srmatrix(n,n));
  D2 = Id(srmatrix(n,n));

  std::vector<real>& D1_val = D1.values();
  std::vector<real>& D2_val = D2.values();
  
  for(int k=0 ; k<3 ; k++) {

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=0 ; i<n ; i++) {
    D1_val[i] = 1.0 / Norm2(T[Col(i+1)]);
  }

  T = D1 * A;

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=0 ; i<n ; i++) {
    D2_val[i] = 1.0 / Norm2(T[Col(i+1)]);
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++) {
    D1_val[i-1] = ::pow(2, static_cast<int>(_double(q_log2(D1_val[i-1]))) );
    D2_val[i-1] = ::pow(2, static_cast<int>(_double(q_log2(D2_val[i-1]))) );
  }
  
  T = T * D2;
  
  }
}


/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(rvector& def, const rvector& b, const srmatrix& A, const rmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);  
  
  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * xapp[Col(i)];
    opdotprec = 1;
    return;
  }

  std::vector<dotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif

  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<real>& val = A.values();
  
  for(int i=1 ; i<=n ; i++) {
    res.push_back(dotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
        #ifdef _OPENMP
	omp_set_lock(&(locks[ri[k]]));
        #endif
        for(int c=1 ; c<=cols ; c++) {
          accumulate(res[ri[k]], -val[k], xapp[j+1][c-1+lb]);
	}
        #ifdef _OPENMP
        omp_unset_lock(&(locks[ri[k]]));
        #endif
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++)
    def[i] = rnd(res[i-1],RND_UP);
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(ivector& def, const rvector& b, const srmatrix& A, const rmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);

  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * ivector(xapp[Col(i)]);
    opdotprec = 1;
    return;
  }
    
  std::vector<dotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif  

  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<real>& val = A.values();

  for(int i=1 ; i<=n ; i++) {
    res.push_back(dotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif    
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
      #ifdef _OPENMP
      omp_set_lock(&(locks[ri[k]]));
      #endif	
      for(int c=1 ; c<=cols ; c++) {
        accumulate(res[ri[k]], -val[k], xapp[j+1][c-1+lb]);
      }
      #ifdef _OPENMP
      omp_unset_lock(&(locks[ri[k]]));
      #endif	
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++)
    rnd(res[i-1],def[i]);
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(ivector& def, const ivector& b, const srmatrix& A, const rmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);

  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * ivector(xapp[Col(i)]);
    opdotprec = 1;
    return;
  }  
  
  std::vector<idotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif  
  
  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<real>& val = A.values();

  for(int i=1 ; i<=n ; i++) {
    res.push_back(idotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif     
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
      #ifdef _OPENMP
      omp_set_lock(&(locks[ri[k]]));
      #endif	
      for(int c=1 ; c<=cols ; c++) {
        accumulate(res[ri[k]], -val[k], xapp[j+1][c-1+lb]);
      }
      #ifdef _OPENMP
      omp_unset_lock(&(locks[ri[k]]));
      #endif
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++) {
    rnd(res[i-1],def[i]);
  }
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(ivector& def, const ivector& b, const simatrix& A, const rmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  
  int lb = Lb(xapp,COL);

  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * xapp[Col(i)];
    opdotprec = 1;
    return;
  }  
  
  std::vector<idotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif    
  ivector x(n);
  dotprecision dot(0.0);

  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<interval>& val = A.values();


  for(int i=1 ; i<=n ; i++) {
    res.push_back(idotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif     
  }
  
  #ifdef _OPENMP
  #pragma omp parallel for private(dot) 
  #endif     
  for(int i=1 ; i<=n ; i++) {
    dot = 0.0;
    for(int c=1 ; c<=cols ; c++)
      dot += xapp[i][c-1+lb];
    rnd(dot, x[i]);
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
      #ifdef _OPENMP
      omp_set_lock(&(locks[ri[k]]));
      #endif      
      accumulate(res[ri[k]], -val[k], x[j+1]);
      #ifdef _OPENMP
      omp_unset_lock(&(locks[ri[k]]));
      #endif
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++) {
    rnd(res[i-1],def[i]);
  }
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(cvector& def, const cvector& b, const scmatrix& A, const cmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);
  
  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * xapp[Col(i)];
    opdotprec = 1;
    return;
  }  
  
  std::vector<cdotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif  
  
  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<complex>& val = A.values();


  for(int i=1 ; i<=n ; i++) {
    res.push_back(cdotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif     
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
        #ifdef _OPENMP
        omp_set_lock(&(locks[ri[k]]));
        #endif	
        for(int c=1 ; c<=cols ; c++) {
          accumulate(res[ri[k]], -val[k], xapp[j+1][c-1+lb]);
	}
        #ifdef _OPENMP
        omp_unset_lock(&(locks[ri[k]]));
        #endif	
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++)
    def[i] = rnd(res[i-1],RND_UP);
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(civector& def, const cvector& b, const scmatrix& A, const cmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);
  
  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * civector(xapp[Col(i)]);
    opdotprec = 1;
    return;
  }  
  
  std::vector<cdotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif    

  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<complex>& val = A.values();


  for(int i=1 ; i<=n ; i++) {
    res.push_back(cdotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif     
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
      #ifdef _OPENMP
      omp_set_lock(&(locks[ri[k]]));
      #endif	
      for(int c=1 ; c<=cols ; c++) {
        accumulate(res[ri[k]], -val[k], xapp[j+1][c-1+lb]);
      }
      #ifdef _OPENMP
      omp_unset_lock(&(locks[ri[k]]));
      #endif      
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++)
    rnd(res[i-1],def[i]);
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(civector& def, const civector& b, const scmatrix& A, const cmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);

  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * civector(xapp[Col(i)]);
    opdotprec = 1;
    return;
  }  
  
  std::vector<cidotprecision> res;
  res.reserve(n);
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif  
  
  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<complex>& val = A.values();


  for(int i=1 ; i<=n ; i++) {
    res.push_back(cidotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif     
  }

  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t)
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
      #ifdef _OPENMP
      omp_set_lock(&(locks[ri[k]]));
      #endif	      
      for(int c=1 ; c<=cols ; c++) {
        accumulate(res[ri[k]], -val[k], xapp[j+1][c-1+lb]);
      }
      #ifdef _OPENMP
      omp_unset_lock(&(locks[ri[k]]));
      #endif
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++) {
    rnd(res[i-1],def[i]);
  }
}

/*! Computes the residuum def=b-A*xapp of an approximate solution represented by the sum of the 
 *  first cols columns of xapp. Uses the long accumulator, results are only one rounding away from
 *  the exact result. */
void residual(civector& def, const civector& b, const scimatrix& A, const cmatrix& xapp, int cols, bool dotp) {
  int n = VecLen(b);
  int lb = Lb(xapp,COL);
  
  if(!dotp) {
    opdotprec = 3;
    def = b;
    for(int i=1 ; i<=cols ; i++)
      def -= A * xapp[Col(i)];
    opdotprec = 1;
    return;
  }  
  
  std::vector<cidotprecision> res;
  res.reserve(n);
  civector x(n);
  cdotprecision dot;
  #ifdef _OPENMP
  std::vector<omp_lock_t> locks(n);
  #endif  
  
  const std::vector<int>& cp   = A.column_pointers();
  const std::vector<int>& ri   = A.row_indices();
  const std::vector<cinterval>& val = A.values();


  for(int i=1 ; i<=n ; i++) {
    res.push_back(cidotprecision(b[i]));
    #ifdef _OPENMP
    omp_init_lock(&(locks[i-1]));
    #endif     
  }
  
  #ifdef _OPENMP
  int t = omp_get_max_threads();
  #pragma omp parallel for schedule(static,n/t) private(dot)
  #endif     
  for(int i=1 ; i<=n ; i++) {
    dot = 0.0;
    for(int c=1 ; c<=cols ; c++)
      dot += xapp[i][c-1+lb];
    rnd(dot, x[i]);
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif    
  for(int j=0 ; j<n ; j++) {
    for(int k=cp[j] ; k<cp[j+1] ; k++) {
     #ifdef _OPENMP
     omp_set_lock(&(locks[ri[k]]));
     #endif	
     accumulate(res[ri[k]], -val[k], x[j+1]);
     #ifdef _OPENMP
     omp_unset_lock(&(locks[ri[k]]));
     #endif
    }
  }

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++) {
    rnd(res[i-1],def[i]);
  }
}

//! Computes D=L*U-A. Function is used to prevent compiler optimizations eliminating the rounding mode setting
void LUMinusA(srmatrix& D, const srmatrix& L, const srmatrix& U, const srmatrix& A) {
  D = L*U - A;
}

//! Computes D=L*U-A. Function is used to prevent compiler optimizations eliminating the rounding mode setting
void LUMinusA(scmatrix& D, const scmatrix& L, const scmatrix& U, const scmatrix& A) {
  D = L*U - A;
}


//! Compute the columns with indices from start to end of D=L*U-A
inline void LUMinusAChunk(const srmatrix& L, const srmatrix& U, const srmatrix& A, srmatrix& D, int start, int end) {
  int n = RowLen(A);

  D = srmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<real>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<real>& U_x  = U.values();
  const std::vector<int>& A_ind = A.row_indices();
  const std::vector<int>& A_p   = A.column_pointers();
  const std::vector<real>& A_x  = A.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<real>& D_x  = D.values();

  std::vector<real> dot(n,real(0.0));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
          dot[L_ind[l]] = L_x[l] * U_x[k];
          nnz++;
        } else {
          dot[L_ind[l]] += L_x[l] * U_x[k];
        }
      } 
    }
    
    for(int l=A_p[i] ; l<A_p[i+1] ; l++) {
      if(w[A_ind[l]] < i) {
        w[A_ind[l]] = i;
        D_ind.push_back(A_ind[l]);
        dot[A_ind[l]] = -A_x[l];
        nnz++;
      } else {
        dot[A_ind[l]] -= A_x[l];
      }
    }   

    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       D_x.push_back(dot[D_ind[j]]);
    }

    D_p[i-start+1+1] = nnz;
  }
  
}

//! Compute the columns with indices from start to end of D=L*U-A
inline void LUMinusAChunk(const scmatrix& L, const scmatrix& U, const scmatrix& A, scmatrix& D, int start, int end) {
  int n = RowLen(A);

  D = scmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<complex>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<complex>& U_x  = U.values();
  const std::vector<int>& A_ind = A.row_indices();
  const std::vector<int>& A_p   = A.column_pointers();
  const std::vector<complex>& A_x  = A.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<complex>& D_x  = D.values();

  std::vector<complex> dot(n,complex(0.0));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
          //dot[L_ind[l]] = L_x[l] * U_x[k];
          fp_mult(L_x[l],U_x[k],dot[L_ind[l]]);	  
          nnz++;
        } else {
          //dot[L_ind[l]] += L_x[l] * U_x[k];
          fp_multadd(L_x[l],U_x[k],dot[L_ind[l]]);          
        }
      }
    }
    
    for(int l=A_p[i] ; l<A_p[i+1] ; l++) {
      if(w[A_ind[l]] < i) {
        w[A_ind[l]] = i;
        D_ind.push_back(A_ind[l]);
        dot[A_ind[l]] = -A_x[l];
        nnz++;
      } else {
        dot[A_ind[l]] -= A_x[l];
      }
    }    

    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       D_x.push_back(dot[D_ind[j]]);
    }

    D_p[i-start+1+1] = nnz;
  }
  
}


//! Compute the columns with indices from start to end of D=L*U-A in higher precision (rounding mode is not used here)
inline void LUMinusAHighPrecChunk(const srmatrix& L, const srmatrix& U, const srmatrix& A, srmatrix& D, int start, int end) {
  int n = RowLen(A);

  D = srmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<real>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<real>& U_x  = U.values();
  const std::vector<int>& A_ind = A.row_indices();
  const std::vector<int>& A_p   = A.column_pointers();
  const std::vector<real>& A_x  = A.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<real>& D_x  = D.values();

  std::vector<sparse_dot> dot(n,sparse_dot(opdotprec));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
	  dot[L_ind[l]].reset();
          dot[L_ind[l]].add_dot(L_x[l],U_x[k]);
          nnz++;
        } else {
          dot[L_ind[l]].add_dot(L_x[l],U_x[k]);
        }
      } 
    }
    
    for(int l=A_p[i] ; l<A_p[i+1] ; l++) {
      if(w[A_ind[l]] < i) {
        w[A_ind[l]] = i;
        D_ind.push_back(A_ind[l]);
        dot[A_ind[l]].reset();
        dot[A_ind[l]].add_dot(-1,A_x[l]);
        nnz++;
      } else {
        dot[A_ind[l]].add_dot(-1,A_x[l]);
      }
    }   

    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);
    dotprecision accu;

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       dot[D_ind[j]].result(accu);
       interval tmp;
       rnd(accu,tmp);
       D_x.push_back(Sup(abs(tmp)));
    }

    D_p[i-start+1+1] = nnz;
  }
  
}

//! Compute the columns with indices from start to end of D=L*U-A in higher precision (rounding mode is not used here)
inline void LUMinusAHighPrecChunk(const scmatrix& L, const scmatrix& U, const scmatrix& A, scmatrix& D, int start, int end) {
  int n = RowLen(A);

  D = scmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<complex>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<complex>& U_x  = U.values();
  const std::vector<int>& A_ind = A.row_indices();
  const std::vector<int>& A_p   = A.column_pointers();
  const std::vector<complex>& A_x  = A.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<complex>& D_x  = D.values();

  std::vector<sparse_cdot> dot(n,sparse_cdot(opdotprec));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
	  dot[L_ind[l]].reset();
          dot[L_ind[l]].add_dot(L_x[l],U_x[k]);
          nnz++;
        } else {
          dot[L_ind[l]].add_dot(L_x[l],U_x[k]);
        }
      } 
    }
    
    for(int l=A_p[i] ; l<A_p[i+1] ; l++) {
      if(w[A_ind[l]] < i) {
        w[A_ind[l]] = i;
        D_ind.push_back(A_ind[l]);
        dot[A_ind[l]].reset();
        dot[A_ind[l]].add_dot(-1,A_x[l]);
        nnz++;
      } else {
        dot[A_ind[l]].add_dot(-1,A_x[l]);
      }
    }   

    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);
    cdotprecision accu;

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       dot[D_ind[j]].result(accu);
       cinterval tmp;
       rnd(accu,tmp);
       D_x.push_back( complex(Sup(abs(Re(tmp))), Sup(abs(Im(tmp)))) );
    }

    D_p[i-start+1+1] = nnz;
  }
  
}

//! Compute the columns with indices from start to end of D=|L|*|U|
inline void absLUChunk(const srmatrix& L, const srmatrix& U, srmatrix& D, int start, int end) {
  int n = RowLen(L);

  D = srmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<real>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<real>& U_x  = U.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<real>& D_x  = D.values();

  std::vector<real> dot(n,real(0.0));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
          dot[L_ind[l]] = abs(L_x[l]) * abs(U_x[k]);
          nnz++;
        } else {
          dot[L_ind[l]] += abs(L_x[l]) * abs(U_x[k]);
        }
      } 
    }
    
    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       D_x.push_back(dot[D_ind[j]]);
    }

    D_p[i-start+1+1] = nnz;
  }
  
}

//! Compute the columns with indices from start to end of D=|L|*|U|
inline void absLUChunk(const scmatrix& L, const scmatrix& U, srmatrix& D, int start, int end) {
  int n = RowLen(L);

  D = srmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<complex>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<complex>& U_x  = U.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<real>& D_x  = D.values();

  std::vector<real> dot(n,real(0.0));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
          //dot[L_ind[l]] = abs(L_x[l]) * abs(U_x[k]);
          fp_mult(sqrt(Re(L_x[l])*Re(L_x[l])+Im(L_x[l])*Im(L_x[l])),
		  sqrt(Re(U_x[k])*Re(U_x[k])+Im(U_x[k])*Im(U_x[k])),
		  dot[L_ind[l]]);	  
          nnz++;
        } else {
          //dot[L_ind[l]] += abs(L_x[l]) * abs(U_x[k]);
          fp_multadd(sqrt(Re(L_x[l])*Re(L_x[l])+Im(L_x[l])*Im(L_x[l])),
		     sqrt(Re(U_x[k])*Re(U_x[k])+Im(U_x[k])*Im(U_x[k])),
		     dot[L_ind[l]]);	  
        }
      }
    }

    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       D_x.push_back(dot[D_ind[j]]);
    }

    D_p[i-start+1+1] = nnz;
  }
  
}


//! Compute the columns with indices from start to end of D=L*U
inline void LUChunk(const srmatrix& L, const srmatrix& U, srmatrix& D, int start, int end) {
  int n = RowLen(L);

  D = srmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<real>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<real>& U_x  = U.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<real>& D_x  = D.values();

  std::vector<real> dot(n,real(0.0));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
          dot[L_ind[l]] = L_x[l] * U_x[k];
          nnz++;
        } else {
          dot[L_ind[l]] += L_x[l] * U_x[k];
        }
      } 
    }
    
    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       D_x.push_back(dot[D_ind[j]]);
    }

    D_p[i-start+1+1] = nnz;
  }
  
}

//! Compute the columns with indices from start to end of D=L*U
inline void LUChunk(const scmatrix& L, const scmatrix& U, scmatrix& D, int start, int end) {
  int n = RowLen(L);

  D = scmatrix(n,end-start+1);

  const std::vector<int>& L_ind = L.row_indices();
  const std::vector<int>& L_p   = L.column_pointers();
  const std::vector<complex>& L_x  = L.values();
  const std::vector<int>& U_ind = U.row_indices();
  const std::vector<int>& U_p   = U.column_pointers();
  const std::vector<complex>& U_x  = U.values();
  std::vector<int>& D_ind = D.row_indices();
  std::vector<int>& D_p   = D.column_pointers();
  std::vector<complex>& D_x  = D.values();

  std::vector<complex> dot(n,complex(0.0));
  std::vector<int> w(n,-1);
  int nnz = 0;

  for(int i=start-1 ; i<end ; i++) {

    for(int k=U_p[i] ; k<U_p[i+1] ; k++) {
      for(int l=L_p[U_ind[k]] ; l<L_p[U_ind[k]+1] ; l++) {
        if(w[L_ind[l]] < i) {
          w[L_ind[l]] = i;
          D_ind.push_back(L_ind[l]);
          fp_mult(L_x[l],U_x[k],dot[L_ind[l]]);
          nnz++;
        } else {
          fp_multadd(L_x[l],U_x[k],dot[L_ind[l]]);  
        }
      }
    }

    sort(D_ind.begin()+D_p[i-start+1], D_ind.begin()+nnz);

    for(int j=D_p[i-start+1] ; j<nnz ; j++) {
       D_x.push_back(dot[D_ind[j]]);
    }

    D_p[i-start+1+1] = nnz;
  }
  
}

    
/*! Determines the size of the chunks for the matrix matrix products computed by
 *  the *chunk functions above. Uses a heuristic based on the number of threads
 *  and the distribution of non-zeros in U */
int determine_chunk_size(int n, int threads, const std::vector<int>& p) {
    if(threads == 1) return n;
    
    int nnz = p[n];
    int size = nnz / (10*threads);
    int chunk = n;
    int cs = 0;
    int tmp = 0;
    
    
    for(int i=0 ; i<n ; i++) {
        cs += p[i+1]-p[i];
        tmp++;
        if(cs >= size) {
            if(tmp < chunk)
                chunk = tmp;
            tmp = 0;
            cs = 0;
        }
    }
    
    if(chunk < 50)
        chunk = (50 < n) ? 50 : n;
    
    return chunk;
}


//! Main computing routine for computation of upper bound for C=|L||U|
template<class TF>
void absLUMain(const TF& L, const TF& U, srmatrix& C) {
   int n = RowLen(L);
   
   #ifdef _OPENMP
      //srmatrix Labs = abs(L);
      C = srmatrix(n,n);
      std::vector<int>& C_ind = C.row_indices();
      std::vector<int>& C_p   = C.column_pointers();
      std::vector<real>& C_x  = C.values();      
      
      const int chunk = determine_chunk_size(n,omp_get_max_threads(),U.column_pointers());
     
      std::vector<srmatrix> D_p((n/chunk)+1);
      std::vector<int> startpos;
      
      #pragma omp parallel default(shared)
      {

	fesetround(FE_UPWARD);
	
        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = Labs * abs(TF(U(1,n,i*chunk+1,(i+1)*chunk)));
          absLUChunk(L,U,D_p[i],i*chunk+1,(i+1)*chunk);
	}
	
        #pragma omp barrier

        #pragma omp single
        {
	  if(n%chunk != 0) {
            //D_p[n/chunk] = Labs * abs(TF(U(1,n,(n/chunk)*chunk+1,n)));
            absLUChunk(L,U,D_p[n/chunk],(n/chunk)*chunk+1,n);
	  }
	  
	  startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  C_ind = std::vector<int>(nnzsum);
	  C_x = std::vector<real>(nnzsum);
	  C_p[0] = 0;
	  C_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<real>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            C_ind[startpos[j]+i] = ind[i];
            C_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            C_p[j*chunk+i] = p[i] + startpos[j];
          }
        }

        fesetround(FE_TONEAREST);
      }

      D_p.clear();

    #else 

      fesetround(FE_UPWARD);
      C = abs(L)*abs(U);
      fesetround(FE_TONEAREST);

   #endif         
      
}

void absLU(const srmatrix& L, const srmatrix& U, srmatrix& C) {
  absLUMain(L,U,C);
}

void absLU(const scmatrix& L, const scmatrix& U, srmatrix& C) {
  absLUMain(L,U,C);
}

//! Main computing routine for computation of upper bound for C=LU-A
template<class TF, class TPoint>
void LUMinusAMain(const TF& L, const TF& U, const TF& A, TF& D_inf, TF& D_sup) {
    int n = RowLen(L);
    D_inf = TF(n,n);
    D_sup = TF(n,n);

    #ifdef _OPENMP
      std::vector<int>& Di_ind = D_inf.row_indices();
      std::vector<int>& Di_p   = D_inf.column_pointers();
      std::vector<TPoint>& Di_x  = D_inf.values();      
      std::vector<int>& Ds_ind = D_sup.row_indices();
      std::vector<int>& Ds_p   = D_sup.column_pointers();
      std::vector<TPoint>& Ds_x  = D_sup.values();      
      
      const int chunk = determine_chunk_size(n,omp_get_max_threads(),U.column_pointers());      
      std::vector<TF> D_p((n/chunk)+1);
      std::vector<int> startpos;      

      #pragma omp parallel default(shared)
      {

	fesetround(FE_UPWARD);

        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = L * (U(1,n,i*chunk+1,(i+1)*chunk)) - (A(1,n,i*chunk+1,(i+1)*chunk));
          LUMinusAChunk(L,U,A,D_p[i],i*chunk+1,(i+1)*chunk);
	  D_p[i].dropzeros();
	}

        #pragma omp barrier

        #pragma omp single
        {
	  if(n%chunk != 0) {
            //D_p[n/chunk] = L * (U(1,n,(n/chunk)*chunk+1,n)) - (A(1,n,(n/chunk)*chunk+1,n));
            LUMinusAChunk(L,U,A,D_p[n/chunk],(n/chunk)*chunk+1,n);
	    D_p[n/chunk].dropzeros();
	  }

	  
          startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  Ds_ind = std::vector<int>(nnzsum);
	  Ds_x = std::vector<TPoint>(nnzsum);
	  Ds_p[0] = 0;
	  Ds_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<TPoint>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            Ds_ind[startpos[j]+i] = ind[i];
            Ds_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            Ds_p[j*chunk+i] = p[i] + startpos[j];
          }
          
        }
        
        #pragma omp single
        {
          D_p.clear();
          D_p.resize((n/chunk)+1,TF(0,0,0));
          startpos.clear();
        }

        #pragma omp barrier

	fesetround(FE_DOWNWARD);
	
        #pragma omp for schedule(dynamic,1) 
        for(int i=0 ; i<n/chunk ; i++) {
            //D_p[i] = L * (U(1,n,i*chunk+1,(i+1)*chunk)) - (A(1,n,i*chunk+1,(i+1)*chunk));
            LUMinusAChunk(L,U,A,D_p[i],i*chunk+1,(i+1)*chunk);          
	    D_p[i].dropzeros();
	}
	
        #pragma omp barrier

        #pragma omp single
        {
          if(n%chunk != 0) {
	    //D_p[n/chunk] = L * (U(1,n,(n/chunk)*chunk+1,n)) - (A(1,n,(n/chunk)*chunk+1,n));
            LUMinusAChunk(L,U,A,D_p[n/chunk],(n/chunk)*chunk+1,n);
	    D_p[n/chunk].dropzeros();
	  }
	  
	  startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  Di_ind = std::vector<int>(nnzsum);
	  Di_x = std::vector<TPoint>(nnzsum);
	  Di_p[0] = 0;
	  Di_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<TPoint>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            Di_ind[startpos[j]+i] = ind[i];
            Di_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            Di_p[j*chunk+i] = p[i] + startpos[j];
          }
          
          D_p[j] = TF(0,0,0);
        }

        fesetround(FE_TONEAREST);
      }


    #else 

      fesetround(FE_DOWNWARD);
      LUMinusA(D_inf,L,U,A);
      D_inf.dropzeros();
      fesetround(FE_UPWARD);
      LUMinusA(D_sup,L,U,A);
      D_sup.dropzeros();
      fesetround(FE_TONEAREST);

   #endif  
}

template<typename TMat>
void LUMinusA(const TMat& L, const TMat& U, const TMat& A, TMat& C) {
   C = L*U - A;
}

//! Main computing routine for computation of upper bound for C=LU-[A]
template<class TF, class TI, class TPoint, class TInterval>
void LUMinusAIntervalMain(const TF& L, const TF& U, const TI& A, TI& D) {
    int n = RowLen(L);

    #ifdef _OPENMP
     std::vector<int>& D_ind = D.row_indices();
     std::vector<int>& DD_p   = D.column_pointers();
     std::vector<TInterval>& D_x = D.values();
      
      const int chunk = determine_chunk_size(n,omp_get_max_threads(),U.column_pointers());
      std::vector<TF> D_p((n/chunk)+1);
      std::vector<int> startpos;      

      #pragma omp parallel default(shared)
      {

	fesetround(FE_UPWARD);
	
        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = L * TF(U(1,n,i*chunk+1,(i+1)*chunk));
          LUChunk(L,U,D_p[i],i*chunk+1,(i+1)*chunk);
	}
	
        #pragma omp barrier

        #pragma omp single
        {
	  if(n%chunk != 0)  {
            //D_p[n/chunk] = L * TF(U(1,n,(n/chunk)*chunk+1,n));
            LUChunk(L,U,D_p[n/chunk],(n/chunk)*chunk+1,n);
	  }
	  
	  startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  D_ind = std::vector<int>(nnzsum);
	  D_x = std::vector<TInterval>(nnzsum);
	  DD_p[0] = 0;
	  DD_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<TPoint>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            D_ind[startpos[j]+i] = ind[i];
            D_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            DD_p[j*chunk+i] = p[i] + startpos[j];
          }
            
            D_p[j] = TF(0,0,0);
        }
        

          fesetround(FE_DOWNWARD);
	
        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = L * TF(U(1,n,i*chunk+1,(i+1)*chunk));
          LUChunk(L,U,D_p[i],i*chunk+1,(i+1)*chunk);
	}
	
        #pragma omp barrier

        #pragma omp single
        {
          if(n%chunk != 0) {
	    //D_p[n/chunk] = L * TF(U(1,n,(n/chunk)*chunk+1,n));
	    LUChunk(L,U,D_p[n/chunk],(n/chunk)*chunk+1,n);
	  }
	  
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<TPoint>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            D_x[startpos[j]+i]   |= x[i];
          }
            
          D_p[j] = TF(0,0,0);
        }

        fesetround(FE_TONEAREST);
      }

      D_p.clear();
      D -= A;
      D.dropzeros();

    #else 

      TF D_inf, D_sup;
      fesetround(FE_DOWNWARD);
      LUMinusA(D_inf,L,U,Sup(A));
      fesetround(FE_UPWARD);
      LUMinusA(D_sup,L,U,Inf(A));
      fesetround(FE_TONEAREST);

   #endif  
}

void LUMinusA(const srmatrix& L, const srmatrix& U, const srmatrix& A, srmatrix& D_inf, srmatrix& D_sup) {
  LUMinusAMain<srmatrix,real>(L,U,A,D_inf,D_sup);
}

void LUMinusA(const scmatrix& L, const scmatrix& U, const scmatrix& A, scmatrix& D_inf, scmatrix& D_sup) {
  LUMinusAMain<scmatrix,complex>(L,U,A,D_inf,D_sup);
}


void LUMinusA(const srmatrix& L, const srmatrix& U, const simatrix& A, simatrix& D) {
  LUMinusAIntervalMain<srmatrix,simatrix,real,interval>(L,U,A,D);
}

void LUMinusA(const scmatrix& L, const scmatrix& U, const scimatrix& A, scimatrix& D) {
  LUMinusAIntervalMain<scmatrix,scimatrix,complex,cinterval>(L,U,A,D);
}


//! Main computing routine for computation of upper bound for C=LU-A in higher precision
template<class TF, class TI, class TPoint, class TInterval>
void LUMinusAHighPrecMain(const TF& L, const TF& U, const TF& A, TF& C) {
   int n = RowLen(L);
   
   #ifdef _OPENMP
      C = TF(n,n);
      std::vector<int>& C_ind = C.row_indices();
      std::vector<int>& C_p   = C.column_pointers();
      std::vector<TPoint>& C_x  = C.values();      
      
      const int chunk = determine_chunk_size(n,omp_get_max_threads(),U.column_pointers());
      std::vector<TF> D_p((n/chunk)+1);
      std::vector<int> startpos;
      
      #pragma omp parallel default(shared)
      {
	
        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = Sup( L * TI(U(1,n,i*chunk+1,(i+1)*chunk)) - A(1,n,i*chunk+1,(i+1)*chunk) );
          LUMinusAHighPrecChunk(L,U,A,D_p[i],i*chunk+1,(i+1)*chunk);
	  D_p[i].dropzeros();
	}
	
        #pragma omp barrier

        #pragma omp single
        {
          if(n%chunk != 0) {
	    //D_p[n/chunk] = Sup( L * TI(U(1,n,(n/chunk)*chunk+1,n)) - A(1,n,(n/chunk)*chunk+1,n) );
	    LUMinusAHighPrecChunk(L,U,A,D_p[n/chunk],(n/chunk)*chunk+1,n);
	    D_p[n/chunk].dropzeros();
	  }
	  
	  startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  C_ind = std::vector<int>(nnzsum);
	  C_x = std::vector<TPoint>(nnzsum);
	  C_p[0] = 0;
	  C_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<TPoint>& x  = D_p[j].values();

	  
          for(unsigned int i=0 ; i<ind.size() ; i++) {
            C_ind[startpos[j]+i] = ind[i];
            C_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            C_p[j*chunk+i] = p[i] + startpos[j];
          }
        }
        
      }

      D_p.clear();

    #else 

      C = Sup(L * TI(U) - A);

   #endif         
   
}

void LUMinusAHighPrec(const srmatrix& L, const srmatrix& U, const srmatrix& A, srmatrix& C) {
  LUMinusAHighPrecMain<srmatrix,simatrix,real,interval>(L,U,A,C);
}

void LUMinusAHighPrec(const scmatrix& L, const scmatrix& U, const scmatrix& A, scmatrix& C) {
  LUMinusAHighPrecMain<scmatrix,scimatrix,complex,cinterval>(L,U,A,C);
}

//! Parallel computation of general sparse matrix product [C,D]=[A*B]
void spMatMulPar(const srmatrix& A, const srmatrix& B, srmatrix& C, srmatrix& D) {
   int n = RowLen(A);
   srmatrix D_inf(n,n);
   srmatrix D_sup(n,n);

    #ifdef _OPENMP
      std::vector<int>& Di_ind = D_inf.row_indices();
      std::vector<int>& Di_p   = D_inf.column_pointers();
      std::vector<real>& Di_x  = D_inf.values();      
      std::vector<int>& Ds_ind = D_sup.row_indices();
      std::vector<int>& Ds_p   = D_sup.column_pointers();
      std::vector<real>& Ds_x  = D_sup.values();      
      
      const int chunk = determine_chunk_size(n,omp_get_max_threads(),B.column_pointers());
      std::vector<srmatrix> D_p((n/chunk)+1);
      std::vector<int> startpos;      

      #pragma omp parallel default(shared)
      {

	fesetround(FE_DOWNWARD);
	
        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = A * srmatrix(B(1,n,i*chunk+1,(i+1)*chunk));
          LUChunk(A,B,D_p[i],i*chunk+1,(i+1)*chunk);
	}
	
        #pragma omp barrier

        #pragma omp single
        {
	  if(n%chunk != 0) {
            //D_p[n/chunk] = A * srmatrix(B(1,n,(n/chunk)*chunk+1,n));
            LUChunk(A,B,D_p[n/chunk],(n/chunk)*chunk+1,n);
	  }
	  
	  startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  Ds_ind = std::vector<int>(nnzsum);
	  Ds_x = std::vector<real>(nnzsum);
	  Ds_p[0] = 0;
	  Ds_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<real>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            Ds_ind[startpos[j]+i] = ind[i];
            Ds_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            Ds_p[j*chunk+i] = p[i] + startpos[j];
          }
        }

	fesetround(FE_UPWARD);
	
        #pragma omp for schedule(dynamic,1)
        for(int i=0 ; i<n/chunk ; i++) {
          //D_p[i] = A * srmatrix(B(1,n,i*chunk+1,(i+1)*chunk));
          LUChunk(A,B,D_p[i],i*chunk+1,(i+1)*chunk);
	}
	
        #pragma omp barrier

        #pragma omp single
        {
	  if(n%chunk != 0) {
            //D_p[n/chunk] = A * srmatrix(B(1,n,(n/chunk)*chunk+1,n));
            LUChunk(A,B,D_p[n/chunk],(n/chunk)*chunk+1,n);
	  }
	  
	  startpos = std::vector<int>(D_p.size());
	  int nnzsum=0;
	  for(unsigned int i=0 ; i<D_p.size() ; i++) {
	    startpos[i] = nnzsum;
	    nnzsum += D_p[i].get_nnz();
	  }
	  
	  Di_ind = std::vector<int>(nnzsum);
	  Di_x = std::vector<real>(nnzsum);
	  Di_p[0] = 0;
	  Di_p[n] = nnzsum;
	}

        #pragma omp barrier

        #pragma omp for schedule(dynamic)
        for(unsigned int j=0 ; j<D_p.size() ; j++) {
          std::vector<int>& ind = D_p[j].row_indices();
          std::vector<int>& p   = D_p[j].column_pointers();
          std::vector<real>& x  = D_p[j].values();

          for(unsigned int i=0 ; i<ind.size() ; i++) {
            Di_ind[startpos[j]+i] = ind[i];
            Di_x[startpos[j]+i]   = x[i];
          }

          for(unsigned int i=1 ; i<p.size() ; i++) {
            Di_p[j*chunk+i] = p[i] + startpos[j];
          }
        }

        #pragma omp single
        {
	  D = 0.5 * (D_sup - D_inf);
	  Resize(D_sup);
	  C = D_inf + D;
	}
	
        fesetround(FE_TONEAREST);
      }


    #else 

      simatrix T = A * simatrix(B);
      D = 0.5 * diam(T);
      C = mid(T);
      
   #endif    
}

//! Parallel computation of general sparse matrix product [C,D]=[A*B]
void spMatMulPar(const scmatrix& A, const scmatrix& B, scmatrix& C, scmatrix& D) {
  srmatrix tmp1,tmp2;
  
/*  vector<int>& p_t1 = tmp1.column_pointers();
  vector<int>& ind_t1 = tmp1.row_indices();*/
  std::vector<real>& x_t1 = tmp1.values();
/*  vector<int>& p_t2 = tmp2.column_pointers();
  vector<int>& ind_t2 = tmp2.row_indices();*/
  std::vector<real>& x_t2 = tmp2.values();
/*  vector<int>& p_C = C.column_pointers();
  vector<int>& ind_C = C.row_indices();
  vector<complex>& x_C = C.values();
  vector<int>& p_D = D.column_pointers();
  vector<int>& ind_D = D.row_indices();
  vector<complex>& x_D = D.values();*/
  
  
  
  spMatMulPar(Re(A),Re(B),tmp1,tmp2);
  C = tmp1;
  D = tmp2;
  
  spMatMulPar(-Im(A),Im(B),tmp1,tmp2);
  fesetround(FE_UPWARD);
  C += tmp1;
  D += tmp2;
  fesetround(FE_TONEAREST);  

  spMatMulPar(Re(A),Im(B),tmp1,tmp2);
  scmatrix tmp(tmp1);
  std::vector<complex>& x_tmp = tmp.values();
  for(int i=0 ; i<x_tmp.size() ; i++)
    x_tmp[i] = complex(0,x_t1[i]);
  C += tmp;
  tmp = tmp2;
  for(int i=0 ; i<x_tmp.size() ; i++)
    x_tmp[i] = complex(0,x_t2[i]);
  D += tmp;
  //C += complex(0,1)*tmp1;
  //D += complex(0,1)*tmp2;
  
  spMatMulPar(Im(A),Re(B),tmp1,tmp2);
  fesetround(FE_UPWARD);
  tmp = tmp1;
  for(int i=0 ; i<x_tmp.size() ; i++)
    x_tmp[i] =  complex(0,Re(x_tmp[i]));
  C += tmp;
  tmp = tmp2;
  for(int i=0 ; i<x_tmp.size() ; i++)
    x_tmp[i] =  complex(0,Re(x_tmp[i]));
  D += tmp;  
  //C += complex(0,1)*tmp1;
  //D += complex(0,1)*tmp2;
  fesetround(FE_TONEAREST);  
}

} //namespace cxsc
