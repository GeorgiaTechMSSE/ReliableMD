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

/* Implementation of the sparse linear system solver slss based on a verified normwise
 * error bound of an approximate solution. */

#include <sparselinsys.hpp>
#include "cxsc_cholmod.hpp"
#include "cxsc_umfpack.hpp"
#include "utility.hpp"
#include "trisolve.hpp"
#include <intvector.hpp>
#include <scimatrix.hpp>
#include <fenv.h>
#include <algorithm>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "timer.hpp"

using namespace cxsc;
using namespace std;


namespace cxsc {
  
//! Error codes used in this module.
static const int
  NoError      = 0,   // No error occurred
  NotSquare    = 1,   // System to be solved is not square
  DimensionErr = 2,   // Dimensions of A and b are not compatible
  DecompFailed = 3,   // System is probably singular
  VerivFailed  = 4,   // Verification failed, system is probably
                      // ill-conditioned
  NotPosDef    = 5;   // Matrix not positive definite
  

//! Compute  D = max(|D_inf|,|D_sup|)
void AbsMax(srmatrix& D, const srmatrix& D_inf, const srmatrix& D_sup) {
  srmatrix E,F;
  int n = RowLen(D_inf);

  D = srmatrix(n,n);

  E = abs(D_inf);
  F = abs(D_sup);

  vector<int>& ind_D   = D.row_indices();
  vector<int>& p_D     = D.column_pointers();
  vector<real>& val_D  = D.values();

  vector<int>& ind_E   = E.row_indices();
  vector<int>& p_E     = E.column_pointers();
  vector<real>& val_E  = E.values();

  vector<int>& ind_F   = F.row_indices();
  vector<int>& p_F     = F.column_pointers();
  vector<real>& val_F  = F.values();

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {
    int k = p_E[j];
    int l = p_F[j];

    while(k < p_E[j+1]  &&  l < p_F[j+1]) {
        if(ind_E[k] < ind_F[l]) {
          ind_D.push_back(ind_E[k]);
          val_D.push_back(val_E[k]);
          k++;
        } else if(ind_E[k] == ind_F[l]) {
          ind_D.push_back(ind_E[k]);
          val_D.push_back((val_E[k] > val_F[l]) ? val_E[k] : val_F[l]);
          k++; l++;
        } else {
          ind_D.push_back(ind_F[l]);
          val_D.push_back(val_F[l]);
          l++;
        }
        nnz++;
    }

    while(k < p_E[j+1]) {
      ind_D.push_back(ind_E[k]);
      val_D.push_back(val_E[k]);
      k++;
      nnz++;
    }

    while(l < p_F[j+1]) {
      ind_D.push_back(ind_F[l]);
      val_D.push_back(val_F[l]);
      l++;
      nnz++;
    }

    p_D[j+1] = nnz;

  }
  
}


//! Compute  D = max(|D_inf|,|D_sup|)
void AbsMax(scmatrix& D, const scmatrix& D_inf, const scmatrix& D_sup) {
  srmatrix tmp;
  AbsMax(tmp, Re(D_inf), Re(D_sup));
  D = tmp;
  AbsMax(tmp, Im(D_inf), Im(D_sup));
  D += complex(0,1)*tmp;  
}

inline real Re(const real& r) {
  return r;
}

/*! This function computes a verified lower bound for the least singular value
    of the matrix Aorig, based on an approximation lambda. */
template<typename TMat>
bool boundSingMin(const TMat& Aorig, real& lambda, real& tau, slssnconfig& cfg) {
  int n = RowLen(Aorig);
  TMat A(Aorig), L(n,n), D(n,n);
  intvector q(n);
  real factor = 0.9;
  bool success = false;
  int err = NoError;

  while(factor > 0.0 && !success) {
    if(cfg.msg) cout << " Trying factor=" << factor << endl;
    real lambda_fac = factor * lambda ;
    A = Aorig - lambda_fac*Id(A);   
    if(cfg.msg) { cout << " Cholesky decomposition of A-factor*lambda*I..."; cout.flush(); }
    
    chol(A,L,q,err);
    if(err == 0) {
      success = true;
      lambda = lambda_fac;
      if(cfg.msg) cout << endl;
    } else {
      if(cfg.msg)  cout << "failed" << endl;
      success = false;
      factor -= 0.1;
    }
  }

  if(!success) {
    return false;
  }

  if(cfg.msg) { cout << " Trying to bound decomposition error, method 1..."; cout.flush(); }
    
  intvector t(n);
  vector<int>& p = L.column_pointers();
  for(unsigned int i=0 ; i<p.size()-1 ; i++) 
    t[i+1] = p[i+1] - p[i] - 1;
    
  A = A(q,q);
  rvector d(n);
  real Eta = power(2,-1074);
  real maxdiag = 0.0, da;
  for(int i=1 ; i<=n ; i++) {
    da = Re(A(i,i)); //abs to treat complex case (since A is hermitian, this is in principle a cast to real)
    if(da > maxdiag) maxdiag = da;
    real gamma = ((t[i]+2)*Epsilon / (1-(t[i]+2)*Epsilon)) / (1-2*Epsilon);
    d[i] = sqrt( (gamma/(1-gamma)) * da ) / (1-4*Epsilon);
  }
    
  real M = ( (3*(2*n+maxdiag)) / (1-3*Epsilon) );
    
  tau = ( ((d*d)/(1-(2*n-1)*Epsilon)) + ((n*M*Eta)/(1-Epsilon)) ) / (1-Epsilon);
  
  if(tau >= lambda) {
    if(cfg.msg) {
      cout << "failed" << endl;
      if(cfg.msg) { cout << " Trying to bound decomposition error, method 2..."; cout.flush(); }
    }
    TMat U = transpherm(L);
    TMat D_inf(n,n), D_sup(n,n);
       
    opdotprec = cfg.precDecompError;
    if(opdotprec == 1) {
      LUMinusA(L,U,A,D_inf,D_sup);
      AbsMax(D,D_inf,D_sup); 
    } else {
      //D = Sup(L * simatrix(U) - A(q,q));           
      LUMinusAHighPrec(L,U,A,D);
      opdotprec = 1;
    }

    real n1 = Norm1(D);
    real n0 = Norm00(D);
    tau = (sqrt(n1) * sqrt(n0)) / (1 - 3*Epsilon); 
   
  }

  return tau < lambda;
}


//! Main function of the solver for spd systems.
template<typename TMat, typename TRhs, typename TSolution, typename TPointVec, typename TDot>
void SparseLinSolveMain(TMat& Ao, TRhs& bo, TSolution& xx, real& lambda, real& tau, srmatrix& S, TRhs& xappstag, TMat* rA, int& err, slssnconfig cfg) {
  feclearexcept(FE_UNDERFLOW);
  int n = RowLen(Ao);
  
  if(n != ColLen(Ao)) {
    err = NotSquare;
    return;
  }
  
  if(ColLen(Ao) != ColLen(bo)) {
    err = DimensionErr;
    return;
  }  
  
  if(Ao != transpherm(Ao)) {
    err = NotPosDef;
    return;
  }

  TMat       A, L, U;
  TRhs       b(n);
  TPointVec  xapp(n), x(n), y(n), h(n);
  TSolution  def(n);
  real       delta;
  intvector  p(n);
  int        oldprec = opdotprec;
  double     start, end;

  opdotprec = 1;
  err = NoError;
  
  #ifdef _OPENMP
  if(cfg.threads > 0) {
     omp_set_num_threads(cfg.threads);
  } else {
     cfg.threads = omp_get_max_threads(); 
     omp_set_num_threads(cfg.threads);
  }
  if(cfg.msg) std::cout << "Using " << cfg.threads << " thread(s) for computations" << std::endl;    
  #endif  
  
  start = GetTime();
  S = Id(srmatrix(n,n));
  if(cfg.scaleMatrix) {
    if(cfg.msg) cout << "Scaling matrix..." << endl;
    vector<real>& S_val = S.values();
//     #ifdef _OPENMP
//     #pragma omp parallel for
//     #endif
    for(int i=1 ; i<=n && !err ; i++){
      if(!checkPosDef(Ao(i,i))) err = NotPosDef;
      S_val[i-1] = ::pow( 2.0, static_cast<int>(_double(-0.5*q_log2(abs(Ao(i,i))))) );
    }
    
    if(err != NoError) return;
  
    A = S*Ao*S;
    b = S*bo;

    //Check for underflow during scaling
    if(fetestexcept(FE_UNDERFLOW) & FE_UNDERFLOW) {
      if(cfg.msg) cout << " Underflow during scaling, no scaling used" << endl;
      S = Id(srmatrix(n,n));
      A = Ao;
      b = bo;
      feclearexcept(FE_UNDERFLOW);
    }

    if(rA != NULL) {
      fesetround(FE_UPWARD);
      *rA = S * (*rA) * S;
      fesetround(FE_TONEAREST);
    }
    
  } else {
    A = Ao;
    b = bo;
  }
  end = GetTime();
  if(cfg.msg) cout << "Time for scaling: " << end-start << endl;

  start = GetTime();
  if(cfg.msg) cout << "Cholesky decomposition..." << endl;
  chol(A, L, p, err);
  if(err != 0) {
    if(cfg.msg) cout << " Cholesky decomposition failed, aborting" << endl;
    err = DecompFailed;
    opdotprec = oldprec;
    return;
  }
  U = transpherm(L);
  end = GetTime();
  if(cfg.msg) cout << "Time for Cholesky decomposition: " << end-start << endl;
  
  //Apply permutation
  A = A(p,p);
  b = b(p);

  
  start = GetTime();
  if(cfg.msg) cout << "Computing approximation for smallest singular value..." << endl;
  y = 1;
  y = y / Norm2(y); 
  real lambdaold = 0;
  lambda = 1;
  for( int k=1 ; k<=cfg.maxIterPower && (abs(lambda-lambdaold)>1e-6) ; k++ ) {
    x = y / Norm2(y);
    lowtrisolve(L, x, h);
    uptrisolve(U, h, y);
    lambdaold = lambda;
    lambda = Norm2(x) / Norm2(y) ;
  } 
  end = GetTime();
  if(cfg.msg) {
    cout << "Computed approximation: lambda=" << lambda << endl;
    cout << "Time for approximation: " << end-start << endl;
  }

  
  start = GetTime();
  if(cfg.msg) cout << "Computing approximate solution..." << endl;
  int prec = cfg.apprxSolutionStagPrec; 

  int rhs = RowLen(b);
  xappstag = TRhs(n,prec*rhs);
  
  for(int r=1 ; r<=rhs ; r++) {
    h = b[Col(r)];      
  
    for(int p=1 ; p<prec ; p++) {
      lowtrisolve(L,h,y);
      uptrisolve(U,y,xapp);    
      xappstag[Col(p+(r-1)*prec)] = xapp;
      residual(h,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),p,cfg.highPrecApprxSolution);
    }
    
    lowtrisolve(L,h,y);
    uptrisolve(U,y,xapp);
    xappstag[Col(prec+(r-1)*prec)] = xapp;
  }
  end = GetTime();
  if(cfg.msg) cout << "Time for approximate solution: " << end-start << endl;

  start = GetTime();
  if(cfg.msg) cout << "Computing lower bound for least singular value..." << endl;
  if( !boundSingMin(A,lambda,tau,cfg) ) {
    if(cfg.msg) {
      cout << "failed" << endl;
      cout << "Could not verify lower bound for least singular value, aborting" << endl;
    }
    err = VerivFailed;
    opdotprec = oldprec;    
    return;
  }
  end = GetTime();
  if(cfg.msg) cout << endl << "Time for computing lower bound: " << end-start << endl;
    
  TDot dot;

  #ifdef _OPENMP
  #pragma omp parallel for private(dot)
  #endif
  for(int i=1 ; i<=n ; i++) {
    for(int r=1 ; r<=rhs ; r++) {
      dot = 0.0;
      for(int j=1 ; j<=prec ; j++)
        dot += xappstag[i][j+(r-1)*prec];
      rnd(dot,xx[i][r]);
    }
  }
  
  xx = xx(perminv(p));
  xappstag = xappstag(perminv(p));
    
  opdotprec = oldprec;

  return;
} 

//! Main function of the solver for general systems.
template<typename TMat, typename TRhs, typename TSolution, typename TPointVec, typename TIVec, typename TDot>
void SparseLinSolveGeneralMain(TMat& Ao, TRhs& bo, TSolution& xx, real& lambda, real& tau, srmatrix& T, TRhs& xappstag, TMat* rA, int& err, slssnconfig cfg) {
  int n = RowLen(Ao);
  
  if(n != ColLen(Ao)) {
    err = NotSquare;
    return;
  }
  
  if(ColLen(Ao) != ColLen(bo)) {
    err = DimensionErr;
    return;
  }

  TMat       L, U, LT, UT, AT, R;
  TPointVec  xapp(n), x(n), y(n), h(n);
  real       delta;
  intvector  p,q,r;
  TMat       A;
  TRhs       b;
  int        oldprec = opdotprec;
  double     start, end;

  opdotprec = 1;

  feclearexcept(FE_UNDERFLOW);
  
  #ifdef _OPENMP
  if(cfg.threads > 0) {
     omp_set_num_threads(cfg.threads);
  } else {
     cfg.threads = omp_get_max_threads(); 
     omp_set_num_threads(cfg.threads);
  }
  if(cfg.msg) std::cout << "Using " << cfg.threads << " thread(s) for computations" << std::endl;    
  #endif  

  start = GetTime();
  srmatrix S(Id(srmatrix(n,n)));
  T = S;
  if(cfg.scaleMatrix) {
    if(cfg.msg) cout << "Scaling matrix..." << endl;
    //snbin(Ao,10,S,T);
    normbin(Ao,S,T);
    A = S*Ao*T;
    b = S*bo;

    //Check for underflow during scaling
    if(fetestexcept(FE_UNDERFLOW) & FE_UNDERFLOW) {
      if(cfg.msg) cout << "Underflow during scaling, no scaling used" << endl;
      S = T = Id(srmatrix(n,n));
      A = Ao;
      b = bo;
      feclearexcept(FE_UNDERFLOW);
    }

    if(rA != NULL) {
      fesetround(FE_UPWARD);
      *rA = S * (*rA) * T;
      fesetround(FE_TONEAREST);
    }
    
  } else {
    A = Ao;
    b = bo;
  }
  end = GetTime();
  if(cfg.msg) cout << "Time for scaling: " << end-start << endl;

  start = GetTime();
  if(cfg.msg) cout << "LU decomposition..." << endl;
  
  lu_decomp(A,L,U,p,q,err);
  if(err) {
    if(cfg.msg) cout << "Error during LU decomposition, aborting..." << endl;
    err = DecompFailed;
    opdotprec = oldprec;
    return;
  }
  end = GetTime();
  if(cfg.msg) cout << "Time for decomposition: " << end-start << endl;
  
  //Apply permutation
  A = A(p,q);
  b = b(p);
  if(rA != NULL) {
    *rA = (*rA)(p,q);
  }  
  
  opdotprec = 1;

  start = GetTime();
  if(cfg.msg) cout << "Computing AT*A..." << endl;
  AT = transpherm(A);
  TMat AT_A, AT_A_D;
  spMatMulPar(AT,A,AT_A,AT_A_D);
  AT = TMat();
  end = GetTime();
  if(cfg.msg) cout << "Time for computation: " << end-start << endl;

  feclearexcept(FE_UNDERFLOW);

  start = GetTime();
  if(cfg.msg) cout << "Scaling AT*A..." << endl;
  srmatrix S2(Id(srmatrix(n,n)));
  vector<real>& S2_val = S2.values();
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(int i=1 ; i<=n ; i++){
    S2_val[i-1] = ::pow(2,static_cast<int>(_double(-0.5*q_log2(abs(AT_A(i,i))))) );
  }
  R = AT_A;
  AT_A = S2*AT_A*S2;
  AT_A_D = S2*AT_A_D*S2;

  //Check for underflow during scaling
  if(fetestexcept(FE_UNDERFLOW) & FE_UNDERFLOW) {
    if(cfg.msg) cout << "Underflow during scaling, no scaling used" << endl;
    S2 = Id(srmatrix(n,n));
    AT_A = R;
    feclearexcept(FE_UNDERFLOW);
  }

  if(rA != NULL) {
    fesetround(FE_UPWARD);
    *rA = S2 * (*rA);
    fesetround(FE_TONEAREST);
  }

  //Free memory
  R = TMat();
  
  end = GetTime();
  if(cfg.msg) cout << "Time for scaling: " << end-start << endl;

  start = GetTime();
  if(cfg.msg) cout << "Computing approximation for smallest singular value of AT*A..." << endl;
  LT = transpherm(L);
  UT = S2 * transpherm(U);
  TMat US  = U * S2;
  y = 1;
  y = y / Norm2(y); 
  real lambdaold = 0;
  lambda = 1;
  for( int k=1 ; k<=cfg.maxIterPower && (abs(lambdaold-lambda)>1e-6) ; k++ ) {
    x = y / Norm2(y);
    lowtrisolve(UT,x,h);
    uptrisolve(LT,h,y);
    lowtrisolve(L,y,h);
    uptrisolve(US,h,y);
    lambdaold = lambda;
    lambda = Norm2(x) / Norm2(y) ;
  } 
  LT = TMat();
  UT = TMat(); 
  US = TMat();
  if(cfg.msg) cout << "Approximate least singular value of scaled A: " << sqrt(lambda) << endl;
  end = GetTime();
  if(cfg.msg) cout << "Time for approximation: " << end-start << endl;

  feclearexcept(FE_UNDERFLOW);
  
  
  start = GetTime();
  if(cfg.msg) cout << "Computing approximate solution..." << endl;
  int prec = cfg.apprxSolutionStagPrec; 

  int rhs = RowLen(b);
  xappstag = TRhs(n,prec*rhs);
  TIVec def(n);
  A = S2*A;
  b = S2*b;
  L = S2*L;
  
  for(int r=1 ; r<=rhs ; r++) {
    h = b[Col(r)];
  
    for(int p=1 ; p<prec ; p++) {
      lowtrisolve(L,h,y);
      uptrisolve(U,y,xapp);
      xappstag[Col(p+(r-1)*prec)] = xapp;
      residual(h,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),p,cfg.highPrecApprxSolution);
    }
    
    lowtrisolve(L,h,y);
    uptrisolve(U,y,xapp);
    xappstag[Col(prec+(r-1)*prec)] = xapp;
  }
  
  L = TMat();
  U = TMat();
  
  end = GetTime();
  if(cfg.msg) cout << "Time for approximate solution: " << end-start << endl;

  start = GetTime();
  if(cfg.msg) cout << "Computing lower bound for least singular value..." << endl;  
  if( !boundSingMin(AT_A,lambda,tau,cfg) ) {
    if(cfg.msg) {
      cout << "failed" << endl;
      cout << "Could not verify lower bound for least singular value, aborting" << endl;
    }
    err = VerivFailed;
    opdotprec = oldprec;    
    return;    
  } else {
    real nrmAAT = sqrt(Norm1(AT_A_D) * Norm00(AT_A_D)) / (1-2*Epsilon);
    tau = (tau + nrmAAT) / (1-Epsilon);
    if(tau >= lambda) {
      if(cfg.msg) {
        cout << "failed" << endl;
        cout << "Could not verify lower bound for least singular value, aborting" << endl;
      }
      err = VerivFailed;
      opdotprec = oldprec;    
      return;        
    }
  }
  lambda = sqrt(lambda) / (1-Epsilon) ;  
  end = GetTime();
  if(cfg.msg) cout << endl << "Time for computing bound: " << end-start << endl;

  y = 1;
  TDot dot;
  #ifdef _OPENMP
  #pragma omp parallel for private(dot)
  #endif
  for(int i=1 ; i<=n ; i++) {
    for(int r=1 ; r<=rhs ; r++) {
      dot = 0.0;
      for(int j=1 ; j<=prec ; j++)
        dot += xappstag[i][j+(r-1)*prec];
      rnd(dot,xx[i][r]);
    }
  }
  
  xx = xx(perminv(q));
  xappstag = xappstag(perminv(q));

  opdotprec = oldprec;

  return;  
}

//!Start function: Starts main function and computes final enclosure
void slssn_start(srmatrix& A, rmatrix& b, imatrix& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  int oldprec = opdotprec;
  rmatrix xappstag;
  srmatrix S;
  real lambda, tau;
  
  if(cfg.spd)
    SparseLinSolveMain<srmatrix,rmatrix,imatrix,rvector,dotprecision>(A,b,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  else
    SparseLinSolveGeneralMain<srmatrix,rmatrix,imatrix,rvector,ivector,dotprecision>(A,b,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  
  if(err!=0) {
    return;
  }

  int rhs = RowLen(b);
  int prec = RowLen(xappstag) / rhs;
  ivector d(n);
  imatrix def(n,rhs);
  rvector y(n);
  xappstag = S*xappstag;

  for(int r=1 ; r<=rhs ; r++) {
    residual(d,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),prec,cfg.highPrecApprxSolution);    
    def[Col(r)] = d;
  }
  
  rvector nrm(rhs);
  rvector delta(rhs);
 
  for(int r=1 ; r<=rhs ; r++) {
    nrm[r] = Norm2(absmax(def[Col(r)]));
    delta[r] = ( ( 1.0 / (lambda - tau) ) * nrm[r] ) / (1-3*Epsilon);
  }

  if(cfg.msg) {
    cout << "Verified lower bound for least singular value: " << subd(lambda, tau) << endl;
    if(rhs == 1) {
      cout << "Norm(b-A*xapp): " << nrm;
      cout << "Delta: " << delta;
    }
    cout << "Computing solution enclosure..." << endl;
  } 

  y = 1.0;
  
  for(int r=1 ; r<=rhs ; r++)
    xx[Col(r)] = S * (xx[Col(r)] + interval(-delta[r],delta[r]) * y);

  opdotprec = oldprec;
}

//!Start function: Starts main function and computes final enclosure
void slssn_start(srmatrix& A, imatrix& b, imatrix& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  int oldprec = opdotprec;
  rmatrix xappstag;
  srmatrix S;
  real lambda, tau;
  
  rmatrix mb = mid(b);
  
  if(cfg.spd)
    SparseLinSolveMain<srmatrix,rmatrix,imatrix,rvector,dotprecision>(A,mb,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  else
    SparseLinSolveGeneralMain<srmatrix,rmatrix,imatrix,rvector,ivector,dotprecision>(A,mb,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  
  if(err!=0) {
    return;
  }

  int rhs = RowLen(b);
  int prec = RowLen(xappstag) / rhs;
  ivector d(n);
  imatrix def(n,rhs);
  rvector y(n);
  xappstag = S*xappstag;

  for(int r=1 ; r<=rhs ; r++) {
    residual(d,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),prec,cfg.highPrecApprxSolution);    
    def[Col(r)] = d;
  }
  
  rvector nrm(rhs);
  rvector delta(rhs);
 
  for(int r=1 ; r<=rhs ; r++) {
    nrm[r] = Norm2(absmax(def[Col(r)]));
    delta[r] = ( ( 1.0 / (lambda - tau) ) * nrm[r] ) / (1-3*Epsilon);
  }

  if(cfg.msg) {
    cout << "Verified lower bound for least singular value: " << subd(lambda, tau) << endl;
    if(rhs == 1) {
      cout << "Norm(b-A*xapp): " << nrm;
      cout << "Delta: " << delta;
    }
    cout << "Computing solution enclosure..." << endl;
  } 

  y = 1.0;
  
  for(int r=1 ; r<=rhs ; r++)
    xx[Col(r)] = S * (xx[Col(r)] + interval(-delta[r],delta[r]) * y);

  opdotprec = oldprec;
}


//!Start function: Starts main function and computes final enclosure
void slssn_start(simatrix& A, imatrix& b, imatrix& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  int oldprec = opdotprec;
  rmatrix xappstag;
  srmatrix S;
  real lambda, tau;
  
  srmatrix mA = mid(A);
  rmatrix mb = mid(b);
  srmatrix rA = 0.5*diam(A);
  
  if(cfg.spd) {
    SparseLinSolveMain<srmatrix,rmatrix,imatrix,rvector,dotprecision>(mA,mb,xx,lambda,tau,S,xappstag,&rA,err,cfg);    
  } else {
    SparseLinSolveGeneralMain<srmatrix,rmatrix,imatrix,rvector,ivector,dotprecision>(mA,mb,xx,lambda,tau,S,xappstag,&rA,err,cfg);
  }
  
  tau = (tau + sqrt(Norm1(rA) * Norm00(rA))) / (1-3*Epsilon);
  
  if(err!=0 || tau >= lambda) {
    if(cfg.msg) cout << "Could not verify lower bound for least singular value, aborting" << endl;
    return;
  }
  
  int rhs = RowLen(b);
  int prec = RowLen(xappstag) / rhs;
  ivector d(n);
  imatrix def(n,rhs);
  rvector y(n);
  xappstag = S*xappstag;
  
  for(int r=1 ; r<=rhs ; r++) {
    residual(d,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),prec,cfg.highPrecApprxSolution);    
    def[Col(r)] = d;
  }
  
  rvector nrm(rhs);
  rvector delta(rhs);
 
  for(int r=1 ; r<=rhs ; r++) {
    nrm[r] = Norm2(absmax(def[Col(r)]));
    delta[r] = ( ( 1.0 / (lambda - tau) ) * nrm[r] ) / (1-3*Epsilon);
  }
  
  if(cfg.msg) {
    cout << "Verified lower bound for least singular value: " << subd(lambda, tau) << endl;
    if(rhs == 1) {
      cout << "Norm(b-A*xapp): " << nrm;
      cout << "Delta: " << delta;
    }
    cout << "Computing solution enclosure..." << endl;
  } 

  y = 1.0;
  
  for(int r=1 ; r<=rhs ; r++)
    xx[Col(r)] = S * (xx[Col(r)] + interval(-delta[r],delta[r]) * y);

  opdotprec = oldprec;
  
}


//!Start function: Starts main function and computes final enclosure
void slssn_start(scmatrix& A, cmatrix& b, cimatrix& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  int oldprec = opdotprec;
  cmatrix xappstag;
  srmatrix S;
  real lambda, tau;
  
  if(cfg.spd)
    SparseLinSolveMain<scmatrix,cmatrix,cimatrix,cvector,cdotprecision>(A,b,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  else
    SparseLinSolveGeneralMain<scmatrix,cmatrix,cimatrix,cvector,civector,cdotprecision>(A,b,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  
  if(err!=0) {
    return;
  }

  int rhs = RowLen(b);
  int prec = RowLen(xappstag) / rhs;
  civector d(n);
  cimatrix def(n,rhs);
  cvector y(n);
  xappstag = S*xappstag;

  for(int r=1 ; r<=rhs ; r++) {
    residual(d,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),prec,cfg.highPrecApprxSolution);    
    def[Col(r)] = d;
  }

  rvector nrm(rhs);
  rvector delta(rhs);
 
  for(int r=1 ; r<=rhs ; r++) {
    nrm[r] = Norm2(Sup(abs(def[Col(r)])));
    delta[r] = ( ( 1.0 / (lambda - tau) ) * nrm[r] ) / (1-3*Epsilon);
  }

  if(cfg.msg) {
    cout << "Verified lower bound for least singular value: " << subd(lambda, tau) << endl;
    if(rhs == 1) {
      cout << "Norm(b-A*xapp): " << nrm;
      cout << "Delta: " << delta;
    }
    cout << "Computing solution enclosure..." << endl;
  } 

  y = complex(1.0,1.0);
  
  for(int r=1 ; r<=rhs ; r++)
    xx[Col(r)] = S * (civector(xx[Col(r)]) + interval(-delta[r],delta[r]) * scvector(y));

  opdotprec = oldprec;
}

//!Start function: Starts main function and computes final enclosure
void slssn_start(scmatrix& A, cimatrix& b, cimatrix& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  int oldprec = opdotprec;
  cmatrix xappstag;
  srmatrix S;
  real lambda, tau;
  
  cmatrix mb = mid(b);
  
  if(cfg.spd)
    SparseLinSolveMain<scmatrix,cmatrix,cimatrix,cvector,cdotprecision>(A,mb,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  else
    SparseLinSolveGeneralMain<scmatrix,cmatrix,cimatrix,cvector,civector,cdotprecision>(A,mb,xx,lambda,tau,S,xappstag,NULL,err,cfg);
  
  if(err!=0) {
    return;
  }

  int rhs = RowLen(b);
  int prec = RowLen(xappstag) / rhs;
  civector d(n);
  cimatrix def(n,rhs);
  cvector y(n);
  xappstag = S*xappstag;

  for(int r=1 ; r<=rhs ; r++) {
    residual(d,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),prec,cfg.highPrecApprxSolution);    
    def[Col(r)] = d;
  }

  rvector nrm(rhs);
  rvector delta(rhs);
 
  for(int r=1 ; r<=rhs ; r++) {
    nrm[r] = Norm2(Sup(abs(def[Col(r)])));
    delta[r] = ( ( 1.0 / (lambda - tau) ) * nrm[r] ) / (1-3*Epsilon);
  }

  if(cfg.msg) {
    cout << "Verified lower bound for least singular value: " << subd(lambda, tau) << endl;
    if(rhs == 1) {
      cout << "Norm(b-A*xapp): " << nrm;
      cout << "Delta: " << delta;
    }
    cout << "Computing solution enclosure..." << endl;
  } 

  y = complex(1.0,1.0);
  
  for(int r=1 ; r<=rhs ; r++)
    xx[Col(r)] = S * (civector(xx[Col(r)]) + interval(-delta[r],delta[r]) * scvector(y));

  opdotprec = oldprec;
}

//!Start function: Starts main function and computes final enclosure
void slssn_start(scimatrix& A, cimatrix& b, cimatrix& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  int oldprec = opdotprec;
  cmatrix xappstag;
  srmatrix S;
  real lambda, tau;
  
  scmatrix mA = mid(A);
  cmatrix mb = mid(b);
  scmatrix rA = 0.5*diam(A);
  
  if(cfg.spd) {
    SparseLinSolveMain<scmatrix,cmatrix,cimatrix,cvector,cdotprecision>(mA,mb,xx,lambda,tau,S,xappstag,&rA,err,cfg);    
  } else {
    SparseLinSolveGeneralMain<scmatrix,cmatrix,cimatrix,cvector,civector,cdotprecision>(mA,mb,xx,lambda,tau,S,xappstag,&rA,err,cfg);
  }
  
  tau = (tau + sqrt(Norm1(rA) * Norm00(rA))) / (1-3*Epsilon);
  
  if(err!=0 || tau >= lambda) {
    if(cfg.msg) cout << "Could not verify lower bound for least singular value, aborting" << endl;
    return;
  }
  
  int rhs = RowLen(b);
  int prec = RowLen(xappstag) / rhs;
  civector d(n);
  cimatrix def(n,rhs);
  cvector y(n);
  xappstag = S*xappstag;
  
  for(int r=1 ; r<=rhs ; r++) {
    residual(d,b[Col(r)],A,xappstag(1,n,(r-1)*prec+1,r*prec),prec,cfg.highPrecApprxSolution);    
    def[Col(r)] = d;
  }
  
  rvector nrm(rhs);
  rvector delta(rhs);
 
  for(int r=1 ; r<=rhs ; r++) {
    nrm[r] = Norm2(Sup(abs(def[Col(r)])));
    delta[r] = ( ( 1.0 / (lambda - tau) ) * nrm[r] ) / (1-3*Epsilon);
  }
  
  if(cfg.msg) {
    cout << "Verified lower bound for least singular value: " << subd(lambda, tau) << endl;
    if(rhs == 1) {
      cout << "Norm(b-A*xapp): " << nrm;
      cout << "Delta: " << delta;
    }
    cout << "Computing solution enclosure..." << endl;
  } 

  y = complex(1.0,1.0);
  
  for(int r=1 ; r<=rhs ; r++)
    xx[Col(r)] = S * (civector(xx[Col(r)]) + interval(-delta[r],delta[r]) * scvector(y));

  opdotprec = oldprec;
  
}


/*!Checks if system is square. If system id under- or overdetermined, it is transformed
 * into a corresponding square system. */
template<typename TSysMat, typename TRhs, typename TSolution>
inline void slssn_check ( TSysMat& A, TRhs& b, TSolution& xx, int& Err, slssnconfig cfg) {
  int m   = Ub(A,ROW) - Lb(A,ROW) + 1;
  int n   = Ub(A,COL) - Lb(A,COL) + 1;
  int rhs = Ub(b,COL) - Lb(b,COL) + 1;
  int dim = m+n;
  int olddotprec = opdotprec;
  
  if(m == n) {

    //Square
    slssn_start(A, b, xx, Err, cfg);

  } else if (m > n) {

    //Overdetermined
    if(cfg.msg) std::cout << "Overdetermined system, solving least squares problem" << std::endl; 

    TSysMat BIG_A(dim,dim);
    TRhs BIG_b(dim,rhs);
    TSolution  BIG_xx(dim,rhs);

    BIG_A( 1,m, 1,n ) = A;
    BIG_A( 1,m, n+1,n+m ) = -Id(srmatrix(m,m));
    BIG_A( m+1,m+n, n+1,n+m ) = transpherm(A);
    
    BIG_b( m+1,m+n,1,rhs ) = 0.0;    
    BIG_b( 1,m,1,rhs ) = b;
    
    cfg.spd = false;
    
    slssn_start( BIG_A, BIG_b, BIG_xx, Err, cfg); 
 
    xx = BIG_xx( 1,n,1,rhs );    

  } else if (m < n) {

    //Underdetermined
    if(cfg.msg) std::cout << "Underdetermined system, computing minimal norm solution" << std::endl; 

    TSysMat BIG_A(dim,dim);
    TRhs BIG_b(dim,rhs);
    TSolution  BIG_xx(dim,rhs);

    BIG_A( 1,n, 1,m ) = transpherm(A);
    BIG_A( 1,n, m+1,m+n ) = -Id(srmatrix(n,n));
    BIG_A( n+1,n+m, m+1,m+n ) = A;
    
    BIG_b( 1,n,1,rhs ) = 0.0;    
    BIG_b( n+1,n+m,1,rhs ) = b;
    
    cfg.spd = false;
    
    slssn_start( BIG_A, BIG_b, BIG_xx, Err, cfg);    

    xx = BIG_xx( m+1,m+n,1,rhs );    

  }
  
  opdotprec = olddotprec;
}

//Below: Starter functions for different data types called by the user

void slssn(srmatrix& Ao, rmatrix& bo, imatrix& xx, int& err, slssnconfig cfg) {
  slssn_check(Ao,bo,xx,err,cfg);
}

void slssn(srmatrix& Ao, imatrix& bo, imatrix& xx, int& err, slssnconfig cfg) {
  slssn_check(Ao,bo,xx,err,cfg);
}

void slssn(simatrix& A, imatrix& b, imatrix& xx, int& err, slssnconfig cfg) {
  slssn_check(A,b,xx,err,cfg);
}

void slssn(srmatrix& A, rvector& b, ivector& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  rmatrix B(n,1);
  imatrix X(n,1);
  B[Col(1)] = b;
  
  slssn_check(A,B,X,err,cfg);
  
  xx = X[Col(1)];
}

void slssn(srmatrix& A, ivector& b, ivector& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  imatrix B(n,1);
  imatrix X(n,1);
  B[Col(1)] = b;
  
  slssn_check(A,B,X,err,cfg);
  
  xx = X[Col(1)];
}

void slssn(simatrix& A, ivector& b, ivector& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  imatrix B(n,1);
  imatrix X(n,1);
  B[Col(1)] = b;
  
  slssn_check(A,B,X,err,cfg);
  
  xx = X[Col(1)];
}


void slssn(scmatrix& Ao, cmatrix& bo, cimatrix& xx, int& err, slssnconfig cfg) {
  slssn_check(Ao,bo,xx,err,cfg);
}

void slssn(scmatrix& Ao, cimatrix& bo, cimatrix& xx, int& err, slssnconfig cfg) {
  slssn_check(Ao,bo,xx,err,cfg);
}

void slssn(scimatrix& A, cimatrix& b, cimatrix& xx, int& err, slssnconfig cfg) {
  slssn_check(A,b,xx,err,cfg);
}

void slssn(scmatrix& A, cvector& b, civector& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  cmatrix B(n,1);
  cimatrix X(n,1);
  B[Col(1)] = b;
  
  slssn_check(A,B,X,err,cfg);
  
  xx = X[Col(1)];
}

void slssn(scmatrix& A, civector& b, civector& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  cimatrix B(n,1);
  cimatrix X(n,1);
  B[Col(1)] = b;
  
  slssn_check(A,B,X,err,cfg);
  
  xx = X[Col(1)];
}

void slssn(scimatrix& A, civector& b, civector& xx, int& err, slssnconfig cfg) {
  int n = RowLen(A);
  cimatrix B(n,1);
  cimatrix X(n,1);
  B[Col(1)] = b;
  
  slssn_check(A,B,X,err,cfg);
  
  xx = X[Col(1)];
}


} //namespace cxsc