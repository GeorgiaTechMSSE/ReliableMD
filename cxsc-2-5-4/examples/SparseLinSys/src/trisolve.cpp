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

/* Implementation of functions for the solution of triangular point systems
   or triangular systems with interval right hand side */

#include "trisolve.hpp"
#include "QRdecomp.hpp"

#include <vector>
#include <fenv.h>
#include <scivector.hpp>

using namespace std;

namespace cxsc {  
  
//! Compute i/r
inline interval intpointdiv(const interval& i, const real& r) {
  return i/r;
}

//! Compute i/r, faster than operator if Im(r)==0
inline cinterval intpointdiv(const cinterval& i, const complex& r) {
  if(Im(r) == 0.0)
    return i/Re(r);
  else
    return i/r;
}

//! Compute i/r, faster than operator if Im(r)==0
inline complex intpointdiv(const complex& i, const complex& r) {
  if(Im(r) == 0.0)
    return i/Re(r);
  else
    return i/r;
}

//! Solve lower triangular system Lx=b with simple forward substitution
template<typename TRhs, typename Tx>
void lowtrisolve(const srmatrix& L, const TRhs& b, Tx& x) {
    int n = RowLen(L);
    const vector<int>& ind   = L.row_indices();
    const vector<int>& p     = L.column_pointers();
    const vector<real>& val  = L.values();
    
    x = b;
    int lb = Lb(x);

    for(int j=0 ; j<n ; j++) {
      x[j+lb] /= val[p[j]];

     #ifdef _OPENMP
     int d = p[j+1]-p[j];
     int t = (d/50 > omp_get_max_threads()) ? omp_get_max_threads() : d/50;
     #pragma omp parallel for if(d >= 100) num_threads(t)
     #endif
     for(int k=p[j]+1 ; k<p[j+1] ; k++) {
       x[ind[k]+lb] -= val[k] * x[j+lb];
     }
    }  
    
}

//! Solve upper triangular system Ux=b with simple backward substitution
template<typename TRhs, typename Tx>
void uptrisolve(const srmatrix& U, const TRhs& b, Tx& x) {
    int n = RowLen(U);
    const vector<int>& ind   = U.row_indices();
    const vector<int>& p     = U.column_pointers();
    const vector<real>& val  = U.values();

    x = b;
    int lb = Lb(x);
    
    for(int j=n-1 ; j>=0 ; j--) {
      x[lb+j] /= val[p[j+1]-1];

      #ifdef _OPENMP
      int d = p[j+1]-p[j];
      int t = (d/50 > omp_get_max_threads()) ? omp_get_max_threads() : d/50;
      #pragma omp parallel for if(d >= 100) num_threads(t)
      #endif
      for(int k=p[j] ; k<p[j+1]-1 ; k++) {
        x[ind[k]+lb] -= val[k] * x[j+lb];
      }
    }  
}

template void lowtrisolve<rvector,rvector>(const srmatrix& L, const rvector& b, rvector& x);
template void uptrisolve<rvector,rvector>(const srmatrix& U, const rvector& b, rvector& x);
template void lowtrisolve<ivector,ivector>(const srmatrix& L, const ivector& b, ivector& x);
template void uptrisolve<ivector,ivector>(const srmatrix& U, const ivector& b, ivector& x);
template void lowtrisolve<rvector,ivector>(const srmatrix& L, const rvector& b, ivector& x);
template void uptrisolve<rvector,ivector>(const srmatrix& U, const rvector& b, ivector& x);



//! Solve lower triangular system Lx=b with simple forward substitution
template<typename TRhs, typename Tx>
void lowtrisolve(const scmatrix& L, const TRhs& b, Tx& x) {
    int n = RowLen(L);
    const vector<int>& ind      = L.row_indices();
    const vector<int>& p        = L.column_pointers();
    const vector<complex>& val  = L.values();

    x = b;
    int lb = Lb(x);

    for(int j=0 ; j<n ; j++) {
      //x[j+lb] /= val[p[j]];
      x[j+lb] = intpointdiv(x[j+lb],val[p[j]]);

      #ifdef _OPENMP
      int d = p[j+1]-p[j];
      int t = (d/50 > omp_get_max_threads()) ? omp_get_max_threads() : d/50;
      #pragma omp parallel for if(d >= 100) num_threads(t)
      #endif
      for(int k=p[j]+1 ; k<p[j+1] ; k++) {
        //x[ind[k]+lb] -= val[k] * x[j+lb];
        SetRe(x[ind[k]+lb], Re(x[ind[k]+lb]) - (Re(val[k]) * Re(x[j+lb]) - Im(val[k]) * Im(x[j+lb])) );
        SetIm(x[ind[k]+lb], Im(x[ind[k]+lb]) - (Im(val[k]) * Re(x[j+lb]) + Re(val[k]) * Im(x[j+lb])) );
      }
    }  
}

//! Solve upper triangular system Ux=b with simple backward substitution
template<typename TRhs, typename Tx>
void uptrisolve(const scmatrix& U, const TRhs& b, Tx& x) {
    int n = RowLen(U);
    
    const vector<int>& ind      = U.row_indices();
    const vector<int>& p        = U.column_pointers();
    const vector<complex>& val  = U.values();

    x = b;
    int lb = Lb(x);

    for(int j=n-1 ; j>=0 ; j--) {
      //x[j+lb] /= val[p[j+1]-1];
      x[j+lb] = intpointdiv(x[j+lb],val[p[j+1]-1]);

      #ifdef _OPENMP
      int d = p[j+1]-p[j];
      int t = (d/50 > omp_get_max_threads()) ? omp_get_max_threads() : d/50;
      #pragma omp parallel for if(d >= 100) num_threads(t)
      #endif
      for(int k=p[j] ; k<p[j+1]-1 ; k++) {
        //x[ind[k]+lb] -= val[k] * x[j+lb];
        SetRe(x[ind[k]+lb], Re(x[ind[k]+lb]) - (Re(val[k]) * Re(x[j+lb]) - Im(val[k]) * Im(x[j+lb])) );
        SetIm(x[ind[k]+lb], Im(x[ind[k]+lb]) - (Im(val[k]) * Re(x[j+lb]) + Re(val[k]) * Im(x[j+lb])) );	
      }
    }  
}

template void lowtrisolve<cvector,cvector>(const scmatrix& L, const cvector& b, cvector& x);
template void uptrisolve<cvector,cvector>(const scmatrix& U, const cvector& b, cvector& x);
template void lowtrisolve<civector,civector>(const scmatrix& L, const civector& b, civector& x);
template void uptrisolve<civector,civector>(const scmatrix& U, const civector& b, civector& x);
template void lowtrisolve<cvector,civector>(const scmatrix& L, const cvector& b, civector& x);
template void uptrisolve<cvector,civector>(const scmatrix& U, const cvector& b, civector& x);



//! Compute v=max(|inf|,|sup|)
void AbsMax(const srvector& inf, const srvector& sup, srvector& v) {
  v = inf;
  const vector<real>& Sval = sup.values();
  vector<real>& Vval = v.values();
  real tmp;

  for(unsigned int i=0 ; i<Vval.size() ; i++) {
    Vval[i] = abs(Vval[i]);
    tmp = abs(Sval[i]);
    if(Vval[i] < tmp) Vval[i] = tmp;
  }
  
}

//! Compute v=max(|inf|,|sup|)
void AbsMax(const scvector& inf, const scvector& sup, srvector& v) {
  const vector<complex>& Ival = inf.values();
  const vector<complex>& Sval = sup.values();
  vector<real>& Vval = v.values();
  vector<int>& Vind = v.row_indices();
  real tmp;

  v = srvector(VecLen(inf));
  Vind = inf.row_indices();
  Vval.resize(Ival.size());
  
  fesetround(FE_UPWARD);
  
  
  for(unsigned int i=0 ; i<Vval.size() ; i++) {
    Vval[i] = sqrt( Re(Ival[i])*Re(Ival[i]) + Im(Ival[i])*Im(Ival[i]) );
    tmp = sqrt( Re(Sval[i])*Re(Sval[i]) + Im(Sval[i])*Im(Sval[i]) );
    if(Vval[i] < tmp) Vval[i] = tmp;
  }
  
  fesetround(FE_TONEAREST);  
}


//! Solve real lower triangular system by implicitly computing the inverse
void lowtrisolve_inv(const cxsc::srmatrix& L, const cxsc::ivector& bo, cxsc::ivector& x) {
    int n = RowLen(L);

    rvector xapp(n);
    ivector b(n);
    real abssup, absinf;

    const vector<int>& ind   = L.row_indices();
    const vector<int>& p     = L.column_pointers();
    const vector<real>& val  = L.values();

    vector<sivector> rhs;

    lowtrisolve(L,mid(bo),xapp);
    b = bo - L*xapp;

    for(int i=0 ; i<n ; i++) {
      abssup = abs(Sup(b[i+1]));
      absinf = abs(Inf(b[i+1]));
      if(abssup > absinf)
        SetInf(b[i+1], -abssup);
      else if(absinf > abssup)
        SetSup(b[i+1], absinf);

      rhs.push_back(sivector(n));
      rhs[i][i+1] = 1.0;
    }

    for(int j=0 ; j<n ; j++) {
      rhs[j] /= val[p[j]];
      x[j+1] = xapp[j+1] + abs(rhs[j]) * b;

      #ifdef _OPENMP
      #pragma omp parallel 
      {
      #endif
	
        #ifdef _OPENMP
        #pragma omp for schedule(dynamic,1)
        #endif
        for(int k=p[j]+1 ; k<p[j+1] ; k++) 
          rhs[ind[k]] -= val[k] * rhs[j];
      
      #ifdef _OPENMP
      }
      #endif

      Resize(rhs[j]);
    }

    rhs.clear();
}


//! Solve complex lower triangular system by implicitly computing the inverse
void lowtrisolve_inv(const cxsc::scmatrix& L, const cxsc::civector& bo, cxsc::civector& x) {
    int n = RowLen(L);

    cvector xapp(n);
    civector b(n);
    real abssupre, absinfre, absinfim, abssupim;

    const vector<int>& ind      = L.row_indices();
    const vector<int>& p        = L.column_pointers();
    const vector<complex>& val  = L.values();

    vector<scivector> rhs;

    lowtrisolve(L,mid(bo),xapp);
    b = bo - L*xapp;

    for(int i=0 ; i<n ; i++) {
      abssupre = abs(SupRe(b[i+1]));
      absinfre = abs(InfRe(b[i+1]));
      if(abssupre > absinfre)
        SetInf(Re(b[i+1]), -abssupre);
      else if(absinfre > abssupre)
        SetSup(Re(b[i+1]), absinfre);

      abssupim = abs(SupIm(b[i+1]));
      absinfim = abs(InfIm(b[i+1]));
      if(abssupim > absinfim)
        SetInf(Im(b[i+1]), -abssupim);
      else if(absinfim > abssupim)
        SetSup(Im(b[i+1]), absinfim);
     
      rhs.push_back(scivector(n));
      rhs[i][i+1] = 1.0;
    }

    for(int j=0 ; j<n ; j++) {
      rhs[j] /= val[p[j]];
      x[j+1] = xapp[j+1] + abs(rhs[j]) * b;

      #ifdef _OPENMP
      #pragma omp parallel 
      {
      #endif
      
        #ifdef _OPENMP
        #pragma omp for schedule(dynamic,1)
        #endif
        for(int k=p[j]+1 ; k<p[j+1] ; k++) {
	  if(Im(val[k]) == 0.0)
            rhs[ind[k]] += Re(-val[k]) * rhs[j];
	  else {
            rhs[ind[k]] += (Re(-val[k]) * Re(rhs[j]) + Im(val[k]) * Im(rhs[j])) 
	                   + (Re(-val[k]) * Im(rhs[j]) + Im(-val[k]) * Re(rhs[j])) * complex(0,1);
	  }
	}

	
      #ifdef _OPENMP
      }
      #endif

      Resize(rhs[j]);
    }

    rhs.clear();
}


//! Solve real upper triangular system by implicitly computing the inverse
void uptrisolve_inv(const cxsc::srmatrix& U, const cxsc::ivector& bo, cxsc::ivector& x) {
    int n = RowLen(U);

    rvector xapp(n);
    ivector b(n);
    real abssup, absinf;

    const vector<int>& ind   = U.row_indices();
    const vector<int>& p     = U.column_pointers();
    const vector<real>& val  = U.values();

    vector<sivector> rhs;

    uptrisolve(U,mid(bo),xapp);  
    b = bo - U*xapp;
    
    int lbx = Lb(x);
    int lbb = Lb(b);
    int lbxa = Lb(xapp);


    for(int i=0 ; i<n ; i++) {
      abssup = abs(Sup(b[i+lbb]));
      absinf = abs(Inf(b[i+lbb]));
      if(abssup > absinf)
        SetInf(b[i+lbb], -abssup);
      else if(absinf > abssup)
        SetSup(b[i+lbb], absinf);

      rhs.push_back(sivector(n));
      rhs[i][i+1] = 1.0;
    }

    for(int j=n-1 ; j>=0 ; j--) {
      rhs[j] /= val[p[j+1]-1];
      x[j+lbx] = xapp[j+lbxa] + abs(rhs[j]) * b;

      #ifdef _OPENMP
      #pragma omp parallel
      {
      #endif
	
        #ifdef _OPENMP
        #pragma omp for schedule(dynamic,1)
        #endif
        for(int k=p[j] ; k<p[j+1]-1 ; k++) {
          rhs[ind[k]] -= val[k] * rhs[j];
        }
      
      #ifdef _OPENMP
      }
      #endif

      Resize(rhs[j]);
    }
     
}


//! Solve complex upper triangular system by implicitly computing the inverse
void uptrisolve_inv(const cxsc::scmatrix& U, const cxsc::civector& bo, cxsc::civector& x) {

    int n = RowLen(U);

    cvector xapp(n);
    civector b(n);
    real abssupre, absinfre, absinfim, abssupim;

    const vector<int>& ind      = U.row_indices();
    const vector<int>& p        = U.column_pointers();
    const vector<complex>& val  = U.values();

    vector<scivector> rhs;

    uptrisolve(U,mid(bo),xapp);
    b = bo - U*xapp;
    
    int lbx = Lb(x);
    int lbb = Lb(b);
    int lbxa = Lb(xapp);    

    for(int i=0 ; i<n ; i++) {
      abssupre = abs(SupRe(b[i+lbb]));
      absinfre = abs(InfRe(b[i+lbb]));
      if(abssupre > absinfre)
        SetInf(Re(b[i+lbb]), -abssupre);
      else if(absinfre > abssupre)
        SetSup(Re(b[i+lbb]), absinfre);

      abssupim = abs(SupIm(b[i+lbb]));
      absinfim = abs(InfIm(b[i+lbb]));
      if(abssupim > absinfim)
        SetInf(Im(b[i+lbb]), -abssupim);
      else if(absinfim > abssupim)
        SetSup(Im(b[i+lbb]), absinfim);

      rhs.push_back(scivector(n));
      rhs[i][i+1] = 1.0;
    }


    for(int j=n-1 ; j>=0 ; j--) {
      rhs[j] /= val[p[j+1]-1];
      x[j+lbx] = xapp[j+lbxa] + abs(rhs[j]) * b;
      
      #ifdef _OPENMP
      #pragma omp parallel
      {
      #endif
      
        #ifdef _OPENMP
        #pragma omp for schedule(dynamic,1) nowait
        #endif
        for(int k=p[j] ; k<p[j+1]-1 ; k++) {
	  if(Im(val[k]) == 0.0)
            rhs[ind[k]] += Re(-val[k]) * rhs[j];
	  else
	    rhs[ind[k]] += (Re(-val[k]) * Re(rhs[j]) + Im(val[k]) * Im(rhs[j])) 
	                   + (Re(-val[k]) * Im(rhs[j]) + Im(-val[k]) * Re(rhs[j])) * complex(0,1);
        }

      #ifdef _OPENMP
      }
      #endif

      Resize(rhs[j]);
    }

}


//! Forward substituition using coordinate transformations (parallel epipeds)
template<typename TMat, typename TVec, typename TIVec, typename TSVec, typename TSIVec, typename TDenseMat, typename TDenseIMat, typename TPoint, typename TInterval>
void lowtrisolve_band_main(  TMat& A,  TMat& At,  TIVec& b, TIVec& x, int l, TDenseIMat& Basis, bool baseStore, int& err) {
  int n = RowLen(A);
  TIVec xtmp(l);

  x = 0.0;
  xtmp = 0.0;

  if(l == 0) {
    for(int i=1 ; i<=n ; i++)
      x[i] = b[i] / A(i,i);
    return;
  } else if (l == 1) {
    lowtrisolve(A,b,x);
    return;
  }

  //Starting values: Compute x[1] to x[l] by forward substitution
  TIVec xstart(l);
  lowtrisolve_inv(A(1,l,1,l),b(1,l),xstart);
  x(1,l) = xstart;
  xtmp = xstart;  
  
  //Compute remaining x by forward substitution using coordinate transformation
  TDenseMat B0(l,l), B1(l,l);
  TDenseIMat B1_inv(l,l);
  TSIVec bi(l);
 
  B0 = Id(B0);
  //bi = 0.0;
  bool computeBase = RowLen(Basis) == 0 || ColLen(Basis) == 0;
  
  if(baseStore) {
    if(computeBase) {

      Basis = TDenseIMat(3*l,(n-l)*l);
    
      for(int i=l+1 ; i<=n ; i++) {
     
        TDenseMat Ai_mid(l,l);
        TDenseIMat Ai(l,l);
    
        TSVec row(At[Col(i)]);
        vector<int>& rp = row.row_indices();
        vector<TPoint>& rx = row.values();

        Ai = 0.0;
        for(int k=1 ; k<l ; k++)
          Ai[k][k+1] = 1.0;
        for(unsigned int j=0 ; j<rp.size()-1 ; j++) {
          Ai[l][rp[j]+1-(i-(l+1))] = intpointdiv(TInterval(-rx[j]), rx[rp.size()-1]);
        }
    
        Ai = Ai * B0;
        Ai_mid = Inf(Ai) + 0.5 * (Sup(Ai) - Inf(Ai));
    
        QR_sorted(Ai_mid, xtmp, B1, err);
	if(err != 0) return;
	
        QR_inv(B1, B1_inv, err);
        if(err != 0) return;
        
        B0 = B1;
      
        Basis(1,l,(i-l-1)*l+1,(i-l)*l) = (B1_inv * Ai);
        Basis(l+1,2*l,(i-l-1)*l+1,(i-l)*l) = B1_inv;
        Basis(2*l+1,3*l,(i-l-1)*l+1,(i-l)*l) = B1;         
	
        bi[l] = intpointdiv(b[i], rx[rp.size()-1]);
      
        xtmp = Basis(1,l,(i-l-1)*l+1,(i-l)*l) * xtmp + Basis(l+1,2*l,(i-l-1)*l+1,(i-l)*l) * bi;
        TDenseIMat B = Basis(2*l+1,3*l,(i-l-1)*l+1,(i-l)*l);
        x[i] = B[Ub(B,ROW)] * xtmp;      
	
      } //for(i=...
      
    } else { //if(computeBase)
    
      const vector<int>& p = A.column_pointers();
      const vector<TPoint>& rx = A.values();
 
      for(int i=l+1 ; i<=n ; i++) { 
        bi[l] = intpointdiv(b[i],rx[p[i-1]]);
      
        xtmp = Basis(1,l,(i-l-1)*l+1,(i-l)*l) * xtmp + Basis(l+1,2*l,(i-l-1)*l+1,(i-l)*l) * bi;
        TDenseIMat B = Basis(2*l+1,3*l,(i-l-1)*l+1,(i-l)*l);
        x[i] = B[Ub(B,ROW)] * xtmp;      
      }    
    
    }
    
  } else {
    
    for(int i=l+1 ; i<=n ; i++) {
      TDenseMat Ai_mid(l,l);
      TDenseIMat Ai(l,l);
    
      TSVec row(At[Col(i)]);
      vector<int>& rp = row.row_indices();
      vector<TPoint>& rx = row.values();

      Ai = 0.0;
      for(int k=1 ; k<l ; k++)
        Ai[k][k+1] = 1.0;
      for(unsigned int j=0 ; j<rp.size()-1 ; j++) {
        Ai[l][rp[j]+1-(i-(l+1))] = intpointdiv(TInterval(-rx[j]), rx[rp.size()-1]);
      }
      
      bi[l] = intpointdiv(b[i], rx[rp.size()-1]);
	
      Ai = Ai * B0;
      Ai_mid = Inf(Ai) + 0.5 * (Sup(Ai) - Inf(Ai));

      QR_sorted(Ai_mid, xtmp, B1, err);
      if(err != 0) return;
      
      QR_inv(B1, B1_inv, err);
      if(err != 0) return;
        
      xtmp = (B1_inv * Ai) * xtmp + B1_inv * bi;
      x[i] = B1[l] * xtmp;
        
      B0 = B1;
    }
    
  }

}

void lowtrisolve_band(  srmatrix& A,  srmatrix& At,  ivector& b, ivector& x, int l, imatrix& Basis, bool baseStore, int& err) {
  lowtrisolve_band_main<srmatrix,rvector,ivector,srvector,sivector,rmatrix,imatrix,real,interval>(A,At,b,x,l,Basis,baseStore,err);
}

void lowtrisolve_band(  scmatrix& A,  scmatrix& At,  civector& b, civector& x, int l, cimatrix& Basis, bool baseStore, int& err) {
  lowtrisolve_band_main<scmatrix,cvector,civector,scvector,scivector,cmatrix,cimatrix,complex,cinterval>(A,At,b,x,l,Basis,baseStore,err);
}


//! Backward substituition using coordinate transformations (parallel epipeds)
template<typename TMat, typename TVec, typename TIVec, typename TSVec, typename TSIVec, typename TDenseMat, typename TDenseIMat, typename TPoint, typename TInterval>
void uptrisolve_band_main(  TMat& A,  TMat& At,  TIVec& b, TIVec& x, int l, TDenseIMat& Basis, bool baseStore, int& err) {
  int n = RowLen(A);
  TIVec xtmp(l);

  x = 0.0;
  xtmp = 0.0;

  if(l == 0) {
    for(int i=1 ; i<=n ; i++)
      x[i] = b[i] / A(i,i);
    return;
  } else if (l == 1) {
    uptrisolve(A,b,x);
    return;
  }

  //Starting values: Compute x[n] to x[n-l] by backward substitution
  TIVec xstart(l);
  uptrisolve_inv(A(n-l+1,n,n-l+1,n), b(n-l+1,n), xstart);
  x(n-l+1,n) = xstart;
  for(int i=1 ; i<=l ; i++)
    xtmp[i] = xstart[Ub(xstart)-i+1];  


  //Compute remaining x by backward substitution using coordinate transformation
  TDenseMat B0(l,l), B1(l,l);
  TDenseIMat B1_inv(l,l);
  TSIVec bi(l);
 
  B0 = Id(B0);
  bool computeBase = RowLen(Basis) == 0 || ColLen(Basis) == 0;

  if(baseStore) {
    
    if(computeBase) {      
      Basis = TDenseIMat(3*l,(n-l)*l);
    
      for(int i=n-l ; i>0 ; i--) {
        TDenseMat Ai_mid(l,l);
        TDenseIMat Ai(l,l);
    
        TSVec row(At[Col(i)]);
        vector<int>& rp = row.row_indices();
        vector<TPoint>& rx = row.values();
	
        Ai = 0.0;
        for(int k=1 ; k<l ; k++)
          Ai[k][k+1] = 1.0;
        for(unsigned int j=rp.size()-1 ; j>0 ; j--)  {
          Ai[l][n-(rp[j]+1)-(n-l-i)+1] = intpointdiv(TInterval(-rx[j]), rx[0]);
        }
    
        Ai = Ai * B0; 
        Ai_mid = Inf(Ai) + 0.5 * (Sup(Ai) - Inf(Ai));    
	
        QR_sorted(Ai_mid, xtmp, B1, err);
	if(err != 0) return;
	
        QR_inv(B1, B1_inv, err);
	if(err != 0) return;
	
        B0 = B1;
	
        Basis(1,l,(i-1)*l+1,i*l) = (B1_inv * Ai);
        Basis(l+1,2*l,(i-1)*l+1,i*l) = B1_inv;
        Basis(2*l+1,3*l,(i-1)*l+1,i*l) = B1;   	
	
        bi[l] = intpointdiv(b[i], rx[0]);
	
        xtmp = Basis(1,l,(i-1)*l+1,i*l) * xtmp + Basis(l+1,2*l,(i-1)*l+1,i*l) * bi;
        TDenseIMat B = Basis(2*l+1,3*l,(i-1)*l+1,i*l);
        x[i] = B[Ub(B,ROW)] * xtmp;	
      }
      
    } else {   

      const vector<int>& p = A.column_pointers();
      const vector<TPoint>& rx = A.values();
      for(int i=n-l ; i>0 ; i--) {    
        bi[l] = intpointdiv(b[i],rx[p[i]-1]);
	
        xtmp = Basis(1,l,(i-1)*l+1,i*l) * xtmp + Basis(l+1,2*l,(i-1)*l+1,i*l) * bi;
        TDenseIMat B = Basis(2*l+1,3*l,(i-1)*l+1,i*l);
        x[i] = B[Ub(B,ROW)] * xtmp;
      }

    }

  } else {
   
    for(int i=n-l ; i>0 ; i--) {
      TDenseMat Ai_mid(l,l);
      TDenseIMat Ai(l,l);
    
      TSVec row(At[Col(i)]);
      vector<int>& rp = row.row_indices();
      vector<TPoint>& rx = row.values();
      
      Ai = 0.0;
      for(int k=1 ; k<l ; k++)
        Ai[k][k+1] = 1.0;
      for(unsigned int j=rp.size()-1 ; j>0 ; j--)  {
        Ai[l][n-(rp[j]+1)-(n-l-i)+1] = intpointdiv(TInterval(-rx[j]), rx[0]);
      }
  
      bi[l] = intpointdiv(b[i], rx[0]);
      
      Ai = Ai * B0; 
      Ai_mid = Inf(Ai) + 0.5 * (Sup(Ai) - Inf(Ai));    
      
      QR_sorted(Ai_mid, xtmp, B1, err);
      if(err != 0) return;
      
      QR_inv(B1, B1_inv, err);
      if(err != 0) return;
      
      xtmp = (B1_inv * Ai) * xtmp + B1_inv * bi;
      x[i] = B1[l] * xtmp;
      B0 = B1;
    }
    
  }

} // end lss_upper


void uptrisolve_band(  srmatrix& A,  srmatrix& At,  ivector& b, ivector& x, int l, imatrix& Basis, bool baseStore, int& err) {
  uptrisolve_band_main<srmatrix,rvector,ivector,srvector,sivector,rmatrix,imatrix,real,interval>(A,At,b,x,l,Basis,baseStore,err);
}

void uptrisolve_band(  scmatrix& A,  scmatrix& At,  civector& b, civector& x, int l, cimatrix& Basis, bool baseStore, int& err) {
  uptrisolve_band_main<scmatrix,cvector,civector,scvector,scivector,cmatrix,cimatrix,complex,cinterval>(A,At,b,x,l,Basis,baseStore,err);
}

} //namespace cxsc