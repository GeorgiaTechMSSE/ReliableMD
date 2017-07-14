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
 * This is a small test program checking the correctness of th dot product
 * algorithms and if the BLAS library used interferes with a previously changed
 * rounding mode
 */

#include <iostream>
#include <dot.hpp>
#include <idot.hpp>
#include <cdot.hpp>
#include <cidot.hpp>
#include <testgen.hpp>

using namespace cxsc;
using namespace std;

int main() {
  dotprecision dot;
  idotprecision idot;
  cdotprecision cdot;
  cidotprecision cidot;
  real cond=1e50;
  int n = 10000;
  real res;
  complex resc;
  bool pass_dot1 = true;
    
  cout << "Testing DotK computations..." << endl;
  cout << "  RDotK...";
  cout.flush();  
  rvector x,y;  
  GenDot(n,cond,x,y,res);
  for(int k=2 ; k<=5 ; k++) {
    dot.set_k(k);
    dot = 0.0;
    accumulate(dot,x,y);
    interval tmp;
    rnd(dot,tmp);
    if(Inf(tmp) > res  ||  Sup(tmp) < res) {
      pass_dot1 = false;
      break;
    }
  }
  
  if(pass_dot1) {
    cout << "ok" << endl;
  } else {
    cout << "failed" << endl;
  }

  cout << "  IDotK...";
  cout.flush();
  bool pass_dot2 = true;
  ivector xi,yi;  
  GenIDot(n,cond,xi,yi,res);
  for(int k=2 ; k<=5 ; k++) {
    idot.set_k(k);
    idot = 0.0;
    accumulate(idot,xi,yi);
    interval tmp = rnd(idot);
    if(Inf(tmp) > res  ||  Sup(tmp) < res) {
      pass_dot2 = false;
      break;
    }
  }
  
  if(pass_dot2) {
    cout << "ok" << endl;
  } else {
    cout << "failed" << endl;
  }

  cout << "  CDotK...";
  cout.flush();
  bool pass_dot3 = true;
  cvector xc,yc;  
  GenCDot(n,cond,xc,yc,resc);
  for(int k=2 ; k<=5 ; k++) {
    cdot.set_k(k);
    cdot = 0.0;
    accumulate(cdot,xc,yc);
    cinterval tmp = cinterval(rnd(cdot,RND_DOWN),rnd(cdot,RND_UP));
    if(InfRe(tmp) > Re(resc)  ||  SupRe(tmp) < Re(resc) ||
       InfIm(tmp) > Im(resc)  ||  SupIm(tmp) < Im(resc)) {
      pass_dot3 = false;
      break;
    }
  }
  
  if(pass_dot3) {
    cout << "ok" << endl;
  } else {
    cout << "failed" << endl;
  }

  cout << "  CIDotK...";
  cout.flush();
  bool pass_dot4 = true;
  civector xci,yci;  
  GenCIDot(n,cond,xci,yci,resc);
  for(int k=2 ; k<=5 ; k++) {
    cidot.set_k(k);
    cidot = 0.0;
    accumulate(cidot,xci,yci);
    cinterval tmp = rnd(cidot);
    if(InfRe(tmp) > Re(resc)  ||  SupRe(tmp) < Re(resc) ||
       InfIm(tmp) > Im(resc)  ||  SupIm(tmp) < Im(resc)) {
      pass_dot4 = false;
      break;
    }
  }
  
  if(pass_dot4) {
    cout << "ok" << endl;
  } else {
    cout << "failed" << endl;
  }
  
  
  cout << "Testing BLAS computations...";
  cout.flush();
  cond = 1e10;
  n = 100;
  rmatrix A(n,n);
  imatrix B(n,n);
  imatrix C1(n,n),C2(n,n);
  GenMat(n,A,cond);
  GenIMat(n,B,cond);  
  idotprecision accu(0);
  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=n ; j++) {
      accu = 0.0;
      accumulate(accu,A[i],B[Col(j)]);
      C1[i][j] = rnd(accu); 
    }
  }

  opdotprec = 1;
  C2 = A*B;
  bool pass_blas = true;
  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=n ; j++) {
      if(Inf(C2[i][j]) > Inf(C1[i][j])  ||  Sup(C2[i][j]) < Sup(C1[i][j])) {
        pass_blas = false;
        break;
      }
    }
    if(!pass_blas) break;
  }
  
  if(pass_blas)
    cout << "ok" << endl;
  else
    cout << "failed" << endl;
  
  cout << endl;  
  if(!pass_dot1 || !pass_dot2 || !pass_dot3 || !pass_dot4)
    cout << "DotK algorithm does not compute correct results. "
         << "You must set appropriate compiler flags to guarantee fully IEEE compliant "
         << "computations when compiling C-XSC, for example by using SSE floating point registers if available." 
         << "If this is not possible, set the dot product precision to K=0 (maximum precision) or K=1 (pure floating point)."<< endl << endl;
  
  if(!pass_blas)
    cout << "Your BLAS library seems not to take a changed floating point environment "
         << "into account (changed rounding mode) leading to wrong enclosures. "
         << "Please use a different BLAS library (ATLAS BLAS should work and provides "
         << "an optimized BLAS-Version for every system)" << endl << endl;
  
  if(pass_blas && pass_dot1 && pass_dot2 && pass_dot3 && pass_dot4)
    cout << "All tests passed" << endl;
  
  return 0;
}

