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
 * This is a test program for performing tests with the dot product algorithms of C-XSC.
 */ 

#include <iostream>
#include <rvector.hpp>
#include <ivector.hpp>
#include <cvector.hpp>
#include <civector.hpp>
#include <idot.hpp>
#include <cidot.hpp>
#include "testgen.hpp"
#include "sys/time.h"
#include <fenv.h>

using namespace cxsc;
using namespace std;

double GetTime() {
   struct timeval _tp;

   gettimeofday(&_tp,0);
   
   return _tp.tv_sec + _tp.tv_usec / 1000000.0;
}

void spreal(int n, real cond, int kmax, int rep) {
    rvector x,y;
    double start,end;
    real exact;
    
    cout << SetPrecision(20,16);   
    
    cout << "Generating vectors..." << endl;
    GenDot(n, cond, x, y, exact);
    cout << "Done." << endl << endl;
    
    cout << "Exact result: " << exact << endl << endl;
    
    dotprecision accu(0);
    interval res_enc;
    for(int k=0 ; k<=kmax ; k++) {
      cout << "k=" << k << ":" << endl;
      accu.set_k(k);
      accu = 0.0;

      start = GetTime();
      for(int r=1 ; r<=rep ; r++) {
        accu = 0.0;
        accumulate(accu,x,y);
        rnd(accu,res_enc);
      }
      end = GetTime();

      cout << res_enc << "\t\t\t" << "diam = " << diam(res_enc) << endl;
      if(in(exact,res_enc))
        cout << "Is an enclosure of the correct result!" << endl;
      else
        cout << "Is NOT an enclosure of the correct result!" << endl;
      cout << "Time used: " << end-start << "s" << endl;
      cout << endl;
    }
   
}


void spinterval(int n, real cond, int kmax, int rep) {
    ivector x,y;
    double start,end;
    real exact;
    
    cout << SetPrecision(20,16);   
    
    cout << "Generating vectors..." << endl;
    GenIDot(n, cond, x, y, exact);
    cout << "Done." << endl << endl;
    
    cout << "Exact result of midpoint dot product: " << exact << endl << endl;
    idotprecision accu(0);
    interval res_enc;
    for(int k=0 ; k<=kmax ; k++) {
      cout << "k=" << k << ":" << endl;
      accu.set_k(k);
      accu = 0.0;

      start = GetTime();
      for(int r=1 ; r<=rep ; r++) {
        accu = 0.0;
        accumulate(accu,x,y);
        res_enc = rnd(accu);
      }
      end = GetTime();

      cout << res_enc << "\t\t\t" << "diam = " << diam(res_enc) << endl;
      if(in(exact,res_enc))
        cout << "Is an enclosure of the correct result!" << endl;
      else
        cout << "Is NOT an enclosure of the correct result!" << endl;
      cout << "Time used: " << end-start << "s" << endl;
      cout << endl;
    }         
}

void spcomplex(int n, real cond, int kmax, int rep) {
    cvector x,y;
    double start,end;
    complex exact;
    
    cout << SetPrecision(20,16);   
    
    cout << "Generating vectors..." << endl;
    GenCDot(n, cond, x, y, exact);
    cout << "Done." << endl << endl;
    
    cout << "Exaxt result: " << exact << endl << endl;
    
    cdotprecision accu(0);
    cinterval res_enc;
    for(int k=0 ; k<=kmax ; k++) {
      cout << "k=" << k << ":" << endl;
      accu.set_k(k);
      accu = 0.0;

      start = GetTime();
      for(int r=1 ; r<=rep ; r++) {
        accu = 0.0;
        accumulate(accu,x,y);
        UncheckedSetInf(res_enc,rnd(accu,RND_DOWN));
	UncheckedSetSup(res_enc,rnd(accu,RND_UP));
      }
      end = GetTime();

      cout << res_enc << "\t\t\t" << "diam = " << diam(res_enc) << endl;
      if(in(Re(exact),Re(res_enc))  &&  in(Im(exact),Im(res_enc)))
        cout << "Is an enclosure of the correct result!" << endl;
      else
        cout << "Is NOT an enclosure of the correct result!" << endl;
      cout << "Time used: " << end-start << "s" << endl;
      cout << endl;
    }         
}

void spcinterval(int n, real cond, int kmax, int rep) {
    civector x,y;
    double start,end;
    complex exact;
    
    cout << SetPrecision(20,16);   
    
    cout << "Generating vectors..." << endl;
    GenCIDot(n, cond, x, y, exact);
    cout << "Done." << endl << endl;
    
    cout << "Exact result of midpoint dot product: " << exact << endl << endl;
    
    cidotprecision accu;
    cinterval res_enc;
    for(int k=0 ; k<=kmax ; k++) {
      cout << "k=" << k << ":" << endl;
      accu.set_k(k);
      accu = 0.0;

      start = GetTime();
      for(int r=1 ; r<=rep ; r++) {
        accu = 0.0;
        accumulate(accu,x,y);
        res_enc = rnd(accu);
      }
      end = GetTime();

      cout << res_enc << "\t\t\t" << "diam = " << diam(res_enc) << endl;
      if(in(Re(exact),Re(res_enc))  &&  in(Im(exact),Im(res_enc)))
        cout << "Is an enclosure of the correct result!" << endl;
      else
        cout << "Is NOT an enclosure of the correct result!" << endl;
      cout << "Time used: " << end-start << "s" << endl;
      cout << endl;
    }              
}


int main() {
    int type, n, kmax, rep;
    real cond;
    
    cout << endl;
    cout << "Dot product test tool" << endl;
    cout << "==========================" << endl;
    cout << endl;
    
    cout << "Please choose a type:" << endl;
    cout << "(1) real" << endl;
    cout << "(2) interval" << endl;
    cout << "(3) complex" << endl;
    cout << "(4) cinterval" << endl;
    cout << "> "; cin  >> type; cout << endl;
    
    cout << "Dimension n: "; 
    cin >> n; cout << endl;
    
    cout << "Desired condition of dot product: ";
    cin >> cond; cout << endl;
    
    cout << "Maximum accuracy to use for DotK: ";
    cin >> kmax; cout << endl;
    
    cout << "How often should the calculations be repeated? ";
    cin >> rep; cout << endl;

    if(rep <= 0) rep = 1;

    cout << endl << endl;    
    
    switch(type) {
      case 1:
           spreal(n,cond,kmax,rep);
           break;
      case 2:
           spinterval(n,cond,kmax,rep);
           break;
      case 3:
           spcomplex(n,cond,kmax,rep);
           break;
      case 4:
           spcinterval(n,cond,kmax,rep);
           break;
      default:
           spreal(n,cond,kmax,rep);
    }
    
    cout << endl << "Done." << endl;
    
    return 0;
}
