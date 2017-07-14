/*
**  FastPLSS: A library of (parallel) verified linear (interval) system 
**  solvers using C-XSC (V 0.3)
**
**  CXSC is a C++ library for eXtended Scientific Computing
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2008 Wiss. Rechnen/Softwaretechnologie
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

#ifndef TESTGEN_HEADER_INCLUDED
#define TESTGEN_HEADER_INCLUDED

#include <real.hpp>
#include <interval.hpp>
#include <cinterval.hpp>
#include <complex.hpp>
#include <rvector.hpp>
#include <cvector.hpp>
#include <ivector.hpp>
#include <civector.hpp>
#include <rmatrix.hpp>
#include <cmatrix.hpp>
#include <imatrix.hpp>
#include <cimatrix.hpp>

using namespace std;
using namespace cxsc;

//! Generates a real matrix choosen by type
void GenMat(int n, rmatrix &A, int type);
//! Generates an interval matrix choosen by type
void GenIMat(int n, imatrix &A, int type);
//! Generates a complex matrix choosen by type
void GenCMat(int n, cmatrix &A, int type);
//! Generates a complex interval matrix choosen by type
void GenCIMat(int n, cimatrix &A, int type);

//! Generates a real matrix with chosen condition
void GenMat(int n, rmatrix &A, real cond);
//! Generates an interval matrix tightly enclosing a real matrix with chosen condition
void GenIMat(int n, imatrix &A, real cond);
//! Generates a complex matrix with chosen condition
void GenCMat(int n, cmatrix &A, real cond);
//! Generates a complex interval matrix tightly enclosing a complex matrix with chosen condition
void GenCIMat(int n, cimatrix &A, real cond);

void GenMat2(int n, rmatrix &A, real cond);

//! Generate a real dot product with given condition
void GenDot(int n, real c, rvector &x, rvector &y, real &res);
//! Generate a complex dot product with given condition
void GenCDot(int n, real c, cvector &x, cvector &y, complex &res);
//! Generate an interval dot product enclosing real dot product with given condition
void GenIDot(int n, real c, ivector &x, ivector &y, real &res);
//! Generate a complex interval dot product enclosing complex dot product with given condition
void GenCIDot(int n, real c, civector &x, civector &y, complex &res);

//! Generates an real matrix with entries randomly distributed in [-1,1]
rmatrix randn(int m, int n);

#endif
