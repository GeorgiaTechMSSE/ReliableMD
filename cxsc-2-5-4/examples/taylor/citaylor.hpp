/*
**  Copyright (C) 1999-2008 F. Blomquist, M. Braeuer, 
**                          W. Hofschuster, W. Kraemer
**                          Wiss. Rechnen/Softwaretechnologie
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
*/

////////////////////////////////////////////////////////////////////////
//
// Headerfile citaylor.hpp for onedimensional complex Taylor-arithmetic
//
////////////////////////////////////////////////////////////////////////

/*----------------------------------------------------------------------

Compiler: g++/gcc version 3.3

Definition of the class citaylor:

class citaylor for calculating all Taylor coefficients up to the maximal 
order p.

  Elements:
       int p .................... maximal order of the Taylor coefficients
       civector tayl ............ complex interval vector; Lb=0, Ub=p;

       static ivector faks ...... storing n! with n<=170
       static int initialized ... switcher for initialization of faks.
       static void initialize().. performs initialization of faks.

  Element functions, methodes:

            look the implementation

Implementation file: citaylor.cpp

-----------------------------------------------------------------------*/

#ifndef _CITAYLOR_H
#define _CITAYLOR_H

// cxsc headers
#include "cimath.hpp"
#include <cinterval.hpp>
#include <civector.hpp>
#include <cidot.hpp>

// C++ standard headers
#include <iostream>

using namespace cxsc;

namespace taylor{
enum{
    _i_ln,

    _i_tan,
    _i_cot,

    _i_asin,
    _i_acos,
    _i_atan,
    _i_acot,

    _i_tanh,
    _i_coth,

    _i_asinh,
    _i_acosh,
    _i_atanh,
    _i_acoth,
};

///////////////////////////////////////////////////////////////
//
// Class citaylor
//
///////////////////////////////////////////////////////////////

// Class citaylor for calculating all Taylor-coefficients up to order p
// of functions with one independent complex variable.

class citaylor{

 private:
  int p;         // max. Taylor-order of the object;
  civector tayl; // Complex interval vector with Taylor-koefficients

  static ivector faks;
  static int initialized;
  static void initialize();
  explicit citaylor(int order);

 public:
  // Constructors and Destructors:
  citaylor();
  citaylor(const citaylor& );
  citaylor(int order, const real& value);    //(x,1,0,..,0) Conversion error!
  citaylor(int order, const complex& value); //(x,1,0,..,0) Conversion error!
  citaylor(int order, const interval& value);  // (x,1,0,...,0)
  citaylor(int order, const cinterval& value); // (x,1,0,...,0)
  ~citaylor(){;};

  // Initialization functions for independent variables (x,1,0,...,0):
  // Caution: (x,1,0,...,0)  for complex x conversion errors are possible!
  friend citaylor var_citaylor(int ord, const real& c);
  friend citaylor var_citaylor(int ord, const complex& c);
  friend citaylor var_citaylor(int ord, const interval& c);
  friend citaylor var_citaylor(int ord, const cinterval& c);

  // Initialization functions for constants (c,0,0,...,0):
  // Caution: (c,0,0,...,0)  for complex c conversion errors are possible!
  friend citaylor const_citaylor(int ord, const real& c);
  friend citaylor const_citaylor(int ord, const complex& c);
  friend citaylor const_citaylor(int ord, const interval& c);
  friend citaylor const_citaylor(int ord, const cinterval& c);

  // assignment operators
  citaylor operator=(const citaylor& );
  citaylor operator=(int);
  citaylor operator=(const real& );
  citaylor operator=(const complex& );
  citaylor operator=(const interval& );
  citaylor operator=(const cinterval& );

  // access to the components:
  friend int get_order(const citaylor& x);
  friend civector get_all_coef(const citaylor& x);
  friend cinterval get_j_coef(const citaylor& x, int j);

  // access to the derivative of order j:
  friend cinterval get_j_derive(const citaylor& x, int j); 
  friend cinterval get_j_derivative(const citaylor& x, int j); 
 
  // Output:
  friend void print_citaylor(const citaylor& x);

  // Overloading the operators for elements of the class citaylor:

  // operator - :
  friend citaylor operator-(const citaylor& x);

  // operators +,-,*,/  for (citaylor, citaylor):
  friend citaylor operator-(const citaylor& x, const citaylor& y);
  friend citaylor operator+(const citaylor& x, const citaylor& y);
  friend citaylor operator*(const citaylor& x, const citaylor& y);
  friend citaylor operator/(const citaylor& x, const citaylor& y);

  // operators +,-,*,/ for (cinterval, citaylor):
  friend citaylor operator-(const cinterval& x, const citaylor& y);
  friend citaylor operator+(const cinterval& x, const citaylor& y);
  friend citaylor operator*(const cinterval& x, const citaylor& y);
  friend citaylor operator/(const cinterval& x, const citaylor& y);

  // operators +,-,*,/ for (citaylor, cinterval):
  friend citaylor operator-(const citaylor& x, const cinterval& y);
  friend citaylor operator+(const citaylor& x, const cinterval& y);
  friend citaylor operator*(const citaylor& x, const cinterval& y);
  friend citaylor operator/(const citaylor& x, const cinterval& y);

  // operators +,-,*,/ for (interval, citaylor):
  friend citaylor operator-(const interval& x, const citaylor& y);
  friend citaylor operator+(const interval& x, const citaylor& y);
  friend citaylor operator*(const interval& x, const citaylor& y);
  friend citaylor operator/(const interval& x, const citaylor& y);

  // operators +,-,*,/ for (citaylor, interval):
  friend citaylor operator-(const citaylor& x, const interval& y);
  friend citaylor operator+(const citaylor& x, const interval& y);
  friend citaylor operator*(const citaylor& x, const interval& y);
  friend citaylor operator/(const citaylor& x, const interval& y);


  // operators +,-,*,/ for (complex, citaylor):
  friend citaylor operator-(const complex& x, const citaylor& y);
  friend citaylor operator+(const complex& x, const citaylor& y);
  friend citaylor operator*(const complex& x, const citaylor& y);
  friend citaylor operator/(const complex& x, const citaylor& y);

  // operators +,-,*,/ for (itaylor, complex):
  friend citaylor operator-(const citaylor& x, const complex& y);
  friend citaylor operator+(const citaylor& x, const complex& y);
  friend citaylor operator*(const citaylor& x, const complex& y);
  friend citaylor operator/(const citaylor& x, const complex& y);

  // operators +,-,*,/ for (real, citaylor):
  friend citaylor operator-(const real& x, const citaylor& y);
  friend citaylor operator+(const real& x, const citaylor& y);
  friend citaylor operator*(const real& x, const citaylor& y);
  friend citaylor operator/(const real& x, const citaylor& y);

  // operators +,-,*,/ for (citaylor, real):
  friend citaylor operator-(const citaylor& x, const real& y);
  friend citaylor operator+(const citaylor& x, const real& y);
  friend citaylor operator*(const citaylor& x, const real& y);
  friend citaylor operator/(const citaylor& x, const real& y);

  // operators +,-,*,/ for (int, citaylor):
  friend citaylor operator-(int x, const citaylor& y);
  friend citaylor operator+(int x, const citaylor& y);
  friend citaylor operator*(int x, const citaylor& y);
  friend citaylor operator/(int x, const citaylor& y);

  // operators +,-,*,/ for (citaylor, int):
  friend citaylor operator-(const citaylor& x, int y);
  friend citaylor operator+(const citaylor& x, int y);
  friend citaylor operator*(const citaylor& x, int y);
  friend citaylor operator/(const citaylor& x, int y);


  // Overloading the standard functions:
  friend citaylor sqr (const citaylor& x);
  friend citaylor sqrt(const citaylor& x);
  friend citaylor sqrt(const citaylor& x, int n);
  friend citaylor sqrt1px2(const citaylor& x);
  friend citaylor sqrt1mx2(const citaylor& x);
  friend citaylor sqrtx2m1(const citaylor& x);
  friend citaylor pow(const citaylor& x, const interval& alpha);
  friend citaylor pow(const citaylor& x, const cinterval& alpha);
  friend citaylor power(const citaylor& x, int n);
  friend citaylor exp(const citaylor& x);
  friend citaylor expm1(const citaylor& x);

  // Help function
  friend void f_g_u(const citaylor& f, const citaylor& g, const citaylor& u,
                    int nb_function);

  friend citaylor ln(const citaylor& x);
  friend citaylor sin(const citaylor& x);
  friend citaylor cos(const citaylor& x);
  friend citaylor tan(const citaylor& x);
  friend citaylor cot(const citaylor& x);

  friend citaylor sinh(const citaylor& x);
  friend citaylor cosh(const citaylor& x);
  friend citaylor tanh(const citaylor& x);
  friend citaylor coth(const citaylor& x);

  friend citaylor asin(const citaylor& x);
  friend citaylor acos(const citaylor& x);
  friend citaylor atan(const citaylor& x);
  friend citaylor acot(const citaylor& x);

  friend citaylor asinh(const citaylor& x);
  friend citaylor acosh(const citaylor& x);
  friend citaylor atanh(const citaylor& x);
  friend citaylor acoth(const citaylor& x);

};



citaylor var_citaylor(int ord, const real& c);
citaylor var_citaylor(int ord, const complex& c);
citaylor var_citaylor(int ord, const interval& c);
citaylor var_citaylor(int ord, const cinterval& c);

citaylor const_citaylor(int ord, const real& c);
citaylor const_citaylor(int ord, const complex& c);
citaylor const_citaylor(int ord, const interval& c);
citaylor const_citaylor(int ord, const cinterval& c);

int get_order(const citaylor& x);
civector get_all_coef(const citaylor& x);
cinterval get_j_coef(const citaylor& x, int j);

cinterval get_j_derive(const citaylor& x, int j); 
cinterval get_j_derivative(const citaylor& x, int j); 
 
void print_citaylor(const citaylor& x);

citaylor operator-(const citaylor& x);
citaylor operator-(const citaylor& x, const citaylor& y);
citaylor operator+(const citaylor& x, const citaylor& y);
citaylor operator*(const citaylor& x, const citaylor& y);
citaylor operator/(const citaylor& x, const citaylor& y);

citaylor operator-(const cinterval& x, const citaylor& y);
citaylor operator+(const cinterval& x, const citaylor& y);
citaylor operator*(const cinterval& x, const citaylor& y);
citaylor operator/(const cinterval& x, const citaylor& y);

citaylor operator-(const citaylor& x, const cinterval& y);
citaylor operator+(const citaylor& x, const cinterval& y);
citaylor operator*(const citaylor& x, const cinterval& y);
citaylor operator/(const citaylor& x, const cinterval& y);

citaylor operator-(const interval& x, const citaylor& y);
citaylor operator+(const interval& x, const citaylor& y);
citaylor operator*(const interval& x, const citaylor& y);
citaylor operator/(const interval& x, const citaylor& y);

citaylor operator-(const citaylor& x, const interval& y);
citaylor operator+(const citaylor& x, const interval& y);
citaylor operator*(const citaylor& x, const interval& y);
citaylor operator/(const citaylor& x, const interval& y);


citaylor operator-(const complex& x, const citaylor& y);
citaylor operator+(const complex& x, const citaylor& y);
citaylor operator*(const complex& x, const citaylor& y);
citaylor operator/(const complex& x, const citaylor& y);

citaylor operator-(const citaylor& x, const complex& y);
citaylor operator+(const citaylor& x, const complex& y);
citaylor operator*(const citaylor& x, const complex& y);
citaylor operator/(const citaylor& x, const complex& y);

citaylor operator-(const real& x, const citaylor& y);
citaylor operator+(const real& x, const citaylor& y);
citaylor operator*(const real& x, const citaylor& y);
citaylor operator/(const real& x, const citaylor& y);

citaylor operator-(const citaylor& x, const real& y);
citaylor operator+(const citaylor& x, const real& y);
citaylor operator*(const citaylor& x, const real& y);
citaylor operator/(const citaylor& x, const real& y);

citaylor operator-(int x, const citaylor& y);
citaylor operator+(int x, const citaylor& y);
citaylor operator*(int x, const citaylor& y);
citaylor operator/(int x, const citaylor& y);

citaylor operator-(const citaylor& x, int y);
citaylor operator+(const citaylor& x, int y);
citaylor operator*(const citaylor& x, int y);
citaylor operator/(const citaylor& x, int y);


citaylor sqr (const citaylor& x);
citaylor sqrt(const citaylor& x);
citaylor sqrt(const citaylor& x, int n);
citaylor sqrt1px2(const citaylor& x);
citaylor sqrt1mx2(const citaylor& x);
citaylor sqrtx2m1(const citaylor& x);
citaylor pow(const citaylor& x, const interval& alpha);
citaylor pow(const citaylor& x, const cinterval& alpha);
citaylor power(const citaylor& x, int n);
citaylor exp(const citaylor& x);
citaylor expm1(const citaylor& x);

void f_g_u(const citaylor& f, const citaylor& g, const citaylor& u,
                    int nb_function);

citaylor ln(const citaylor& x);
citaylor sin(const citaylor& x);
citaylor cos(const citaylor& x);
citaylor tan(const citaylor& x);
citaylor cot(const citaylor& x);

citaylor sinh(const citaylor& x);
citaylor cosh(const citaylor& x);
citaylor tanh(const citaylor& x);
citaylor coth(const citaylor& x);

citaylor asin(const citaylor& x);
citaylor acos(const citaylor& x);
citaylor atan(const citaylor& x);
citaylor acot(const citaylor& x);

citaylor asinh(const citaylor& x);
citaylor acosh(const citaylor& x);
citaylor atanh(const citaylor& x);
citaylor acoth(const citaylor& x);



} // End of namespace taylor

#endif
