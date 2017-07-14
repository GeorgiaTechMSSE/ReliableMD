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

////////////////////////////////////////////////////////////////////////////
//
//  Headerfile litaylor.hpp for onedimensional staggered Taylor-arithmetic,
//  Definition of class l_itaylor;
//
////////////////////////////////////////////////////////////////////////////

/*--------------------------------------------------------------------------

Compiler: g++/gcc version 3.3

Definition of the class l_itaylor:

class l_itaylor for calculating all Taylor coefficients up to the 
maximal order p.

  Elements:
            int p .................... maximal order of the Taylor coefficients
            l_ivector tayl ............. l_interval vector; Lb=0, Ub=p;

	    static l_ivector faks ........ storing n! with n<=170
            static int initialized ....... switcher for initialization of faks.
            static void initialize()...... performs initialization of faks.
	    explicit l_itaylor(int order). no automatic type matching

  Element functions, methodes:

            look the implementation

Implementation file: litaylor.cpp

-----------------------------------------------------------------------*/

#ifndef LITAYLOR_H
#define LITAYLOR_H

// cxsc headers
#include <l_imath.hpp>
#include <l_interval.hpp>
#include <l_ivector.hpp>
#include <idot.hpp>

//C++ standard headers
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
//                   Class l_itaylor
//
///////////////////////////////////////////////////////////////

// Class l_itaylor for calculating all Taylor-coefficients up to order p
// of functions with one independent variable.

class l_itaylor{

 private:
  int p;         // max. Taylor-order of the object;
  l_ivector tayl; // L_interval vector with Taylor-coefficients

  static l_ivector faks;
  static int initialized;
  static void initialize();
  explicit l_itaylor(int order);

 public:
  // Constructors and Destructors:
  l_itaylor();
  l_itaylor(const l_itaylor& );
  l_itaylor(int order, const real& value);       // (x,1,0,...,0)
  l_itaylor(int order, const l_real& value);     // (x,1,0,...,0)
  l_itaylor(int order, const interval& value);   // (x,1,0,...,0)
  l_itaylor(int order, const l_interval& value); // (x,1,0,...,0)
  ~l_itaylor(){;};

  // Initialization functions for independent variables (x,1,0,...,0):
  // Caution: (x,1,0,...,0)  for real x conversion errors are possible!
  friend l_itaylor var_l_itaylor(int ord, const real& c);
  friend l_itaylor var_l_itaylor(int ord, const l_real& c);
  friend l_itaylor var_l_itaylor(int ord, const interval& c);
  friend l_itaylor var_l_itaylor(int ord, const l_interval& c);

  // Initialization functions for constants (c,0,0,...,0):
  // Caution: (c,0,0,...,0)  for real c conversion errors are possible!
  friend l_itaylor const_l_itaylor(int ord, const real& c);
  friend l_itaylor const_l_itaylor(int ord, const l_real& c);
  friend l_itaylor const_l_itaylor(int ord, const interval& c);
  friend l_itaylor const_l_itaylor(int ord, const l_interval& c);

  // assignment operators:
  l_itaylor operator=(const l_itaylor& );
  l_itaylor operator=(int);
  l_itaylor operator=(const real& );
  l_itaylor operator=(const l_real& ); 
  l_itaylor operator=(const interval& );
  l_itaylor operator=(const l_interval& ); 

  // access to the components:
  friend int get_order(const l_itaylor& x);
  friend l_ivector get_all_coef(const l_itaylor& x);
  friend l_interval get_j_coef(const l_itaylor& x, int j);

  // access to the derivative of order j:
  friend l_interval get_j_derive(const l_itaylor& x, int j);
  friend l_interval get_j_derivative(const l_itaylor& x, int j);

  // Output:
  friend void print_l_itaylor(const l_itaylor& x);

  // Overloading the operators for elements of the class l_itaylor:

  // operator - :
  friend l_itaylor operator-(const l_itaylor& x);

  // operators +,-,*,/  for (l_itaylor, l_itaylor):
  friend l_itaylor operator-(const l_itaylor& x, const l_itaylor& y);
  friend l_itaylor operator+(const l_itaylor& x, const l_itaylor& y);
  friend l_itaylor operator*(const l_itaylor& x, const l_itaylor& y);
  friend l_itaylor operator/(const l_itaylor& x, const l_itaylor& y);

  // operators +,-,*,/ for (interval, l_itaylor):
  friend l_itaylor operator-(const interval& x, const l_itaylor& y);
  friend l_itaylor operator+(const interval& x, const l_itaylor& y);
  friend l_itaylor operator*(const interval& x, const l_itaylor& y);
  friend l_itaylor operator/(const interval& x, const l_itaylor& y);

  // operators +,-,*,/ for (l_interval, l_itaylor):
  friend l_itaylor operator-(const l_interval& x, const l_itaylor& y);
  friend l_itaylor operator+(const l_interval& x, const l_itaylor& y);
  friend l_itaylor operator*(const l_interval& x, const l_itaylor& y);
  friend l_itaylor operator/(const l_interval& x, const l_itaylor& y);

  // operators +,-,*,/ for (l_itaylor, interval):
  friend l_itaylor operator-(const l_itaylor& x, const interval& y);
  friend l_itaylor operator+(const l_itaylor& x, const interval& y);
  friend l_itaylor operator*(const l_itaylor& x, const interval& y);
  friend l_itaylor operator/(const l_itaylor& x, const interval& y);

  // operators +,-,*,/ for (l_itaylor, l_interval):
  friend l_itaylor operator-(const l_itaylor& x, const l_interval& y);
  friend l_itaylor operator+(const l_itaylor& x, const l_interval& y);
  friend l_itaylor operator*(const l_itaylor& x, const l_interval& y);
  friend l_itaylor operator/(const l_itaylor& x, const l_interval& y);

  // operators +,-,*,/ for (real, l_itaylor):
  friend l_itaylor operator-(const real& x, const l_itaylor& y);
  friend l_itaylor operator+(const real& x, const l_itaylor& y);
  friend l_itaylor operator*(const real& x, const l_itaylor& y);
  friend l_itaylor operator/(const real& x, const l_itaylor& y);

  // operators +,-,*,/ for (int, l_itaylor):
  friend l_itaylor operator-(int x, const l_itaylor& y);
  friend l_itaylor operator+(int x, const l_itaylor& y);
  friend l_itaylor operator*(int x, const l_itaylor& y);
  friend l_itaylor operator/(int x, const l_itaylor& y);

  // operators +,-,*,/ for (l_real, l_itaylor):
  friend l_itaylor operator-(const l_real& x, const l_itaylor& y);
  friend l_itaylor operator+(const l_real& x, const l_itaylor& y);
  friend l_itaylor operator*(const l_real& x, const l_itaylor& y);
  friend l_itaylor operator/(const l_real& x, const l_itaylor& y);

  // operators +,-,*,/ for (l_itaylor, real):
  friend l_itaylor operator-(const l_itaylor& x, const real& y);
  friend l_itaylor operator+(const l_itaylor& x, const real& y);
  friend l_itaylor operator*(const l_itaylor& x, const real& y);
  friend l_itaylor operator/(const l_itaylor& x, const real& y);

  // operators +,-,*,/ for (l_itaylor, int):
  friend l_itaylor operator-(const l_itaylor& x, int y);
  friend l_itaylor operator+(const l_itaylor& x, int y);
  friend l_itaylor operator*(const l_itaylor& x, int y);
  friend l_itaylor operator/(const l_itaylor& x, int y);

  // operators +,-,*,/ for (l_itaylor, l_real):
  friend l_itaylor operator-(const l_itaylor& x, const l_real& y);
  friend l_itaylor operator+(const l_itaylor& x, const l_real& y);
  friend l_itaylor operator*(const l_itaylor& x, const l_real& y);
  friend l_itaylor operator/(const l_itaylor& x, const l_real& y);

  // Overloading the standard functions:
  friend l_itaylor sqr(const l_itaylor& x);
  friend l_itaylor sqrt(const l_itaylor& x);
  friend l_itaylor sqrt(const l_itaylor& x, int n);
  friend l_itaylor sqrt1px2(const l_itaylor& x);
  friend l_itaylor sqrtp1m1(const l_itaylor& x);
  friend l_itaylor sqrt1mx2(const l_itaylor& x);
  friend l_itaylor sqrtx2m1(const l_itaylor& x);

  friend l_itaylor pow(const l_itaylor& x, const l_interval& alpha);

  friend l_itaylor exp(const l_itaylor& x);
  //friend l_itaylor expm1(const l_itaylor& x);

  // Help function
  friend void f_g_u(const l_itaylor& f, const l_itaylor& g, 
                    const l_itaylor& u, int nb_function);


  friend l_itaylor ln(const l_itaylor& x);
  friend l_itaylor lnp1(const l_itaylor& x);

  friend l_itaylor sin(const l_itaylor& x);
  friend l_itaylor cos(const l_itaylor& x);
  friend l_itaylor tan(const l_itaylor& x);
  friend l_itaylor cot(const l_itaylor& x);

  friend l_itaylor sinh(const l_itaylor& x);
  friend l_itaylor cosh(const l_itaylor& x);
  friend l_itaylor tanh(const l_itaylor& x);
  friend l_itaylor coth(const l_itaylor& x);

  friend l_itaylor asin(const l_itaylor& x);
  friend l_itaylor acos(const l_itaylor& x);
  friend l_itaylor atan(const l_itaylor& x);
  friend l_itaylor acot(const l_itaylor& x);

  friend l_itaylor asinh(const l_itaylor& x);
  friend l_itaylor acosh(const l_itaylor& x);
  friend l_itaylor atanh(const l_itaylor& x);
  friend l_itaylor acoth(const l_itaylor& x);

};


l_itaylor var_l_itaylor(int ord, const real& c);
l_itaylor var_l_itaylor(int ord, const l_real& c);
l_itaylor var_l_itaylor(int ord, const interval& c);
l_itaylor var_l_itaylor(int ord, const l_interval& c);
l_itaylor const_l_itaylor(int ord, const real& c);
l_itaylor const_l_itaylor(int ord, const l_real& c);
l_itaylor const_l_itaylor(int ord, const interval& c);
l_itaylor const_l_itaylor(int ord, const l_interval& c);

int get_order(const l_itaylor& x);
l_ivector get_all_coef(const l_itaylor& x);
l_interval get_j_coef(const l_itaylor& x, int j);
l_interval get_j_derive(const l_itaylor& x, int j);
l_interval get_j_derivative(const l_itaylor& x, int j);

void print_l_itaylor(const l_itaylor& x);

l_itaylor operator-(const l_itaylor& x);
l_itaylor operator-(const l_itaylor& x, const l_itaylor& y);
l_itaylor operator+(const l_itaylor& x, const l_itaylor& y);
l_itaylor operator*(const l_itaylor& x, const l_itaylor& y);
l_itaylor operator/(const l_itaylor& x, const l_itaylor& y);

l_itaylor operator-(const interval& x, const l_itaylor& y);
l_itaylor operator+(const interval& x, const l_itaylor& y);
l_itaylor operator*(const interval& x, const l_itaylor& y);
l_itaylor operator/(const interval& x, const l_itaylor& y);

l_itaylor operator-(const l_interval& x, const l_itaylor& y);
l_itaylor operator+(const l_interval& x, const l_itaylor& y);
l_itaylor operator*(const l_interval& x, const l_itaylor& y);
l_itaylor operator/(const l_interval& x, const l_itaylor& y);

l_itaylor operator-(const l_itaylor& x, const interval& y);
l_itaylor operator+(const l_itaylor& x, const interval& y);
l_itaylor operator*(const l_itaylor& x, const interval& y);
l_itaylor operator/(const l_itaylor& x, const interval& y);

l_itaylor operator-(const l_itaylor& x, const l_interval& y);
l_itaylor operator+(const l_itaylor& x, const l_interval& y);
l_itaylor operator*(const l_itaylor& x, const l_interval& y);
l_itaylor operator/(const l_itaylor& x, const l_interval& y);

l_itaylor operator-(const real& x, const l_itaylor& y);
l_itaylor operator+(const real& x, const l_itaylor& y);
l_itaylor operator*(const real& x, const l_itaylor& y);
l_itaylor operator/(const real& x, const l_itaylor& y);

l_itaylor operator-(int x, const l_itaylor& y);
l_itaylor operator+(int x, const l_itaylor& y);
l_itaylor operator*(int x, const l_itaylor& y);
l_itaylor operator/(int x, const l_itaylor& y);

l_itaylor operator-(const l_real& x, const l_itaylor& y);
l_itaylor operator+(const l_real& x, const l_itaylor& y);
l_itaylor operator*(const l_real& x, const l_itaylor& y);
l_itaylor operator/(const l_real& x, const l_itaylor& y);

l_itaylor operator-(const l_itaylor& x, const real& y);
l_itaylor operator+(const l_itaylor& x, const real& y);
l_itaylor operator*(const l_itaylor& x, const real& y);
l_itaylor operator/(const l_itaylor& x, const real& y);

l_itaylor operator-(const l_itaylor& x, int y);
l_itaylor operator+(const l_itaylor& x, int y);
l_itaylor operator*(const l_itaylor& x, int y);
l_itaylor operator/(const l_itaylor& x, int y);

l_itaylor operator-(const l_itaylor& x, const l_real& y);
l_itaylor operator+(const l_itaylor& x, const l_real& y);
l_itaylor operator*(const l_itaylor& x, const l_real& y);
l_itaylor operator/(const l_itaylor& x, const l_real& y);

l_itaylor sqr(const l_itaylor& x);
l_itaylor sqrt(const l_itaylor& x);
l_itaylor sqrt(const l_itaylor& x, int n);
l_itaylor sqrt1px2(const l_itaylor& x);
l_itaylor sqrtp1m1(const l_itaylor& x);
l_itaylor sqrt1mx2(const l_itaylor& x);
l_itaylor sqrtx2m1(const l_itaylor& x);

l_itaylor pow(const l_itaylor& x, const l_interval& alpha);
l_itaylor exp(const l_itaylor& x);

void f_g_u(const l_itaylor& f, const l_itaylor& g, 
                    const l_itaylor& u, int nb_function);


l_itaylor ln(const l_itaylor& x);
l_itaylor lnp1(const l_itaylor& x);

l_itaylor sin(const l_itaylor& x);
l_itaylor cos(const l_itaylor& x);
l_itaylor tan(const l_itaylor& x);
l_itaylor cot(const l_itaylor& x);

l_itaylor sinh(const l_itaylor& x);
l_itaylor cosh(const l_itaylor& x);
l_itaylor tanh(const l_itaylor& x);
l_itaylor coth(const l_itaylor& x);

l_itaylor asin(const l_itaylor& x);
l_itaylor acos(const l_itaylor& x);
l_itaylor atan(const l_itaylor& x);
l_itaylor acot(const l_itaylor& x);

l_itaylor asinh(const l_itaylor& x);
l_itaylor acosh(const l_itaylor& x);
l_itaylor atanh(const l_itaylor& x);
l_itaylor acoth(const l_itaylor& x);

} // End of namespace taylor

#endif
