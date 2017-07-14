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

#ifndef _LDIM2TAYLOR_H
#define _LDIM2TAYLOR_H

#include <iostream>

#include "l_ivector.hpp"
#include "l_real.hpp"
#include "l_imath.hpp"
#include "idot.hpp"

using namespace cxsc;
namespace taylor{

enum{_ln,_lnp1,_tan,_cot,_asin,_acos,_atan,_acot,_tanh,_coth,
       _asinh,_acosh,_atanh,_acoth,_sqrtp1m1};

/*--------------------------------------------------------------------------

Klasse ldim2taylor: Taylorarithmetik fuer Funktionen in zwei Veraenderlichen
                    in staggered Arithmetic
         int p.....maximale Taylorordnung des Objekts
                   
  l_ivector* dat.....dynamischer Vektor von staggered Intervallvektoren
                   zur Speicherplatzersparnis als Dreiecks-
		   matrix implementiert

--------------------------------------------------------------------------*/




class  ldim2taylor{  
 private:
  
  int p; // Order of Taylor expansion
  l_ivector* dat; // pointer to a block (array) of elements of type
		  // l_ivector, the block is realized as a triangle matrix.
  
 public:

  ldim2taylor();     // Default constructor
  ldim2taylor(int);  // Constructor for special order (int)
  ldim2taylor(const ldim2taylor& ); // Copy constructor
  
  ~ldim2taylor();

  ldim2taylor& operator=(const ldim2taylor& );
  l_ivector& operator[](int n) const; 

  friend ldim2taylor init_var(int, int, const l_interval& );
  friend ldim2taylor init_const(int, const l_interval& );

  int get_p() const {return p;};

  void print_ldim2taylor(); // debug


// Overloadung the arithmetic operators:

  friend ldim2taylor operator-(const ldim2taylor& s);

  friend ldim2taylor operator-(const ldim2taylor&, const ldim2taylor& );
  friend ldim2taylor operator+(const ldim2taylor&, const ldim2taylor& );
  friend ldim2taylor operator*(const ldim2taylor&, const ldim2taylor& );
  friend ldim2taylor operator/(const ldim2taylor&, const ldim2taylor& );

  friend ldim2taylor operator-(const l_interval&, const ldim2taylor& );
  friend ldim2taylor operator+(const l_interval&, const ldim2taylor& );
  friend ldim2taylor operator*(const l_interval&, const ldim2taylor& );
  friend ldim2taylor operator/(const l_interval&, const ldim2taylor& );

  friend ldim2taylor operator-(const ldim2taylor&, const l_interval& );
  friend ldim2taylor operator+(const ldim2taylor&, const l_interval& );
  friend ldim2taylor operator*(const ldim2taylor&, const l_interval& );
  friend ldim2taylor operator/(const ldim2taylor&, const l_interval& );

  friend ldim2taylor operator-(const interval&, const ldim2taylor& );
  friend ldim2taylor operator+(const interval&, const ldim2taylor& );
  friend ldim2taylor operator*(const interval&, const ldim2taylor& );
  friend ldim2taylor operator/(const interval&, const ldim2taylor& );

  friend ldim2taylor operator-(const ldim2taylor&, const interval& );
  friend ldim2taylor operator+(const ldim2taylor&, const interval& );
  friend ldim2taylor operator*(const ldim2taylor&, const interval& );
  friend ldim2taylor operator/(const ldim2taylor&, const interval& );

  // Caution: possible conversion errors
  friend ldim2taylor operator-(const real&, const ldim2taylor& );
  friend ldim2taylor operator+(const real&, const ldim2taylor& );
  friend ldim2taylor operator*(const real&, const ldim2taylor& );
  friend ldim2taylor operator/(const real&, const ldim2taylor& );

  // Caution: possible conversion errors
  friend ldim2taylor operator-(const ldim2taylor&, const real& );
  friend ldim2taylor operator+(const ldim2taylor&, const real& );
  friend ldim2taylor operator*(const ldim2taylor&, const real& );
  friend ldim2taylor operator/(const ldim2taylor&, const real& );

  // Caution: possible conversion errors
  friend ldim2taylor operator-(const l_real&, const ldim2taylor& );
  friend ldim2taylor operator+(const l_real&, const ldim2taylor& );
  friend ldim2taylor operator*(const l_real&, const ldim2taylor& );
  friend ldim2taylor operator/(const l_real&, const ldim2taylor& );

  // Caution: possible conversion errors
  friend ldim2taylor operator-(const ldim2taylor&, const l_real& );
  friend ldim2taylor operator+(const ldim2taylor&, const l_real& );
  friend ldim2taylor operator*(const ldim2taylor&, const l_real& );
  friend ldim2taylor operator/(const ldim2taylor&, const l_real& );

  friend ldim2taylor operator-(int, const ldim2taylor& );
  friend ldim2taylor operator+(int, const ldim2taylor& );
  friend ldim2taylor operator*(int, const ldim2taylor& );
  friend ldim2taylor operator/(int, const ldim2taylor& );

  friend ldim2taylor operator-(const ldim2taylor&, int);
  friend ldim2taylor operator+(const ldim2taylor&, int);
  friend ldim2taylor operator*(const ldim2taylor&, int);
  friend ldim2taylor operator/(const ldim2taylor&, int);

  friend ldim2taylor sqr(const ldim2taylor& );
  friend ldim2taylor sqrt(const ldim2taylor& );
  friend ldim2taylor sqrt1px2(const ldim2taylor& );
  friend ldim2taylor sqrtx2m1(const ldim2taylor& );
  friend ldim2taylor sqrtp1m1(const ldim2taylor& );
  friend ldim2taylor power(const ldim2taylor& , const l_interval& );
  friend ldim2taylor power(const ldim2taylor& , int );

  friend ldim2taylor exp(const ldim2taylor& );
  friend ldim2taylor ln(const ldim2taylor& );
  friend ldim2taylor lnp1(const ldim2taylor& );

  friend ldim2taylor sin(const ldim2taylor& );
  friend ldim2taylor cos(const ldim2taylor& );
  friend ldim2taylor tan(const ldim2taylor& );
  friend ldim2taylor cot(const ldim2taylor& );

  friend ldim2taylor sinh(const ldim2taylor& );
  friend ldim2taylor cosh(const ldim2taylor& );
  friend ldim2taylor tanh(const ldim2taylor& );
  friend ldim2taylor coth(const ldim2taylor& );

  friend ldim2taylor asin(const ldim2taylor& );
  friend ldim2taylor acos(const ldim2taylor& );
  friend ldim2taylor atan(const ldim2taylor& );
  friend ldim2taylor acot(const ldim2taylor& );

  friend ldim2taylor asinh(const ldim2taylor& );
  friend ldim2taylor acosh(const ldim2taylor& );
  friend ldim2taylor atanh(const ldim2taylor& );
  friend ldim2taylor acoth(const ldim2taylor& );

  friend void f_g_u(const ldim2taylor& , const ldim2taylor& , const ldim2taylor& ,
			  int ); //war: int&, mg2005
};

ldim2taylor init_var(int, int, const l_interval& );
ldim2taylor init_const(int, const l_interval& );

ldim2taylor operator-(const ldim2taylor& s);
ldim2taylor operator-(const ldim2taylor&, const ldim2taylor& );
ldim2taylor operator+(const ldim2taylor&, const ldim2taylor& );
ldim2taylor operator*(const ldim2taylor&, const ldim2taylor& );
ldim2taylor operator/(const ldim2taylor&, const ldim2taylor& );

ldim2taylor operator-(const l_interval&, const ldim2taylor& );
ldim2taylor operator+(const l_interval&, const ldim2taylor& );
ldim2taylor operator*(const l_interval&, const ldim2taylor& );
ldim2taylor operator/(const l_interval&, const ldim2taylor& );

ldim2taylor operator-(const ldim2taylor&, const l_interval& );
ldim2taylor operator+(const ldim2taylor&, const l_interval& );
ldim2taylor operator*(const ldim2taylor&, const l_interval& );
ldim2taylor operator/(const ldim2taylor&, const l_interval& );

ldim2taylor operator-(const interval&, const ldim2taylor& );
ldim2taylor operator+(const interval&, const ldim2taylor& );
ldim2taylor operator*(const interval&, const ldim2taylor& );
ldim2taylor operator/(const interval&, const ldim2taylor& );

ldim2taylor operator-(const ldim2taylor&, const interval& );
ldim2taylor operator+(const ldim2taylor&, const interval& );
ldim2taylor operator*(const ldim2taylor&, const interval& );
ldim2taylor operator/(const ldim2taylor&, const interval& );

ldim2taylor operator-(const real&, const ldim2taylor& );
ldim2taylor operator+(const real&, const ldim2taylor& );
ldim2taylor operator*(const real&, const ldim2taylor& );
ldim2taylor operator/(const real&, const ldim2taylor& );

ldim2taylor operator-(const ldim2taylor&, const real& );
ldim2taylor operator+(const ldim2taylor&, const real& );
ldim2taylor operator*(const ldim2taylor&, const real& );
ldim2taylor operator/(const ldim2taylor&, const real& );

ldim2taylor operator-(const l_real&, const ldim2taylor& );
ldim2taylor operator+(const l_real&, const ldim2taylor& );
ldim2taylor operator*(const l_real&, const ldim2taylor& );
ldim2taylor operator/(const l_real&, const ldim2taylor& );

ldim2taylor operator-(const ldim2taylor&, const l_real& );
ldim2taylor operator+(const ldim2taylor&, const l_real& );
ldim2taylor operator*(const ldim2taylor&, const l_real& );
ldim2taylor operator/(const ldim2taylor&, const l_real& );

ldim2taylor operator-(int, const ldim2taylor& );
ldim2taylor operator+(int, const ldim2taylor& );
ldim2taylor operator*(int, const ldim2taylor& );
ldim2taylor operator/(int, const ldim2taylor& );

ldim2taylor operator-(const ldim2taylor&, int);
ldim2taylor operator+(const ldim2taylor&, int);
ldim2taylor operator*(const ldim2taylor&, int);
ldim2taylor operator/(const ldim2taylor&, int);

ldim2taylor sqr(const ldim2taylor& );
ldim2taylor sqrt(const ldim2taylor& );
ldim2taylor sqrt1px2(const ldim2taylor& );
ldim2taylor sqrtx2m1(const ldim2taylor& );
ldim2taylor sqrtp1m1(const ldim2taylor& );
ldim2taylor power(const ldim2taylor& , const l_interval& );
ldim2taylor power(const ldim2taylor& , int );

ldim2taylor exp(const ldim2taylor& );
ldim2taylor ln(const ldim2taylor& );
ldim2taylor lnp1(const ldim2taylor& );

ldim2taylor sin(const ldim2taylor& );
ldim2taylor cos(const ldim2taylor& );
ldim2taylor tan(const ldim2taylor& );
ldim2taylor cot(const ldim2taylor& );

ldim2taylor sinh(const ldim2taylor& );
ldim2taylor cosh(const ldim2taylor& );
ldim2taylor tanh(const ldim2taylor& );
ldim2taylor coth(const ldim2taylor& );

ldim2taylor asin(const ldim2taylor& );
ldim2taylor acos(const ldim2taylor& );
ldim2taylor atan(const ldim2taylor& );
ldim2taylor acot(const ldim2taylor& );

ldim2taylor asinh(const ldim2taylor& );
ldim2taylor acosh(const ldim2taylor& );
ldim2taylor atanh(const ldim2taylor& );
ldim2taylor acoth(const ldim2taylor& );

void f_g_u(const ldim2taylor& , const ldim2taylor& , const ldim2taylor& , int ); 


/*--------------------------------------------------------------------------

Class ldim2taylor_vector: For definition and initialization of objects
                         of the class ldim2taylor respective to the two
                         independent variables x and y.

--------------------------------------------------------------------------*/
class ldim2taylor_vector{  // Vector for independent variables
 private:
  
  int dim; // Dimension of the vector
  int lb;
  int ub;

  int p_el;
  ldim2taylor* comp; // Pointer to a dyn. array of elements of type ldim2taylor
  
 public:
    // Constructors:
    ldim2taylor_vector(); // Default constructor
    ldim2taylor_vector(int, int, int);  // order, lb, ub;
    ldim2taylor_vector(const ldim2taylor_vector& );
  
    ~ldim2taylor_vector();
    // Member functions
    int get_dim()  const {return dim;};   // Number of elements 
                                          // of type ldim2taylor
    int get_lb ()  const {return lb;};    // lb index of the array
    int get_ub ()  const {return ub;};    // ub index of the array
    int get_p_el() const {return p_el;};  // p_el: Order of the
                                          // Taylor expansion; 
    // Operators
    ldim2taylor_vector& operator=(const ldim2taylor_vector& );
    ldim2taylor& operator[](int n) const;  

    // Function for initialization of objects of type ldim2taylor_vector:
    friend ldim2taylor_vector init_var(int order, l_ivector& values);

};

ldim2taylor_vector init_var(int order, l_ivector& values);

} // end of namespace taylor

#endif
