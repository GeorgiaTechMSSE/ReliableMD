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

//////////////////////////////////////////////////////////////
//
//       Implementation of class l_itaylor in litaylor.cpp
//
//////////////////////////////////////////////////////////////

#include "litaylor.hpp"

///////////////////////////////////////////////////////////////
//
//                     class l_itaylor
//
///////////////////////////////////////////////////////////////

namespace taylor {

l_ivector l_itaylor::faks(1);
int l_itaylor::initialized=0;

void l_itaylor::initialize()
{
 Resize(l_itaylor::faks,0,170);
 l_itaylor::faks[0]=l_interval(1.0);
 l_itaylor::faks[1]=l_interval(1.0);

 for(int i=2; i<=170; i++)
   l_itaylor::faks[i]=l_itaylor::faks[i-1]*l_interval(i);
}

//-----------------------------------------------------------------------

// Constructors:

l_itaylor::l_itaylor()
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
}

//-----------------------------------------------------------------------

l_itaylor::l_itaylor(const l_itaylor& s)
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
 p=s.p;
 Resize(tayl,0,p);
 tayl=s.tayl;
}

//-----------------------------------------------------------------------

l_itaylor::l_itaylor(int order)
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
 if(order<0 || order>170)
  {
    std::cerr << "l_itaylor::l_itaylor: incorrect order! 0<=order<=170"
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
}


//-----------------------------------------------------------------------

l_itaylor::l_itaylor(int order, const real& value)
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
 if(order<0 || order>170)
  {
    std::cerr << "l_itaylor::l_itaylor: incorrect order! 0<=order<=170"
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 tayl[0] = l_interval(value);
 if(p>0)
   {
     tayl[1]=interval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=interval(0.0);
   }
}

//-----------------------------------------------------------------------

l_itaylor::l_itaylor(int order, const l_real& value)
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
 if(order<0 || order>170)
  {
    std::cerr << "l_itaylor::l_itaylor:  incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 tayl[0] = l_interval(value);
 if(p>0)
   {
     tayl[1]=interval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=interval(0.0);
   }
}

//-----------------------------------------------------------------------

l_itaylor::l_itaylor(int order, const interval& value)
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "l_itaylor::l_itaylor: incorrect order! 0<=order<=170"
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 tayl[0]=l_interval(value);
 if(p>0)
   { 
     tayl[1]=interval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=interval(0.0);
   }
}

//-----------------------------------------------------------------------

l_itaylor::l_itaylor(int order, const l_interval& value)
{
 if(!l_itaylor::initialized){l_itaylor::initialize();l_itaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "l_itaylor::l_itaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 tayl[0]=value;
 if(p>0)
   { 
     tayl[1]=interval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=interval(0.0);
   }
}

// Functions for initialization of independent variables:
l_itaylor var_l_itaylor(int ord, const real& x)
{
 l_itaylor erg(ord,x);
 return erg;
}

//-----------------------------------------------------------------------

l_itaylor var_l_itaylor(int ord, const l_real& x)
{
 l_itaylor erg(ord,x);
 return erg;
}

//-----------------------------------------------------------------------

l_itaylor var_l_itaylor(int ord, const interval& x)
{
 l_itaylor erg(ord,x);
 return erg;
}

//-----------------------------------------------------------------------

l_itaylor var_l_itaylor(int ord, const l_interval& x)
{
 l_itaylor erg(ord,x);
 return erg;
}

//-----------------------------------------------------------------------

// Functions for initialization of constants:
l_itaylor const_l_itaylor(int ord, const real& c)
{
 l_itaylor erg(ord);
 erg.tayl[0]=l_interval(c);
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=interval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

l_itaylor const_l_itaylor(int ord, const l_real& c)
{
 l_itaylor erg(ord);
 erg.tayl[0]=l_interval(c);
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=interval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

l_itaylor const_l_itaylor(int ord, const interval& c)
{
 l_itaylor erg(ord);
 erg.tayl[0]=l_interval(c);
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=interval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

l_itaylor const_l_itaylor(int ord, const l_interval& c)
{
 l_itaylor erg(ord);
 erg.tayl[0]=c;
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=interval(0.0);
 return erg;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

// assignment operators:

l_itaylor l_itaylor::operator=(const l_itaylor& s)
{
    p=s.p;  // p: order of the left part; s.p: order of the right part; 
    Resize(tayl,0,s.p);
    tayl=s.tayl;
    return *this; // The left operand gets the order of the right operand!
}

//-----------------------------------------------------------------------

l_itaylor l_itaylor::operator=(int c)
{
    tayl[0] = l_interval(c);
    for (int j=1; j<=p; j++) tayl[j]=0.0;
    return *this;
}

//-----------------------------------------------------------------------

l_itaylor l_itaylor::operator=(const real& c)
{
    tayl[0] = l_interval(c);
    for (int j=1; j<=p; j++) tayl[j]=0.0;
    return *this;
}

//-----------------------------------------------------------------------

l_itaylor l_itaylor::operator=(const l_real& c)
{
    tayl[0] = l_interval(c);
    for (int j=1; j<=p; j++) tayl[j]=0.0;
    return *this;
}

//-----------------------------------------------------------------------

l_itaylor l_itaylor::operator=(const interval& c)
{
    tayl[0] = l_interval(c);
    for (int j=1; j<=p; j++) tayl[j]=0.0;
    return *this;
}

//-----------------------------------------------------------------------

l_itaylor l_itaylor::operator=(const l_interval& c)
{
    tayl[0] = c;
    for (int j=1; j<=p; j++) tayl[j]=0.0;
    return *this;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

// class components:

// returning the maximal order
int get_order(const l_itaylor& x)
{
 return x.p;
}

//-----------------------------------------------------------------------

// returning all Taylor-coefficients by an l_interval vector
l_ivector get_all_coef(const l_itaylor& x)
{
 return x.tayl;
}

//-----------------------------------------------------------------------

// returning Taylor-coefficient of order j
l_interval get_j_coef(const l_itaylor& x, int j)
{
 return x.tayl[j];
}

//-----------------------------------------------------------------------

// returning derivative of order j,  j <= 170;
l_interval get_j_derive(const l_itaylor& x, int j)
{
  return x.tayl[j]*l_itaylor::faks[j];
}

l_interval get_j_derivative(const l_itaylor& x, int j)
{
  return x.tayl[j]*l_itaylor::faks[j];
}

//-----------------------------------------------------------------------

// Output of all taylor coefficients 
void print_l_itaylor(const l_itaylor& x)
{
 std::cerr <<"Ausgabe l_itaylor der Ordnung " << x.p << " " << std::endl;
 for(int i=Lb(x.tayl); i<=Ub(x.tayl);i++) 
  {
   std::cerr << "i  " << i << "  component: " << x.tayl[i] << std::endl;
  };
 std::cerr << std::endl;
}

//-----------------------------------------------------------------------

// Overloading of operators:

//-----------------------------------------------------------------------

// - operator:
l_itaylor operator-(const l_itaylor& x)
{
 int order=get_order(x);
 l_itaylor erg(order);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= -x.tayl[j];
 return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:  +,-,*,/  for (l_itaylor, l_itaylor):
// All operands are independent variables.
l_itaylor operator-(const l_itaylor& x, const l_itaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in l_itaylor, operator - : different orders " 
             << std::endl;
   exit(1);
  };

 l_itaylor erg(order1);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= x.tayl[j]-y.tayl[j];
 return erg; 
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_itaylor& x, const l_itaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in l_itaylor, operator + : different orders " 
             << std::endl;
   exit(1);
  };

 l_itaylor erg(order1);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= x.tayl[j]+y.tayl[j];
 return erg; 
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_itaylor& x, const l_itaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in l_itaylor, operator * : different orders " 
             << std::endl;
   exit(1);
  };

 l_itaylor erg(order1);
 l_interval sum; 
 idotprecision sum_idot; // for accumulate(...), scalar product

 for(int j=0; j<=order1; j++) 
 {
  sum_idot=interval(0);
  for(int i=0; i<=j; i++)
   {
    accumulate(sum_idot, x.tayl[i],y.tayl[j-i]);
   }
  sum = sum_idot;
  erg.tayl[j]= sum;
 }
 return erg; 
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_itaylor& x, const l_itaylor& y)
{
 int order1(get_order(x));
 int order2(get_order(y));
 if(order1 != order2) 
  {
   std::cerr << "Error in l_itaylor, operator / : different orders " 
             << std::endl;
   exit(1);
  };
 
 if(0 <= y.tayl[0]) 
  {
   std::cerr << "Error in l_itaylor, operator / : 0 in interval " << std::endl;
   exit(1);
  };

 l_itaylor erg(order1);
 l_interval sum; 
 idotprecision sum_idot; // for accumulate(...), scalar product

 for(int j=0; j<=order1; j++) 
 {
  sum_idot=interval(0);
  for(int i=1; i<=j; i++)
   {
    accumulate(sum_idot, y.tayl[i],erg.tayl[j-i]);
   }
  sum = sum_idot;
  erg.tayl[j]= (x.tayl[j]-sum)/y.tayl[0];
 }
 return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (interval, l_itaylor):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const interval& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const interval& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const interval& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const interval& x, const l_itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in l_itaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    idotprecision idot;
    l_interval sum;
    int order(get_order(y));
    l_itaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	sum = idot;
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}


// Operators with two operands:   +,-,*,/  for (l_interval, l_itaylor):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_interval& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_interval& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_interval& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_interval& x, const l_itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in l_itaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    idotprecision idot;
    l_interval sum;
    int order(get_order(y));
    l_itaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	sum = idot;
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (l_itaylor, interval):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_itaylor& x, const interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_itaylor& x, const interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_itaylor& x, const interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_itaylor& x, const interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}


// Operators with two operands:   +,-,*,/  for (l_itaylor, l_interval):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_itaylor& x, const l_interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_itaylor& x, const l_interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_itaylor& x, const l_interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_itaylor& x, const l_interval& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (real, l_itaylor):
// The operand of type real is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const real& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const real& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const real& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/( const real& x, const l_itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in l_itaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    idotprecision idot;
    l_interval sum;
    int order(get_order(y));
    l_itaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	sum = idot;  // rnd(idot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

// Operators with two operands:   +,-,*,/  for (int, l_itaylor):
// The operand of type real is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(int x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = -y;
    erg.tayl[0] = interval(x) - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(int x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = y;
    erg.tayl[0] = interval(x) + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(int x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = interval(x)*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(int x, const l_itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in l_itaylor, operator / : 0 in interval "
                  << std::endl;
	exit(1);
    };
    idotprecision idot;
    l_interval sum;
    int order(get_order(y));
    l_itaylor w(order);
    w.tayl[0] = interval(x) / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++)
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	sum = idot;  // rnd(idot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}


//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (l_real, l_itaylor):
// The operand of type real is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_real& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_real& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_real& x, const l_itaylor& y)
{
    int order(get_order(y));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_real& x, const l_itaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in l_itaylor, operator / : 0 in interval "
                  << std::endl;
	exit(1);
    };
    idotprecision idot;
    l_interval sum;
    int order(get_order(y));
    l_itaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	idot=0;
	for (int j=1; j<=k; j++)
	    accumulate(idot,y.tayl[j],w.tayl[k-j]);
	sum = idot;  // rnd(idot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}


//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (l_itaylor, real):
// The operand of type real is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_itaylor& x, const real& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_itaylor& x, const real& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_itaylor& x, const real& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_itaylor& x, const real& y)
{
    if (y==0)
    {
	std::cerr << "Error in l_itaylor: division by 0"
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

// Operators with two operands:   +,-,*,/  for (l_itaylor, int):
// The operand of type real is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_itaylor& x, int y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - interval(y);
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_itaylor& x, int y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + interval(y);
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_itaylor& x, int y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*interval(y);
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_itaylor& x, int y)
{
    if (y==0)
    {
	std::cerr << "Error in l_itaylor: division by 0"
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/interval(y);
    return erg;
}


//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (l_itaylor, l_real):
// The operand of type real is assumed to be a constant and not
// an independent variable!

l_itaylor operator-(const l_itaylor& x, const l_real& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator+(const l_itaylor& x, const l_real& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator*(const l_itaylor& x, const l_real& y)
{
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

l_itaylor operator/(const l_itaylor& x, const l_real& y)
{
    if (y==0)
    {
	std::cerr << "Error in l_itaylor: division by 0"
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    l_itaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

//-----------------------------------------------------------------------

// Overloading the standard functions:

//-----------------------------------------------------------------------

// sqr-function
l_itaylor sqr(const l_itaylor& x)
{
    idotprecision idot;
    l_interval sum;
    int order(get_order(x)),m;
    l_itaylor erg(order);

    erg.tayl[0] = sqr(x.tayl[0]);
    for (int k=1; k<=order; k++)
    {
	m = (k+1) / 2;
	idot = 0;
	for (int j=0; j<=m-1; j++)
	    accumulate(idot,x.tayl[j],x.tayl[k-j]);
	sum = idot;
	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(x.tayl[m]); // k even
    }
    return erg;
}

//-----------------------------------------------------------------------

// Square-root

l_itaylor sqrt(const l_itaylor& x)
{
    idotprecision idot;
    l_interval sum,h;
    int order(get_order(x)),m;
    l_itaylor erg(order);
    if (0<=x.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, sqrt: 0 in interval" << std::endl;
	exit(1);
    };
    erg.tayl[0] = sqrt(x.tayl[0]);
    h = erg.tayl[0];
    times2pown(h,1);
    for (int k=1; k<=order; k++)
    {
	m = (k+1) / 2;
	idot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	sum = idot;
	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]) / h;
    }
    return erg;
}


//-----------------------------------------------------------------------

// sqrt(x,n)

l_itaylor sqrt(const l_itaylor& x, int n)
{
    int order(get_order(x));
    l_itaylor erg(order);
    if (0<=x.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, sqrt(x,n): 0 in interval" 
                  << std::endl;
	exit(1);
    };
    erg.tayl[0] = sqrt(x.tayl[0],n); // element No. 0
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=0; j<=k-1; j++)
	    erg.tayl[k] += (l_interval(k-j)/real(n)-interval(j))
		           * erg.tayl[j] * x.tayl[k-j];
	erg.tayl[k] /= (interval(k)*x.tayl[0]);
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(1+x^2)
l_itaylor sqrt1px2(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg(order);
    const real c = 500.0;

    if (Inf(x.tayl[0]) > c) erg = x*sqrt(real(1)+real(1)/sqr(x));
    else if (Sup(x.tayl[0]) < -c) erg = -x*sqrt(real(1)+real(1)/sqr(x)); 
    else erg = sqrt(real(1)+sqr(x));  
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(1+x)-1
l_itaylor sqrtp1m1(const l_itaylor& x)
{
    int order(get_order(x)),m;
    l_itaylor erg(order);
    idotprecision idot;
    l_interval h,Ne;

    if (Inf(x.tayl[0])<=-1)
    {
	std::cerr << "Error in l_itaylor, sqrtp1m1: wrong argument" 
                  << std::endl;
	exit(1);
    };

    erg.tayl[0] = sqrtp1m1(x.tayl[0]);
    Ne = real(1.0)+erg.tayl[0];
    times2pown(Ne,1);

    for (int k=1; k<=order; k++)
    {
	m = (k+1)/2;
	idot = 0;
	for (int j=1; j<=m-1; j++)
	    accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	h = idot;
	times2pown(h,1);
	erg.tayl[k] = h;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]);
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]) / Ne; 
    }
    return erg;
} // sqrtp1m1


//-----------------------------------------------------------------------

// sqrt(1-x^2):
l_itaylor sqrt1mx2(const l_itaylor& x)
{
    idotprecision idot;
    l_interval sum,h;
    int order(get_order(x)),m;
    l_itaylor erg(order), g(order);

    if (Inf(x.tayl[0])<=-1 || Sup(x.tayl[0])>=1)
    {
	std::cerr << "Error in l_itaylor, sqrt1mx2: wrong argument" 
                  << std::endl;
	exit(1);
    };
    erg.tayl[0]=sqrt1mx2(x.tayl[0]);  
    h = real(-1)/erg.tayl[0];
    times2pown(h,-1);
    g = sqr(x);
    for (int k=1; k<=order; k++)
    {
	m = (k+1)/2;
	idot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	sum = idot;
	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (g.tayl[k]+erg.tayl[k])*h;
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(x^2-1):
l_itaylor sqrtx2m1(const l_itaylor& x)
{
    const real c = 30.0; 
    idotprecision idot;
    l_interval sum,h;
    interval xi;
    int order(get_order(x)),m;
    l_itaylor erg(order), g(order);

    xi = x.tayl[0];
    if (Disjoint(xi,interval(-1,1))==0)
    {
	std::cerr << "Error in l_itaylor, sqrtx2m1: wrong argument" 
                  << std::endl;
	exit(1);
    };

    if (Inf(x.tayl[0])>c) erg = x*sqrt1mx2(real(1)/x);
    else if (Sup(x.tayl[0])<-c) erg = -x*sqrt1mx2(real(1)/x);
    else {
	erg.tayl[0]=sqrtx2m1(x.tayl[0]); 
	g = sqr(x);
	h = real(1)/erg.tayl[0];
	times2pown(h,-1);
	for (int k=1; k<=order; k++)
	{
	    m = (k+1)/2;
	    idot = 0;
	    for (int j=1; j<=m-1; j++) 
		accumulate(idot,erg.tayl[j],erg.tayl[k-j]);
	    sum = idot;
	    times2pown(sum,1);  // Multiplication with 2
	    erg.tayl[k] = sum;
	    if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	    erg.tayl[k] = (g.tayl[k]-erg.tayl[k])*h;
	}
    }
    return erg;
}

//-----------------------------------------------------------------------

// power-function
l_itaylor pow(const l_itaylor& x, const l_interval& alpha)
{
    int order(get_order(x));
    l_itaylor erg(order);

    if (0<=x.tayl[0])
    {
	std::cerr << "Error in l_itaylor, pow(x,a): 0 in interval x" 
                  << std::endl;
	exit(1);
    };

    erg.tayl[0] = pow(x.tayl[0],alpha); // element No. 0
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=0; j<=k-1; j++)
	    erg.tayl[k] += (interval(k-j)*alpha-interval(j))
		           * erg.tayl[j] * x.tayl[k-j];
	erg.tayl[k] /= (interval(k)*x.tayl[0]);
    }
    return erg;
}


//-----------------------------------------------------------------------

// Exponential-function
l_itaylor exp(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg(order);

    erg.tayl[0] = exp(x.tayl[0]); // element No. 0; function value
    for(int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for(int j=0; j<=k-1; j++)
	    erg.tayl[k] += interval(k-j)*erg.tayl[j]*x.tayl[k-j]; 
	erg.tayl[k] /= interval(k);
    }
    return erg; 
}

//-----------------------------------------------------------------------

// exp(x)-1;
/*
l_itaylor expm1(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg(order);

    erg.tayl[0] = exp(x.tayl[0]); // element No. 0; function value
    for(int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for(int j=0; j<=k-1; j++)
	    erg.tayl[k] += interval(k-j)*erg.tayl[j]*x.tayl[k-j]; 
	erg.tayl[k] /= interval(k);
    }
    erg.tayl[0] = expm1(x.tayl[0]); 
    return erg; 
}
*/

//-----------------------------------------------------------------------

// Help function
void f_g_u(const l_itaylor& f, const l_itaylor& g, const l_itaylor& u, int nb_function)
{
 int order1=get_order(f);
 int order2=get_order(g);
 int order3=get_order(u);

 // The following errors should be caught before 
 // but for security here again:
 if(order1 != order2) 
  {
   std::cerr << "Error1 in f_g_u: different orders " 
             << std::endl;
   exit(1);
  };

 if(order3 != order2) 
  {
   std::cerr << "Error2 in f_g_u: different orders " << std::endl;
   exit(1);
  };

 if(0 <= g.tayl[0]) 
  {
   std::cerr << "Error in f_g_u : wrong argument " << std::endl;
   exit(1);
  }; 

 switch(nb_function) // element No. 0
   {
    case _i_ln:f.tayl[0]=ln(u.tayl[0]); break;

    case _i_tan:f.tayl[0]=tan(u.tayl[0]); break;
    case _i_cot:f.tayl[0]=cot(u.tayl[0]); break;

    case _i_asin:f.tayl[0]=asin(u.tayl[0]); break;
    case _i_acos:f.tayl[0]=acos(u.tayl[0]); break;
    case _i_atan:f.tayl[0]=atan(u.tayl[0]); break;
    case _i_acot:f.tayl[0]=acot(u.tayl[0]); break;

    case _i_tanh:f.tayl[0]=tanh(u.tayl[0]); break;
    case _i_coth:f.tayl[0]=coth(u.tayl[0]); break; 

    case _i_asinh:f.tayl[0]=asinh(u.tayl[0]); break;
    case _i_acosh:f.tayl[0]=acosh(u.tayl[0]); break;
    case _i_atanh:f.tayl[0]=atanh(u.tayl[0]); break;

    case _i_acoth:f.tayl[0]=acoth(u.tayl[0]); break;

   }
 
 // remaining elements:
 l_interval sum; 
 for(int j=1; j<=Ub(f.tayl); j++) 
 {
     sum = l_interval(0.0);
     for(int i=1; i<=j-1; i++)
	 sum += interval(i)*f.tayl[i]*g.tayl[j-i];
     f.tayl[j] = (u.tayl[j]-sum/interval(j)) / g.tayl[0];
 }
} // f_g_u

//-----------------------------------------------------------------------

// Logarithm-function
l_itaylor ln(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);

    f_g_u(f,x,x,_i_ln);
  
    return f;
}

//-----------------------------------------------------------------------

// ln(1+x)
l_itaylor lnp1(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg(order), g(order);

    g = interval(1) + x;
    erg.tayl[0] = lnp1(x.tayl[0]);
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=1; j<=k-1; j++)
	    erg.tayl[k] += interval(j) * erg.tayl[j] * g.tayl[k-j];
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]/interval(k)) / g.tayl[0];
    }
    return erg;
}

//-----------------------------------------------------------------------

// Sinus-function
l_itaylor sin(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg1(order);   // sin
    l_itaylor erg2(order);   // cos
    l_interval s1,s2;

    erg1.tayl[0]=sin(x.tayl[0]); // Element No. 0:  erg1 (sin)
    erg2.tayl[0]=cos(x.tayl[0]); // Element No. 0:  erg2 (cos)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= -s2/interval(j);
    }
    return erg1; 
}

//-----------------------------------------------------------------------

// Cosinus-function
l_itaylor cos(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg1(order);   // sin
    l_itaylor erg2(order);   // cos
    l_interval s1,s2;

    erg1.tayl[0]=sin(x.tayl[0]); // Element No. 0:  erg1 (sin)
    erg2.tayl[0]=cos(x.tayl[0]); // Element No. 0:  erg2 (cos)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0.0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= -s2/interval(j);
    }
 return erg2; 
}

//-----------------------------------------------------------------------

// Tangens-function
l_itaylor tan(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=sqr(cos(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, tan : wrong argument" << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_tan);
  
    return f;
}

//-----------------------------------------------------------------------

// Cotangens-function
l_itaylor cot(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=-sqr(sin(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, cot : wrong argument" << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_cot);
  
    return f;
}

//-----------------------------------------------------------------------

// Sinushyperbolicus-function
l_itaylor sinh(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg1(order);  // sinh
    l_itaylor erg2(order);  // cosh
    l_interval s1,s2;

    erg1.tayl[0]=sinh(x.tayl[0]); // element No. 0:  erg1 (sinh)
    erg2.tayl[0]=cosh(x.tayl[0]); // element No. 0:  erg2 (cosh)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= s2/interval(j);
    }
    return erg1; 
}

//-----------------------------------------------------------------------

// Cosinushyperbolicus-function
l_itaylor cosh(const l_itaylor& x)
{
    int order(get_order(x));
    l_itaylor erg1(order); // sinh
    l_itaylor erg2(order); // cosh
    l_interval s1,s2;

    erg1.tayl[0]=sinh(x.tayl[0]); // element No. 0:  erg1 (sinh)
    erg2.tayl[0]=cosh(x.tayl[0]); // element No. 0:  erg2 (cosh)

    // remaining elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=interval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= s2/interval(j);
    }
    return erg2; 
}

//-----------------------------------------------------------------------

//Tangenshyperbolicus-function
l_itaylor tanh(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=sqr(cosh(x)); 
 
    f_g_u(f,g,x,_i_tanh);
  
    return f;
}

//-----------------------------------------------------------------------
 
// Cotangenshyperbolicus-function
l_itaylor coth(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=-sqr(sinh(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, coth : wrong argument " 
                  << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_coth);
  
    return f;
}

//-----------------------------------------------------------------------

//Arcsinusfunktion
l_itaylor asin(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=sqrt1mx2(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, asin : wrong argument " 
                  << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_asin);
  
    return f;
}

//-----------------------------------------------------------------------

//Arccosinusfunktion
l_itaylor acos(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g = -sqrt1mx2(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, acos : wrong argument " 
                  << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_acos);
  
    return f;
}

//-----------------------------------------------------------------------

//Arctan-function
l_itaylor atan(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=interval(1.0)+sqr(x);

    f_g_u(f,g,x,_i_atan);
  
    return f;
}

//-----------------------------------------------------------------------

//Arccotan-function
l_itaylor acot(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=-(interval(1.0)+sqr(x));

    f_g_u(f,g,x,_i_acot);
  
    return f;
}

//-----------------------------------------------------------------------

//Areasinh-function
l_itaylor asinh(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=sqrt1px2(x);

    f_g_u(f,g,x,_i_asinh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areacosh-function
l_itaylor acosh(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=sqrtx2m1(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, acosh : wrong argument " 
                  << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_acosh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areatanh-function
l_itaylor atanh(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=interval(1.0)-sqr(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, atanh : wrong argument " 
                  << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_atanh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areacotanh-function
l_itaylor acoth(const l_itaylor& x)
{
    int order=get_order(x);
    l_itaylor f(order);
    l_itaylor g(order);

    g=interval(1.0)-sqr(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in l_itaylor, acoth : wrong argument " 
                  << std::endl;
	exit(1);
    };  

    f_g_u(f,g,x,_i_acoth);
  
    return f;
}

} // end namespace taylor
