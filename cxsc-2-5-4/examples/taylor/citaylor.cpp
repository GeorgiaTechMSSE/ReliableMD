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
//       Implementation of class citaylor in citaylor.cpp     
//
//////////////////////////////////////////////////////////////

#include "citaylor.hpp"

///////////////////////////////////////////////////////////////
//
//                     class citaylor
//
///////////////////////////////////////////////////////////////

namespace taylor {

ivector citaylor::faks(1);
int citaylor::initialized=0;

void citaylor::initialize()
{
 Resize(citaylor::faks,0,170);
 citaylor::faks[0]=interval(1.0);
 citaylor::faks[1]=interval(1.0);

 for(int i=2; i<=170; i++)
   citaylor::faks[i]=citaylor::faks[i-1]*interval(i);
}

//-----------------------------------------------------------------------

// Constructors:

citaylor::citaylor()
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
}

//-----------------------------------------------------------------------

citaylor::citaylor(int order)
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "citaylor::citaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
}

//-----------------------------------------------------------------------

citaylor::citaylor(const citaylor& s)
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
 p=s.p;
 Resize(tayl,0,p);
 tayl=s.tayl;
}

//-----------------------------------------------------------------------

citaylor::citaylor(int order, const real& value)
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "citaylor::citaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 
 cinterval cinterval_value = cinterval(value);
 tayl[0] = cinterval_value;
 if(p>0)
   {
     tayl[1] = cinterval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i] = cinterval(0.0);
   }
}

//-----------------------------------------------------------------------

citaylor::citaylor(int order, const complex& value)
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "citaylor::citaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 
 cinterval cinterval_value=cinterval(value);
 tayl[0]=cinterval_value;
 if(p>0)
   {
     tayl[1]=cinterval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=cinterval(0.0);
   }
}

//-----------------------------------------------------------------------

citaylor::citaylor(int order, const interval& value)
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "citaylor::citaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 tayl[0] = value;
 if(p>0)
   { 
     tayl[1] = cinterval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i] = cinterval(0.0);
   }
}

//-----------------------------------------------------------------------

citaylor::citaylor(int order, const cinterval& value)
{
 if(!citaylor::initialized){citaylor::initialize();citaylor::initialized=1;};
 if(order<0 || order>170) 
  {
    std::cerr << "citaylor::citaylor: incorrect order! 0<=order<=170" 
              << std::endl;
    exit(1);
  }
 p=order;
 Resize(tayl,0,p);
 tayl[0]=value;
 if(p>0)
   { 
     tayl[1]=cinterval(1.0);
     for(int i=2; i<=Ub(tayl);i++) tayl[i]=cinterval(0.0);
   }
}

//-----------------------------------------------------------------------

// Functions for initialization of independent variables:

citaylor var_citaylor(int ord, const real& x)
{
    citaylor erg(ord,x);
    return erg;
}

//-----------------------------------------------------------------------
citaylor var_citaylor(int ord, const complex& x)
{
    citaylor erg(ord,x);
    return erg;
}

//-----------------------------------------------------------------------

citaylor var_citaylor(int ord, const cinterval& x)
{
    citaylor erg(ord,x);
    return erg;
}

//-----------------------------------------------------------------------

citaylor var_citaylor(int ord, const interval& x)
{
    citaylor erg(ord,x);
    return erg;
}
//-----------------------------------------------------------------------


// Functions for initialization of constants:
citaylor const_citaylor(int ord, const real& c)
{
 citaylor erg(ord);
 erg.tayl[0] = cinterval(c);
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i] = cinterval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

citaylor const_citaylor(int ord, const complex& c)
{
 citaylor erg(ord);
 erg.tayl[0]=cinterval(c);
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=cinterval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

citaylor const_citaylor(int ord, const interval& c)
{
 citaylor erg(ord);
 erg.tayl[0] = c;
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i] = cinterval(0.0);
 return erg;
}

//-----------------------------------------------------------------------

citaylor const_citaylor(int ord, const cinterval& c)
{
 citaylor erg(ord);
 erg.tayl[0]=c;
 for(int i=1; i<=Ub(erg.tayl);i++) erg.tayl[i]=cinterval(0.0);
 return erg;
}
//-----------------------------------------------------------------------



//-----------------------------------------------------------------------
// assignment operators
//-----------------------------------------------------------------------

citaylor citaylor::operator=(const citaylor& s)
{
 p=s.p;
 Resize(tayl,0,s.p);
 tayl=s.tayl;
 return *this;
}

citaylor citaylor::operator=(int n)
{
 tayl[0] = cinterval(n);
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

citaylor citaylor::operator=(const real& r)
{
 tayl[0] = cinterval(r);
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

citaylor citaylor::operator=(const complex& x)
{
 tayl[0] = cinterval(x);
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

citaylor citaylor::operator=(const interval& x)
{
 tayl[0] = x;
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

citaylor citaylor::operator=(const cinterval& x)
{
 tayl[0] = x;
 for (int j=1; j<=p; j++) tayl[j]=0.0;
 return *this;
}

//-----------------------------------------------------------------------

// class components:

// returning the maximal order
int get_order(const citaylor& x)
{
 return x.p;
}

//-----------------------------------------------------------------------

// returning all Taylor-coefficients by an interval vector
civector get_all_coef(const citaylor& x)
{
 return x.tayl;
}

//-----------------------------------------------------------------------

// returning Taylor-coefficient of order j
cinterval get_j_coef(const citaylor& x, int j)
{
 return x.tayl[j];
}

//-----------------------------------------------------------------------

// returning derivative of order j,  j <= 170;
cinterval get_j_derive(const citaylor& x, int j)
{
  return x.tayl[j]*citaylor::faks[j];
}

// returning derivative of order j,  j <= 170;
cinterval get_j_derivative(const citaylor& x, int j)
{
  return x.tayl[j]*citaylor::faks[j];
}

//-----------------------------------------------------------------------

// Output of all taylor coefficients 
void print_citaylor(const citaylor& x)
{
 std::cerr <<"Output citaylor of order " << x.p << " " << std::endl;
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
citaylor operator-(const citaylor& x)
{
 int order=get_order(x);
 citaylor erg(order);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= -x.tayl[j];
 return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:  +,-,*,/  for (citaylor, citaylor):
// All operands are independent variables.
citaylor operator-(const citaylor& x, const citaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in citaylor, operator - : different orders " 
             << std::endl;
   exit(1);
  };

 citaylor erg(order1);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= x.tayl[j]-y.tayl[j];
 return erg; 
}

//-----------------------------------------------------------------------

citaylor operator+(const citaylor& x, const citaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in citaylor, operator + : different orders " 
             << std::endl;
   exit(1);
  };

 citaylor erg(order1);
 for(int j=Lb(x.tayl); j<=Ub(x.tayl); j++) erg.tayl[j]= x.tayl[j]+y.tayl[j];
 return erg; 
}

//-----------------------------------------------------------------------

citaylor operator*(const citaylor& x, const citaylor& y)
{
 int order1=get_order(x);
 int order2=get_order(y);
 if(order1 != order2) 
  {
   std::cerr << "Error in citaylor, operator * : different orders " 
             << std::endl;
   exit(1);
  };

 citaylor erg(order1);
 cinterval sum; 
 cidotprecision sum_cidot; // for accumulate(...), scalar product

 for(int j=0; j<=order1; j++) 
 {
  sum_cidot=cinterval(0);
  for(int i=0; i<=j; i++)
   {
    accumulate(sum_cidot, x.tayl[i],y.tayl[j-i]);
   }
  rnd(sum_cidot,sum);
  erg.tayl[j]= sum;
 }
 return erg; 
}

//-----------------------------------------------------------------------

citaylor operator/(const citaylor& x, const citaylor& y)
{
 int order1(get_order(x));
 int order2(get_order(y));
 if(order1 != order2) 
  {
   std::cerr << "Error in citaylor, operator / : different orders " 
             << std::endl;
   exit(1);
  };
 
 if(0 <= y.tayl[0]) 
  {
   std::cerr << "Error in citaylor, operator / : 0 in interval " << std::endl;
   exit(1);
  };

 citaylor erg(order1);
 cinterval sum; 
 cidotprecision sum_cidot; // for accumulate(...), scalar product

 for(int j=0; j<=order1; j++) 
 {
  sum_cidot=cinterval(0);
  for(int i=1; i<=j; i++)
   {
    accumulate(sum_cidot, y.tayl[i],erg.tayl[j-i]);
   }
  rnd(sum_cidot,sum);
  erg.tayl[j]= (x.tayl[j]-sum)/y.tayl[0];
 }
 return erg;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (cinterval, citaylor):
// The operand of type cinterval is assumed to be a constant and not
// an independent variable!

citaylor operator-(const cinterval& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const cinterval& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const cinterval& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const cinterval& x, const citaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in citaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    cidotprecision cidot;
    cinterval sum;
    int order(get_order(y));
    citaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	cidot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(cidot,y.tayl[j],w.tayl[k-j]);
	rnd(cidot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (citaylor, cinterval):
// The operand of type cinterval is assumed to be a constant and not
// an independent variable!

citaylor operator-(const citaylor& x, const cinterval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const citaylor& x, const cinterval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const citaylor& x, const cinterval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const citaylor& x, const cinterval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (interval, citaylor):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

citaylor operator-(const interval& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const interval& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const interval& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x * y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const interval& x, const citaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in citaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    cidotprecision cidot;
    cinterval sum;
    int order(get_order(y));
    citaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	cidot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(cidot,y.tayl[j],w.tayl[k-j]);
	rnd(cidot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (citaylor, interval):
// The operand of type interval is assumed to be a constant and not
// an independent variable!

citaylor operator-(const citaylor& x, const interval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const citaylor& x, const interval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const citaylor& x, const interval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const citaylor& x, const interval& y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

//-----------------------------------------------------------------------


// Operators with two operands:   +,-,*,/  for (complex, citaylor):
// The operand of type complex is assumed to be a constant and not
// an independent variable!

citaylor operator-(const complex& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const complex& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const complex& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const complex& x, const citaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in citaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    cidotprecision cidot;
    cinterval sum;
    int order(get_order(y));
    citaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	cidot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(cidot,y.tayl[j],w.tayl[k-j]);
	rnd(cidot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (citaylor, complex):
// The operand of type complex is assumed to be a constant and not
// an independent variable!

citaylor operator-(const citaylor& x, const complex& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const citaylor& x, const complex& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const citaylor& x, const complex& y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const citaylor& x, const complex& y)
{
    if (y==0)
    {
	std::cerr << "Error in citaylor: division by 0" 
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}

// Operators with two operands:   +,-,*,/  for (real, citaylor):
// The operand of type real is assumed to be a constant and not
// an independent variable!

citaylor operator-(const real& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = -y;
    erg.tayl[0] = x - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const real& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = y;
    erg.tayl[0] = x + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const real& x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const real& x, const citaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in citaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    cidotprecision cidot;
    cinterval sum;
    int order(get_order(y));
    citaylor w(order);
    w.tayl[0] = x / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	cidot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(cidot,y.tayl[j],w.tayl[k-j]);
	rnd(cidot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (citaylor, real):
// The operand of type real is assumed to be a constant and not
// an independent variable!

citaylor operator-(const citaylor& x, const real& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const citaylor& x, const real& y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const citaylor& x, const real& y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*y;
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const citaylor& x, const real& y)
{
    if (y==0)
    {
	std::cerr << "Error in citaylor: division by 0" 
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/y;
    return erg;
}


// Operators with two operands:   +,-,*,/  for (int, citaylor):
// The operand of type int is assumed to be a constant and not
// an independent variable!

citaylor operator-(int x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = -y;
    erg.tayl[0] = cinterval(x) - y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(int x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    erg = y;
    erg.tayl[0] = cinterval(x) + y.tayl[0];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(int x, const citaylor& y)
{
    int order(get_order(y));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = cinterval(x)*y.tayl[j];
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(int x, const citaylor& y)
{
    if (0<=y.tayl[0])
    {
	std::cerr << "Error in citaylor, operator / : 0 in interval " 
                  << std::endl;
	exit(1);
    }; 
    cidotprecision cidot;
    cinterval sum;
    int order(get_order(y));
    citaylor w(order);
    w.tayl[0] = cinterval(x) / y.tayl[0];
    for (int k=1; k<=order; k++)
    {
	cidot=0;
	for (int j=1; j<=k; j++) 
	    accumulate(cidot,y.tayl[j],w.tayl[k-j]);
	rnd(cidot,sum);
	w.tayl[k] = -sum / y.tayl[0];
    }
    return w;
}

//-----------------------------------------------------------------------

// Operators with two operands:   +,-,*,/  for (citaylor, int):
// The operand of type int is assumed to be a constant and not
// an independent variable!

citaylor operator-(const citaylor& x, int y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] - cinterval(y);
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator+(const citaylor& x, int y)
{
    int order(get_order(x));
    citaylor erg(order);
    erg = x;
    erg.tayl[0] = x.tayl[0] + cinterval(y);
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator*(const citaylor& x, int y)
{
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]*cinterval(y);
    return erg;
}

//-----------------------------------------------------------------------

citaylor operator/(const citaylor& x, int y)
{
    if (y==0)
    {
	std::cerr << "Error in citaylor: division by 0" 
                  << std::endl;
	exit(1);
    };
    int order(get_order(x));
    citaylor erg(order);
    for (int j=0; j<=order; j++) erg.tayl[j] = x.tayl[j]/cinterval(y);
    return erg;
}

//-----------------------------------------------------------------------

// Overloading the standard functions:

//-----------------------------------------------------------------------

// sqr-function
citaylor sqr(const citaylor& x)
{
    cidotprecision cidot;
    cinterval sum;
    int order(get_order(x)),m;
    citaylor erg(order);

    erg.tayl[0] = sqr(x.tayl[0]);
    for (int k=1; k<=order; k++)
    {
	m = (k+1) / 2;
	cidot = 0;
	for (int j=0; j<=m-1; j++) 
	    accumulate(cidot,x.tayl[j],x.tayl[k-j]);
	rnd(cidot,sum);
//	times2pown(sum,1);  // Multiplication with 2
	erg.tayl[k] = 2*sum; // bei 'times2pown(sum,1)' entfaellt der Faktor 2.
	if (k%2==0) erg.tayl[k] += sqr(x.tayl[m]); // k even 
    }
    return erg; 
}

//-----------------------------------------------------------------------

// Square-root

citaylor sqrt(const citaylor& x)
{
    cidotprecision cidot;
    cinterval sum,h;
    int order(get_order(x)),m;
    citaylor erg(order);
    if (0<=x.tayl[0]) 
    {
	std::cerr << "Error in citaylor, sqrt: 0 in interval" << std::endl;
	exit(1);
    };
    erg.tayl[0] = sqrt(x.tayl[0]);
    h = erg.tayl[0];
//    times2pown(h,1); 
    h = 2*h;
    for (int k=1; k<=order; k++)
    {
	m = (k+1) / 2;
	cidot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(cidot,erg.tayl[j],erg.tayl[k-j]);
	rnd(cidot,sum);
//	times2pown(sum,1);  // Fast multiplication with 2
	erg.tayl[k] = 2*sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (x.tayl[k]-erg.tayl[k]) / h;
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(x,n)

citaylor sqrt(const citaylor& x, int n)
{
    int order(get_order(x));
    citaylor erg(order);
    if (0<=x.tayl[0]) 
    {
      std::cerr << "Error in citaylor, sqrt(x,n): 0 in interval" << std::endl;
      exit(1);
    };
    erg.tayl[0] = sqrt(x.tayl[0],n); // element No. 0
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=0; j<=k-1; j++)
	    erg.tayl[k] += (interval(k-j)/real(n)-interval(j))
		           * erg.tayl[j] * x.tayl[k-j];
	erg.tayl[k] /= (interval(k)*x.tayl[0]);
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(1+x^2)
citaylor sqrt1px2(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg(order);
    const real c = 500.0;
    cinterval f0(x.tayl[0]);
    interval R(Re(f0)), I(Im(f0));

    if (0<=R) erg = sqrt(1+sqr(x));
    else 
	if (Sup(abs(R))>c or Sup(abs(I))>c) 
	    if (Inf(R)>0) erg = x*sqrt(1+sqr(1/x));
	    else erg = -x*sqrt(1+sqr(1/x));
	else erg = sqrt(1+sqr(x));
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(1-x^2):
citaylor sqrt1mx2(const citaylor& x)
{
    cidotprecision cidot;
    cinterval sum,h;
    int order(get_order(x)),m;
    citaylor erg(order), g(order);

    erg.tayl[0] = sqrt(1-sqr(x.tayl[0]));  
    h = -1/erg.tayl[0];
    h = h/2;
    g = sqr(x);
    for (int k=1; k<=order; k++)
    {
	m = (k+1)/2;
	cidot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(cidot,erg.tayl[j],erg.tayl[k-j]);
	rnd(cidot,sum);
	sum = 2*sum; 
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (g.tayl[k]+erg.tayl[k])*h;
    }
    return erg;
}

//-----------------------------------------------------------------------

// sqrt(x^2-1):
citaylor sqrtx2m1(const citaylor& x)
{
    cidotprecision cidot;
    cinterval sum,h;
    int order(get_order(x)),m;
    citaylor erg(order), g(order);

    erg.tayl[0] = sqrt(sqr(x.tayl[0])-1);  
    g = sqr(x);
    h = 1/erg.tayl[0];
    h = h/2;
    for (int k=1; k<=order; k++)
    {
	m = (k+1)/2;
	cidot = 0;
	for (int j=1; j<=m-1; j++) 
	    accumulate(cidot,erg.tayl[j],erg.tayl[k-j]);
	rnd(cidot,sum);
	sum = 2*sum; 
	erg.tayl[k] = sum;
	if (k%2==0) erg.tayl[k] += sqr(erg.tayl[m]); // k even 
	erg.tayl[k] = (g.tayl[k]-erg.tayl[k])*h;
    }
    return erg;
}


//-----------------------------------------------------------------------

// power-function
citaylor pow(const citaylor& x, const interval& alpha)
{
    int order(get_order(x));
    citaylor erg(order);

    if (0<=x.tayl[0])
    {
	std::cerr << "Error in citaylor, pow(x,a): 0 in x" 
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

// power-function
citaylor pow(const citaylor& x, const cinterval& alpha)
{
    int order(get_order(x));
    citaylor erg(order);

    if (0<=x.tayl[0])
    {
	std::cerr << "Error in citaylor, pow(x,a): 0 in x" 
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

// power-function
citaylor power(const citaylor& x, int n)
{
    int order(get_order(x));
    citaylor erg(order);

    erg.tayl[0] = power(x.tayl[0],n); // element No. 0
    for (int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for (int j=0; j<=k-1; j++)
	    erg.tayl[k] += (interval(k-j)*n-interval(j))
		           * erg.tayl[j] * x.tayl[k-j];
	erg.tayl[k] /= (interval(k)*x.tayl[0]);
    }
    return erg;
}


//-----------------------------------------------------------------------

// Exponential-function
citaylor exp(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg(order);

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
citaylor expm1(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg(order);

    erg.tayl[0] = exp(x.tayl[0]); // element No. 0; function value
    for(int k=1; k<=order; k++)
    {
	erg.tayl[k] = 0;
	for(int j=0; j<=k-1; j++)
	    erg.tayl[k] += interval(k-j)*erg.tayl[j]*x.tayl[k-j]; 
	erg.tayl[k] /= interval(k);
    }
    erg.tayl[0] = exp(x.tayl[0])-real(1); // = expm1(x.tayl[0]); Blomquist
    return erg; 
}

//-----------------------------------------------------------------------


// Help function
void f_g_u(const citaylor& f, const citaylor& g, const citaylor& u, 
           int nb_function)
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
    case _i_ln:f.tayl[0]=Ln(u.tayl[0]); break;

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
 cinterval sum; 
 for(int j=1; j<=Ub(f.tayl); j++) 
 {
     sum = cinterval(0);
     for(int i=1; i<=j-1; i++)
	 sum += interval(i)*f.tayl[i]*g.tayl[j-i];
     f.tayl[j] = (u.tayl[j]-sum/interval(j)) / g.tayl[0];
 }
} // f_g_u


// Logarithm-function
citaylor ln(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);

    f_g_u(f,x,x,_i_ln);
  
    return f;
}

//-----------------------------------------------------------------------


// Sinus-function
citaylor sin(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg1(order);   // sin
    citaylor erg2(order);   // cos
    cinterval s1,s2;

    erg1.tayl[0]=sin(x.tayl[0]); // Element No. 0:  erg1 (sin)
    erg2.tayl[0]=cos(x.tayl[0]); // Element No. 0:  erg2 (cos)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=cinterval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= real(-1.0)/interval(j)*s2;
    }
    return erg1; 
}


//-----------------------------------------------------------------------

// Cosinus-function
citaylor cos(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg1(order);   // sin
    citaylor erg2(order);   // cos
    cinterval s1,s2;

    erg1.tayl[0]=sin(x.tayl[0]); // Element No. 0:  erg1 (sin)
    erg2.tayl[0]=cos(x.tayl[0]); // Element No. 0:  erg2 (cos)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=cinterval(0);
	for(int i=0; i<=j-1; i++)
	{
	    s1 += interval(j-i) * erg2.tayl[i] * x.tayl[j-i];
	    s2 += interval(j-i) * erg1.tayl[i] * x.tayl[j-i];
	}
	erg1.tayl[j]= s1/interval(j);
	erg2.tayl[j]= real(-1.0)/interval(j)*s2;
    }
 return erg2; 
}

//-----------------------------------------------------------------------

// Tangens-function
citaylor tan(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g=sqr(cos(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, tan : wrong argument" << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_tan);
  
    return f;
}

//-----------------------------------------------------------------------


// Cotangens-function
citaylor cot(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g=-sqr(sin(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, cot : wrong argument" << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_cot);
  
    return f;
}

//-----------------------------------------------------------------------


// Sinushyperbolicus-function
citaylor sinh(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg1(order);  // sinh
    citaylor erg2(order);  // cosh
    cinterval s1,s2;

    erg1.tayl[0]=sinh(x.tayl[0]); // element No. 0:  erg1 (sinh)
    erg2.tayl[0]=cosh(x.tayl[0]); // element No. 0:  erg2 (cosh)

    // remainig elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=cinterval(0);
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
citaylor cosh(const citaylor& x)
{
    int order(get_order(x));
    citaylor erg1(order); // sinh
    citaylor erg2(order); // cosh
    cinterval s1,s2;

    erg1.tayl[0]=sinh(x.tayl[0]); // element No. 0:  erg1 (sinh)
    erg2.tayl[0]=cosh(x.tayl[0]); // element No. 0:  erg2 (cosh)

    // remaining elements: 
    for(int j=1; j<=Ub(x.tayl); j++) 
    {
	s1=s2=cinterval(0);
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
citaylor tanh(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = sqr(cosh(x)); 
 
    f_g_u(f,g,x,_i_tanh);
  
    return f;
}

//-----------------------------------------------------------------------

// Cotangenshyperbolicus-function
citaylor coth(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g=-sqr(sinh(x));

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, coth : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_coth);
  
    return f;
}

//-----------------------------------------------------------------------

//Arcsinusfunktion
citaylor asin(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = sqrt1mx2(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, asin : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_asin);
  
    return f;
}

//-----------------------------------------------------------------------

//Arccosinusfunktion
citaylor acos(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = -sqrt1mx2(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, acos : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_acos);
  
    return f;
}

//-----------------------------------------------------------------------

//Arctan-function
citaylor atan(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = interval(1.0) + sqr(x);

    f_g_u(f,g,x,_i_atan);
  
    return f;
}

//-----------------------------------------------------------------------

//Arccotan-function
citaylor acot(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g=-(interval(1.0)+sqr(x));

    f_g_u(f,g,x,_i_acot);
  
    return f;
}

//-----------------------------------------------------------------------

//Areasinh-function
citaylor asinh(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = sqrt1px2(x);

    f_g_u(f,g,x,_i_asinh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areacosh-function
citaylor acosh(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = sqrtx2m1(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, acosh : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_acosh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areatanh-function
citaylor atanh(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = interval(1.0)-sqr(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, atanh : wrong argument " << std::endl;
	exit(1);
    };  
 
    f_g_u(f,g,x,_i_atanh);
  
    return f;
}

//-----------------------------------------------------------------------

//Areacotanh-function
citaylor acoth(const citaylor& x)
{
    int order=get_order(x);
    citaylor f(order);
    citaylor g(order);

    g = interval(1.0)-sqr(x);

    if(0 <= g.tayl[0]) 
    {
	std::cerr << "Error in citaylor, acoth : wrong argument " << std::endl;
	exit(1);
    };  

    f_g_u(f,g,x,_i_acoth);
  
    return f;
}

//-----------------------------------------------------------------------

} // end namespace taylor
