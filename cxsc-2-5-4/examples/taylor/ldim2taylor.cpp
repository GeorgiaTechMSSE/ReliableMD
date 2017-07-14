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

#include "ldim2taylor.hpp"

using namespace std;
using namespace cxsc;

//------------------------------------------------------------------------

// class ldim2taylor: 2-dim. Staggered Taylor arithmetic 

//------------------------------------------------------------------------

namespace taylor{
// Constructors:

ldim2taylor::ldim2taylor()
{
 p=0;
 dat=new l_ivector[p+1];
 for(int i=0; i<=p; i++) Resize(dat[i], 0, p+1-i);
}


//-----------------------------------------------------------------------

ldim2taylor::ldim2taylor(int order)
{
 p=order;
 dat=new l_ivector[p+1];
 for(int i=0; i<=p; i++) Resize(dat[i], 0, p+1-i);
}

//------------------------------------------------------------------------

ldim2taylor::ldim2taylor(const ldim2taylor& s)
{
 p=s.p;
 dat=new l_ivector[p+1];
 for(int i=0; i<=p ;i++) dat[i]=s.dat[i];

}
//------------------------------------------------------------------------

ldim2taylor::~ldim2taylor()
{
 delete[] dat;
 dat=NULL;
}

//------------------------------------------------------------------------
// Assignment operator:

ldim2taylor& ldim2taylor::operator=(const ldim2taylor& s)
{
 if(this != &s)
 {
  delete[] dat;
  dat=NULL;
  p=s.p;
  dat=new l_ivector[p+1];

  for(int i=0; i<=p ; i++) dat[i]=s.dat[i]; 

 }
 return *this;
}

//------------------------------------------------------------------------

l_ivector& ldim2taylor::operator[](int n) const
{
    return dat[n];
} 

//------------------------------------------------------------------------
  
ldim2taylor init_var(int order, int nr, const l_interval& value)
{

 ldim2taylor t(order);

 if( (nr<1) && (nr>2) ) std::cerr << "Error in ldim2taylor::init_var" 
                                  << std::endl;

 t[0][0]=value;

 if(t.p>0)
  {
   if(nr==1) t[1][0]=interval(1.0);
   else t[1][0]=interval(0.0);

   if(nr==2) t[0][1]=interval(1.0);
   else t[0][1]=interval(0.0);
   
   for(int i=2; i<=t.p; i++) {t[0][i]=interval(0.0);}   // Rest 0. line
   for(int i=1; i<=t.p-1; i++) {t[1][i]=interval(0.0);} // Rest 1. line  

   // remaining elements
   for(int j=2; j<=t.p ; j++) 
     for(int i=0; i<=t.p-j ; i++) t[j][i]=interval(0.0);       
  }

 return t;
}

//------------------------------------------------------------------------

ldim2taylor init_const(int order, const l_interval& value)
{

 ldim2taylor t(order);

 if(t.p>0)
  {
   for(int j=0; j<=t.p ; j++) 
     for(int i=0; i<=t.p-j ; i++) t[j][i]=interval(0.0);       
  }

 t[0][0]=value;

 return t;
}

//------------------------------------------------------------------------

void ldim2taylor::print_ldim2taylor() // screen output
{
    std::cout << SetDotPrecision(16*stagprec,16*stagprec) << Scientific;
 for(int j=0; j<=p ; j++) 
  {
   for(int i=0; i<=p-j ; i++) std::cout << dat[j][i]<< " ";
   std::cout << std::endl;
  }
 std::cout << std::endl;
}

//------------------------------------------------------------------------

// Overloading of the elementary operators 

ldim2taylor operator-(const ldim2taylor& s)
{
 int order=s.p;
 ldim2taylor erg(order);

 for(int j=0; j<=erg.p ; j++) 
  {
   for(int i=0; i<=erg.p-j ; i++) erg[j][i]=-s[j][i];
  }
 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const ldim2taylor& s, const ldim2taylor& t) 
{
 int order1=s.p;
 int order2=t.p;

 if(order1 != order2)
   {
     std::cerr << "ldim2taylor operator - : Operands with different orders";
     std::cerr << std::endl;
     exit(1);
   }

 ldim2taylor erg(order1);

 for(int j=0; j<=erg.p ; j++) 
  {
   for(int i=0; i<=erg.p-j ; i++) erg[j][i]=s[j][i]-t[j][i];
  }
 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator+(const ldim2taylor& s, const ldim2taylor& t)
{
int order1=s.p;
int order2=t.p;

 if(order1 != order2)
   {
     cerr << "ldim2taylor operator + : Operands with different orders";
     cerr << endl;
     exit(1);
   }

 ldim2taylor erg(order1);

 for(int j=0; j<=erg.p ; j++) 
  {
   for(int i=0; i<=erg.p-j ; i++) erg[j][i]=s[j][i]+t[j][i];
  }

 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator*(const ldim2taylor& s, const ldim2taylor& t)
{
 int order1=s.p;
 int order2=t.p;

 if(order1 != order2)
   {
     cerr << "ldim2taylor operator* : Operands with different orders";
     cerr << endl;
     exit(1);
   }

 ldim2taylor erg(order1);

 idotprecision sum_idot;

 for(int k=0; k<=erg.p; k++)
  {
    for(int i=0; i<=k; i++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i; l++) // calculating erg(i,k-i)
	{
	 for(int m=0; m<=k-i; m++)
	  {
	    accumulate(sum_idot, s[l][m], t[i-l][k-i-m]);
	  } 
	}  

//      rnd(sum_idot,erg[i][k-i]);  
      erg[i][k-i] = l_interval(sum_idot);    
     }
  }

 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator/(const ldim2taylor& s, const ldim2taylor& t)
{
 int order1=s.p;
 int order2=t.p;

 if(order1 != order2)
   {
     cerr << "ldim2taylor operator / : Operands with different orders";
     cerr << endl;
     exit(1);
   }

 ldim2taylor erg(order1);

 idotprecision sum_idot;

 l_interval h, sum;

 if(0.0 <= t[0][0]) 
  {
    cerr << "ldim2taylor operator/ : 0 in denominator" << endl;
    exit(1);
  }
 h = interval(1.0)/t[0][0];
 

 for(int k=0; k<=erg.p; k++) // calculating erg(i,k-i)
  {
    for(int i=0; i<=k; i++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i; l++) // Without Coeff. (l,m)=(0,0)
	{
	 for(int m=1; m<=k-i; m++) // m >= 1
	  {
	    accumulate(sum_idot, t[l][m], erg[i-l][k-i-m]);
	  }  //for m
	}    //for l

      for(int l=1; l<=i; l++) //  m=0, l!=0
	{
	 accumulate(sum_idot, t[l][0], erg[i-l][k-i]);
	}  // for l

//      rnd(sum_idot, sum);
      sum = l_interval(sum_idot);  
      erg[i][k-i] = h*(s[i][k-i]-sum);
     }
  }
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const l_interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 ldim2taylor s_ty = init_const(t.p, s);

 erg=s_ty-t;
 return erg;
}

ldim2taylor operator+(const l_interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 ldim2taylor s_ty=init_const(t.p, s);

 erg=s_ty+t;
 return erg;
}

ldim2taylor operator*(const l_interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 ldim2taylor s_ty=init_const(t.p, s);

 erg=s_ty*t;
 return erg;
}

ldim2taylor operator/(const l_interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 ldim2taylor s_ty=init_const(t.p, s);

 erg=s_ty/t;
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const ldim2taylor& s, const l_interval& t)
{
 ldim2taylor erg(s.p);
 ldim2taylor t_ty=init_const(s.p, t);

 erg=s-t_ty;
 return erg;
}

ldim2taylor operator+(const ldim2taylor& s, const l_interval& t)
{
 ldim2taylor erg(s.p);
 ldim2taylor t_ty=init_const(s.p, t);

 erg=s+t_ty;
 return erg;
}

ldim2taylor operator*(const ldim2taylor& s, const l_interval& t)
{
 ldim2taylor erg(s.p);
 ldim2taylor t_ty=init_const(s.p, t);

 erg=s*t_ty;
 return erg;
}

ldim2taylor operator/(const ldim2taylor& s, const l_interval& t)
{
 ldim2taylor erg(s.p);

 ldim2taylor t_ty=init_const(s.p, t);

 erg=s/t_ty;
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg = s_i-t;
 return erg;
}

ldim2taylor operator+(const real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i+t; 
 return erg;
}

ldim2taylor operator*(const real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i*t; 
 return erg;
}

ldim2taylor operator/(const real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i/t; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const ldim2taylor& s, const real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s-t_i; 
 return erg;
}

ldim2taylor operator+(const ldim2taylor& s, const real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s+t_i; 
 return erg;
}

ldim2taylor operator*(const ldim2taylor& s, const real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s*t_i; 
 return erg;
}

ldim2taylor operator/(const ldim2taylor& s, const real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s/t_i; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg = s_i-t;
 return erg;
}

ldim2taylor operator+(const interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i+t; 
 return erg;
}

ldim2taylor operator*(const interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i*t; 
 return erg;
}

ldim2taylor operator/(const interval& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i/t; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const ldim2taylor& s, const interval& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s-t_i; 
 return erg;
}

ldim2taylor operator+(const ldim2taylor& s, const interval& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s+t_i; 
 return erg;
}

ldim2taylor operator*(const ldim2taylor& s, const interval& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s*t_i; 
 return erg;
}

ldim2taylor operator/(const ldim2taylor& s, const interval& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s/t_i; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(int n, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(n); 
 
 erg=s_i-t; 
 return erg;
}

ldim2taylor operator+(int n, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(n);
 
 erg=s_i+t; 
 return erg;
}

ldim2taylor operator*(int n, const ldim2taylor& t) 
{
 ldim2taylor erg(t.p);
 l_interval s_i(n);
 
 erg = s_i*t; 
 return erg;
}

ldim2taylor operator/(int n, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(n);
 
 erg=s_i/t; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const ldim2taylor& s, int n)
{
 ldim2taylor erg(s.p);
 l_interval t_i(n);
 
 erg=s-t_i; 
 return erg;
}

ldim2taylor operator+(const ldim2taylor& s, int n)
{
 ldim2taylor erg(s.p);
 l_interval t_i(n);
 
 erg=s+t_i; 
 return erg;
}

ldim2taylor operator*(const ldim2taylor& s, int n)
{
 ldim2taylor erg(s.p);
 l_interval t_i(n);
 
 erg=s*t_i; 
 return erg;
}

ldim2taylor operator/(const ldim2taylor& s, int n)
{
 ldim2taylor erg(s.p);
 l_interval t_i(n);
 
 erg=s/t_i; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const l_real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg = s_i-t;
 return erg;
}

ldim2taylor operator+(const l_real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i+t; 
 return erg;
}

ldim2taylor operator*(const l_real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i*t; 
 return erg;
}

ldim2taylor operator/(const l_real& s, const ldim2taylor& t)
{
 ldim2taylor erg(t.p);
 l_interval s_i(s); 
 
 erg=s_i/t; 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor operator-(const ldim2taylor& s, const l_real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s-t_i; 
 return erg;
}

ldim2taylor operator+(const ldim2taylor& s, const l_real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s+t_i; 
 return erg;
}

ldim2taylor operator*(const ldim2taylor& s, const l_real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s*t_i; 
 return erg;
}

ldim2taylor operator/(const ldim2taylor& s, const l_real& t)
{
 ldim2taylor erg(s.p);
 l_interval t_i(t);
 
 erg=s/t_i; 
 return erg;
}

//------------------------------------------------------------------------

// Overloading the elementary functions

//------------------------------------------------------------------------


//------------------------------------------------------------------------

ldim2taylor sqr(const ldim2taylor& s)
{
 ldim2taylor erg(s.p);
 
 idotprecision sum_idot;
 l_interval sum1(0.0);
 l_interval sum2(0.0);

 erg[0][0]=sqr(s[0][0]); // Koeff. (0,0)

 if(erg.p > 0)
  {
    
    for(int k=0; k<=erg.p; k++) // calculating the remaining coefficients
      {
	for(int i=0; i<=k; i++)
	  {
	    if(i%2==1) // i: odd
	      {
		sum_idot=interval(0.0);
		for(int l=0; l<=(i-1)/2; l++) // calculating erg(i,k-i)
		  {
		    for(int m=0; m<=k-i; m++)
		      {
			accumulate(sum_idot, s[l][m], s[i-l][k-i-m]);
		      } //for m
		  }     //for l

		sum1 = l_interval(sum_idot);
		erg[i][k-i] = real(2.0)*sum1;	  	
	      }
	    else //i%2==0: i: even 
	      {
		if( (k-i)%2==1 ) //k-i: odd
		  {
		    sum_idot=interval(0.0);	
		    for(int l=0; l<=(i-2)/2; l++) // calculating erg(i,k-i)
		      {
			for(int m=0; m<=k-i; m++)
			  {
			    accumulate(sum_idot, s[l][m], s[i-l][k-i-m]);
			  } //for m
		      }     //for l
                    sum1 = l_interval(sum_idot);
		    sum_idot = interval(0.0);
		    for(int m=0; m<=(k-i-1)/2; m++)
		      {
			accumulate(sum_idot, s[i/2][m], s[i/2][k-i-m]);
		      }//for m
                    sum2 = l_interval(sum_idot);
		    erg[i][k-i]=interval(2.0)*(sum1+sum2);  
             
		  }
		else //(k-i)%2==0: k-i even
		  {
		    sum_idot=interval(0.0);	
		    for(int l=0; l<=(i-2)/2; l++) // calculating erg(i,k-i)
		      {
			for(int m=0; m<=k-i; m++)
			  {
			    accumulate(sum_idot, s[l][m], s[i-l][k-i-m]);
			  } // for m
		      }     // for l
                    sum1 = l_interval(sum_idot);
		    sum_idot=interval(0.0);
		    for(int m=0; m<=(k-i-2)/2; m++)
		      {
			accumulate(sum_idot, s[i/2][m], s[i/2][k-i-m]);
		      } // for m
                    sum2 = l_interval(sum_idot);
		    erg[i][k-i] = real(2.0)*(sum1+sum2);
		    erg[i][k-i] = erg[i][k-i]+sqr(s[i/2][(k-i)/2]);   
		  }
	      }
	  } // for i
      }     // for k
  }
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor sqrt(const ldim2taylor& s)
// New fast version, Blomquist 17.10.2005;
// In the sum (5.25) of Braeuers thesis by twos summands are equal,
// so the evaluation can be simplified. The first and the last
// summands are equal and must not be added. After the simplification
// only the first summand must not be added, see:
// if (z==1) continue; e.g. skipping the first summand.
{ 
 ldim2taylor erg(s.p);
 idotprecision sum_idot;
 l_interval sum(0.0);

 erg[0][0] = sqrt(s[0][0]);  // Coeff. (0,0) --> function value

 if(erg.p > 0)
  { 
    if(0.0 <= erg[0][0])
      {
	cerr << "Error in ldim2taylor sqrt: 0 in argument interval";
	cerr << endl;
	exit(1);
      }

    l_interval h = erg[0][0]; 
    times2pown(h,1); // fast multiplication with 2;
    int ki,il,z;
    for(int k=1; k<=erg.p; k++) // calculating all other coefficients
    {   
	for(int i=0; i<=k; i++)
	{   
	    sum_idot=interval(0.0);
    // do not sum (l,m)=(0,0) and (l,m)=(i,k-i), see continue assignments
	    ki = k-i;   z = 0; // z: Numbering the following summands 
	    if (i%2==1) // i: odd
	    { 
		for(int l=0; l<=(i-1)/2; l++) 
		{
		    il = i-l;   
		    for(int m=0; m<=ki; m++)
		    {
			z++; // Numbering the summands, z = 1,2,3,...
			// (l,m) = (0,0) is the first summand (z=1),
			if (z==1) continue; // skipping accumulate(...)
			accumulate(sum_idot, erg[l][m], erg[il][ki-m]);
		    }
		}
		sum = l_interval(sum_idot);
		times2pown(sum,1);
		erg[i][ki] = (s[i][ki]-sum)/h;
	    }
	    else // i: even
	    {
		if (ki%2==1) // ki=k-i: odd
		{
		    for(int l=0; l<=(i-2)/2; l++) 
		    {
			il = i-l;
			for(int m=0; m<=ki; m++)
			{
			    z++;
			    if (z==1) continue;
			    accumulate(sum_idot,erg[l][m],erg[il][ki-m]);
			}
		    }
		    for(int m=0; m<=(ki-1)/2; m++)
		    {
			z++;
			if (z==1) continue;
			accumulate(sum_idot,erg[i/2][m],erg[i/2][ki-m]);
		    }
		    sum = l_interval(sum_idot);
		    times2pown(sum,1);
		    erg[i][ki] = (s[i][ki]-sum)/h;
		}
		else // ki=k-i: even
		{
		    for(int l=0; l<=(i-2)/2; l++)
		    {
			il = i-l;
			for(int m=0; m<=ki; m++)
			{
			    z++;
			    if (z==1) continue;
			    accumulate(sum_idot,erg[l][m],erg[il][ki-m]);
			}
		    }
		    for(int m=0; m<=(ki-2)/2; m++)
		    {
			z++;
			if (z==1) continue;
			accumulate(sum_idot,erg[i/2][m],erg[i/2][ki-m]);
		    }
		    sum = l_interval(sum_idot);
		    times2pown(sum,1);
		    sum += sqr(erg[i/2][ki/2]);
		    erg[i][ki] = (s[i][ki]-sum)/h;
		}
	    }
	} // for i
      }   // for k
  }
 return erg;
 }

//------------------------------------------------------------------------

ldim2taylor sqrt1px2(const ldim2taylor& x){
    // sqrt(1+x^2);
    ldim2taylor erg;
    if (Inf(x[0][0])>1) erg = x * sqrt(1+sqr(1/x));
    else 
	if (Sup(x[0][0])<-1)
	    erg = -x * sqrt(1+sqr(1/x));
	else erg = sqrt(1+sqr(x));
    return erg; 
}

//------------------------------------------------------------------------

ldim2taylor sqrtx2m1(const ldim2taylor& x){
    // sqrt(x^2-1);
    ldim2taylor erg;
    if (Inf(x[0][0])>4) erg = x * sqrt( 1-sqr(1/x) );
    else 
	if (Sup(x[0][0])<-4)
	    erg = -x * sqrt( 1-sqr(1/x) );
	else erg = sqrt(sqr(x)-1);
    return erg; 
    }

//------------------------------------------------------------------------

ldim2taylor sqrtp1m1(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);
 g = 2*sqrt(1+s);;

 if(0<=interval(1)+s[0][0])
   {
     cerr << "Error in ldim2taylor sqrtp1m1 : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _sqrtp1m1);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor power(const ldim2taylor& s, const l_interval& alpha)
{
 ldim2taylor erg(s.p);

 idotprecision sum_idot;
 l_interval sum1, sum2, h;

 erg[0][0] = pow(s[0][0], alpha); 
 
 for(int j=1; j<=erg.p; j++)
   {
     sum_idot=interval(0.0);
     for(int i=0; i<=j-1; i++)
       {
	 h=alpha*(interval(j)-interval(i))-interval(i);
	 accumulate(sum_idot, h*erg[0][i], s[0][j-i]);
       }
//     rnd(sum_idot, sum1);
     sum1 = l_interval(sum_idot);
     erg[0][j]=sum1/interval(j)/s[0][0];
   }

 for(int i=1; i<=erg.p; i++) // now calculating the remaining erg(i,k)
  {
    for(int k=0; k<=erg.p-i; k++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i-1; l++) // Do not sum the coeff. (l,m)=(0,0) 
	{
	 h=alpha*(_interval(i)-_interval(l))-_interval(l);

	 for(int m=0; m<=k; m++)
	  {
	    accumulate(sum_idot, h*erg[l][m], s[i-l][k-m]);
	  }
	}
      sum1 = l_interval(sum_idot);  
      sum_idot=interval(0.0);
      for(int m=1; m<=k; m++) 
	{
	 accumulate(sum_idot, s[0][m], erg[i][k-m]);
	}

      sum2 = l_interval(sum_idot);
      erg[i][k]=(sum1/interval(i)-sum2)/s[0][0];
     }
  } 
 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor power(const ldim2taylor& s, int n)
{
 ldim2taylor erg(s.p);

 idotprecision sum_idot;
 l_interval sum1, sum2, h;

 erg[0][0]=power(s[0][0], n); 
 
 for(int j=1; j<=erg.p; j++)
   {
     sum_idot=interval(0.0);
     for(int i=0; i<=j-1; i++)
       {
	 h = l_interval(n)*(interval(j)-interval(i))-interval(i);
	 accumulate(sum_idot, h*erg[0][i], s[0][j-i]);
       }
     sum1 = l_interval(sum_idot);
     erg[0][j]=sum1/_interval(j)/s[0][0];
   }

 for(int i=1; i<=erg.p; i++) // now calculating the remaining erg(i,k)
  {
    for(int k=0; k<=erg.p-i; k++)
     {
      sum_idot=interval(0.0);

      for(int l=0; l<=i-1; l++) // Don't sum coeff. (l,m)=(0,0)
	{
	 h=l_interval(n)*(interval(i)-interval(l))-interval(l);

	 for(int m=0; m<=k; m++)
	  {
	    accumulate(sum_idot, h*erg[l][m], s[i-l][k-m]);
	  }
	}
      sum1 = l_interval(sum_idot);  
      sum_idot=_interval(0.0);
      for(int m=1; m<=k; m++) 
	{
	 accumulate(sum_idot, s[0][m], erg[i][k-m]);
	}
      sum2 = l_interval(sum_idot);
      erg[i][k]=(sum1/interval(i)-sum2)/s[0][0];
     }
  } 
 
 return erg;
}

//------------------------------------------------------------------------

ldim2taylor exp(const ldim2taylor& s)
{
 ldim2taylor erg(s.p);

 idotprecision sum_idot;
 l_interval sum;

 erg[0][0]=exp(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating erg(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h = interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot, h*erg[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

             sum = l_interval(sum_idot); 
	     erg[i][k-i]=sum/interval(k);
	   }
       }
   }

 return erg;
}

//------------------------------------------------------------------------

ldim2taylor ln(const ldim2taylor& s)
{
 ldim2taylor f(s.p);

 if(0<=s[0][0])
   {
     cerr << "Error in ldim2taylor ln : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, s, s, _ln);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor lnp1(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);
 g = interval(1) + s;

 if(0<=interval(1)+s[0][0])
   {
     cerr << "Error in ldim2taylor lnp1 : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _lnp1);
 return f;
}
 
//------------------------------------------------------------------------

ldim2taylor sin(const ldim2taylor& s)
{
 ldim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 l_interval sum1, sum2;

 erg1[0][0]=sin(s[0][0]);
 erg2[0][0]=cos(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1=interval(0.0);
	     sum_idot2=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     sum1 = l_interval(sum_idot1);  
	     sum2 = l_interval(sum_idot2);

	     erg1[i][k-i]=sum1/interval(k);
	     erg2[i][k-i]=-sum2/interval(k);
	   }
       }
   }

 return erg1;
}

//------------------------------------------------------------------------

ldim2taylor cos(const ldim2taylor& s)
{
 ldim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 l_interval sum1, sum2;

 erg1[0][0] = sin(s[0][0]);
 erg2[0][0] = cos(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1 = interval(0.0);
	     sum_idot2 = interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h=interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     sum1 = l_interval(sum_idot1);
	     sum2 = l_interval(sum_idot2);

	     erg1[i][k-i]=sum1/_interval(k);
	     erg2[i][k-i]=-sum2/_interval(k);
	   }
       }
   }

 return erg2;
}

//------------------------------------------------------------------------

ldim2taylor tan(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g = sqr(cos(s));

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor tan : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _tan);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor cot(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g=-sqr(sin(s));

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor cot : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _cot);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor sinh(const ldim2taylor& s)
{
 ldim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 l_interval sum1, sum2;

 erg1[0][0]=sinh(s[0][0]);
 erg2[0][0]=cosh(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1 = interval(0.0);
	     sum_idot2 = interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h = interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   }  // for m
	       }      // for l

	     sum1 = l_interval(sum_idot1);
	     sum2 = l_interval(sum_idot2);

	     erg1[i][k-i]=sum1/_interval(k);
	     erg2[i][k-i]=sum2/_interval(k);
	   }
       }
   }

 return erg1;
}

//------------------------------------------------------------------------

ldim2taylor cosh(const ldim2taylor& s)
{
 ldim2taylor erg1(s.p), erg2(s.p); 

 idotprecision sum_idot1, sum_idot2;
 l_interval sum1, sum2;

 erg1[0][0]=sinh(s[0][0]);
 erg2[0][0]=cosh(s[0][0]);

 if(s.p>0)
   {
     for(int k=1; k<=erg1.p; k++)
       {
	 for(int i=0; i<=k; i++)
	   {
	     sum_idot1=interval(0.0);
	     sum_idot2=interval(0.0);

	     for(int l=0; l<=i; l++) // now calculating ergs(i,k-i)
	       {
		 for(int m=0; m<=k-i; m++)
		   { 
		     interval h = interval(k)-interval(l)-interval(m);
		     accumulate(sum_idot1, h*erg2[l][m], s[i-l][k-i-m]);
		     accumulate(sum_idot2, h*erg1[l][m], s[i-l][k-i-m]);
		   } // for m
	       }     // for l

	     sum1 = l_interval(sum_idot1); 
	     sum2 = l_interval(sum_idot2); 

	     erg1[i][k-i]=sum1/_interval(k);
	     erg2[i][k-i]=sum2/_interval(k);
	   }
       }
   }

 return erg2;
}

//------------------------------------------------------------------------

ldim2taylor tanh(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g=sqr(cosh(s));

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor tanh : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _tanh);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor coth(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g=-sqr(sinh(s));

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor coth : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _coth);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor asin(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g=sqrt( l_interval(1.0)-sqr(s) );

 if(0<=g[0][0])
   {
     cerr << "Errorr in ldim2taylor asin : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _asin);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor acos(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g=-sqrt( l_interval(1.0)-sqr(s) );

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor acos : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acos);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor atan(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g = l_interval(1.0)+sqr(s);

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor atan : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _atan);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor acot(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g=-(l_interval(1.0)+sqr(s));

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor acot : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acot);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor asinh(const ldim2taylor& s)
{
    ldim2taylor f(s.p), g(s.p);

//  g=sqrt(l_interval(1.0)+sqr(s));
    g = sqrt1px2(s);  // Blomquist, 08.10.2005;

    if(0<=g[0][0])
    {
	cerr << "Error in ldim2taylor asinh : 0 in interval" << endl;
	exit(1);
    }
    f_g_u(f, g, s, _asinh);
    return f;
}

//------------------------------------------------------------------------

ldim2taylor acosh(const ldim2taylor& s)
{
    ldim2taylor f(s.p), g(s.p);

//  g=sqrt(sqr(s)-l_interval(1.0));
    g = sqrtx2m1(s);  // Blomquist, 08.10.2005;
    if(0<=g[0][0])
    {
	cerr << "Error in ldim2taylor acosh : 0 in interval" << endl;
	exit(1);
    }
    f_g_u(f, g, s, _acosh);
    return f;
}

//------------------------------------------------------------------------

ldim2taylor atanh(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g = l_interval(1.0)-sqr(s);

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor atanh : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _atanh);
 return f;
}

//------------------------------------------------------------------------

ldim2taylor acoth(const ldim2taylor& s)
{
 ldim2taylor f(s.p), g(s.p);

 g = l_interval(1.0)-sqr(s);

 if(0<=g[0][0])
   {
     cerr << "Error in ldim2taylor acoth : 0 in interval" << endl;
     exit(1);
   }
 f_g_u(f, g, s, _acoth);
 return f;
}

//------------------------------------------------------------------------

void f_g_u(const ldim2taylor& f, const ldim2taylor& g, const ldim2taylor& u,
		 int _fkt) 
{
 idotprecision sum_idot;

 l_interval sum1, sum2;
 
 switch(_fkt) 
   { 
   case _ln:   {f[0][0]=ln(u[0][0]);   break;}
   case _tan:  {f[0][0]=tan(u[0][0]);  break;}
   case _cot:  {f[0][0]=cot(u[0][0]);  break;}

   case _asin:  {f[0][0]=asin(u[0][0]);  break;}
   case _acos:  {f[0][0]=acos(u[0][0]);  break;}
   case _atan:  {f[0][0]=atan(u[0][0]);  break;}
   case _acot:  {f[0][0]=acot(u[0][0]);  break;}

   case _tanh:  {f[0][0]=tanh(u[0][0]);  break;}
   case _coth:  {f[0][0]=coth(u[0][0]);  break;}

   case _asinh:  {f[0][0]=asinh(u[0][0]);  break;}
   case _acosh:  {f[0][0]=acosh(u[0][0]);  break;}
   case _atanh:  {f[0][0]=atanh(u[0][0]);  break;}
   case _acoth:  {f[0][0]=acoth(u[0][0]);  break;}
   } 

 if(f.p>0)
   {
     for(int k=1; k<=f.p; k++)
       {
	 sum_idot=interval(0.0);

	 for(int j=1; j<=k-1; j++)
	   {
	     accumulate(sum_idot, _interval(j)*f[0][j], g[0][k-j]);
	   }
	 sum1 = l_interval(sum_idot);
	 f[0][k]=(u[0][k]-sum1/_interval(k))/g[0][0];
       }

     for(int i=1; i<=f.p; i++) // now calculating the remainding f(i,j)
       {
	 for(int j=0; j<=f.p-i; j++)
	   {
	     sum_idot=interval(0.0);
	     
	     for(int l=1; l<=i-1; l++) 
	       {
		 for(int m=0; m<=j; m++)
		  {
		  accumulate(sum_idot, interval(l)*f[l][m], g[i-l][j-m]);
		  }
	       }
	     sum1 = l_interval(sum_idot);	     
	     sum_idot = interval(0.0);
	     for(int m=1; m<=j; m++) 
	       {
		 accumulate(sum_idot, g[0][m], f[i][j-m]);
	       }
	     sum2 = l_interval(sum_idot);
	     f[i][j]=(u[i][j]-sum1/_interval(i)-sum2)/g[0][0];
	   }
       } 
   }
}

//------------------------------------------------------------------------

// Class ldim2taylor_vector

//------------------------------------------------------------------------
ldim2taylor_vector::ldim2taylor_vector() // Default constructor
{
    dim=2;
    lb=1;
    ub=2;
    p_el=1;
    comp = new ldim2taylor[dim]; 
    for (int i=0; i<dim; i++)
      comp[i]=ldim2taylor(p_el);
}


ldim2taylor_vector::ldim2taylor_vector(int ordnung, int Lb, int Ub)
{
    lb=Lb;
    ub=Ub;
    dim=ub-lb+1;
    p_el=ordnung;
    comp = new ldim2taylor[dim]; 
    for (int i=0; i<dim; i++)
      comp[i]=ldim2taylor(p_el);
}

//------------------------------------------------------------------------

ldim2taylor_vector::ldim2taylor_vector(const ldim2taylor_vector& s)
{
 dim=s.dim;
 lb=s.lb;
 ub=s.ub;
 p_el=s.p_el;
 comp=new ldim2taylor[dim];
          //Taylor order p_el and memory are handled by ldim2taylor.operator=
 for(int i=0; i<dim; i++) comp[i]=s.comp[i];
}

//------------------------------------------------------------------------

ldim2taylor_vector::~ldim2taylor_vector()
{
 delete[] comp;
 comp = NULL;
}

//------------------------------------------------------------------------

ldim2taylor_vector& ldim2taylor_vector::operator=(const ldim2taylor_vector& s)
{
 if(this!=&s)
 {
  delete[] comp;

  p_el=s.p_el;
  dim=s.dim;
  lb=s.lb;
  ub=s.ub;

  comp=new ldim2taylor[dim];
          //Taylor order p_el and memory are handled by dim2taylor.operator=
  for(int i=0; i<dim ; i++) comp[i]=s.comp[i];
 }

 return *this;
}

//------------------------------------------------------------------------

ldim2taylor& ldim2taylor_vector::operator[](int n) const
{
    return comp[n-lb];  // allowed indices: n = lb,lb+1,...,ub;
} 

//------------------------------------------------------------------------

ldim2taylor_vector init_var(int order, l_ivector& values)
{
 int dimension = Ub(values)-Lb(values)+1;

 if(dimension != 2) 
  {
   cerr << "Error in ldim2taylor_vector::init_var" << endl;
   cerr << "! 2-dimensional Taylor arithmetic !" << endl; 
   exit(1);
  } 

 ldim2taylor_vector erg(order, Lb(values), Ub(values));

 erg[Lb(values)] = init_var(order,1,values[Lb(values)]);
 erg[Ub(values)] = init_var(order,2,values[Ub(values)]);

 return erg;
}

} // end of namespace taylor
