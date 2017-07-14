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

#include "testgen.hpp"
#include <real.hpp>
#include <rvector.hpp>
#include <rmatrix.hpp>
#include <rmath.hpp>
#include <dot.hpp>
#include <cvector.hpp>
#include <complex.hpp>
#include <ivector.hpp>
#include <civector.hpp>
#include <interval.hpp>
#include <cinterval.hpp>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <matinv.hpp>

using namespace std;
using namespace cxsc;

//! Logarithm of base 10 
real log(real x) {
  return ln(x) / ln(real(10.0));
}

/*!
 * Generates an mxn real matrix, whose elements are randomly distributed in
 * [0,1]
 *
 * \param m First dimension of matrix
 * \param n second dimension of matrix
 * \return Random mxn matrix
 */ 
rmatrix rando(int m, int n) {
  rmatrix x(m, n);
  int i, j;
  real rand_max = real(RAND_MAX);
  
  for(i=0; i<m; i++) {
     for(j=0; j<n; j++) {
        x[i+1][j+1] = real(std::rand()) / rand_max;
     }
  } 
  
  return x;
}

/*!
 * Generates an mxn real matrix, whose elements are randomly distributed in
 * [-1,1]
 *
 * \param m First dimension of matrix
 * \param n second dimension of matrix
 * \return Random mxn matrix
 */ 
rmatrix randn(int m, int n) {  
  rmatrix x(m, n);
  int k;
  double pi  = 4.0 * atan(1.0);
  double square, amp=0, angle=0;

  rmatrix u = rando(m * n + 1 , 1);
  
  k = 0;
  for(int i = 0; i < m; i++) {   
     for(int j = 0; j < n; j++) {    
        if( k % 2 == 0 ) {    
           square = - 2.0 * log( _double(u[k+1][1]) );
           if( square < 0.0 ) square = 0.0;
           amp = sqrt(square);
           angle = 2.0 * pi * _double(u[k+2][1]);
           x[i+1][j+1] = amp * sin( angle );
        } else { 
           x[i+1][j+1] = amp * cos( angle );
        }
        k++;
     }
  }
  
  return x;
}


/*!
 * Generates two real vectors of length n, whose dot product has condition
 * cond and whose exact result is res.
 *
 * \param n Length of vectors
 * \param cnd Desired condition
 * \param x First vector of dot product
 * \param y Second vector of dot product
 * \param res Exact result of dot product
 *
 */
void GenDot(int n, real cnd, rvector &x, rvector &y, real &res) {
  x = rvector(n);
  y = rvector(n);
  real eps = pow(2.0,-53);
  int m = n/2;
  real L_real = log(cnd) / (-log(eps));
  int L = (int)floor(_double(L_real));

  if(L == 0) L=1;

  if(n % 2 == 0) {
    rvector r(m-2), c(m-2), b(m-2);
    
    for(int i=1 ; i<=m-2 ; i++) {
      r[i] = i % L;
    }
    
    rmatrix tmp = randn(1,m-2);
    for(int i=1 ; i<=m-2 ; i++) {
      c[i] = tmp[1][i] * pow(eps,r[i]);
    }
    
    for(int i=1 ; i<=n ; i++) {
      if(i == 1) {
        x[i] = 1;
      } else if(i > 1  &&  i < m) {
        x[i] = c[i-1];
      } else if(i == m) {
        x[i] = 0.5 * pow(cnd, -1.0);
      } else if(i == m+1) {
        x[i] = -1;
      } else if(i > m+1  &&  i < 2*m) {
        x[i] = -c[i-(m+1)];
      } else {
        x[i] = 0.5 * pow(cnd, -1.0);
      }
    }
    
    tmp = randn(1,m-2);
    b = tmp[1];


    for(int i=1 ; i<=n ; i++) {
      if(i == 1  ||  i == m  ||  i == m+1  ||  i == n) {
        y[i] = 1;
      } else if(i > 1  &&  i < m) {
        y[i] = b[i-1];
      } else if(i > m+1  &&  i < 2*m) {
        y[i] = b[i-(m+1)];
      }
    }


  } else {


    rvector r(m-1), c(m-1), b(m-1);
    
    for(int i=1 ; i<=m-1 ; i++) {
      r[i] = i % L;
    }
    
    rmatrix tmp = randn(1,m-1);
    for(int i=1 ; i<=m-1 ; i++) {
      c[i] = tmp[1][i] * pow(eps,r[i]);
    }
    
    for(int i=1 ; i<=n ; i++) {
      if(i == 1) {
        x[i] = 1;
      } else if(i > 1  &&  i < m+1) {
        x[i] = c[i-1];
      } else if(i == m+1) {
        x[i] = pow(cnd, -1.0);
      } else if(i == m+2) {
        x[i] = -1;
      } else if(i > m+2) {
        x[i] = -c[i-(m+2)];
      }
    }
    
    tmp = randn(1,m-1);
    b = tmp[1];

    for(int i=1 ; i<=n ; i++) {
      if(i == 1  ||  i == m+1  ||  i == m+2) {
        y[i] = 1;
      } else if(i > 1  &&  i < m+1) {
        y[i] = b[i-1];
      } else if(i > m+2) {
        y[i] = b[i-(m+2)];
      }
    }
  }
  
  res = power(cnd, -1);
}


/*!
 * Generates two interval vectors of length n, tightly enclosing two real
 * vectors whose dot product is of condition cond and whose result is res
 *
 * \param n Length of vectors
 * \param cnd Desired condition
 * \param x First vector of dot product
 * \param y Second vector of dot product
 * \param res Exact result of dot product
 *
 */
void GenIDot(int n, real cnd, ivector &x, ivector &y, real &res) {
  rvector a(n),b(n);
  real eps = power(2,-53);
  real inv = power(cnd,-1);
  real delta = abs(inv * eps);
  x = ivector(n);
  y = ivector(n);
  
  //Generate real dot product
  GenDot(n, cnd, a, b, res);
  
  //Tightly enclose dot product
  for(int i=1 ; i<=n ; i++) {
    x[i] = interval(a[i]-delta, a[i]+delta);
    y[i] = interval(b[i]-delta, b[i]+delta);
  }
} 


/*!
 * Generates two complex vectors of length n, whose dot product has condition
 * cond and whose exact result is res.
 *
 * \param n Length of vectors
 * \param c Desired condition
 * \param x First vector of dot product
 * \param y Second vector of dot product
 * \param res Exact result of dot product
 *
 */
void GenCDot(int n, real c, cvector &x, cvector &y, complex &res) {
  rvector a,b;
  x = cvector(n);
  y = cvector(n);
  real res_re, res_im=0;
  
  //Generate real dot product of length 2*n
  GenDot(2*n, c, a, b, res_re);
  
  //Set vectors such that the real part of the result ist condition ^-1
  for(int i=1 ; i<=n ; i++) {
    SetRe(x[i], a[2*i-1]);
    SetIm(x[i], -a[2*i]);
    SetRe(y[i], b[2*i-1]);
    SetIm(y[i], b[2*i]);
  }
  
  //Determine imaginary part of the result
  if(n%2 == 0)
    res_im = res_re * (-Re(y[n]));
  else {
    dotprecision accu(0);
    for(int i=1 ; i<=n ; i++) {
       accumulate(accu, Re(x[i]),Im(y[i]));
       accumulate(accu, Im(x[i]),Re(y[i]));       
    }
    res_im = rnd(accu);
  }
  res = complex(res_re,res_im);
}

/*!
 * Generates two complex interval vectors of length n, tightly enclosing two 
 * complex vectors whose dot product is of condition cond and whose result is 
 * res
 *
 * \param n Length of vectors
 * \param c Desired condition
 * \param x First vector of dot product
 * \param y Second vector of dot product
 * \param res Exact result of dot product
 *
 */
void GenCIDot(int n, real c, civector &x, civector &y, complex &res) {
  ivector a,b;
  x = civector(n);
  y = civector(n);
  real res_re, res_im=0;
  
  GenIDot(2*n, c, a, b, res_re);
  
  for(int i=1 ; i<=n ; i++) {
    SetRe(x[i], a[2*i-1]);
    SetIm(x[i], -a[2*i]);
    SetRe(y[i], b[2*i-1]);
    SetIm(y[i], b[2*i]);
  }
  
  if(n%2 == 0)
    res_im = res_re * (-Re(mid(y[n])));
  else {
    dotprecision accu(0);
    for(int i=1 ; i<=n ; i++) {
       accumulate(accu, Re(mid(x[i])),Im(mid(y[i])));
       accumulate(accu, Im(mid(x[i])),Re(mid(y[i])));       
    }
    res_im = rnd(accu);
  } 
  
  res = complex(res_re, res_im);
}



//------------------------------------------------------------------------------


/*!
 * Generates a selection of real matrices of dimesion n. Type chooses which 
 * matrix to generate
 *
 * \param n Dimension
 * \param A Generated Matrix
 * \param type Type of generated matrix
 */
void GenMat(int n, rmatrix &A, int type) {
     A = rmatrix(n,n);
     int i,j;
     real x,y;
     rvector a;
     
     switch(type) {
        case 1:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/(n+1)) * mid(sin((i*j*Pi())/(n+1)));
                }
             break;
        case 2:             
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = n - cxsc::abs(i-j);
                }
             }
             break;
        case 3:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = (double)i/j;
                   else
                     A[i][j] = (double)j/i;
                }
             }
             break;             
        case 4:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = n+1-j;
                   else
                     A[i][j] = n+1-i;
                }
             }
             break;             
        case 5:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = j;
                   else
                     A[i][j] = i;
                }
             }
             break;              
        case 6:
             a = 0.5*randn(1,n)[Row(1)];
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = 1;
                   else
                     A[i][j] = a[j];
                }
             }
             break;           
        case 7:
             x = randn(1,1)[1][1] + 1.0;
             y = randn(1,1)[1][1] + 1.0;
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i==j)
                     A[i][j] = x+y;
                   else
                     A[i][j] = y;
                }
             }
             break;                                                     
        default:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/(n+1)) * mid(sin((i*j*Pi())/(n+1)));
                }
             break;
        }
}

/*!
 * Generates a selection of tight interval matrices of dimesion n. 
 * Type chooses which matrix to generate
 *
 * \param n Dimension
 * \param A Generated Matrix
 * \param type Type of generated matrix
 */
void GenIMat(int n, imatrix &A, int type) {
     A = imatrix(n,n);
     int i,j;
     real x,y;
     rvector a;
          
     switch(type) {
        case 1:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/n+1) * mid(sin((i*j*Pi())/(n+1)));
                }  
             break;
        case 2:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = n - cxsc::abs(i-j);
                }
             }
             break;

        case 3:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i <= j)
                      A[i][j] = interval(i)/interval(j);
                   else 
                      A[i][j] = interval(j)/interval(i);
                }             
             }
             break;
        case 4:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = n+1-j;
                   else
                     A[i][j] = n+1-i;
                }
             }
             break;             
        case 5:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = j;
                   else
                     A[i][j] = i;
                }
             }
             break;              
        case 6:
             a = 0.5*randn(1,n)[Row(1)];
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = 1;
                   else
                     A[i][j] = a[j];
                }
             }
             break;           
        case 7:
             x = randn(1,1)[1][1] + 1.0;
             y = randn(1,1)[1][1] + 1.0;
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i==j)
                     A[i][j] = x+y;
                   else
                     A[i][j] = y;
                }
             }
             break;                 
        default:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/n+1) * mid(sin((i*j*Pi())/(n+1)));
                }                           
             break;
     }
}

/*!
 * Generates a selection of complex matrices of dimesion n. Type chooses which 
 * matrix to generate
 *
 * \param n Dimension
 * \param A Generated Matrix
 * \param type Type of generated matrix
 */
void GenCMat(int n, cmatrix &A, int type) {
     A = cmatrix(n,n);
     int i,j;
     complex x,y;
     cvector a;
     
     switch(type) {
        case 1:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/(n+1)) * mid(sin((i*j*Pi())/(n+1)));
                }
             break;
        case 2:             
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = n - cxsc::abs(i-j);
                }
             }
             break;
        case 3:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = (double)i/j;
                   else
                     A[i][j] = (double)j/i;
                }
             }
             break;             
        case 4:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = n+1-j;
                   else
                     A[i][j] = n+1-i;
                }
             }
             break;             
        case 5:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = j;
                   else
                     A[i][j] = i;
                }
             }
             break;              
        case 6:
             a = cvector(n);
             SetRe(a, 0.5*randn(1,n)[Row(1)]);
             SetIm(a, 0.5*randn(1,n)[Row(1)]);
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = 1;
                   else
                     A[i][j] = a[j];
                }
             }
             break;           
        case 7:
             x = complex(randn(1,1)[1][1],randn(1,1)[1][1]) + 1.0;
             y = complex(randn(1,1)[1][1],randn(1,1)[1][1]) + 1.0;
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i==j)
                     A[i][j] = x+y;
                   else
                     A[i][j] = y;
                }
             }
             break;                                                     
        default:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/(n+1)) * mid(sin((i*j*Pi())/(n+1)));
                }
             break;
        }
}

/*!
 * Generates a selection of tight interval matrices of dimesion n. 
 * Type chooses which matrix to generate
 *
 * \param n Dimension
 * \param A Generated Matrix
 * \param type Type of generated matrix
 */
void GenCIMat(int n, cimatrix &A, int type) {
     A = cimatrix(n,n);
     int i,j;
     complex x,y;
     cvector a;
          
     switch(type) {
        case 1:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/n+1) * mid(sin((i*j*Pi())/(n+1)));
                }  
             break;
        case 2:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = n - cxsc::abs(i-j);
                }
             }
             break;

        case 3:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i <= j)
                      A[i][j] = interval(i)/interval(j);
                   else 
                      A[i][j] = interval(j)/interval(i);
                }             
             }
             break;
        case 4:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = n+1-j;
                   else
                     A[i][j] = n+1-i;
                }
             }
             break;             
        case 5:
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = j;
                   else
                     A[i][j] = i;
                }
             }
             break;              
        case 6:
             a = cvector(n);             
             SetRe(a, 0.5*randn(1,n)[Row(1)]);
             SetIm(a, 0.5*randn(1,n)[Row(1)]);
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i<=j)
                     A[i][j] = 1;
                   else
                     A[i][j] = a[j];
                }
             }
             break;           
        case 7:
             x = complex(randn(1,1)[1][1],randn(1,1)[1][1]) + 1.0;
             y = complex(randn(1,1)[1][1],randn(1,1)[1][1]) + 1.0;
             for(i=1 ; i<=n ; i++) {
                for(j=1 ; j<=n ; j++) {
                   if(i==j)
                     A[i][j] = x+y;
                   else
                     A[i][j] = y;
                }
             }
             break;                
        default:
             for(i=1 ; i<=n ; i++) 
                for(j=1 ; j<=n ; j++) {
                   A[i][j] = sqrt(2.0/n+1) * mid(sin((i*j*Pi())/(n+1)));
                }                           
             break;
     }
}

/*!
 * Generates a real matrix of dimension n with choosen condition
 *
 * \param n Dimension of matrix
 * \param A The generated matrix
 * \param cond The desired condition of the generated matrix
 */
void GenMat(int n, rmatrix &A, real cond) {
  A = rmatrix(n,n);
  rmatrix BB(n,n);
  dotprecision accu(0);
  real norm = 0;

  //Generate random 1xn matrix
  rmatrix B = randn(1,n);

  //Compute norm of B
  for(int i=1 ; i<=n ; i++) {
    accu += sqr(B[1][i]);
  }
  norm = sqrt(rnd(accu));
  
  //Normalize B
  B = (1.0/norm)*B;
  
  //Compute matrix A
  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=n ; j++) {      
      BB[i][j] = B[1][i]*B[1][j];
      A[i][j] = cond*BB[i][j];
      if(i==j) A[i][j] += 1;               
    }
  }
  
}

/*!
 * Generates an interval matrix of dimension n, tightly enclosing a real
 * matrix with choosen condition
 *
 * \param n Dimension of matrix
 * \param A The generated matrix
 * \param cond The desired condition of the midpoint matrix
 */
void GenIMat(int n, imatrix &A, real cond) {
  A = imatrix(n,n);
  rmatrix B = rando(1,n);
  rmatrix BB(n,n);
  dotprecision accu(0);
  real norm = 0;
  const real eps = power(2,-53);
  
  for(int i=1 ; i<=n ; i++) {
    accu += sqr(B[1][i]);
  }
  norm = power(sqrt(rnd(accu)),-1);
  B = B*norm;
  
  
  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=n ; j++) {      
      BB[i][j] = B[1][i]*B[1][j];
      A[i][j] = BB[i][j]*cond;
      if(i==j) A[i][j] += 1;               
      A[i][j] += A[i][j] * interval(-eps,eps);
    }
  }

}


/*!
 * Generates a complex matrix of dimension n with choosen condition
 *
 * \param n Dimension of matrix
 * \param A The generated matrix
 * \param cond The desired condition of the generated matrix
 */
void GenCMat(int n, cmatrix &A, real cond) {
  A = cmatrix(n,n);
  cmatrix B(1,n);
  cmatrix BB(n,n);
  cdotprecision accu(0);
  real norm = 0;
  
  SetRe(B,rando(1,n));
  SetIm(B,rando(1,n));

  for(int i=1 ; i<=n ; i++)
    accu += sqr(abs(B[1][i]));    
  norm = power(sqrt(abs(rnd(accu))),-1);

  B = B*norm;

  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=n ; j++) {      
      complex tmp(Re(B[1][j]), -Im(B[1][j]));
      BB[i][j] = B[1][i]*tmp;
      A[i][j] = BB[i][j]*cond;
      if(i==j) 
        A[i][j] += 1;               
    }
  }
}

/*!
 * Generates a complex interval matrix of dimension n, tightly enclosing a 
 * complex matrix with choosen condition
 *
 * \param n Dimension of matrix
 * \param A The generated matrix
 * \param cond The desired condition of the midpoint matrix
 */
void GenCIMat(int n, cimatrix &A, real cond) {
  A = cmatrix(n,n);
  cimatrix B(1,n);
  cimatrix BB(n,n);
  cidotprecision accu;
  interval norm(0,0);

  SetRe(B,imatrix(rando(1,n)));
  SetIm(B,imatrix(rando(1,n)));

  accu = 0.0;
  for(int i=1 ; i<=n ; i++) {
    accu += sqr(abs(B[1][i]));
  }
  norm = power(sqrt(abs(rnd(accu))),-1);

  B = B*norm;

  for(int i=1 ; i<=n ; i++) {
    for(int j=1 ; j<=n ; j++) {      
      cinterval tmp(Re(B[1][j]), -Im(B[1][j]));
      BB[i][j] = B[1][i]*tmp;
      A[i][j] = BB[i][j]*cond;
      if(i==j) 
        A[i][j] += 1;               
    }
  }

}

//!Computes the transpose of a real matrix 
// static rmatrix transp ( const rmatrix& A )
// {                                            
//   int      i,j;
//   rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));
// 
//   for (i = Lb(A,1); i <= Ub(A,1); i++) 
//     for (j = Lb(A,2); j <= Ub(A,2); j++) 
//       res[j][i] = A[i][j];
//       
//   return res;
// }

//!Computes the maximum norm of a real matrix
static real MaxNorm ( const rmatrix& M )
{
  int          i, j;
  real         Max, Tmp;
  dotprecision Accu;

  Max = 0.0;
  for (i = Lb(M,1); i <= Ub(M,1); i++) {
    // Compute Tmp = #>( sum(j=1(1)n,|M_ij|) )
    //----------------------------------------
    Accu = 0.0;
    for (j = Lb(M,2); j <= Ub(M,2); j++)
      Accu += abs(M[i][j]);
    Tmp = rnd(Accu,RND_UP);
    if (Tmp > Max) Max = Tmp;
  }
  return Max;
} // MaxNorm

/*!
 * Alternative algorithm for generation of ill conditioned matrices.
 * This algorithm is very slow, condition number may be worse than requested.
 * All elements of the generated matrix are integers and thus can be represented
 * as floting point numbers.
 *
 * \param n Dimension of the matrix (n x n)
 * \param A The generated matrix
 * \param cond The desired condition number
 *
 */
void GenMat2(int n, rmatrix &A, real cond) {
  real log10cond = log(cond);
  int k = 0;
  rmatrix R, tmp; int Err;
  A = rmatrix(n,n);
  rmatrix L(n,n);

  while (log10cond > 14) {
    k = k+1;
    log10cond = log10cond/2;
  }
  
  real log10c = 0;
  real e = 1;
  int index = 0;

  while(abs(log10cond-log10c)>1) {
    index = index+1;
    if(index>5) {
      index = 0;
      if(log10c<log10cond) 
        e = 0.5 + (e-0.5)*1.1;
      else
        e = 0.5 + (e-0.5)*0.9;
    }

    tmp = randn(n,n)*e;
    for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=n ; j++) {
        if(j>i)
          A[i][j] = 0.0;
        else if(i==j)
          A[i][j] = 1.0;
        else
          A[i][j] = round(_double(tmp[i][j]));
      }
    }

    tmp = randn(n,n)*e;
    for(int i=1 ; i<=n ; i++) {
      for(int j=1 ; j<=n ; j++) {
        if(j>i)
          L[i][j] = 0.0;
        else if(i==j)
          L[i][j] = 1.0;
        else
          L[i][j] = round(_double(tmp[i][j]));
      }
    }

    A = transp(L)*A;
    MatInv(A,R,Err);
    log10c = log(MaxNorm(A)*MaxNorm(R));
  }
  
  for(int i=1 ; i<=k ; i++) {
    A = A*transp(A);
  }

  A = L*A;

  MatInv(A,R,Err);
  cout << "Kondition:" << MaxNorm(A)*MaxNorm(R) << endl;  
}


