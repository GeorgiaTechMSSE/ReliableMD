//=======================================================================
// Example: Parametric Linear Systems - input data from a file
//=======================================================================
//-----------------------------------------------------------------------
// Structure of the data file:
//    'n'        : integer > 0 specifying the dimension of the system.
//    'k'          integer > 0 specifying the number of parameters. 
//    'SharpC'   : 0 or 1; 1 for computing sharp enclosure of the
//                 iteration  matrix C; 0 - otherwise.
//    'Eps'	 : constant for the epsilon inflation.
//    'Inner'	 : 0 or 1; 1 for computing inner estimation.
//    'Ap'    	 : matrix of the system in factorized 
//                 representation Ap = (A0, A1, . . ., Ak)'.     
//    'bp'       : right-hand side of the system in a factorized
//                 matrix representation bp = (b0, b1, ... bk).
//    'p'       : vector of the interval values for the k parameters. 
//                
//------------------------------------------------------------------------

#include <iostream>
#include <parlinsys.hpp>  // Parametric Linear System Solver
#include <ivector.hpp>    // Interval vector arithmetic

using namespace std;
using namespace cxsc;


static real Sharpness(interval& xx, interval& yy)
{
   real sh;

   if (diam(xx) == 0.0) sh = 1.0;
   else if (Inf(yy) > Sup(yy)) sh = 0.0;
        else 
          if (!(Inf(xx) <= Inf(yy) && Sup(yy) <= Sup(xx)))
            sh = -100;
          else sh = diam(yy)/diam(xx);

   return sh;
}


int main()
{
	int n, k, i, j, Err, SharpC, Inner;	
	real Eps;
	
	cout << SetPrecision(10, 12) << Scientific << endl;   //Output format 
	   
	cin >> n; 
	cin >> k;
	cin >> SharpC;
	cin >> Eps;
	cin >> Inner;
	
	rmatrix Ap(n*(k+1), n), bp(n, k+1), A(n, n);
	ivector p(k+1), xx, yy;
	
	for(i = 1, j = 0; i <= k+1; ++i,j = n*(i-1))
	{
	    cin >> A;
	    Ap(j+1, i*n, 1, n) = A;
	}
	
	cin >> bp;
	
	p[1] = interval(1);
	for(i = 2; i <= k+1; i++)  cin >> p[i];

	cout << "n = " << n << "  ";
	if (SharpC) cout << "Sharp  eps =  " << Eps << endl;
	 else cout << "Rough  eps = " << Eps << endl;


	ParLinSolve(Ap, bp, p, SharpC, Eps, Inner, xx, yy, Err);
	
	if(!Err)
	 { cout << "Verified solution enclosure found in:" << endl << xx <<  endl;
           if (Inner)
           { 
	    rvector sharpness(n);
	   
	    cout << "Inner Estimation:" << endl;
            for (i=1; i<=n; i++)
              if (Inf(yy[i]) > Sup(yy[i]))
                {
		   cout << "In[" << i << "] = Empty" << endl;
		   sharpness[i] = Sharpness(xx[i], yy[i]);
		}
	      else 
		{
		  cout << yy[i] << endl;  
		  sharpness[i] = Sharpness(xx[i], yy[i]);
		};
	    cout << endl << "Sharpness of the solution enclosure:" << endl
		<< SetPrecision(5, 2) << sharpness << endl;
	   } // end Inner
	 }
	else
	  cout << ParLinSolveErrMsg(Err) << endl;
	
	cout << endl;	
	return 0;
}

