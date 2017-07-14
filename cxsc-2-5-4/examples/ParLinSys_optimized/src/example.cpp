/*
 *  Example program for the usage of the parametric linear system solver
 *  based on the Krawczyk Operator
 */

#include <iostream>
#include <string>
#include <parlinsys.hpp>  // Parametric Linear System
#include <ivector.hpp>      
#include <timer.hpp>
#include <vector>

using namespace std;
using namespace cxsc;

/*
 * Example:
 * 
 * Code for solving the system
 * 
 * (3 p p)     (1)
 * (p 3 p) x = (0)
 * (p p 3)     (0)
 * 
 * with p in [0,2]
 */

int main() {
  //STL-vector holding coefficient matrices for each parameter
  vector<CoeffMatrix> Ap;
  
  //Matrix A used to temporarily store each coefficient matrix
  //matrix can be either dense or sparse (data type srmatrix)
  rmatrix A(3,3);
  
  //Define non-parametric part of system matrix
  A = 0.0; //set all elements to 0.0
  A[1][1] = 3.0;
  A[2][2] = 3.0;
  A[3][3] = 3.0;
  //Add non-parametric part to Ap by wrapping A inside a CoeffMatrix Object
  Ap.push_back(CoeffMatrix(A));
  
  //Part of the system matrix depending on the parameter (here only one parameter)
  A = 1.0; //set all elements to 1.0
  A[1][1] = 0.0;
  A[2][2] = 0.0;
  A[3][3] = 0.0;
  //Add coefficient matrix for parameter to Ap by wrapping A inside a CoeffMatrix Object
  Ap.push_back(CoeffMatrix(A));
   
  //Define right hand side
  //One column per parameter + one column for non-parameter dependent part
  //Here we need two columns
  rmatrix bp(3,2); //bp could also be sparse
  bp=0.0; //Set all elements to 0.0
  bp[1][1]=1.0;

  //Define vector of parameter intervals
  ivector p(2);
  p[1] = 1.0; //Entry for non-parameter dependent part, usually 1.0
  p[2] = interval(0,2); //Interval for parameter p
	
  //Vector to store the outer enclosure
  ivector x(3);
  //Variable for Error code
  int Err;

  //Create an instance of configuration struct to override some default values
  parlinsysconfig cfg;
  cfg.refinement = true; //Activate iterative refinement
  cfg.msg = true;        //Activate status message output
  
  //Call parametric solver
  cout << "Starting parametric solver..." << endl << endl;
  ParLinSolve(Ap, bp, p, x, Err, cfg);	
  cout << endl;
  
  if(!Err) {
    //No error occured
    cout << "Verified outer enclosure: " << endl << x << endl;
  } else {
    //An error occured, put out corresponding error message
    cout << ParLinSolveErrMsg(Err) << endl;
  }

  return 0;
}

