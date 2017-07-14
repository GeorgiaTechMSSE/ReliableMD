/*
 *  Example program for the usage of the parametric linear system solver
 *  based on alternative algorithm by Neumaier/Pownuk
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
 * (p1+p2 p1-p2) (x1)  (6)
 *                   =
 * (p1-p2 p1+p2) (x2)  (6)
 * 
 * with p1,p2 in [1-delta,1+delta], 0<=delta<1
 * 
 * System can be transformed into the form (A^T D A)x = Fb
 * with
 * A = (1  1)   D=(p1  0)  F = I   b = (6)
 *     (1  -1)    (0  p2)              (6)
 */

int main() {
  //Define matrices A,F of system

  rmatrix A(2,2); //A could also be sparse
  A = 1.0; 
  A[2][2] = -1.0;
  
  //F is unity matrix of same size as A
  rmatrix F(Id(A));

  //Define b
  ivector b(2);
  b[1] = b[2] = 6.0;

  //Define D with delta = 0.5
  simatrix D(2,2);
  real delta = 0.5;
  D[1][1] = interval(1-delta,1+delta);
  D[2][2] = interval(1-delta,1+delta);

  //Variable for error code
  int Err;
  //Solution vector
  ivector x(2);

  //Create instance of configuration variable to use non standard configuration
  parlinsysnpconfig cfg;
  cfg.verified = true; //activate verification
  cfg.msg = true;      //activate status messages
  cfg.maxIter = 5;     //perform only 5 iteration steps

  cout << "Computing result using Neumaier/Pownuk algorithm" << endl << endl;
  //Start parametric solver
  ParLinSolve(A, D, x, F, b, Err, cfg);

  if(Err) {
    cout << endl << "An error occured: " << ParLinSolveErrMsg(Err) << endl;
  } else {
    cout << endl << "Verified outer enclosure: " << endl << x << endl;   
  }


  //Now try to solve the system using the Krawczyk operator
  
  //Define coefficient matrices
  vector<CoeffMatrix> Apv;
  //Non-parametric part
  A=0.0;
  Apv.push_back(CoeffMatrix(A));
  //For p1
  A=1.0;
  Apv.push_back(CoeffMatrix(A));
  //For p2
  A[1][2] = -1.0;
  A[2][1] = -1.0;
  Apv.push_back(CoeffMatrix(A));

  //Define right hand side
  rmatrix bp(2,3);
  bp = 0.0;
  bp[Col(1)]=6.0;
      
  //Define parameter intervals
  ivector p(3);
  p[1] = 1.0;
  p[2] = p[3] = interval(1-delta, 1+delta);
      
  //Start solver
  cout << endl << "Trying to solve system with Krawczyk-Operator..." << endl;
  ParLinSolve(Apv, bp, p, x, Err);
  if(!Err) {
    //No error occured
    cout << endl << "Verified outer enclosure: " << endl << x << endl;
  } else {
    //An error occured, put out corresponding error message
    cout << ParLinSolveErrMsg(Err) << endl;
  }
  
  
  return 0;
}
