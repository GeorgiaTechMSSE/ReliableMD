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

/*
** Example program for real interval systems. 
**
** Solves a random system that is stored in process 0 at the beginning and 
** distributed from there. Code for distribution using function pointers
** as an alternative is also provided but commented out.
**/

#include "mpi.h"
#include <cxscplss.hpp>
#include "testgen.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "sys/time.h"


using namespace cxsc;
using namespace std;

double GetTime() {
   struct timeval _tp;

   gettimeofday(&_tp,0);
   
   return _tp.tv_sec + _tp.tv_usec / 1000000.0;
}

static int m = 1000;
static int n = 1000;

//Function defining A
void makeA(int i, int j, interval& r) {
  r = sqrt(2.0/(n+1))*sin((i*j*M_PI)/(n+1));
  r = r + r*interval(-1e-12,1e-12);
}

//Function defining b
void makeb(int i, int j, interval& r) {
  if(i==j) r=1.0; else r=0.0;
}

//Computes relative error of solution
static real RelErr(const ivector &x) {
  dotprecision accu(0);
  
  for(int i=1 ; i<=VecLen(x) ; i++) {
    real d = 0.5 * diam(x[i]);
    real m = abs(mid(x[i]));
    if(m == 0) {
      m = 1;
    }
    
    accu += d / m;
  }

  return rnd(accu) / VecLen(x);
}

int main (int argc, char** argv)
{
  int  Err;
  clock_t t;
  int myid, np;    
  int K=2;
  
  srand(time(NULL));
  
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD, &np);  
  
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  

  ostringstream s;
  s << "output" << myid << ".txt";  

  //Number of right hand sides
  int rhs = 1;

  //Dot product precision K, Dimension n and number of right hand
  //sides rhs can be given as arguments
  if (argc == 4) {
    K = atoi(argv[1]);
    n = atoi(argv[2]);
    m = n;
    rhs = atoi(argv[3]);
  } 

  int dim = rhs/np + 1;
  if(rhs % np <= myid) dim--;

  imatrix A(m,n);
  imatrix b(m,rhs);
  imatrix x(n,dim);

  //b is unity matrix (blown up a little)
  b = 0.0;
  for(int i=1 ; i<=rhs ; i++) {
     b[i][i] = 1.0;
  }
  b = b + b*interval(-1e-12,1e-12);

  //A is a random matrix
  A = randn(m,n);
  A = A + A*interval(-1e-12,1e-12);

  //Output file stream for status messages (and final result)
  ofstream out(s.str().c_str());  
  cout << SetPrecision(23,15) << Scientific;   // Output format

  MPI_Barrier(MPI_COMM_WORLD);

  double start, end;
 
  //Struct for configuration of solver
  //uncomment to set non default values
  struct plssconfig cfg;
  cfg.K = K;                         //Dot product precision
  cfg.lssparts = LSS_BOTH_PARTS;     //Solver stages to use (LSS_BOTH_PARTS | LSS_ONLY_PART_ONE | LSS_ONLY_PART_TWO)
  //cfg.threads = -1;                //Number of threads for OpenMP (-1 = OMP default)
  //cfg.nb = 256;                    //Blocksize for ScaLAPACK
  //cfg.maxIterResCorr = 5;          //Maximum number of iterations during residual correction
  //cfg.maxIterVer = 5;              //Maximum number of iterations during the verification step
  //cfg.refinement = false;          //Perform an iterative refinement?
  //cfg.maxIterRef = 5;              //Maximum number of iterations during the refinement step
  //cfg.epsVer = 0.1;                //Epsilon for the verification step
  //cfg.epsRef = 1e-5;               //Epsilon fot the refinement step
  //cfg.matrixMode = false;          //Eenable Matrix Mode (forces K=1)

  start = GetTime();
  //Version 1 using matrices A and b locally stored in process 0
  pilss(A,b,x,m,n,np,myid,Err,out,cfg);
  //Version 2 using function pointers to define A and b (allowing larger systems)
  //pilss(makeA,makeb,x,m,n,rhs,np,myid,Err,out,cfg);
  end = GetTime();

  printf("Process %d:  %lf s\n", myid, end-start);
  
  if (!Err && K!=1) {
    //When not using matrix mode, result is column cyclically distributed among processes
    for(int i=1 ; i<=dim ; i++) {
      out << endl << "Result for right hand side " << (myid+1)+(i-1)*np << ": " << endl;
      out << "Relative error: "  << RelErr(x[Col(i)]) << endl;
      //uncomment to put out complete result vector to output file
      //out << "Result: " << endl << x[Col(i)] << endl << endl;
    }
  } else if(!Err && K==1) {
    //When using matrix mode and data was stored in process 0, result is stored completely in process 0
    //Attention: when using function pointer method, result will instead be distributed in two-dimensional block cyclic
    //distribution!
    for(int i=1 ; i<=rhs ; i++) {
      out << endl << "Result for right hand side " << i << ": " << endl;
      out << "Relative error: "  << RelErr(x[Col(i)]) << endl;
      //uncomment to put out complete result vector to output file
      //out << "Result: " << endl << x[Col(i)] << endl << endl;
    }
  } else out << "Error: " << LinSolveErrMsg(Err) << endl;
  
  MPI_Finalize();
  
  return 0;
}

