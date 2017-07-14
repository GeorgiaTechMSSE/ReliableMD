/*
**  CXSC is a C++ library for eXtended Scientific Computing 
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2006 Wiss. Rechnen/Softwaretechnologie
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

#include <iostream>
#include <imatrix.hpp>
#include "cxsc_mpicomm.hpp"

using namespace std;
using namespace cxsc;

int main(int argc, char *argv[]) {

  int mypid;  // process ID
  MPI_Status status;

  // MPI Init
  MPI_Init(&argc, &argv);

  //Determine own rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mypid);

  //initialize communication buffer  
  //This function MUST be called before the first communication!
  //Alternativly you can also call init_CXSC_MPI(int n), where n
  //is the size of the buffer in byte
  init_CXSC_MPI();

  imatrix A(5,5);
  A = 0.0;
 
  if (mypid==0) {

     //Fill matrix in process 1
     for(int i=1 ; i<=5 ; i++)
       A[i][i] = interval(i,i+1);
  
     cout << "Matrix stored by process 0:" << endl;
     cout << A << endl;

     // Send matrix A from process 0 to process 1
     MPI_Send(A, 1, 0, MPI_COMM_WORLD);

  } else if(mypid==1) {

     //Process 1 receives A from process 0
     MPI_Recv(A, 0, 0, MPI_COMM_WORLD, &status);

     cout << "Matrix received by process 1:" << endl;
     cout << A << endl;
  }
   
  MPI_Finalize();

  //Delete the communication buffer.
  //This function should be called as soon as there are no
  //more C-XSC MPI communications necessary to free the memory used
  //for the buffer
  finish_CXSC_MPI();
     
  return 0;
}

