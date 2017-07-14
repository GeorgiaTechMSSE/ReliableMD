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


// MPI Example for C-XSC-MPI-Communication
// Author: Markus Grimmer
// Date: 5.3.2006
//
// Incorporated Data Types: real
// Incorporated Communication Routines: MPI_Send, MPI_Recv

// Sending C-XSC objects around a ring


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime> //For timing

#include "cxsc_mpicomm.h"

using namespace std;
using namespace cxsc;

/* Catch MPI errors */

void MPI_Fehler(int error, char *mpi_routine, char *subroutine, int rank) {
  char error_message[MPI_MAX_ERROR_STRING];
  int message_length;
  
  cout << "Errorr " << error << " in MPI routine " << mpi_routine  
       << " in subroutine " << subroutine << " in process " << rank << endl;
  MPI_Error_string(error, error_message, &message_length);
  cout << "MPI Error message: " << error_message << endl;
  MPI_Abort(MPI_COMM_WORLD,1);
}

void start_clock(clock_t& t1)
{
   t1= clock();
   if (t1 == clock_t(-1))  //terminate if timer does not work properly
   {
      cerr << "Sorry, no clock\n";
      exit(1);
   }
   return;
}

void print_time_used(clock_t t1, ofstream& out)
{
   // Writes to an output FILE

   clock_t t2= clock();
   if (t2 == clock_t(-1))
   {
     cerr<< "Sorry, clock overflow\n";
     exit(2);
   }

   out << "Time  used: " << 1000*double(t2-t1)/CLOCKS_PER_SEC
        << " msec" << endl;
   return;
} 


int main(int argc, char *argv[]) {

  // Variables for communication

  int procs;  // total number of procs
  int mypid;  // process ID
  int leftpid; // left neighbor
  int rightpid; // right neighbor
  int length = 100; // max. message length
  int error;
  MPI_Status status;

  /* MPI Init */
  if ((error=MPI_Init(&argc, &argv)) != MPI_SUCCESS)
    MPI_Fehler(error,"MPI_Init", "main",-1);
    
  if ((error=MPI_Comm_rank(MPI_COMM_WORLD, &mypid)) != MPI_SUCCESS)
    MPI_Fehler(error,"MPI_Comm_rank", "main",mypid);
    
  if ((error=MPI_Comm_size(MPI_COMM_WORLD, &procs)) != MPI_SUCCESS)
    MPI_Fehler(error,"MPI_Comm_size", "main",procs);
 
  /* Auffinden der Nachbarn im Ring */
  if (procs <= 2){
    cout << "Too few processors" << endl;
    cout << "Abort.\n";
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
 
  rightpid = (mypid+1)%procs;
  leftpid  = (mypid-1+procs)%procs;

  // MPI-CXSC data type definition (see cxsc_mpicomm package):

  MPI_Define_CXSC_Types();

  // ------------------------------------------------------------------------
  // Define and open output file
  // 
  // This example works with output files. If you don't like to have output 
  // from every process or if you don't like the file handling inside the program,
  // remove the corresponding sections (i.e. the following and the one at the very
  // end of this file...
  //
  // The file handling works on appropriate Linux/Un*x systems only.
  // ------------------------------------------------------------------------

      ofstream ausg;
      ostringstream argvns;
      ostringstream ausgdats;
      string ausgabe, ausgabedatei;
      string argvnull;

      // The following is to handle file names with full path specification
      // Check if the practice below works on your system!
      // You can leave this part out if you don't work with file names as input
      // parameters...

      int substr_start=0;
      int substr_len=0;

      argvns << argv[0];
      argvnull=argvns.str();
      substr_start=argvnull.find_last_of("/")+1;
      if ((substr_start>argvnull.size())||(substr_start<0))
        substr_start=0;
      substr_len=argvnull.size()-substr_start;
      ausgabe=argvnull.substr(substr_start,substr_len);

      // *************************************************************************
      // Determine output file directory
      // *************************************************************************

      ostringstream datadirs, resultdirs;
      string datadir, resultdir;

      if (argc>2) 
      {
        datadirs << argv[1];
        resultdirs << argv[2];
        datadir = datadirs.str();
        resultdir = resultdirs.str();
      }
      else        // Default values ************************************************
      {
        datadir="/data/username";         // an example, reset it according 
                                          // to your needs 
        resultdir="/home/username/some_subdirectory";   
      }
      // *************************************************************************

      // The following could be helpful on your system, use it if necessary:
      /*
      string syscmd1= "if test -d "
                   +datadir+"; then echo \"Directory exists\"; else mkdir "
                   +datadir+"; fi;";
      system(syscmd1.c_str());
      */
      //(end couldbe)
 
      ausgdats << datadir << "/" << ausgabe << "." << mypid << ".out";
      ausgabedatei=ausgdats.str();

      // Open output file
 
      ausg.open(ausgabedatei.c_str());
      
  // ------------------------------------------------------------------------
  // (end output file) 
  // ------------------------------------------------------------------------

  // CXSC-Variables for all processes:

  int width=23, digits=18;

  // Data:
 
  string inp, inpcopy;  
  real r1, rlb, rub;
  interval ri1;
 
  //...insert more...

  // Every process adds something:

  real rinc=1.0;         
  interval iinc(0.5,0.5);

for (int dist=1000;dist<=10000;dist*=10)   //(For timing: First go 1000 laps,
                                           // then 10000, see below)
{

 if (mypid==0) {      // process 0

   inp="0.1234567890123456789012345678901234567890";   
   inpcopy=inp;
   inp >> RndUp >> rub;
   inpcopy >> RndDown >> rlb;

   r1=rlb;   
   ri1=_interval(rlb,rub);

   // Test output

   ausg << "r1= " << endl << SetPrecision (width,digits) << r1 << endl;
   ausg << "ri1= " << endl << SetPrecision (width,digits) << ri1 << endl;

 } else {

   // Nothing to do...

 }
 
 clock_t t;        //defined in <ctime> 
 if (mypid==0) {
   start_clock(t); //store starting time in t
 }

 int laps=dist; // see above

 // ***** Communikation ******************************************

 for(int lapcount=0; lapcount<laps; lapcount++)
 {

   if (mypid==0) {

     // Send proc. 0
 
     if ((error=MPI_Send(r1, rightpid, 0, MPI_COMM_WORLD)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Send", "main",mypid);
     if ((error=MPI_Send(ri1, rightpid, 0, MPI_COMM_WORLD)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Send", "main",mypid);
  
     // Receive proc. 0

     if ((error=MPI_Recv(r1, leftpid, 0, MPI_COMM_WORLD, &status)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Recv", "main",mypid);
     if ((error=MPI_Recv(ri1, leftpid, 0, MPI_COMM_WORLD, &status)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Recv", "main",mypid);

   } else {

     // Receive proc. 1..n

     if ((error=MPI_Recv(r1, leftpid, 0, MPI_COMM_WORLD, &status)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Recv", "main",mypid);
     if ((error=MPI_Recv(ri1, leftpid, 0, MPI_COMM_WORLD, &status)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Recv", "main",mypid);

     // Do some computation...

     r1+=rinc;
     ri1+=iinc;

     // Send proc. 1..n

     if ((error=MPI_Send(r1, rightpid, 0, MPI_COMM_WORLD)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Send", "main",mypid);
     if ((error=MPI_Send(ri1, rightpid, 0, MPI_COMM_WORLD)) != MPI_SUCCESS)
       MPI_Fehler(error,"MPI_Send", "main",mypid);
   }

 } // end for
 
 if (mypid==0) {

   print_time_used(t, ausg); // print time difference to output stream "ausg"

   ausg << "After communication:" << endl;

   ausg << "r1= " << endl << SetPrecision (width,digits) << r1 << endl;
   ausg << "ri1= " << endl << SetPrecision (width,digits) << ri1 << endl;

 } else {

 }

} //for dist

 // ------------------------------------------------------------------------
 // For all processes:
 // ------------------------------------------------------------------------


 // ------------------------------------------------------------------------
 // Close and move output file
 // ------------------------------------------------------------------------
 
      ausg.close();
          
      string syscmd_mvres= "mv "+ ausgabedatei + " " + resultdir;
      system(syscmd_mvres.c_str());

 // ------------------------------------------------------------------------
   
  if ((error=MPI_Finalize()))
     MPI_Fehler(error,"MPI_Finalize","main",mypid);
     
  return 0;


}

