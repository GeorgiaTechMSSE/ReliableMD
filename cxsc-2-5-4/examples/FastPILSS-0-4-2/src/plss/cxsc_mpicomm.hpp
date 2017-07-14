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
**  
**  Author: Michael Zimmer
**
**  This software is based on:
**    - Module LinSys of the C-XSC-Toolbox
**      Authors: Rolf Hammer, Matthias Hocks, Dietmar Ratz
**    - Self-verifying solver for a dense system of linear equations
**      Authors: Carlos Holbig, Walter Kraemer, Paulo Sergio Morandi Junior,
**               Bernardo Frederes Kramer Alcalde, 
**    - Parallel Interval Linear System Solver
**      Author: Markus Grimmer
**/


#ifndef _CXSC_MPICOMM_INCLUDE
#define _CXSC_MPICOMM_INCLUDE

#include <mpi.h>


#include <real.hpp>
#include <interval.hpp>
#include <complex.hpp>
#include <cinterval.hpp>
#include <l_real.hpp>
#include <l_interval.hpp>
#include <l_complex.hpp>
#include <iostream>

#include <vector>


namespace cxsc {

class rvector;
class ivector;
class cvector;
class civector;
class srvector;
class sivector;
class scvector;
class scivector;
class l_rvector;
class l_ivector;
class l_cvector;
class l_civector;
class rmatrix;
class imatrix;
class cmatrix;
class cimatrix;
class srmatrix;
class simatrix;
class scmatrix;
class scimatrix;
class l_rmatrix;
class l_imatrix;
class l_cmatrix;
class l_cimatrix;
class dotprecision;
class idotprecision;
class cdotprecision;
class cidotprecision;


extern int MPI_CXSC_BUFFERLEN;    // Buffer length for MPI_PACK
                                  // sets bounds to the size of 
                                  // objects that can be sent!

//static char commbuffer[MPI_CXSC_BUFFERLEN]; 
extern char* commbuffer;  // Allocated by init function


//Initialize MPI communication buffer with default length
void init_CXSC_MPI();
//Initialize MPI communication buffer with length n
void init_CXSC_MPI(int n);
//Delete MPI communication buffer
void finish_CXSC_MPI();


// Definition of new MPI data types:
extern bool MPI_CXSC_TYPES_DEFINED;
extern MPI_Datatype MPI_CXSC_REAL;
extern MPI_Datatype MPI_CXSC_COMPLEX;
extern MPI_Datatype MPI_CXSC_INTERVAL;
extern MPI_Datatype MPI_CXSC_CINTERVAL;

int MPI_Define_CXSC_Types();

// **************************************************************
// Communication functions for STL-vector
// **************************************************************

template<class T> 
inline int MPI_Pack (std::vector<T> &, void *, int, int *,  MPI_Comm);
template<class T> 
inline int MPI_Unpack (void*, int, int *, std::vector<T> &, MPI_Comm);

template<class T>
inline int MPI_Send(std::vector<T> &, int, int, MPI_Comm);
template<class T>
inline int MPI_Bsend(std::vector<T> &, int, int, MPI_Comm);
template<class T>
inline int MPI_Ssend(std::vector<T> &, int, int, MPI_Comm);
template<class T>
inline int MPI_Rsend(std::vector<T> &, int, int, MPI_Comm);

template<class T>
inline int MPI_Isend(std::vector<T> &, int, int, MPI_Comm, MPI_Request*);
template<class T>
inline int MPI_Ibsend(std::vector<T> &, int, int, MPI_Comm, MPI_Request*);
template<class T>
inline int MPI_Issend(std::vector<T> &, int, int, MPI_Comm, MPI_Request*);
template<class T>
inline int MPI_Irsend(std::vector<T> &, int, int, MPI_Comm, MPI_Request*);

template<class T> 
inline int MPI_Recv (std::vector<T> &, int, int, MPI_Comm, MPI_Status*);
//@}

template<class T>
inline int MPI_Bcast(std::vector<T>&, int, MPI_Comm);
template<class T>
inline int MPI_Bcast(std::vector<T>&, int, int, int, int, int, MPI_Comm);

// **************************************************************
// Pack and Unpack for C-XSC data types
// **************************************************************

// Scalar data types: (not necessary since these are defined as 
// MPI data types - only for completeness purposes)

int MPI_Pack (real&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, real&, MPI_Comm);
int MPI_Pack (interval&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, interval&, MPI_Comm);
int MPI_Pack (complex&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, complex&, MPI_Comm);
int MPI_Pack (cinterval&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cinterval&, MPI_Comm);

// Multiple precision scalar data types: 

int MPI_Pack (l_real&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_real&, MPI_Comm);
int MPI_Pack (l_interval&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_interval&, MPI_Comm);
int MPI_Pack (l_complex&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_complex&, MPI_Comm);

// Vector data types:

int MPI_Pack (rvector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, rvector&, MPI_Comm);
int MPI_Pack (ivector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, ivector&, MPI_Comm);
int MPI_Pack (cvector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cvector&, MPI_Comm);
int MPI_Pack (civector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, civector&, MPI_Comm);
int MPI_Pack (rvector&, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, rvector&, int, int, MPI_Comm);
int MPI_Pack (ivector&, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, ivector&, int, int, MPI_Comm);
int MPI_Pack (cvector&, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cvector&, int, int, MPI_Comm);
int MPI_Pack (civector&,  int, int,void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, civector&, int, int, MPI_Comm);


//Sparse
int MPI_Pack (srvector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, srvector&, MPI_Comm);
int MPI_Pack (sivector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, sivector&, MPI_Comm);
int MPI_Pack (scvector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, scvector&, MPI_Comm);
int MPI_Pack (scivector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, scivector&, MPI_Comm);

// Multiple precision Vector data types: 

int MPI_Pack (l_rvector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_rvector&, MPI_Comm);
int MPI_Pack (l_ivector&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_ivector&, MPI_Comm);

// Matrix data types:

int MPI_Pack (rmatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, rmatrix&, MPI_Comm);
int MPI_Pack (imatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, imatrix&, MPI_Comm);
int MPI_Pack (cmatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cmatrix&, MPI_Comm);
int MPI_Pack (cimatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cimatrix&, MPI_Comm);
int MPI_Pack (rmatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, rmatrix&, int, int, int, int, MPI_Comm);
int MPI_Pack (imatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, imatrix&, int, int, int, int, MPI_Comm);
int MPI_Pack (cmatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cmatrix&, int, int, int, int, MPI_Comm);
int MPI_Pack (cimatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cimatrix&, int, int, int, int, MPI_Comm);


//Sparse
int MPI_Pack (srmatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, srmatrix&, MPI_Comm);
int MPI_Pack (simatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, simatrix&, MPI_Comm);
int MPI_Pack (scmatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, scmatrix&, MPI_Comm);
int MPI_Pack (scimatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, scimatrix&, MPI_Comm);
int MPI_Pack (srmatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, srmatrix&, int, int, int, int, MPI_Comm);
int MPI_Pack (simatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, simatrix&, int, int, int, int, MPI_Comm);
int MPI_Pack (scmatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, scmatrix&, int, int, int, int, MPI_Comm);
int MPI_Pack (scimatrix&, int, int, int, int, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, scimatrix&, MPI_Comm);

// Multiple precision Matrix data types: 

int MPI_Pack (l_rmatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_rmatrix&, MPI_Comm);
int MPI_Pack (l_imatrix&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, l_imatrix&, MPI_Comm);

// Dotprecision types

int MPI_Pack (dotprecision&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, dotprecision&, MPI_Comm);
int MPI_Pack (idotprecision&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, idotprecision&, MPI_Comm);
int MPI_Pack (cdotprecision&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cdotprecision&, MPI_Comm);
int MPI_Pack (cidotprecision&, void *, int, int *,  MPI_Comm);
int MPI_Unpack (void*, int, int *, cidotprecision&, MPI_Comm);

// **************************************************************
// Point-to-Point communication functions for C-XSC data types
// **************************************************************

// **************************************************************
// Template functions for Point-to-Point Communication
// **************************************************************

template<class T>
int MPI_Send(T&, int, int, MPI_Comm);
template<class T>
int MPI_Send(T&, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Send(T&, int, int, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Bsend(T&, int, int, MPI_Comm);
template<class T>
int MPI_Bsend(T&, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Bsend(T&, int, int, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Ssend(T&, int, int, MPI_Comm);
template<class T>
int MPI_Ssend(T&, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Ssend(T&, int, int, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Rsend(T&, int, int, MPI_Comm);
template<class T>
int MPI_Rsend(T&, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Rsend(T&, int, int, int, int, int, int, MPI_Comm);

template<class T>
int MPI_Isend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Isend(T&, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Isend(T&, int, int, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Ibsend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Ibsend(T&, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Ibsend(T&, int, int, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Issend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Issend(T&, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Issend(T&, int, int, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Irsend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Irsend(T&, int, int, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Irsend(T&, int, int, int, int, int, int, MPI_Comm, MPI_Request*);

template<class T>
int MPI_Recv(T&, int, int, MPI_Comm, MPI_Status*);
template<class T>
int MPI_Recv(T&, int, int, int, int, MPI_Comm, MPI_Status*);
template<class T>
int MPI_Recv(T&, int, int, int, int, int, int, MPI_Comm, MPI_Status*);

template<class T>
int MPI_Bcast(T&, int, MPI_Comm);
template<class T>
int MPI_Bcast(T&, int, int, int, int, int, MPI_Comm);
template<class T>
int MPI_Bcast(T&, int, int, int, MPI_Comm);




// **************************************************************
// Definitions for STL-vector functions
// **************************************************************


template<class T>
inline int MPI_Pack (std::vector<T>& tv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int len=tv.size();

  if (err=::MPI_Pack(&len,1,MPI_INT,buff,bufflen,pos,MC)!=MPI_SUCCESS)
    return err;

  for (int i=0; i<len; i++)
    {
      if (err=MPI_Pack(tv[i],buff,bufflen,pos,MC)!= MPI_SUCCESS)
      return err;
    }

  return err;
}

template<class T>
inline int MPI_Unpack (void* buff, int bufflen, int* pos, std::vector<T>& tv, MPI_Comm MC)
{
  int err;
  int len;
  
  if (err=::MPI_Unpack(buff, bufflen, pos, &len, 1, MPI_INT, MC)!=MPI_SUCCESS)
    return err;

  tv.resize(len);

  for (int i=0; i<len; i++)
    {
      if (err=MPI_Unpack(buff,bufflen,pos,tv[i],MC)!= MPI_SUCCESS)
      return err;
    }
  return err;
}

// =============================================================================
// Point-to-point Communication
// =============================================================================


template<class T>
inline int MPI_Send(std::vector<T>& tv, int i1, int i2, MPI_Comm MC)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Send((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
inline int MPI_Bsend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Bsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
inline int MPI_Ssend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Ssend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
inline int MPI_Rsend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Rsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
inline int MPI_Isend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Isend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
inline int MPI_Ibsend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Ibsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
inline int MPI_Issend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Issend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
inline int MPI_Irsend(std::vector<T>& tv, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the vector, then send it.

  int pos=0;
  int err;

  // Call overloaded MPI_Pack:

  if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
    return err;
  err=::MPI_Irsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
inline int MPI_Recv(std::vector<T>& tv, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int pos=0;
  int err;

  if (err=MPI_Recv((void*)commbuffer, MPI_CXSC_BUFFERLEN, MPI_PACKED, 
                   i1, i2, MC, MS)!= MPI_SUCCESS)
    return err;

  // Call overloaded MPI_Unpack

  err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, tv, MC);  
  return err;
}


// =============================================================================
// Collective Communication: MPI_Bcast
// =============================================================================


template<class T>
inline int MPI_Bcast(std::vector<T>& tv, int root, MPI_Comm MC)
{
  int pos=0;
  int len=0;
  int err;

  int mypid;
  err=MPI_Comm_rank(MC, &mypid);
  
  if (mypid==root)  
  {
    // Pack the vector, then broadcast it.
    // It is assumed that a suitable
    // MPI_Pack is available!

    // Call overloaded MPI_Pack:

    if (err=MPI_Pack(tv, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC)!=MPI_SUCCESS)
      return err;
    len=pos;
    if (err=::MPI_Bcast(&len, 1, MPI_INT, root, MC)!=MPI_SUCCESS)
      return err;
    if (err=::MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC)!=MPI_SUCCESS)
      return err;
  }
  else
  {

    // Broadcast-receive the vector, then unpack it.
    // It is assumed that a suitable
    // MPI_Unpack is available!

    if (err=::MPI_Bcast(&len, 1, MPI_INT, root, MC)!= MPI_SUCCESS)
      return err;
    if (err=::MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC)!= MPI_SUCCESS)
      return err;

    // Call overloaded MPI_Unpack:

    err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, tv, MC);
  }
  
  return err;
}



// =============================================================================
// Specializations
// =============================================================================


template<>
inline int MPI_Pack<double> (std::vector<double>& tv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int len=tv.size();

  if (err=::MPI_Pack(&len,1,MPI_INT,buff,bufflen,pos,MC)!=MPI_SUCCESS)
    return err;

  for (int i=0; i<len; i++)
    {
      if (err=::MPI_Pack(&tv[i],1,MPI_DOUBLE,buff,bufflen,pos,MC)!= MPI_SUCCESS)
      return err;
    }

  return err;
}

template<>
inline int MPI_Unpack<double> (void* buff, int bufflen, int* pos, std::vector<double>& tv, MPI_Comm MC)
{
  int err=0;
  int len=0;
  
  if (err=::MPI_Unpack(buff, bufflen, pos, &len, 1, MPI_INT, MC)!=MPI_SUCCESS)
    return err;

  tv.resize(len);

  for (int i=0; i<len; i++)
  {
    if (err=::MPI_Unpack(buff,bufflen,pos,&tv[i],1,MPI_DOUBLE,MC)!= MPI_SUCCESS)
    return err;
  }
  return err;
}

template<>
inline int MPI_Pack<int> (std::vector<int>& tv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int len=tv.size();

  if (err=::MPI_Pack(&len,1,MPI_INT,buff,bufflen,pos,MC)!=MPI_SUCCESS)
    return err;

  for (int i=0; i<len; i++)
    {
      if (err=::MPI_Pack(&tv[i],1,MPI_INT,buff,bufflen,pos,MC)!= MPI_SUCCESS)
      return err;
    }

  return err;
}

template<>
inline int MPI_Unpack<int> (void* buff, int bufflen, int* pos, std::vector<int>& tv, MPI_Comm MC)
{
  int err=0;
  int len=0;
  
  if (err=::MPI_Unpack(buff, bufflen, pos, &len, 1, MPI_INT, MC)!=MPI_SUCCESS)
    return err;

  tv.resize(len);

  for (int i=0; i<len; i++)
  {
    if (err=::MPI_Unpack(buff,bufflen,pos,&tv[i],1,MPI_INT,MC)!= MPI_SUCCESS)
    return err;
  }
  return err;
}

} //namespace cxsc

#endif
