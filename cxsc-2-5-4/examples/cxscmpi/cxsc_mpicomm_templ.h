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

// MPI-Routines for C-XSC data types
//
// Author: Markus Grimmer
// Date: 5.3.2006
//
// Incorporated Data Types: real, interval, rvector, ivector, rmatrix, imatrix,
//                          complex, cinterval, cvector, civector, cmatrix, cimatrix,
//                          l_real, l_interval, l_complex,
//                          l_rvector, l_ivector, l_rmatrix, l_imatrix
//                          dotprecision, idotprecision, 
//                          cdotprecision, cidotprecision
// Incorporated Communication Routines: MPI_Send, MPI_Bsend, MPI_Ssend, MPI_Rsend
//                                      MPI_Isend, MPI_Ibsend, MPI_Issend, MPI_Irsend
//                                      MPI_Recv
// Incorporated Utility Routines: MPI_Pack, MPI_Unpack

#ifndef _CXSC_MPICOMM_INCLUDE
#define _CXSC_MPICOMM_INCLUDE

#include <mpi.h>

#include <real.hpp>
#include <interval.hpp>
#include <rvector.hpp>
#include <ivector.hpp>
#include <rmatrix.hpp>
#include <imatrix.hpp>

#include <complex.hpp>
#include <cinterval.hpp>
#include <cvector.hpp>
#include <civector.hpp>
#include <cmatrix.hpp>
#include <cimatrix.hpp>

#include <l_real.hpp>
#include <l_interval.hpp>
#include <l_rvector.hpp>
#include <l_ivector.hpp>
#include <l_rmatrix.hpp>
#include <l_imatrix.hpp>

#include <l_complex.hpp>

#include <dot.hpp>
#include <idot.hpp>
#include <cdot.hpp>
#include <cidot.hpp>

using namespace cxsc;


const int MPI_CXSC_BUFFERLEN=1000000; // Buffer length for MPI_PACK
                                      // sets bounds to the size of 
                                      // objects that can be sent!
                                      // Reset according to your needs!
                                      // You can also turn it into a non-const
                                      // variable.

// Definition of new MPI data types:

extern bool MPI_CXSC_TYPES_DEFINED;
extern MPI_Datatype MPI_CXSC_REAL;
extern MPI_Datatype MPI_CXSC_COMPLEX;
extern MPI_Datatype MPI_CXSC_INTERVAL;
extern MPI_Datatype MPI_CXSC_CINTERVAL;

int MPI_Define_CXSC_Types();

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
int MPI_Bsend(T&, int, int, MPI_Comm);
template<class T>
int MPI_Ssend(T&, int, int, MPI_Comm);
template<class T>
int MPI_Rsend(T&, int, int, MPI_Comm);

template<class T>
int MPI_Isend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Ibsend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Issend(T&, int, int, MPI_Comm, MPI_Request*);
template<class T>
int MPI_Irsend(T&, int, int, MPI_Comm, MPI_Request*);

template<class T>
int MPI_Recv(T&, int, int, MPI_Comm, MPI_Status*);

#endif
