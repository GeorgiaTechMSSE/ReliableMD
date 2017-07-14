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
// "Packed" version: Sends data in Packed format
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
// Incorporated Communication Routines: MPI_Send, MPI_Isend, 
//                                      MPI_Recv

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

// MPI-Send/Recv for real:

int MPI_Send (real&, int, int, MPI_Comm);
int MPI_Isend(real&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (real&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for complex:

int MPI_Send (complex&, int, int, MPI_Comm);
int MPI_Isend(complex&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (complex&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for interval:

int MPI_Send (interval&, int, int, MPI_Comm);
int MPI_Isend(interval&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (interval&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for cinterval:

int MPI_Send (cinterval&, int, int, MPI_Comm);
int MPI_Isend(cinterval&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (cinterval&, int, int, MPI_Comm, MPI_Status* );

// **************************************************************

// MPI-Send/Recv for l_real:

int MPI_Send (l_real&, int, int, MPI_Comm);
int MPI_Isend(l_real&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_real&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for l_interval:

int MPI_Send (l_interval&, int, int, MPI_Comm);
int MPI_Isend(l_interval&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_interval&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for l_complex:

int MPI_Send (l_complex&, int, int, MPI_Comm);
int MPI_Isend(l_complex&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_complex&, int, int, MPI_Comm, MPI_Status* );

// **************************************************************

// MPI-Send/Recv for rvector:

int MPI_Send (rvector&, int, int, MPI_Comm);
int MPI_Isend(rvector&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (rvector&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for cvector:

int MPI_Send (cvector&, int, int, MPI_Comm);
int MPI_Isend(cvector&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (cvector&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for ivector:

int MPI_Send (ivector&, int, int, MPI_Comm);
int MPI_Isend(ivector&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (ivector&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for civector:

int MPI_Send (civector&, int, int, MPI_Comm);
int MPI_Isend(civector&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (civector&, int, int, MPI_Comm, MPI_Status* );

// **************************************************************

// MPI-Send/Recv for rmatrix:

int MPI_Send (rmatrix&, int, int, MPI_Comm);
int MPI_Isend(rmatrix&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (rmatrix&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for cmatrix:

int MPI_Send (cmatrix&, int, int, MPI_Comm);
int MPI_Isend(cmatrix&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (cmatrix&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for imatrix:

int MPI_Send (imatrix&, int, int, MPI_Comm);
int MPI_Isend(imatrix&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (imatrix&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for cimatrix:

int MPI_Send (cimatrix&, int, int, MPI_Comm);
int MPI_Isend(cimatrix&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (cimatrix&, int, int, MPI_Comm, MPI_Status* );

// **************************************************************

// MPI-Send/Recv for l_rvector:

int MPI_Send (l_rvector&, int, int, MPI_Comm);
int MPI_Isend(l_rvector&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_rvector&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for l_ivector:

int MPI_Send (l_ivector&, int, int, MPI_Comm);
int MPI_Isend(l_ivector&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_ivector&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for l_rmatrix:

int MPI_Send (l_rmatrix&, int, int, MPI_Comm);
int MPI_Isend(l_rmatrix&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_rmatrix&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for l_imatrix:

int MPI_Send (l_imatrix&, int, int, MPI_Comm);
int MPI_Isend(l_imatrix&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (l_imatrix&, int, int, MPI_Comm, MPI_Status* );

// **************************************************************

// MPI-Send/Recv for dotprecision:

int MPI_Send (dotprecision&, int, int, MPI_Comm);
int MPI_Isend(dotprecision&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (dotprecision&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for idotprecision:

int MPI_Send (idotprecision&, int, int, MPI_Comm);
int MPI_Isend(idotprecision&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (idotprecision&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for cdotprecision:

int MPI_Send (cdotprecision&, int, int, MPI_Comm);
int MPI_Isend(cdotprecision&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (cdotprecision&, int, int, MPI_Comm, MPI_Status* );

// MPI-Send/Recv for cidotprecision:

int MPI_Send (cidotprecision&, int, int, MPI_Comm);
int MPI_Isend(cidotprecision&, int, int, MPI_Comm, MPI_Request*);
int MPI_Recv (cidotprecision&, int, int, MPI_Comm, MPI_Status* );

#endif
