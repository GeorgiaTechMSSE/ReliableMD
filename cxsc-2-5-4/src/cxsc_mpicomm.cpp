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
// Author: Markus Grimmer, Michael Zimmer
// Date: 15.12.2009
//

#include "cxsc_mpicomm.hpp" 
#include <rvector.hpp>
#include <ivector.hpp>
#include <srvector.hpp>
#include <sivector.hpp>
#include <rmatrix.hpp>
#include <imatrix.hpp>
#include <srmatrix.hpp>
#include <simatrix.hpp>

#include <cvector.hpp>
#include <civector.hpp>
#include <scvector.hpp>
#include <scivector.hpp>
#include <cmatrix.hpp>
#include <cimatrix.hpp>
#include <scmatrix.hpp>
#include <scimatrix.hpp>

#include <l_rvector.hpp>
#include <l_ivector.hpp>
#include <l_rmatrix.hpp>
#include <l_imatrix.hpp>

#include <dot.hpp>
#include <idot.hpp>
#include <cdot.hpp>
#include <cidot.hpp>

using namespace std;

namespace cxsc {

int MPI_CXSC_BUFFERLEN = 100000000; //default buffer length
char* commbuffer = 0;

//Initialize MPI communication buffer with default length
void init_CXSC_MPI() {
  commbuffer = new char[MPI_CXSC_BUFFERLEN];
}

//Initialize MPI communication buffer with length n
void init_CXSC_MPI(int n) {
  MPI_CXSC_BUFFERLEN = n;
  commbuffer = new char[MPI_CXSC_BUFFERLEN];
}

//Delete MPI communication buffer
void finish_CXSC_MPI() {
  delete[] commbuffer;
}


bool MPI_CXSC_TYPES_DEFINED=false;
MPI_Datatype MPI_CXSC_REAL;
MPI_Datatype MPI_CXSC_COMPLEX;
MPI_Datatype MPI_CXSC_INTERVAL;
MPI_Datatype MPI_CXSC_CINTERVAL;

// Definition of new MPI data types:

int MPI_Define_CXSC_Types()
{
  int err;
  if ((err=MPI_Type_contiguous(1,MPI_DOUBLE,&MPI_CXSC_REAL))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_commit(&MPI_CXSC_REAL))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_contiguous(2,MPI_CXSC_REAL,&MPI_CXSC_COMPLEX))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_commit(&MPI_CXSC_COMPLEX))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_contiguous(2,MPI_CXSC_REAL,&MPI_CXSC_INTERVAL))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_commit(&MPI_CXSC_INTERVAL))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_contiguous(2,MPI_CXSC_INTERVAL,&MPI_CXSC_CINTERVAL))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Type_commit(&MPI_CXSC_CINTERVAL))!=MPI_SUCCESS)
    return err;
  MPI_CXSC_TYPES_DEFINED=true;
  return err;
}

// =============================================================================
// =============================================================================
//
// Pack and Unpack functions for C-XSC data types
//
// =============================================================================
// =============================================================================

// =============================================================================
//
// Scalar data types (not necessary since these are defined as 
// MPI data types - only for completeness purposes)
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for real:
// =============================================================================


int MPI_Pack (real& r, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Pack(&r,1,MPI_CXSC_REAL,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, real& r, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Unpack(buff,bufflen,pos,&r,1,MPI_CXSC_REAL,MC);
  return err;
}

// =============================================================================
// MPI-Pack/Unpack for interval:
// =============================================================================

int MPI_Pack (interval& ri, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Pack(&ri,1,MPI_CXSC_INTERVAL,buff,bufflen,pos,MC);
  return err;

}

int MPI_Unpack (void* buff, int bufflen, int* pos, interval& ri, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Unpack(buff,bufflen,pos,&ri,1,MPI_CXSC_INTERVAL,MC);
  return err;
}

// =============================================================================
// MPI-Pack/Unpack for complex:
// =============================================================================

int MPI_Pack (complex& c, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Pack(&c,1,MPI_CXSC_COMPLEX,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, complex& c, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Unpack(buff,bufflen,pos,&c,1,MPI_CXSC_COMPLEX,MC);
  return err;
}

// =============================================================================
// MPI-Pack/Unpack for cinterval:
// =============================================================================

int MPI_Pack (cinterval& ci, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Pack(&ci,1,MPI_CXSC_CINTERVAL,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, cinterval& ci, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Unpack(buff,bufflen,pos,&ci,1,MPI_CXSC_CINTERVAL,MC);
  return err;
}

// =============================================================================
//
// Multiple Precision Scalar data types
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for l_real:
// =============================================================================


int MPI_Pack (l_real& l_r, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int prec=StagPrec(l_r);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&prec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&l_r[1],prec,MPI_CXSC_REAL,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_real& l_r, MPI_Comm MC)
{
  int err;
  int lrprec;
  int stagprectmp;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lrprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
    
  // Temporarily change stagprec to declare a new variable with the received
  // precision:
  
  stagprectmp=stagprec;
  stagprec=lrprec;
  l_real lr_temp; // Create temporary variable with received precision
  stagprec=stagprectmp;

  err=MPI_Unpack(buff, bufflen, pos, &lr_temp[1], lrprec, MPI_CXSC_REAL, MC);  
  l_r=lr_temp;
  return err;

}

// =============================================================================
// MPI-Pack/Unpack for l_interval:
// =============================================================================

int MPI_Pack (l_interval& l_i, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int prec=StagPrec(l_i);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&prec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&l_i[1],prec+1,MPI_CXSC_REAL,buff,bufflen,pos,MC);
  return err;

}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_interval& l_i, MPI_Comm MC)
{
  int err;
  int liprec;
  int stagprectmp;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &liprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  // Temporarily change stagprec to declare a new variable with the received
  // precision:
  
  stagprectmp=stagprec;
  stagprec=liprec;
  l_interval li_temp;  // Create temporary variable with received precision
  stagprec=stagprectmp;

  err=MPI_Unpack(buff, bufflen, pos, &li_temp[1], liprec+1, MPI_CXSC_REAL, MC);  
  l_i=li_temp;
  return err;

}

// =============================================================================
// MPI-Pack/Unpack for l_complex:
// =============================================================================

int MPI_Pack (l_complex& l_c, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int reprec=StagPrec(Re(l_c));
  int imprec=StagPrec(Im(l_c));
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&reprec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&(Re(l_c)[1]),reprec,MPI_CXSC_REAL,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&imprec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&(Im(l_c)[1]),imprec,MPI_CXSC_REAL,buff,bufflen,pos,MC);
  return err;

}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_complex& l_c, MPI_Comm MC)
{
  int err;
  int reprec,imprec;
  int stagprectmp;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &reprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  // Temporarily change stagprec to declare a new variable with the received
  // precision for re:
  
  stagprectmp=stagprec;
  stagprec=reprec;
  l_real re_temp; // Create temporary variable with received precision
  stagprec=stagprectmp;

  if ((err=MPI_Unpack(buff, bufflen, pos, &re_temp[1], reprec, MPI_CXSC_REAL, MC))!=MPI_SUCCESS)  
    return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &imprec, 1, MPI_INT, MC))!=MPI_SUCCESS)  
    return err;

  // Temporarily change stagprec to declare a new variable with the received
  // precision for im:
  
  stagprectmp=stagprec;
  stagprec=imprec;
  l_real im_temp; // Create temporary variable with received precision
  stagprec=stagprectmp;

  err=MPI_Unpack(buff, bufflen, pos, &im_temp[1], imprec, MPI_CXSC_REAL, MC);  

  Re(l_c)=re_temp; Im(l_c)=im_temp;
  return err;
}


// =============================================================================
//
// Vector data types
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for rvector:
// =============================================================================


int MPI_Pack (rvector& rv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int lb=Lb(rv);
  int ub=Ub(rv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types() != MPI_SUCCESS))
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&rv[lb],ub-lb+1,MPI_CXSC_REAL,buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (rvector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  err=MPI_Pack(&rm[lb1], ub1-lb1+1,MPI_CXSC_REAL,buff,
             bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, rvector& rv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  Resize(rv,vlen);
  SetLb(rv,lb);                         
  //rv.SetUb(ub);  // Should not be necessary, you can try it for test purposes
  
  err=MPI_Unpack(buff, bufflen, pos, &rv[lb] , vlen, MPI_CXSC_REAL, MC);  
  return err;
}


int MPI_Unpack (void* buff, int bufflen, int* pos, 
                rvector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     


  err=MPI_Unpack(buff, bufflen, pos, &rm[lb1], ub1-lb1+1, MPI_CXSC_REAL, MC);  
  return err;

}

// =============================================================================
// MPI-Pack/Unpack for ivector:
// =============================================================================

int MPI_Pack (ivector& iv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb=Lb(iv);
  int ub=Ub(iv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&iv[lb],ub-lb+1,MPI_CXSC_INTERVAL,buff,bufflen,pos,MC);
  return err;

}

int MPI_Pack (ivector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  err=MPI_Pack(&rm[lb1], ub1-lb1+1,MPI_CXSC_INTERVAL,buff,
             bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, ivector& iv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  Resize(iv,vlen);
  SetLb(iv,lb);                         
  //iv.SetUb(ub);  // Should not be necessary, you can try it for test purposes
  
  err=MPI_Unpack(buff, bufflen, pos, &iv[lb] , vlen, MPI_CXSC_INTERVAL, MC);  
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                ivector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     


  err=MPI_Unpack(buff, bufflen, pos, &rm[lb1], ub1-lb1+1, MPI_CXSC_INTERVAL, MC);  
  return err;

}

// =============================================================================
// MPI-Pack/Unpack for cvector:
// =============================================================================

int MPI_Pack (cvector& cv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb=Lb(cv);
  int ub=Ub(cv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&cv[lb],ub-lb+1,MPI_CXSC_COMPLEX,buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (cvector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  err=MPI_Pack(&rm[lb1], ub1-lb1+1,MPI_CXSC_COMPLEX,buff,
             bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, cvector& cv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  Resize(cv,vlen);
  SetLb(cv,lb);                         
  //cv.SetUb(ub);  // Should not be necessary, you can try it for test purposes
  
  err=MPI_Unpack(buff, bufflen, pos, &cv[lb] , vlen, MPI_CXSC_COMPLEX, MC);  
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                cvector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     


  err=MPI_Unpack(buff, bufflen, pos, &rm[lb1], ub1-lb1+1, MPI_CXSC_COMPLEX, MC);  
  return err;

}
// =============================================================================
// MPI-Pack/Unpack for civector:
// =============================================================================

int MPI_Pack (civector& civ, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb=Lb(civ);
  int ub=Ub(civ);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&civ[lb],ub-lb+1,MPI_CXSC_CINTERVAL,buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (civector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  err=MPI_Pack(&rm[lb1], ub1-lb1+1,MPI_CXSC_CINTERVAL,buff,
             bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, civector& civ, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  Resize(civ,vlen);
  SetLb(civ,lb);                         
  //civ.SetUb(ub);  // Should not be necessary, you can try it for test purposes
  
  err=MPI_Unpack(buff, bufflen, pos, &civ[lb] , vlen, MPI_CXSC_CINTERVAL, MC);  
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                civector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     


  err=MPI_Unpack(buff, bufflen, pos, &rm[lb1], ub1-lb1+1, MPI_CXSC_CINTERVAL, MC);  
  return err;

}


// =============================================================================
// MPI-Pack/Unpack for srvector:
// =============================================================================


int MPI_Pack (srvector& rv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int lb=Lb(rv);
  int ub=Ub(rv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rv.row_indices(),buff,bufflen,pos,MC)) != MPI_SUCCESS)
    return err;
  err=MPI_Pack(rv.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (srvector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  srvector x = rm(lb1,ub1);
  err=MPI_Pack(x,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, srvector& rv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  SetLb(rv,lb);                         
  
  if((err=MPI_Unpack(buff, bufflen, pos, rv.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, rv.values(), MC);  

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                srvector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     

  srvector x(ub1-lb1+1);
  if((err=MPI_Unpack(buff, bufflen, pos, x, MC)) != MPI_SUCCESS) return err;
  rm(l1,u1) = x;

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for sivector:
// =============================================================================

int MPI_Pack (sivector& iv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int lb=Lb(iv);
  int ub=Ub(iv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(iv.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(iv.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (sivector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  sivector x = rm(lb1,ub1);
  err=MPI_Pack(x,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, sivector& iv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  SetLb(iv,lb);                         
  
  if( (err=MPI_Unpack(buff, bufflen, pos, iv.row_indices(), MC)) != MPI_SUCCESS) 
    return err;
  err=MPI_Unpack(buff, bufflen, pos, iv.values(), MC);  

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                sivector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     

  sivector x(u1-l1+1);
  if((err=MPI_Unpack(buff, bufflen, pos, x, MC)) != MPI_SUCCESS) return err;
  rm(l1,u1) = x;

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for scvector:
// =============================================================================

int MPI_Pack (scvector& cv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int lb=Lb(cv);
  int ub=Ub(cv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(cv.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(cv.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (scvector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  scvector x = rm(lb1,ub1);
  err=MPI_Pack(x,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, scvector& cv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  SetLb(cv,lb);                         
  
  if((err=MPI_Unpack(buff, bufflen, pos, cv.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, cv.values(), MC);  

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                scvector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     

  scvector x(ub1-lb1+1);
  if((err=MPI_Unpack(buff, bufflen, pos, x, MC)) != MPI_SUCCESS) return err;
  rm(l1,u1) = x;

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for scivector:
// =============================================================================

int MPI_Pack (scivector& civ, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int lb=Lb(civ);
  int ub=Ub(civ);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(civ.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(civ.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (scivector& rm, int lb1, int ub1, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm))||(lb1>Ub(rm))
     ||(ub1<Lb(rm))||(ub1>Ub(rm)))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  scivector x = rm(lb1,ub1);
  err=MPI_Pack(x,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, scivector& civ, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  SetLb(civ,lb);                         
  
  if((err=MPI_Unpack(buff, bufflen, pos, civ.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, civ.values(), MC);  

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                scivector& rm, int l1, int u1, MPI_Comm MC)
{
  int err;
  int lb1, ub1;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm))||(l1>Ub(rm))
     ||(u1<Lb(rm))||(u1>Ub(rm)))
    return MPI_ERR_DIMS;     

  scivector x(ub1-lb1+1);
  if((err=MPI_Unpack(buff, bufflen, pos, x, MC)) != MPI_SUCCESS) return err;
  rm(l1,u1) = x;

  return err;
}

// =============================================================================
//
// Multiple Precision Vector data types
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for l_rvector:
// =============================================================================


int MPI_Pack (l_rvector& l_rv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  int lb=Lb(l_rv);
  int ub=Ub(l_rv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  for (int i=lb;i<=ub;i++)
    {
      int prec=StagPrec(l_rv[i]);

      if ((err=MPI_Pack(&prec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
        return err;
      if ((err=MPI_Pack(&l_rv[i][1],prec,MPI_CXSC_REAL,buff,bufflen,pos,MC))!=MPI_SUCCESS)
        return err;
    }

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_rvector& l_rv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  int stagprectmp;

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  Resize(l_rv,vlen);
  SetLb(l_rv,lb);                         
  //l_rv.SetUb(ub);  // Should not be necessary, you can try it for test purposes
  
  for (int i=lb; i<=ub; i++)
    {
      int lrprec;
      if ((err=MPI_Unpack(buff, bufflen, pos, &lrprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
        return err;
      stagprectmp=stagprec;
      stagprec=lrprec;
      l_real lr_temp; // Create temporary variable with received precision
      stagprec=stagprectmp;
    
      if ((err=MPI_Unpack(buff, bufflen, pos, &lr_temp[1], lrprec, MPI_CXSC_REAL, MC))!=MPI_SUCCESS)
        return err;
      l_rv[i]=lr_temp;

    }

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for l_ivector:
// =============================================================================

int MPI_Pack (l_ivector& l_iv, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb=Lb(l_iv);
  int ub=Ub(l_iv);
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  if ((err=MPI_Pack(&lb,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  for (int i=lb;i<=ub;i++)
    {
      int prec=StagPrec(l_iv[i]);

      if ((err=MPI_Pack(&prec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
        return err;
      if ((err=MPI_Pack(&l_iv[i][1],prec+1,MPI_CXSC_REAL,buff,bufflen,pos,MC))!=MPI_SUCCESS)
        return err;
    }

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_ivector& l_iv, MPI_Comm MC)
{
  int err;
  int vlen, lb, ub;
  int stagprectmp;

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  vlen=ub-lb+1;
  Resize(l_iv,vlen);
  SetLb(l_iv,lb);                         
  //l_iv.SetUb(ub);  // Should not be necessary, you can try it for test purposes
  
  for (int i=lb; i<=ub; i++)
    {
      int liprec;
      if ((err=MPI_Unpack(buff, bufflen, pos, &liprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
        return err;
    
      stagprectmp=stagprec;
      stagprec=liprec;
      l_interval li_temp; // Create temporary variable with received precision
      stagprec=stagprectmp;
    
      if ((err=MPI_Unpack(buff, bufflen, pos, &li_temp[1], liprec+1, MPI_CXSC_REAL, MC))!=MPI_SUCCESS)
        return err;
      l_iv[i]=li_temp;
    }
  return err;
}

// =============================================================================
//
// Matrix data types
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for rmatrix:
// =============================================================================


int MPI_Pack (rmatrix& rm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(rm,1);
  int lb2=Lb(rm,2);
  int ub1=Ub(rm,1);
  int ub2=Ub(rm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&rm[lb1][lb2],(ub1-lb1+1)*(ub2-lb2+1),MPI_CXSC_REAL,buff,
           bufflen,pos,MC);

  return err;
}

// -------------
// Submatrix:
// -------------

int MPI_Pack (rmatrix& rm, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm,1))||(lb1>Ub(rm,1))
     ||(ub1<Lb(rm,1))||(ub1>Ub(rm,1))
     ||(lb2<Lb(rm,2))||(lb2>Ub(rm,2))
     ||(ub2<Lb(rm,2))||(ub2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  for (int i1=lb1;i1<=ub1; i1++)
    err=MPI_Pack(&rm[i1][lb2],(ub2-lb2+1),MPI_CXSC_REAL,buff,
             bufflen,pos,MC);
  return err;
}

// ------------
// Full matrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, rmatrix& rm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(rm,mlen1,mlen2);
  SetLb(rm,1,lb1);
  SetLb(rm,2,lb2);
  //SetUb(rm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(rm,2,ub2); // Should not be necessary - try it for test purposes
   
  err=MPI_Unpack(buff, bufflen, pos, &rm[lb1][lb2], 
             mlen1*mlen2, MPI_CXSC_REAL, MC);  
  return err;

}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                rmatrix& rm, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm,1))||(l1>Ub(rm,1))
     ||(u1<Lb(rm,1))||(u1>Ub(rm,1))
     ||(l2<Lb(rm,2))||(l2>Ub(rm,2))
     ||(u2<Lb(rm,2))||(u2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen2=u2-l2+1;
  
  // In constrast to full matrices, the matrix is NOT resized here!
     // Resize(rm,mlen1,mlen2);
     // SetLb(rm,1,lb1);
     // SetLb(rm,2,lb2);

  for (int i1=l1;i1<=u1; i1++)
    err=MPI_Unpack(buff, bufflen, pos, &rm[i1][l2], 
                   mlen2, MPI_CXSC_REAL, MC);  
  return err;

}


// =============================================================================
// MPI-Pack/Unpack for imatrix:
// =============================================================================

// -------------
// Full matrix:
// -------------

int MPI_Pack (imatrix& im, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(im,1);
  int lb2=Lb(im,2);
  int ub1=Ub(im,1);
  int ub2=Ub(im,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&im[lb1][lb2],(ub1-lb1+1)*(ub2-lb2+1),MPI_CXSC_INTERVAL,buff,
           bufflen,pos,MC);
  return err;
}

// -------------
// Submatrix:
// -------------

int MPI_Pack (imatrix& im, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(im,1))||(lb1>Ub(im,1))
     ||(ub1<Lb(im,1))||(ub1>Ub(im,1))
     ||(lb2<Lb(im,2))||(lb2>Ub(im,2))
     ||(ub2<Lb(im,2))||(ub2>Ub(im,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  for (int i1=lb1;i1<=ub1; i1++)
    err=MPI_Pack(&im[i1][lb2],(ub2-lb2+1),MPI_CXSC_INTERVAL,buff,
             bufflen,pos,MC);
  return err;
}

// -------------
// Full matrix:
// -------------

int MPI_Unpack (void* buff, int bufflen, int* pos, imatrix& im, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(im,mlen1,mlen2);
  SetLb(im,1,lb1);
  SetLb(im,2,lb2);
  //SetUb(im,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(im,2,ub2); // Should not be necessary - try it for test purposes
   
  err=MPI_Unpack(buff, bufflen, pos, &im[lb1][lb2], 
             mlen1*mlen2, MPI_CXSC_INTERVAL, MC);  
  return err;
}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                imatrix& im, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(im,1))||(l1>Ub(im,1))
     ||(u1<Lb(im,1))||(u1>Ub(im,1))
     ||(l2<Lb(im,2))||(l2>Ub(im,2))
     ||(u2<Lb(im,2))||(u2>Ub(im,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen2=u2-l2+1;
  
  // In constrast to full matrices, the matrix is NOT resized here!
     // Resize(im,mlen1,mlen2);
     // SetLb(im,1,lb1);
     // SetLb(im,2,lb2);

  for (int i1=l1;i1<=u1; i1++)
    err=MPI_Unpack(buff, bufflen, pos, &im[i1][l2], 
                   mlen2, MPI_CXSC_INTERVAL, MC);  
  return err;

}


// =============================================================================
// MPI-Pack/Unpack for cmatrix:
// =============================================================================

// -------------
// Full matrix:
// -------------

int MPI_Pack (cmatrix& cm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(cm,1);
  int lb2=Lb(cm,2);
  int ub1=Ub(cm,1);
  int ub2=Ub(cm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&cm[lb1][lb2],(ub1-lb1+1)*(ub2-lb2+1),MPI_CXSC_COMPLEX,buff,
           bufflen,pos,MC);
  return err;
}

// -------------
// Submatrix:
// -------------

int MPI_Pack (cmatrix& cm, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(cm,1))||(lb1>Ub(cm,1))
     ||(ub1<Lb(cm,1))||(ub1>Ub(cm,1))
     ||(lb2<Lb(cm,2))||(lb2>Ub(cm,2))
     ||(ub2<Lb(cm,2))||(ub2>Ub(cm,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  for (int i1=lb1;i1<=ub1; i1++)
    err=MPI_Pack(&cm[i1][lb2],(ub2-lb2+1),MPI_CXSC_COMPLEX,buff,
             bufflen,pos,MC);
  return err;
}

// -------------
// Full matrix:
// -------------

int MPI_Unpack (void* buff, int bufflen, int* pos, cmatrix& cm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(cm,mlen1,mlen2);
  SetLb(cm,1,lb1);
  SetLb(cm,2,lb2);
  //SetUb(cm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(cm,2,ub2); // Should not be necessary - try it for test purposes
   
  err=MPI_Unpack(buff, bufflen, pos, &cm[lb1][lb2], 
             mlen1*mlen2, MPI_CXSC_COMPLEX, MC);  
  return err;
}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                cmatrix& cm, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(cm,1))||(l1>Ub(cm,1))
     ||(u1<Lb(cm,1))||(u1>Ub(cm,1))
     ||(l2<Lb(cm,2))||(l2>Ub(cm,2))
     ||(u2<Lb(cm,2))||(u2>Ub(cm,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen2=u2-l2+1;
  
  // In constrast to full matrices, the matrix is NOT resized here!
     // Resize(cm,mlen1,mlen2);
     // SetLb(cm,1,lb1);
     // SetLb(cm,2,lb2);

  for (int i1=l1;i1<=u1; i1++)
    err=MPI_Unpack(buff, bufflen, pos, &cm[i1][l2], 
                   mlen2, MPI_CXSC_COMPLEX, MC);  
  return err;

}

// =============================================================================
// MPI-Pack/Unpack for cimatrix:
// =============================================================================

// -------------
// Full matrix:
// -------------

int MPI_Pack (cimatrix& cim, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(cim,1);
  int lb2=Lb(cim,2);
  int ub1=Ub(cim,1);
  int ub2=Ub(cim,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(&cim[lb1][lb2],(ub1-lb1+1)*(ub2-lb2+1),MPI_CXSC_CINTERVAL,buff,
           bufflen,pos,MC);
  return err;
}

// -------------
// Submatrix:
// -------------

int MPI_Pack (cimatrix& cim, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(cim,1))||(lb1>Ub(cim,1))
     ||(ub1<Lb(cim,1))||(ub1>Ub(cim,1))
     ||(lb2<Lb(cim,2))||(lb2>Ub(cim,2))
     ||(ub2<Lb(cim,2))||(ub2>Ub(cim,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  for (int i1=lb1;i1<=ub1; i1++)
    err=MPI_Pack(&cim[i1][lb2],(ub2-lb2+1),MPI_CXSC_CINTERVAL,buff,
             bufflen,pos,MC);
  return err;
}

// -------------
// Full matrix:
// -------------

int MPI_Unpack (void* buff, int bufflen, int* pos, cimatrix& cim, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(cim,mlen1,mlen2);
  SetLb(cim,1,lb1);
  SetLb(cim,2,lb2);
  //SetUb(cim,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(cim,2,ub2); // Should not be necessary - try it for test purposes
   
  err=MPI_Unpack(buff, bufflen, pos, &cim[lb1][lb2], 
             mlen1*mlen2, MPI_CXSC_CINTERVAL, MC);  
  return err;
}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                cimatrix& cim, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(cim,1))||(l1>Ub(cim,1))
     ||(u1<Lb(cim,1))||(u1>Ub(cim,1))
     ||(l2<Lb(cim,2))||(l2>Ub(cim,2))
     ||(u2<Lb(cim,2))||(u2>Ub(cim,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen2=u2-l2+1;
  
  // In constrast to full matrices, the matrix is NOT resized here!
     // Resize(cim,mlen1,mlen2);
     // SetLb(cim,1,lb1);
     // SetLb(cim,2,lb2);

  for (int i1=l1;i1<=u1; i1++)
    err=MPI_Unpack(buff, bufflen, pos, &cim[i1][l2], 
                   mlen2, MPI_CXSC_CINTERVAL, MC);  
  return err;

}


// =============================================================================
// MPI-Pack/Unpack for srmatrix:
// =============================================================================


int MPI_Pack (srmatrix& rm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(rm,1);
  int lb2=Lb(rm,2);
  int ub1=Ub(rm,1);
  int ub2=Ub(rm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.column_pointers(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(rm.values(),buff,bufflen,pos,MC);
  return err;
}

// -------------
// Submatrix:
// -------------

int MPI_Pack (srmatrix& rm, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm,1))||(lb1>Ub(rm,1))
     ||(ub1<Lb(rm,1))||(ub1>Ub(rm,1))
     ||(lb2<Lb(rm,2))||(lb2>Ub(rm,2))
     ||(ub2<Lb(rm,2))||(ub2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  srmatrix A = rm(lb1,ub1,lb2,ub2);
  err=MPI_Pack(A,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, srmatrix& rm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(rm,mlen1,mlen2);
  SetLb(rm,1,lb1);
  SetLb(rm,2,lb2);
  //SetUb(rm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(rm,2,ub2); // Should not be necessary - try it for test purposes
   
  if((err=MPI_Unpack(buff, bufflen, pos, rm.column_pointers(), MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Unpack(buff, bufflen, pos, rm.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, rm.values(), MC);
  
  return err;

}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                srmatrix& rm, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen1,mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm,1))||(l1>Ub(rm,1))
     ||(u1<Lb(rm,1))||(u1>Ub(rm,1))
     ||(l2<Lb(rm,2))||(l2>Ub(rm,2))
     ||(u2<Lb(rm,2))||(u2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen1=u1-l1+1;
  mlen2=u2-l2+1;
  
  srmatrix A(mlen1,mlen2);

  if((err=MPI_Unpack(buff, bufflen, pos, A, MC))!=MPI_SUCCESS)
    return err;

  rm(l1,u1,l2,u2) = A;  

  return err;

}

// =============================================================================
// MPI-Pack/Unpack for simatrix:
// =============================================================================

int MPI_Pack (simatrix& rm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(rm,1);
  int lb2=Lb(rm,2);
  int ub1=Ub(rm,1);
  int ub2=Ub(rm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.column_pointers(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(rm.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (simatrix& rm, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm,1))||(lb1>Ub(rm,1))
     ||(ub1<Lb(rm,1))||(ub1>Ub(rm,1))
     ||(lb2<Lb(rm,2))||(lb2>Ub(rm,2))
     ||(ub2<Lb(rm,2))||(ub2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  simatrix A = rm(lb1,ub1,lb2,ub2);
  err=MPI_Pack(A,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, simatrix& rm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(rm,mlen1,mlen2);
  SetLb(rm,1,lb1);
  SetLb(rm,2,lb2);
  //SetUb(rm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(rm,2,ub2); // Should not be necessary - try it for test purposes
   
  if((err=MPI_Unpack(buff, bufflen, pos, rm.column_pointers(), MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Unpack(buff, bufflen, pos, rm.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, rm.values(), MC);
  
  return err;
}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                simatrix& rm, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen1,mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm,1))||(l1>Ub(rm,1))
     ||(u1<Lb(rm,1))||(u1>Ub(rm,1))
     ||(l2<Lb(rm,2))||(l2>Ub(rm,2))
     ||(u2<Lb(rm,2))||(u2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen1=u1-l1+1;
  mlen2=u2-l2+1;
  
  simatrix A(mlen1,mlen2);

  if((err=MPI_Unpack(buff, bufflen, pos, A, MC))!=MPI_SUCCESS)
    return err;

  rm(l1,u1,l2,u2) = A;  

  return err;

}


// =============================================================================
// MPI-Pack/Unpack for scmatrix:
// =============================================================================

int MPI_Pack (scmatrix& rm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(rm,1);
  int lb2=Lb(rm,2);
  int ub1=Ub(rm,1);
  int ub2=Ub(rm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.column_pointers(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(rm.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (scmatrix& rm, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm,1))||(lb1>Ub(rm,1))
     ||(ub1<Lb(rm,1))||(ub1>Ub(rm,1))
     ||(lb2<Lb(rm,2))||(lb2>Ub(rm,2))
     ||(ub2<Lb(rm,2))||(ub2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  scmatrix A = rm(lb1,ub1,lb2,ub2);
  err=MPI_Pack(A,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, scmatrix& rm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(rm,mlen1,mlen2);
  SetLb(rm,1,lb1);
  SetLb(rm,2,lb2);
  //SetUb(rm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(rm,2,ub2); // Should not be necessary - try it for test purposes
   
  if((err=MPI_Unpack(buff, bufflen, pos, rm.column_pointers(), MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Unpack(buff, bufflen, pos, rm.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, rm.values(), MC);
  
  return err;
}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                scmatrix& rm, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen1,mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm,1))||(l1>Ub(rm,1))
     ||(u1<Lb(rm,1))||(u1>Ub(rm,1))
     ||(l2<Lb(rm,2))||(l2>Ub(rm,2))
     ||(u2<Lb(rm,2))||(u2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen1=u1-l1+1;
  mlen2=u2-l2+1;
  
  scmatrix A(mlen1,mlen2);

  if((err=MPI_Unpack(buff, bufflen, pos, A, MC))!=MPI_SUCCESS)
    return err;

  rm(l1,u1,l2,u2) = A;  

  return err;

}


// =============================================================================
// MPI-Pack/Unpack for scimatrix:
// =============================================================================

int MPI_Pack (scimatrix& rm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(rm,1);
  int lb2=Lb(rm,2);
  int ub1=Ub(rm,1);
  int ub2=Ub(rm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.column_pointers(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Pack(rm.row_indices(),buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Pack(rm.values(),buff,bufflen,pos,MC);
  return err;
}

int MPI_Pack (scimatrix& rm, int lb1, int ub1, int lb2, int ub2, 
              void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;

  if ((lb1<Lb(rm,1))||(lb1>Ub(rm,1))
     ||(ub1<Lb(rm,1))||(ub1>Ub(rm,1))
     ||(lb2<Lb(rm,2))||(lb2>Ub(rm,2))
     ||(ub2<Lb(rm,2))||(ub2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2))
    return MPI_ERR_DIMS;     
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  scimatrix A = rm(lb1,ub1,lb2,ub2);
  err=MPI_Pack(A,buff,bufflen,pos,MC);
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, scimatrix& rm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(rm,mlen1,mlen2);
  SetLb(rm,1,lb1);
  SetLb(rm,2,lb2);
  //SetUb(rm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(rm,2,ub2); // Should not be necessary - try it for test purposes
   
  if((err=MPI_Unpack(buff, bufflen, pos, rm.column_pointers(), MC))!=MPI_SUCCESS)
    return err;
  if((err=MPI_Unpack(buff, bufflen, pos, rm.row_indices(), MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Unpack(buff, bufflen, pos, rm.values(), MC);
  
  return err;
}

// ------------
// Submatrix:
// ------------

int MPI_Unpack (void* buff, int bufflen, int* pos, 
                scimatrix& rm, int l1, int u1, int l2, int u2, MPI_Comm MC)
{
  int err;
  int mlen1,mlen2, lb1, ub1, lb2, ub2;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  if ((l1<Lb(rm,1))||(l1>Ub(rm,1))
     ||(u1<Lb(rm,1))||(u1>Ub(rm,1))
     ||(l2<Lb(rm,2))||(l2>Ub(rm,2))
     ||(u2<Lb(rm,2))||(u2>Ub(rm,2))
     ||(lb1>ub1)||(lb2>ub2)
     ||((ub1-lb1)!=(u1-l1))
     ||((ub2-lb2)!=(u2-l2)))
    return MPI_ERR_DIMS;     

  mlen1=u1-l1+1;
  mlen2=u2-l2+1;
  
  scimatrix A(mlen1,mlen2);

  if((err=MPI_Unpack(buff, bufflen, pos, A, MC))!=MPI_SUCCESS)
    return err;

  rm(l1,u1,l2,u2) = A;  

  return err;

}


// =============================================================================
//
// Multiple Precision Matrix data types
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for l_rmatrix:
// =============================================================================


int MPI_Pack (l_rmatrix& l_rm, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(l_rm,1);
  int lb2=Lb(l_rm,2);
  int ub1=Ub(l_rm,1);
  int ub2=Ub(l_rm,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  for (int i=lb1;i<=ub1;i++)
    for (int j=lb2;j<=ub2;j++)
      {
        int prec=StagPrec(l_rm[i][j]);

        if ((err=MPI_Pack(&prec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
          return err;
        if ((err=MPI_Pack(&l_rm[i][j][1],prec,MPI_CXSC_REAL,buff,
                 bufflen,pos,MC))!=MPI_SUCCESS)
          return err;
      }

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_rmatrix& l_rm, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  int stagprectmp;

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(l_rm,mlen1,mlen2);
  SetLb(l_rm,1,lb1);
  SetLb(l_rm,2,lb2);
  //SetUb(l_rm,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(l_rm,2,ub2); // Should not be necessary - try it for test purposes
  
  for (int i=lb1;i<=ub1;i++)
    for (int j=lb2;j<=ub2;j++)
      {
        int lrprec;
        if ((err=MPI_Unpack(buff, bufflen, pos, &lrprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
          return err;
   
        stagprectmp=stagprec;
        stagprec=lrprec;
        l_real lr_temp; // Create temporary variable with received precision
        stagprec=stagprectmp;
    
        if ((err=MPI_Unpack(buff, bufflen, pos, &lr_temp[1], lrprec, MPI_CXSC_REAL, MC))!=MPI_SUCCESS)
          return err;
        l_rm[i][j]=lr_temp;
    }
  
  return err;
}

// =============================================================================
// MPI-Pack/Unpack for l_imatrix:
// =============================================================================

int MPI_Pack (l_imatrix& l_im, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  
  int lb1=Lb(l_im,1);
  int lb2=Lb(l_im,2);
  int ub1=Ub(l_im,1);
  int ub2=Ub(l_im,2);

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Pack(&lb1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub1,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&lb2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Pack(&ub2,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
    return err;

  for (int i=lb1;i<=ub1;i++)
    for (int j=lb2;j<=ub2;j++)
      {
        int prec=StagPrec(l_im[i][j]);

        if ((err=MPI_Pack(&prec,1,MPI_INT,buff,bufflen,pos,MC))!=MPI_SUCCESS)
          return err;
        if ((err=MPI_Pack(&l_im[i][j][1],prec+1,MPI_CXSC_REAL,buff,
                 bufflen,pos,MC))!=MPI_SUCCESS)
          return err;
      }
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, l_imatrix& l_im, MPI_Comm MC)
{
  int err;
  int mlen1, mlen2, lb1, ub1, lb2, ub2;
  int stagprectmp;

  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  if ((err=MPI_Unpack(buff, bufflen, pos, &lb1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub1, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &lb2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;
  if ((err=MPI_Unpack(buff, bufflen, pos, &ub2, 1, MPI_INT, MC))!=MPI_SUCCESS)
    return err;

  mlen1=ub1-lb1+1;
  mlen2=ub2-lb2+1;
  Resize(l_im,mlen1,mlen2);
  SetLb(l_im,1,lb1);
  SetLb(l_im,2,lb2);
  //SetUb(l_im,1,ub1); // Should not be necessary - try it for test purposes
  //SetUb(l_im,2,ub2); // Should not be necessary - try it for test purposes
  
  for (int i=lb1;i<=ub1;i++)
    for (int j=lb2;j<=ub2;j++)
      {
        int liprec;
        if ((err=MPI_Unpack(buff, bufflen, pos, &liprec, 1, MPI_INT, MC))!=MPI_SUCCESS)
          return err;
    
        stagprectmp=stagprec;
        stagprec=liprec;
        l_interval li_temp; // Create temporary variable with received precision
        stagprec=stagprectmp;
    
        if ((err=MPI_Unpack(buff, bufflen, pos, &li_temp[1], liprec+1, MPI_CXSC_REAL, MC))!=MPI_SUCCESS)
          return err;
        l_im[i][j]=li_temp;
    }
  
  return err;
}

// =============================================================================
//
// Dotprecision data types
//
// =============================================================================

// =============================================================================
// MPI-Pack/Unpack for rdotprecision:
// =============================================================================


int MPI_Pack (dotprecision& dotac, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
        
  a_btyp* akkuP=*(dotac.ptr()); // Points to the same address as akku

  // BUFFERSIZE is an internal C-XSC variable controlling the 
  // length of a dotprecision akku  
  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack K
  int k = dotac.get_k();
  if((err=MPI_Pack(&k,1,MPI_INT,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  real error = dotac.get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;
  
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, dotprecision& dotac, MPI_Comm MC)
{
  int err;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  a_btyp* akkuP=*(dotac.ptr()); // Points to the same address as akku

  // BUFFERSIZE is an internal C-XSC variable controlling the 
  // length of a dotprecision akku
  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  int k;
  if((err=MPI_Unpack(buff, bufflen, pos, &k, 1, MPI_INT, MC)) != MPI_SUCCESS)
	return err;  
  dotac.set_k(k);

  real e;
  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  dotac.set_err(e);

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for idotprecision:
// =============================================================================

int MPI_Pack (idotprecision& idotac, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
        
  a_btyp* akkuP=*(Inf(idotac).ptr()); // Points to the same address as akku

  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  akkuP=*(Sup(idotac).ptr()); // Points to the same address as akku
  
  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  int k = idotac.get_k();
  err=MPI_Pack(&k,1,MPI_INT,buff,bufflen,pos,MC);

  //Pack error value
  real error = Sup(idotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  error = Inf(idotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, idotprecision& idotac, MPI_Comm MC)
{
  int err;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  a_btyp* akkuP=*(Inf(idotac).ptr()); // Points to the same address as akku

  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  akkuP=*(Sup(idotac).ptr()); // Points to the same address as akku

  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  int k;
  if((err=MPI_Unpack(buff, bufflen, pos, &k, 1, MPI_INT, MC)) != MPI_SUCCESS)
	return err;  
  idotac.set_k(k);

  real e;
  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  Inf(idotac).set_err(e);

  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  Sup(idotac).set_err(e);

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for cdotprecision:
// =============================================================================

int MPI_Pack (cdotprecision& cdotac, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
        
  a_btyp* akkuP=*(Re(cdotac).ptr()); // Points to the same address as akku
  
  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  akkuP=*(Im(cdotac).ptr()); // Points to the same address as akku

  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  int k = cdotac.get_k();
  if((err=MPI_Pack(&k,1,MPI_INT,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  real error = Re(cdotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  error = Im(cdotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;
  
  return err;
}

int MPI_Unpack (void* buff, int bufflen, int* pos, cdotprecision& cdotac, MPI_Comm MC)
{
  int err;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  a_btyp* akkuP=*(Re(cdotac).ptr()); // Points to the same address as akku

  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  
                 // BUFFERSIZE is an internal C-XSC variable controlling the 
                 // length of a dotprecision akku

  akkuP=*(Im(cdotac).ptr()); // Points to the same address as akku

  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  int k;
  if((err=MPI_Unpack(buff, bufflen, pos, &k, 1, MPI_INT, MC)) != MPI_SUCCESS)
	return err;  
  cdotac.set_k(k);

  real e;
  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  Re(cdotac).set_err(e);

  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  Im(cdotac).set_err(e); 

  return err;
}

// =============================================================================
// MPI-Pack/Unpack for cidotprecision:
// =============================================================================

int MPI_Pack (cidotprecision& cidotac, void* buff, int bufflen, int* pos,  MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
        
  a_btyp* akkuP=*(InfRe(cidotac).ptr()); // Points to the same address as akku
  
  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;
                 // BUFFERSIZE is an internal C-XSC variable controlling the 
                 // length of a dotprecision akku

  akkuP=*(InfIm(cidotac).ptr()); // Points to the same address as akku
  
  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  akkuP=*(SupRe(cidotac).ptr()); // Points to the same address as akku
  
  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  akkuP=*(SupIm(cidotac).ptr()); // Points to the same address as akku

  if((err=MPI_Pack((void*)akkuP,BUFFERSIZE,MPI_BYTE,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  int k = cidotac.get_k();
  if((err=MPI_Pack(&k,1,MPI_INT,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  real error = InfRe(cidotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  error = InfIm(cidotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  error = SupRe(cidotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;

  //Pack error value
  error = SupIm(cidotac).get_err();
  if((err=MPI_Pack(&error,1,MPI_CXSC_REAL,buff,bufflen,pos,MC)) != MPI_SUCCESS)
	return err;
  
  return err;

}

int MPI_Unpack (void* buff, int bufflen, int* pos, cidotprecision& cidotac, MPI_Comm MC)
{
  int err;
  
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  a_btyp* akkuP=*(InfRe(cidotac).ptr()); // Points to the same address as akku
  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  
                 // BUFFERSIZE is an internal C-XSC variable controlling the 
                 // length of a dotprecision akku

  akkuP=*(InfIm(cidotac).ptr()); // Points to the same address as akku
  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  akkuP=*(SupRe(cidotac).ptr()); // Points to the same address as akku
  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  akkuP=*(SupIm(cidotac).ptr()); // Points to the same address as akku
  if((err=MPI_Unpack(buff, bufflen, pos, akkuP, BUFFERSIZE, MPI_BYTE, MC)) != MPI_SUCCESS)
	return err;  

  int k;
  if((err=MPI_Unpack(buff, bufflen, pos, &k, 1, MPI_INT, MC)) != MPI_SUCCESS)
	return err;  
  cidotac.set_k(k);

  real e;
  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  InfRe(cidotac).set_err(e);

  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  InfIm(cidotac).set_err(e); 

  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  SupRe(cidotac).set_err(e);

  if((err=MPI_Unpack(buff, bufflen, pos, &e, 1, MPI_CXSC_REAL, MC)) != MPI_SUCCESS)
	return err;  
  SupIm(cidotac).set_err(e); 

  return err;
}

// *****************************************************************************
// *****************************************************************************
//
// POINT-TO-POINT COMMUNICATION FUNCTIONS
//
// *****************************************************************************
// *****************************************************************************


// =============================================================================
// =============================================================================
//
// TEMPLATE VERSIONS based on Pack/Unpack
//
// =============================================================================
// =============================================================================

// =============================================================================
// Template MPI_Send, MPI_Bsend, MPI_Ssend, MPI_Rsend
// =============================================================================


template<class T>
int MPI_Send(T& Tobj, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Send((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Send(T& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Send((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Send(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2, MPI_Comm MC)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Send((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Bsend(T& Tobj, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Bsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Bsend(T& Tobj, int lb1, int lb2, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, lb2, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Bsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Bsend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2, MPI_Comm MC)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Bsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Ssend(T& Tobj, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Ssend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Ssend(T& Tobj, int lb1, int lb2, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, lb2, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Ssend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Ssend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2, MPI_Comm MC)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Ssend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Rsend(T& Tobj, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Rsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Rsend(T& Tobj, int lb1, int lb2, int i1, int i2, MPI_Comm MC)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, lb2, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Rsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

template<class T>
int MPI_Rsend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2, MPI_Comm MC)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Rsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC);
  return err;
}

// =============================================================================
// Template MPI_Isend, MPI_Ibsend, MPI_Issend, MPI_Irsend
// =============================================================================

template<class T>
int MPI_Isend(T& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Isend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Isend(T& Tobj, int lb1, int lb2, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, lb2, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Isend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Isend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2,
              MPI_Comm MC, MPI_Request* MR)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Isend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Ibsend(T& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Ibsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Ibsend(T& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Ibsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Ibsend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2,
               MPI_Comm MC, MPI_Request* MR)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Ibsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Issend(T& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Issend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Issend(T& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Issend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Issend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2,
               MPI_Comm MC, MPI_Request* MR)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Issend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Irsend(T& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Irsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Irsend(T& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  // Pack the T object, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type or other type parameter:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Irsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

template<class T>
int MPI_Irsend(T& Tobj, int lb1, int ub1, int lb2, int ub2, int i1, int i2,
               MPI_Comm MC, MPI_Request* MR)
{
  // For (short) matrix data types only! 

  // Pack the T subobject, then send it.
  // It is assumed that a suitable 
  // MPI_Pack(T&,void*,int,int,int,int,int,int*,MPI_Comm) is available!

  int pos=0;
  int err;

  // Call overloaded MPI_Pack for C-XSC data type:
  // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, 
                   commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
    return err;
  err=MPI_Irsend((void*)commbuffer, pos, MPI_PACKED, i1, i2, MC, MR);
  return err;
}

// =============================================================================
// Template MPI_Recv
// =============================================================================

template<class T>
int MPI_Recv(T& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int pos=0;
  int err;

  if ((err=MPI_Recv((void*)commbuffer, MPI_CXSC_BUFFERLEN, MPI_PACKED, 
                   i1, i2, MC, MS))!= MPI_SUCCESS)
    return err;

  // Call overloaded MPI_Unpack for C-XSC data type or other type parameter:
  // (MPI_Unpack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, Tobj, MC);  
  return err;
}

template<class T>
int MPI_Recv(T& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int pos=0;
  int err;

  if ((err=MPI_Recv((void*)commbuffer, MPI_CXSC_BUFFERLEN, MPI_PACKED, 
                   i1, i2, MC, MS))!= MPI_SUCCESS)
    return err;

  // Call overloaded MPI_Unpack for C-XSC data type or other type parameter:
  // (MPI_Unpack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, Tobj, lb1, ub1, MC);  
  return err;
}

template<class T>
int MPI_Recv(T& Tobj, int lb1, int ub1, int lb2, int ub2, 
             int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  // For (short) matrix data types only! 

  int pos=0;
  int err;

  if ((err=MPI_Recv((void*)commbuffer, MPI_CXSC_BUFFERLEN, MPI_PACKED, 
                   i1, i2, MC, MS))!= MPI_SUCCESS)
    return err;

  // Call overloaded MPI_Unpack for C-XSC data type or other type parameter:
  // (MPI_Unpack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true) 

  err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, Tobj, lb1, ub1, lb2, ub2, MC);  
  return err;
}

// *****************************************************************************
// *****************************************************************************
//
// COLLECTIVE COMMUNICATION FUNCTIONS
//
// *****************************************************************************
// *****************************************************************************


// =============================================================================
// =============================================================================
//
// TEMPLATE VERSIONS based on Pack/Unpack
//
// =============================================================================
// =============================================================================

// =============================================================================
// Template MPI_Bcast
// =============================================================================


template<class T>
int MPI_Bcast(T& Tobj, int root, MPI_Comm MC)
{
  int pos = 0, len = 0;
  int err;
  
  int mypid;
  err = MPI_Comm_rank(MC, &mypid);
  
  if (mypid == root)  {
    // Pack the T object, then broadcast it.
    // It is assumed that a suitable
    // MPI_Pack is available!

    // Call overloaded MPI_Pack for C-XSC data type:
    // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true)

    if (err = MPI_Pack(Tobj, (void*)commbuffer, MPI_CXSC_BUFFERLEN, &pos, MC) != MPI_SUCCESS)
      return err;
    len = pos;
  }
  
  if ((err = MPI_Bcast(&len, 1, MPI_INT, root, MC)) != MPI_SUCCESS) {
    return err;
  }
  
  if ((err = MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC)) != MPI_SUCCESS) {
    return err;
  }


  if (mypid != root) {
    // Call overloaded MPI_Unpack for C-XSC data type or other type parameter:
    // (MPI_Unpack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true)

    err = MPI_Unpack((void*)commbuffer, MPI_CXSC_BUFFERLEN, &pos, Tobj, MC);
  }
  return err;
}


template<class T>
int MPI_Bcast(T& Tobj, int lb1, int ub1, int lb2, int ub2, int root, MPI_Comm MC)
{
  // For (short) matrix data types only! 

  int pos=0;
  int len=0;
  int err;

  int mypid;
  err=MPI_Comm_rank(MC, &mypid);
  
  if (mypid==root)  
  {
    // Pack the T object, then broadcast it.
    // It is assumed that a suitable
    // MPI_Pack is available!

    // Call overloaded MPI_Pack for C-XSC data type:
    // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true)

    if ((err=MPI_Pack(Tobj, lb1, ub1, lb2, ub2, commbuffer, 
        MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
      return err;
    len=pos;
    if ((err=MPI_Bcast(&len, 1, MPI_INT, root, MC))!=MPI_SUCCESS)
      return err;
    if ((err=MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC))!=MPI_SUCCESS)
      return err;

  }
  else
  {

    // Broadcast-receive the T object, then unpack it.
    // It is assumed that a suitable
    // MPI_Unpack is available!

    if ((err=MPI_Bcast(&len, 1, MPI_INT, root, MC))!= MPI_SUCCESS)
      return err;
    if ((err=MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC))!= MPI_SUCCESS)
      return err;

    // Call overloaded MPI_Unpack for C-XSC data type or other type parameter:
    // (MPI_Unpack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true)

    err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, Tobj,
                   lb1, ub1, lb2, ub2, MC);
  }
  
  return err;
}

template<class T>
int MPI_Bcast(T& Tobj, int lb1, int ub1, int root, MPI_Comm MC)
{
  // For (short) vector data types only! 
  int pos=0;
  int len=0;
  int err;

  int mypid;
  err=MPI_Comm_rank(MC, &mypid);
  
  if (mypid==root)  
  {
    // Pack the T object, then broadcast it.
    // It is assumed that a suitable
    // MPI_Pack is available!

    // Call overloaded MPI_Pack for C-XSC data type:
    // (MPI_Pack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true)

    if ((err=MPI_Pack(Tobj, lb1, ub1, commbuffer, 
        MPI_CXSC_BUFFERLEN, &pos, MC))!=MPI_SUCCESS)
      return err;
    len=pos;
    if ((err=MPI_Bcast(&len, 1, MPI_INT, root, MC))!=MPI_SUCCESS)
      return err;
    if ((err=MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC))!=MPI_SUCCESS)
      return err;

  }
  else
  {

    // Broadcast-receive the T object, then unpack it.
    // It is assumed that a suitable
    // MPI_Unpack is available!

    if ((err=MPI_Bcast(&len, 1, MPI_INT, root, MC))!= MPI_SUCCESS)
      return err;
    if ((err=MPI_Bcast((void*)commbuffer, len, MPI_PACKED, root, MC))!= MPI_SUCCESS)
      return err;

    // Call overloaded MPI_Unpack for C-XSC data type or other type parameter:
    // (MPI_Unpack for C-XSC types checks if MPI_CXSC_TYPES_DEFINED is already true)

    err=MPI_Unpack(commbuffer, MPI_CXSC_BUFFERLEN, &pos, Tobj,
                   lb1, ub1, MC);
  }
  
  return err;
}

// =============================================================================
//
// Specializations for C-XSC base types
//
// C-XSC Base types real, interval, complex, cinterval do not require 
// separate implementations since being available as MPI data types.
//
// Specializations are nevertheless provided to avoid unnecessary packing 
// that would be done in instances of the templates given below.
//
// =============================================================================

// =============================================================================
// MPI-Send/Recv for real:
// =============================================================================
// Versions without packing to avoid template instantiation for base types

template<> int MPI_Send <real> (real& r, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Send(&r,1,MPI_CXSC_REAL,i1,i2,MC);
  return err;
}

template<> int MPI_Bsend <real> (real& r, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Bsend(&r,1,MPI_CXSC_REAL,i1,i2,MC);
  return err;
}

template<> int MPI_Ssend <real> (real& r, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Ssend(&r,1,MPI_CXSC_REAL,i1,i2,MC);
  return err;
}

template<> int MPI_Rsend <real> (real& r, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Rsend(&r,1,MPI_CXSC_REAL,i1,i2,MC);
  return err;
}

template<> int MPI_Isend <real> (real& r, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Isend(&r,1,MPI_CXSC_REAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Ibsend <real> (real& r, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ibsend(&r,1,MPI_CXSC_REAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Issend <real> (real& r, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Issend(&r,1,MPI_CXSC_REAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Irsend <real> (real& r, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Irsend(&r,1,MPI_CXSC_REAL,i1,i2,MC,MR);
  return err;
}


template<> int MPI_Recv <real> (real& r, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Recv(&r,1,MPI_CXSC_REAL,i1,i2,MC,MS);
  return err;
}

template<> int MPI_Bcast <real> (real& r, int root, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Bcast(&r,1,MPI_CXSC_REAL,root,MC);
  return err;
}

// =============================================================================
// MPI-Send/Recv for complex:
// =============================================================================
// Versions without packing to avoid template instantiation for base types

template<> int MPI_Send <complex> (complex& c, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Send(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC);
  return err;
}

template<> int MPI_Bsend <complex> (complex& c, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Bsend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC);
  return err;
}

template<> int MPI_Ssend <complex> (complex& c, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ssend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC);
  return err;
}

template<> int MPI_Rsend <complex> (complex& c, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Rsend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC);
  return err;
}


template<> int MPI_Isend <complex> (complex& c, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Isend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Ibsend <complex> (complex& c, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ibsend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Issend <complex> (complex& c, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Issend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Irsend <complex> (complex& c, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Irsend(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC,MR);
  return err;
}


template<> int MPI_Recv <complex> (complex& c, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Recv(&c,1,MPI_CXSC_COMPLEX,i1,i2,MC,MS);
  return err;
}

template<> int MPI_Bcast <complex> (complex& c, int root, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Bcast(&c,1,MPI_CXSC_COMPLEX,root,MC);
  return err;
}

// =============================================================================
// MPI-Send/Recv for interval:
// =============================================================================
// Versions without packing to avoid template instantiation for base types

template<> int MPI_Send <interval> (interval& ri, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Send(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Bsend <interval> (interval& ri, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Bsend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Ssend <interval> (interval& ri, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ssend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Rsend <interval> (interval& ri, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Rsend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Isend <interval> (interval& ri, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Isend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Ibsend <interval> (interval& ri, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ibsend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Issend <interval> (interval& ri, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Issend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Irsend <interval> (interval& ri, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Irsend(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Recv <interval> (interval& ri, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Recv(&ri,1,MPI_CXSC_INTERVAL,i1,i2,MC,MS);
  return err;
}

template<> int MPI_Bcast <interval> (interval& ri, int root, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Bcast(&ri,1,MPI_CXSC_INTERVAL,root,MC);
  return err;
}

// =============================================================================
// MPI-Send/Recv for cinterval:
// =============================================================================
// Versions without packing to avoid template instantiation for base types

template<> int MPI_Send <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Send(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Bsend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Bsend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Ssend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ssend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC);
  return err;
}

template<> int MPI_Rsend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Rsend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC);
  return err;
}


template<> int MPI_Isend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Isend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Ibsend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Ibsend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Issend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Issend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Irsend <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC, MPI_Request* MR)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Irsend(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC,MR);
  return err;
}

template<> int MPI_Recv <cinterval> (cinterval& ci, int i1, int i2, MPI_Comm MC, MPI_Status* MS)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;
  err=MPI_Recv(&ci,1,MPI_CXSC_CINTERVAL,i1,i2,MC,MS);
  return err;
}

template<> int MPI_Bcast <cinterval> (cinterval& ci, int root, MPI_Comm MC)
{
  int err;
  if (!MPI_CXSC_TYPES_DEFINED) 
     if ((err=MPI_Define_CXSC_Types())!=MPI_SUCCESS)
        return err;

  err=MPI_Bcast(&ci,1,MPI_CXSC_CINTERVAL,root,MC);
  return err;
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Explicit instantiation of templates
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

template int MPI_Send <l_real>(l_real& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <l_interval>(l_interval& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <l_complex>(l_complex& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <rvector>(rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <ivector>(ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <cvector>(cvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <civector>(civector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <srvector>(srvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <sivector>(sivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <scvector>(scvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <scivector>(scivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <rvector>(rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Send <ivector>(ivector& Tobj, int lb1, int ub1,  int i1, int i2, MPI_Comm MC);
template int MPI_Send <cvector>(cvector& Tobj, int lb1, int ub1,  int i1, int i2, MPI_Comm MC);
template int MPI_Send <civector>(civector& Tobj, int lb1, int ub1,  int i1, int i2, MPI_Comm MC);
template int MPI_Send <srvector>(srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Send <sivector>(sivector& Tobj, int lb1, int ub1,  int i1, int i2, MPI_Comm MC);
template int MPI_Send <scvector>(scvector& Tobj, int lb1, int ub1,  int i1, int i2, MPI_Comm MC);
template int MPI_Send <scivector>(scivector& Tobj, int lb1, int ub1,  int i1, int i2, MPI_Comm MC);
template int MPI_Send <l_rvector>(l_rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <l_ivector>(l_ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <rmatrix>(rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <imatrix>(imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <cmatrix>(cmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <cimatrix>(cimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <srmatrix>(srmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <simatrix>(simatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <scmatrix>(scmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <scimatrix>(scimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <rmatrix>(rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <imatrix>(imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <cmatrix>(cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <cimatrix>(cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <srmatrix>(srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <simatrix>(simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <scmatrix>(scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <scimatrix>(scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Send <l_rmatrix>(l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <l_imatrix>(l_imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <dotprecision>(dotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <idotprecision>(idotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <cdotprecision>(cdotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Send <cidotprecision>(cidotprecision& Tobj, int i1, int i2, MPI_Comm MC);

template int MPI_Bsend <l_real>(l_real& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <l_interval>(l_interval& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <l_complex>(l_complex& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <rvector>(rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <ivector>(ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cvector>(cvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <civector>(civector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <srvector>(srvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <sivector>(sivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scvector>(scvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scivector>(scivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <rvector>(rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <ivector>(ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cvector>(cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <civector>(civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <srvector>(srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <sivector>(sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scvector>(scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scivector>(scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <l_rvector>(l_rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <l_ivector>(l_ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <rmatrix>(rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <imatrix>(imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cmatrix>(cmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cimatrix>(cimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <srmatrix>(srmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <simatrix>(simatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scmatrix>(scmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scimatrix>(scimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <rmatrix>(rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <imatrix>(imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cmatrix>(cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cimatrix>(cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <srmatrix>(srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <simatrix>(simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scmatrix>(scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <scimatrix>(scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <l_rmatrix>(l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <l_imatrix>(l_imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <dotprecision>(dotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <idotprecision>(idotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cdotprecision>(cdotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Bsend <cidotprecision>(cidotprecision& Tobj, int i1, int i2, MPI_Comm MC);

template int MPI_Ssend <l_real>(l_real& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <l_interval>(l_interval& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <l_complex>(l_complex& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <rvector>(rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <ivector>(ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cvector>(cvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <civector>(civector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <srvector>(srvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <sivector>(sivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scvector>(scvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scivector>(scivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <rvector>(rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <ivector>(ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cvector>(cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <civector>(civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <srvector>(srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <sivector>(sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scvector>(scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scivector>(scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <l_rvector>(l_rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <l_ivector>(l_ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <rmatrix>(rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <imatrix>(imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cmatrix>(cmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cimatrix>(cimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <srmatrix>(srmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <simatrix>(simatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scmatrix>(scmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scimatrix>(scimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <rmatrix>(rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <imatrix>(imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cmatrix>(cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cimatrix>(cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <srmatrix>(srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <simatrix>(simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scmatrix>(scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <scimatrix>(scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <l_rmatrix>(l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <l_imatrix>(l_imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <dotprecision>(dotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <idotprecision>(idotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cdotprecision>(cdotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Ssend <cidotprecision>(cidotprecision& Tobj, int i1, int i2, MPI_Comm MC);

template int MPI_Rsend <l_real>(l_real& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <l_interval>(l_interval& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <l_complex>(l_complex& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <rvector>(rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <ivector>(ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cvector>(cvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <civector>(civector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <srvector>(srvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <sivector>(sivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scvector>(scvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scivector>(scivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <l_rvector>(l_rvector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <l_ivector>(l_ivector& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <rvector>(rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <ivector>(ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cvector>(cvector& Tobj, int lb1 , int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <civector>(civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <srvector>(srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <sivector>(sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scvector>(scvector& Tobj, int lb1 , int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scivector>(scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <rmatrix>(rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <imatrix>(imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cmatrix>(cmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cimatrix>(cimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <srmatrix>(srmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <simatrix>(simatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scmatrix>(scmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scimatrix>(scimatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <rmatrix>(rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <imatrix>(imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cmatrix>(cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cimatrix>(cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <srmatrix>(srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <simatrix>(simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scmatrix>(scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2, 
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <scimatrix>(scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <l_rmatrix>(l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <l_imatrix>(l_imatrix& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <dotprecision>(dotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <idotprecision>(idotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cdotprecision>(cdotprecision& Tobj, int i1, int i2, MPI_Comm MC);
template int MPI_Rsend <cidotprecision>(cidotprecision& Tobj, int i1, int i2, MPI_Comm MC);

template int MPI_Isend <l_real> (l_real& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <l_interval> (l_interval& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <l_complex> (l_complex& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <rvector> (rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <ivector> (ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cvector> (cvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <civector> (civector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <srvector> (srvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <sivector> (sivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scvector> (scvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scivector> (scivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <rvector> (rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <ivector> (ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cvector> (cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <civector> (civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <srvector> (srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <sivector> (sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scvector> (scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scivector> (scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <l_rvector> (l_rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <l_ivector> (l_ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <rmatrix> (rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <imatrix> (imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cmatrix> (cmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cimatrix> (cimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <srmatrix> (srmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <simatrix> (simatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scmatrix> (scmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scimatrix> (scimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <rmatrix> (rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <imatrix> (imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cmatrix> (cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cimatrix> (cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <srmatrix> (srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <simatrix> (simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scmatrix> (scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <scimatrix> (scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <l_rmatrix> (l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <l_imatrix> (l_imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <dotprecision> (dotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <idotprecision> (idotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cdotprecision> (cdotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Isend <cidotprecision> (cidotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);

template int MPI_Ibsend <l_real> (l_real& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <l_interval> (l_interval& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <l_complex> (l_complex& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <rvector> (rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <ivector> (ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cvector> (cvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <civector> (civector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <srvector> (srvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <sivector> (sivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scvector> (scvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scivector> (scivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <l_rvector> (l_rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <l_ivector> (l_ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <rvector> (rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <ivector> (ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cvector> (cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <civector> (civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <srvector> (srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <sivector> (sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scvector> (scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scivector> (scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <rmatrix> (rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <imatrix> (imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cmatrix> (cmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cimatrix> (cimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <srmatrix> (srmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <simatrix> (simatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scmatrix> (scmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scimatrix> (scimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <rmatrix> (rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <imatrix> (imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cmatrix> (cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cimatrix> (cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <srmatrix> (srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <simatrix> (simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scmatrix> (scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <scimatrix> (scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <l_rmatrix> (l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <l_imatrix> (l_imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <dotprecision> (dotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <idotprecision> (idotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cdotprecision> (cdotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Ibsend <cidotprecision> (cidotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);

template int MPI_Issend <l_real> (l_real& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <l_interval> (l_interval& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <l_complex> (l_complex& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <rvector> (rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <ivector> (ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cvector> (cvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <civector> (civector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <srvector> (srvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <sivector> (sivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scvector> (scvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scivector> (scivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <rvector> (rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <ivector> (ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cvector> (cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <civector> (civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <srvector> (srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <sivector> (sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scvector> (scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scivector> (scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <l_rvector> (l_rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <l_ivector> (l_ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <rmatrix> (rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <imatrix> (imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cmatrix> (cmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cimatrix> (cimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <srmatrix> (srmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <simatrix> (simatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scmatrix> (scmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scimatrix> (scimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <rmatrix> (rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <imatrix> (imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cmatrix> (cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cimatrix> (cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <srmatrix> (srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <simatrix> (simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scmatrix> (scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <scimatrix> (scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <l_rmatrix> (l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <l_imatrix> (l_imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <dotprecision> (dotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <idotprecision> (idotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cdotprecision> (cdotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Issend <cidotprecision> (cidotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);

template int MPI_Irsend <l_real> (l_real& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <l_interval> (l_interval& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <l_complex> (l_complex& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <rvector> (rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <ivector> (ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cvector> (cvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <civector> (civector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <srvector> (srvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <sivector> (sivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scvector> (scvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scivector> (scivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <rvector> (rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <ivector> (ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cvector> (cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <civector> (civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <srvector> (srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <sivector> (sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scvector> (scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scivector> (scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <l_rvector> (l_rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <l_ivector> (l_ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <rmatrix> (rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <imatrix> (imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cmatrix> (cmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cimatrix> (cimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <srmatrix> (srmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <simatrix> (simatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scmatrix> (scmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scimatrix> (scimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <rmatrix> (rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <imatrix> (imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cmatrix> (cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cimatrix> (cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <srmatrix> (srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <simatrix> (simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scmatrix> (scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <scimatrix> (scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <l_rmatrix> (l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <l_imatrix> (l_imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <dotprecision> (dotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <idotprecision> (idotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cdotprecision> (cdotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);
template int MPI_Irsend <cidotprecision> (cidotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Request* MR);

template int MPI_Recv <l_real> (l_real& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <l_interval> (l_interval& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <l_complex> (l_complex& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <rvector> (rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <ivector> (ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cvector> (cvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <civector> (civector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <srvector> (srvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <sivector> (sivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scvector> (scvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scivector> (scivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <rvector> (rvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <ivector> (ivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cvector> (cvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <civector> (civector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <srvector> (srvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <sivector> (sivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scvector> (scvector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scivector> (scivector& Tobj, int lb1, int ub1, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <l_rvector> (l_rvector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <l_ivector> (l_ivector& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <rmatrix> (rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <imatrix> (imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cmatrix> (cmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cimatrix> (cimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <srmatrix> (srmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <simatrix> (simatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scmatrix> (scmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scimatrix> (scimatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <rmatrix> (rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <imatrix> (imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cmatrix> (cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cimatrix> (cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <srmatrix> (srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <simatrix> (simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scmatrix> (scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <scimatrix> (scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                 int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <l_rmatrix> (l_rmatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <l_imatrix> (l_imatrix& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <dotprecision> (dotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <idotprecision> (idotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cdotprecision> (cdotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);
template int MPI_Recv <cidotprecision> (cidotprecision& Tobj, int i1, int i2, MPI_Comm MC, MPI_Status* MS);

template int MPI_Bcast <l_real> (l_real& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <l_interval> (l_interval& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <l_complex> (l_complex& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <rvector> (rvector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <ivector> (ivector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <cvector> (cvector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <civector> (civector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <srvector> (srvector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <sivector> (sivector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <scvector> (scvector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <scivector> (scivector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <rvector> (rvector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <ivector> (ivector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <cvector> (cvector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <civector> (civector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <srvector> (srvector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <sivector> (sivector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <scvector> (scvector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <scivector> (scivector& Tobj, int lb1, int ub1, int root, MPI_Comm MC);
template int MPI_Bcast <l_rvector> (l_rvector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <l_ivector> (l_ivector& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <rmatrix> (rmatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <imatrix> (imatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <cmatrix> (cmatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <cimatrix> (cimatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <srmatrix> (srmatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <simatrix> (simatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <scmatrix> (scmatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <scimatrix> (scimatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <rmatrix> (rmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <imatrix> (imatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <cmatrix> (cmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <cimatrix> (cimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <srmatrix> (srmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <simatrix> (simatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <scmatrix> (scmatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <scimatrix> (scimatrix& Tobj, int lb1, int ub1, int lb2, int ub2,
                                  int root, MPI_Comm MC);
template int MPI_Bcast <l_rmatrix> (l_rmatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <l_imatrix> (l_imatrix& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <dotprecision> (dotprecision& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <idotprecision> (idotprecision& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <cdotprecision> (cdotprecision& Tobj, int root, MPI_Comm MC);
template int MPI_Bcast <cidotprecision> (cidotprecision& Tobj, int root, MPI_Comm MC);

} //namespace cxsc

// End explicit instantiation of templates
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






