/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "stdlib.h"
#include "atom_vec_atomic_rad.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecAtomicRad::AtomVecAtomicRad(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;

  comm_x_only = comm_f_only = 0;
  size_forward = 12;
  size_reverse = 12;
  size_border = 15;
  size_velocity = 12;
  size_data_atom = 14;
  size_data_vel = 13;
  xcol_data = 3;

  atom->xrad_flag = 1;
  atom->vrad_flag = 1;
  atom->frad_flag = 1;
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecAtomicRad::grow(int n)
{
  if (n == 0) grow_nmax();
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  xrad = memory->grow(atom->xrad,nmax,3,"atom:xrad");
  xub = memory->grow(atom->xub,nmax*comm->nthreads,3,"atom:xub");
  xlb = memory->grow(atom->xlb,nmax*comm->nthreads,3,"atom:xlb");

  v = memory->grow(atom->v,nmax,3,"atom:v");
  vrad = memory->grow(atom->vrad,nmax,3,"atom:vrad");
  vub = memory->grow(atom->vub,nmax*comm->nthreads,3,"atom:vub");
  vlb = memory->grow(atom->vlb,nmax*comm->nthreads,3,"atom:vlb");

  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");
  frad = memory->grow(atom->frad,nmax*comm->nthreads,3,"atom:frad");
  fub = memory->grow(atom->fub,nmax*comm->nthreads,3,"atom:fub");
  flb = memory->grow(atom->flb,nmax*comm->nthreads,3,"atom:flb");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecAtomicRad::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  xrad = atom->xrad; 
  vrad = atom->vrad;
  frad = atom->frad;
  xub  = atom->xub;
  xlb  = atom->xlb;
  vub  = atom->vub;
  vlb  = atom->vlb;
  fub  = atom->fub;
  flb  = atom->flb;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecAtomicRad::copy(int i, int j, int delflag)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  xrad[j][0] = xrad[i][0];
  xrad[j][1] = xrad[i][1];
  xrad[j][2] = xrad[i][2];
  xub[j][0] = xub[i][0];
  xub[j][1] = xub[i][1];
  xub[j][2] = xub[i][2];
  xlb[j][0] = xlb[i][0];
  xlb[j][1] = xlb[i][1];
  xlb[j][2] = xlb[i][2];

  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];
  vrad[j][0] = vrad[i][0];
  vrad[j][1] = vrad[i][1];
  vrad[j][2] = vrad[i][2];
  vub[j][0] = vub[i][0];
  vub[j][1] = vub[i][1];
  vub[j][2] = vub[i][2];
  vlb[j][0] = vlb[i][0];
  vlb[j][1] = vlb[i][1];
  vlb[j][2] = vlb[i][2];



  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  // std::cout << "pack_comm() function is used" << std::endl;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = xrad[j][0];
      buf[m++] = xrad[j][1];
      buf[m++] = xrad[j][2];
      buf[m++] = xub[j][0];
      buf[m++] = xub[j][1];
      buf[m++] = xub[j][2];
      buf[m++] = xlb[j][0];
      buf[m++] = xlb[j][1];
      buf[m++] = xlb[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = xrad[j][0];
      buf[m++] = xrad[j][1];
      buf[m++] = xrad[j][2];
      buf[m++] = xub[j][0];
      buf[m++] = xub[j][1];
      buf[m++] = xub[j][2];
      buf[m++] = xlb[j][0];
      buf[m++] = xlb[j][1];
      buf[m++] = xlb[j][2];
    }
  }

  // std::cout << "xrad[j] = [" << xrad[j][0] << ", " 
  //   << xrad[j][1] << ", "
  //   << xrad[j][2] << "]" << std::endl;
  // std::cout << "x[j] = [" << x[j][0] << ", " 
  //   << x[j][1] << ", "
  //   << x[j][2] << "]" << std::endl;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_comm_vel(int n, int *list, double *buf,
                                 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;
  // std::cout << "pack_comm_vel() function is used" << std::endl;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = xrad[j][0];
      buf[m++] = xrad[j][1];
      buf[m++] = xrad[j][2];
      buf[m++] = xub[j][0];
      buf[m++] = xub[j][1];
      buf[m++] = xub[j][2];
      buf[m++] = xlb[j][0];
      buf[m++] = xlb[j][1];
      buf[m++] = xlb[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = vrad[j][0];
      buf[m++] = vrad[j][1];
      buf[m++] = vrad[j][2];
      buf[m++] = vub[j][0];
      buf[m++] = vub[j][1];
      buf[m++] = vub[j][2];
      buf[m++] = vlb[j][0];
      buf[m++] = vlb[j][1];
      buf[m++] = vlb[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = xrad[j][0];
        buf[m++] = xrad[j][1];
        buf[m++] = xrad[j][2];
        buf[m++] = xub[j][0];
        buf[m++] = xub[j][1];
        buf[m++] = xub[j][2];
        buf[m++] = xlb[j][0];
        buf[m++] = xlb[j][1];
        buf[m++] = xlb[j][2];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = vrad[j][0];
        buf[m++] = vrad[j][1];
        buf[m++] = vrad[j][2];
        buf[m++] = vub[j][0];
        buf[m++] = vub[j][1];
        buf[m++] = vub[j][2];
        buf[m++] = vlb[j][0];
        buf[m++] = vlb[j][1];
        buf[m++] = vlb[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = xrad[j][0];
        buf[m++] = xrad[j][1];
        buf[m++] = xrad[j][2];
        buf[m++] = xub[j][0];
        buf[m++] = xub[j][1];
        buf[m++] = xub[j][2];
        buf[m++] = xlb[j][0];
        buf[m++] = xlb[j][1];
        buf[m++] = xlb[j][2];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
          buf[m++] = vrad[j][0];
          buf[m++] = vrad[j][1];
          buf[m++] = vrad[j][2];
          buf[m++] = vub[j][0];
          buf[m++] = vub[j][1];
          buf[m++] = vub[j][2];
          buf[m++] = vlb[j][0];
          buf[m++] = vlb[j][1];
          buf[m++] = vlb[j][2];

        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = vrad[j][0];
          buf[m++] = vrad[j][1];
          buf[m++] = vrad[j][2];
          buf[m++] = vub[j][0];
          buf[m++] = vub[j][1];
          buf[m++] = vub[j][2];
          buf[m++] = vlb[j][0];
          buf[m++] = vlb[j][1];
          buf[m++] = vlb[j][2];
        }
      }
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecAtomicRad::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    xrad[i][0] = buf[m++];
    xrad[i][1] = buf[m++];
    xrad[i][2] = buf[m++];
    xub[i][0] = buf[m++]; 
    xub[i][1] = buf[m++];
    xub[i][2] = buf[m++];
    xlb[i][0] = buf[m++];
    xlb[i][1] = buf[m++];
    xlb[i][2] = buf[m++];
    // std::cout << "x[i][*] = [" << xrad[i][0] << ","
    //   << xrad[i][1] << ","
    //   << xrad[i][2] << "]" << std::endl;
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecAtomicRad::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;
  // std::cout << "unpack_comm_vel() function is used" << std::endl;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    xrad[i][0] = buf[m++];
    xrad[i][1] = buf[m++];
    xrad[i][2] = buf[m++];
    xub[i][0] = buf[m++]; 
    xub[i][1] = buf[m++];
    xub[i][2] = buf[m++];
    xlb[i][0] = buf[m++];
    xlb[i][1] = buf[m++];
    xlb[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    vrad[i][0] = buf[m++];
    vrad[i][1] = buf[m++];
    vrad[i][2] = buf[m++];
    vub[i][0] = buf[m++]; 
    vub[i][1] = buf[m++];
    vub[i][2] = buf[m++];
    vlb[i][0] = buf[m++];
    vlb[i][1] = buf[m++];
    vlb[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
    buf[m++] = frad[i][0];
    buf[m++] = frad[i][1];
    buf[m++] = frad[i][2];
    buf[m++] = fub[i][0];
    buf[m++] = fub[i][1];
    buf[m++] = fub[i][2];
    buf[m++] = flb[i][0];
    buf[m++] = flb[i][1];
    buf[m++] = flb[i][2];
    // test
    // buf[m++] = 0.0;
    // buf[m++] = 0.0;
    // buf[m++] = 0.0;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecAtomicRad::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
    frad[j][0] += buf[m++];
    frad[j][1] += buf[m++];
    frad[j][2] += buf[m++];
    fub[j][0] += buf[m++];
    fub[j][1] += buf[m++];
    fub[j][2] += buf[m++];
    flb[j][0] += buf[m++];
    flb[j][1] += buf[m++];
    flb[j][2] += buf[m++];
    //test 
    // std::cout << "frad = [" << frad[j][0] << ", "
    // << frad[j][1] << "," 
    // << frad[j][1] << "]" << std::endl;
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_border(int n, int *list, double *buf,
                               int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  // std::cout << "pack_boder() function is used" << std::endl;
  // std::cout << "xrad = [" << xrad[i][0] << ", "
  //     << xrad[i][1] << ", "
  //     << xrad[i][2] << "] "
  //     << std::endl;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = xrad[j][0];
      buf[m++] = xrad[j][1];
      buf[m++] = xrad[j][2];
      buf[m++] = xub[j][0];
      buf[m++] = xub[j][1];
      buf[m++] = xub[j][2];
      buf[m++] = xlb[j][0];
      buf[m++] = xlb[j][1];
      buf[m++] = xlb[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = xrad[j][0];
      buf[m++] = xrad[j][1];
      buf[m++] = xrad[j][2];
      buf[m++] = xub[j][0];
      buf[m++] = xub[j][1];
      buf[m++] = xub[j][2];
      buf[m++] = xlb[j][0];
      buf[m++] = xlb[j][1];
      buf[m++] = xlb[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_border_vel(int n, int *list, double *buf,
                                   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = xrad[j][0];
      buf[m++] = xrad[j][1];
      buf[m++] = xrad[j][2];
      buf[m++] = xub[j][0];
      buf[m++] = xub[j][1];
      buf[m++] = xub[j][2];
      buf[m++] = xlb[j][0];
      buf[m++] = xlb[j][1];
      buf[m++] = xlb[j][2];
      buf[m++] = ubuf(tag[j]).d;
      buf[m++] = ubuf(type[j]).d;
      buf[m++] = ubuf(mask[j]).d;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = vrad[j][0];
      buf[m++] = vrad[j][1];
      buf[m++] = vrad[j][2];
      buf[m++] = vub[j][0];
      buf[m++] = vub[j][1];
      buf[m++] = vub[j][2];
      buf[m++] = vlb[j][0];
      buf[m++] = vlb[j][1];
      buf[m++] = vlb[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = xrad[j][0];
        buf[m++] = xrad[j][1];
        buf[m++] = xrad[j][2];
        buf[m++] = xub[j][0];
        buf[m++] = xub[j][1];
        buf[m++] = xub[j][2];
        buf[m++] = xlb[j][0];
        buf[m++] = xlb[j][1];
        buf[m++] = xlb[j][2];
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = vrad[j][0];
        buf[m++] = vrad[j][1];
        buf[m++] = vrad[j][2];
        buf[m++] = vub[j][0];
        buf[m++] = vub[j][1];
        buf[m++] = vub[j][2];
        buf[m++] = vlb[j][0];
        buf[m++] = vlb[j][1];
        buf[m++] = vlb[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = xub[j][0];
        buf[m++] = xub[j][1];
        buf[m++] = xub[j][2];
        buf[m++] = xlb[j][0];
        buf[m++] = xlb[j][1];
        buf[m++] = xlb[j][2];
        buf[m++] = ubuf(tag[j]).d;
        buf[m++] = ubuf(type[j]).d;
        buf[m++] = ubuf(mask[j]).d;
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
          buf[m++] = vrad[j][0];
          buf[m++] = vrad[j][1];
          buf[m++] = vrad[j][2];
          buf[m++] = vub[j][0];
          buf[m++] = vub[j][1];
          buf[m++] = vub[j][2];
          buf[m++] = vlb[j][0];
          buf[m++] = vlb[j][1];
          buf[m++] = vlb[j][2];
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = vrad[j][0];
          buf[m++] = vrad[j][1];
          buf[m++] = vrad[j][2];
          buf[m++] = vub[j][0];
          buf[m++] = vub[j][1];
          buf[m++] = vub[j][2];
          buf[m++] = vlb[j][0];
          buf[m++] = vlb[j][1];
          buf[m++] = vlb[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecAtomicRad::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    xrad[i][0] = buf[m++];
    xrad[i][1] = buf[m++];
    xrad[i][2] = buf[m++];
    xub[i][0] = buf[m++];
    xub[i][1] = buf[m++];
    xub[i][2] = buf[m++];
    xlb[i][0] = buf[m++];
    xlb[i][1] = buf[m++];
    xlb[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    // std::cout << "xrad = [" << xrad[i][0]
    //   << xrad[i][1] << ", "
    //   << xrad[i][2] << "] "
    //   << std::endl;
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecAtomicRad::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    xrad[i][0] = buf[m++];
    xrad[i][1] = buf[m++];
    xrad[i][2] = buf[m++];
    xub[i][0] = buf[m++];
    xub[i][1] = buf[m++];
    xub[i][2] = buf[m++];
    xlb[i][0] = buf[m++];
    xlb[i][1] = buf[m++];
    xlb[i][2] = buf[m++];
    tag[i] = (tagint) ubuf(buf[m++]).i;
    type[i] = (int) ubuf(buf[m++]).i;
    mask[i] = (int) ubuf(buf[m++]).i;
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    vrad[i][0] = buf[m++];
    vrad[i][1] = buf[m++];
    vrad[i][2] = buf[m++];
    vub[i][0] = buf[m++];
    vub[i][1] = buf[m++];
    vub[i][2] = buf[m++];
    vlb[i][0] = buf[m++];
    vlb[i][1] = buf[m++];
    vlb[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = xrad[i][0];
  buf[m++] = xrad[i][1];
  buf[m++] = xrad[i][2];
  buf[m++] = xub[i][0];
  buf[m++] = xub[i][1];
  buf[m++] = xub[i][2];
  buf[m++] = xlb[i][0];
  buf[m++] = xlb[i][1];
  buf[m++] = xlb[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = vrad[i][0];
  buf[m++] = vrad[i][1];
  buf[m++] = vrad[i][2];
  buf[m++] = vub[i][0];
  buf[m++] = vub[i][1];
  buf[m++] = vub[i][2];
  buf[m++] = vlb[i][0];
  buf[m++] = vlb[i][1];
  buf[m++] = vlb[i][2];
  // buf[m++] = 3.0;
  // buf[m++] = 4.0;
  // buf[m++] = 5.0;
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;

  // std::cout << "vrad[i] = [" 
  //   << vrad[i][0] << ", "
  //   << vrad[i][1] << ", "
  //   << vrad[i][2] << "] "
  //   << std::endl;
  // std::cout << "v[i] = [" 
  //   << v[i][0] << ", "
  //   << v[i][1] << ", "
  //   << v[i][2] << "] "
  //   << std::endl;
  // std::cout << "xrad[i] = [" 
  //   << xrad[i][0] << ", "
  //   << xrad[i][1] << ", "
  //   << xrad[i][2] << "] "
  //   << std::endl;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  // begin debug - check if the vub,vlb can be retrieved from the buf[] message
  // for (int k = 0; k < 3; k++ ) {
  //   std::cout << "v["   << i << "][" << k << "] =" << v[i][k]   << std::endl;
  //   std::cout << "vub[" << i << "][" << k << "] =" << vub[i][k] << std::endl;
  //   std::cout << "vlb[" << i << "][" << k << "] =" << vlb[i][k] << std::endl;
  // }
  // end debug
  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecAtomicRad::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  xrad[nlocal][0] = buf[m++];
  xrad[nlocal][1] = buf[m++];
  xrad[nlocal][2] = buf[m++];
  xub[nlocal][0] = buf[m++];
  xub[nlocal][1] = buf[m++];
  xub[nlocal][2] = buf[m++];
  xlb[nlocal][0] = buf[m++];
  xlb[nlocal][1] = buf[m++];
  xlb[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  vrad[nlocal][0] = buf[m++];
  vrad[nlocal][1] = buf[m++];
  vrad[nlocal][2] = buf[m++];
  vub[nlocal][0] = buf[m++];
  vub[nlocal][1] = buf[m++];
  vub[nlocal][2] = buf[m++];
  vlb[nlocal][0] = buf[m++];
  vlb[nlocal][1] = buf[m++];
  vlb[nlocal][2] = buf[m++];
    // std::cout << "vrad[" << nlocal << "] = [" << vrad[nlocal][0] 
    //   << "," << vrad[nlocal][1]
    //   << "," << vrad[nlocal][2] << "]" << std::endl;
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
        unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecAtomicRad::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 17 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
        n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecAtomicRad::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = xrad[i][0];
  buf[m++] = xrad[i][1];
  buf[m++] = xrad[i][2];
  buf[m++] = xub[i][0];
  buf[m++] = xub[i][1];
  buf[m++] = xub[i][2];
  buf[m++] = xlb[i][0];
  buf[m++] = xlb[i][1];
  buf[m++] = xlb[i][2];
  buf[m++] = ubuf(tag[i]).d;
  buf[m++] = ubuf(type[i]).d;
  buf[m++] = ubuf(mask[i]).d;
  buf[m++] = ubuf(image[i]).d;
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = vrad[i][0];
  buf[m++] = vrad[i][1];
  buf[m++] = vrad[i][2];
  buf[m++] = vub[i][0];
  buf[m++] = vub[i][1];
  buf[m++] = vub[i][2];
  buf[m++] = vlb[i][0];
  buf[m++] = vlb[i][1];
  buf[m++] = vlb[i][2];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecAtomicRad::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  xrad[nlocal][0] = buf[m++];
  xrad[nlocal][1] = buf[m++];
  xrad[nlocal][2] = buf[m++];
  xub[nlocal][0] = buf[m++];
  xub[nlocal][1] = buf[m++];
  xub[nlocal][2] = buf[m++];
  xlb[nlocal][0] = buf[m++];
  xlb[nlocal][1] = buf[m++];
  xlb[nlocal][2] = buf[m++];
  tag[nlocal] = (tagint) ubuf(buf[m++]).i;
  type[nlocal] = (int) ubuf(buf[m++]).i;
  mask[nlocal] = (int) ubuf(buf[m++]).i;
  image[nlocal] = (imageint) ubuf(buf[m++]).i;
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  vrad[nlocal][0] = buf[m++];
  vrad[nlocal][1] = buf[m++];
  vrad[nlocal][2] = buf[m++];
  vub[nlocal][0] = buf[m++];
  vub[nlocal][1] = buf[m++];
  vub[nlocal][2] = buf[m++];
  vlb[nlocal][0] = buf[m++];
  vlb[nlocal][1] = buf[m++];
  vlb[nlocal][2] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecAtomicRad::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  xrad[nlocal][0] = 0.0;
  xrad[nlocal][1] = 0.0;
  xrad[nlocal][2] = 0.0;
  xub[nlocal][0] = coord[0];
  xub[nlocal][1] = coord[1];
  xub[nlocal][2] = coord[2];
  xlb[nlocal][0] = coord[0];
  xlb[nlocal][1] = coord[1];
  xlb[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((imageint) IMGMAX << IMG2BITS) |
    ((imageint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  vrad[nlocal][0] = 0.0;
  vrad[nlocal][1] = 0.0;
  vrad[nlocal][2] = 0.0;
  vub[nlocal][0] = 0.0;
  vub[nlocal][1] = 0.0;
  vub[nlocal][2] = 0.0;
  vlb[nlocal][0] = 0.0;
  vlb[nlocal][1] = 0.0;
  vlb[nlocal][2] = 0.0;
  // std::cout << "create_atom() in atom style" << std::endl;
  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecAtomicRad::data_atom(double *coord, imageint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = ATOTAGINT(values[0]);
  type[nlocal] = atoi(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  xrad[nlocal][0] = 0.0;
  xrad[nlocal][1] = 0.0;
  xrad[nlocal][2] = 0.0;
  xub[nlocal][0] = coord[0];
  xub[nlocal][1] = coord[1];
  xub[nlocal][2] = coord[2];
  xlb[nlocal][0] = coord[0];
  xlb[nlocal][1] = coord[1];
  xlb[nlocal][2] = coord[2];
  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  vrad[nlocal][0] = 0.0;
  vrad[nlocal][1] = 0.0;
  vrad[nlocal][2] = 0.0;
  vub[nlocal][0] = 0.0;
  vub[nlocal][1] = 0.0;
  vub[nlocal][2] = 0.0;
  vlb[nlocal][0] = 0.0;
  vlb[nlocal][1] = 0.0;
  vlb[nlocal][2] = 0.0;


  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecAtomicRad::pack_data(double **buf)
{
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0]  = ubuf(tag[i]).d;
    buf[i][1]  = ubuf(type[i]).d;
    buf[i][2]  = x[i][0];
    buf[i][3]  = x[i][1];
    buf[i][4]  = x[i][2];
    buf[i][5]  = xrad[i][0];
    buf[i][6]  = xrad[i][1];
    buf[i][7]  = xrad[i][2];
    buf[i][8]  = xub[i][0];
    buf[i][9]  = xub[i][1];
    buf[i][10] = xub[i][2];
    buf[i][11] = xlb[i][0];
    buf[i][12] = xlb[i][1];
    buf[i][13] = xlb[i][2];
    buf[i][14] = ubuf((image[i] & IMGMASK) - IMGMAX).d;
    buf[i][15] = ubuf((image[i] >> IMGBITS & IMGMASK) - IMGMAX).d;
    buf[i][16] = ubuf((image[i] >> IMG2BITS) - IMGMAX).d;
  }
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecAtomicRad::write_data(FILE *fp, int n, double **buf)
{
  for (int i = 0; i < n; i++)
    fprintf(fp,TAGINT_FORMAT " %d %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %-1.16e %d %d %d\n",
            (tagint) ubuf(buf[i][0]).i,(int) ubuf(buf[i][1]).i,
            buf[i][2],buf[i][3],buf[i][4],
            buf[i][5],buf[i][6],buf[i][7],
            buf[i][8],buf[i][9],buf[i][10],
            buf[i][11],buf[i][12],buf[i][13],
            (int) ubuf(buf[i][14]).i,(int) ubuf(buf[i][15]).i,
            (int) ubuf(buf[i][16]).i);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecAtomicRad::memory_usage()
{
  bigint bytes = 0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);
  if (atom->memcheck("xrad")) bytes += memory->usage(xrad,nmax,3);
  if (atom->memcheck("xub")) bytes += memory->usage(xub,nmax,3);
  if (atom->memcheck("xlb")) bytes += memory->usage(xlb,nmax,3);
  if (atom->memcheck("vrad")) bytes += memory->usage(vrad,nmax,3);
  if (atom->memcheck("vub")) bytes += memory->usage(vub,nmax,3);
  if (atom->memcheck("vlb")) bytes += memory->usage(vlb,nmax,3);
  if (atom->memcheck("frad")) bytes += memory->usage(frad,nmax*comm->nthreads,3);
  if (atom->memcheck("fub")) bytes += memory->usage(fub,nmax*comm->nthreads,3);
  if (atom->memcheck("flb")) bytes += memory->usage(flb,nmax*comm->nthreads,3);

  return bytes;
}
