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

#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "compute_pressure.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePressure::ComputePressure(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal compute pressure command");
  if (igroup) error->all(FLERR,"Compute pressure must use group all");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 0;
  pressflag = 1;
  timeflag = 1;

  // store temperature ID used by pressure computation
  // insure it is valid for temperature computation

  if (strcmp(arg[3],"NULL") == 0) id_temp = NULL;
  else {
    int n = strlen(arg[3]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[3]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute pressure temperature ID");
    if (modify->compute[icompute]->tempflag == 0)
      error->all(FLERR,"Compute pressure temperature ID does not "
                 "compute temperature");
  }

  // process optional args

  if (narg == 4) {
    keflag = 1;
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = fixflag = 1;
  } else {
    keflag = 0;
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = fixflag = 0;
    int iarg = 4;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"ke") == 0) keflag = 1;
      else if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg],"fix") == 0) fixflag = 1;
      else if (strcmp(arg[iarg],"virial") == 0) {
        pairflag = 1;
        bondflag = angleflag = dihedralflag = improperflag = 1;
        kspaceflag = fixflag = 1;
      } else error->all(FLERR,"Illegal compute pressure command");
      iarg++;
    }
  }

  // error check

  if (keflag && id_temp == NULL) 
    error->all(FLERR,"Compute pressure requires temperature ID "
	       "to include kinetic energy");

  vector = new double[6];
  vectorrad = new double[6];  // add rad
  vector_lb = new double[6];   // add rad
  vector_ub = new double[6];   // add rad
  nvirial = 0;
  nvirialrad = 0; // add rad
  nvirial_lb = 0; // add rad
  nvirial_ub = 0; // add rad
  vptr = NULL;
  vptrrad = NULL; // add rad
  vptr_ub = NULL; // add rad
  vptr_lb = NULL; // add rad
}

/* ---------------------------------------------------------------------- */

ComputePressure::~ComputePressure()
{
  delete [] id_temp;
  delete [] vector;
  delete [] vectorrad; // add rad
  delete [] vector_ub; // add rad
  delete [] vector_lb; // add rad
  delete [] vptr;
  delete [] vptrrad; // add rad
  delete [] vptr_ub; // add rad
  delete [] vptr_lb; // add rad
}

/* ---------------------------------------------------------------------- */

void ComputePressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // set temperature compute, must be done in init()
  // fixes could have changed or compute_modify could have changed it

  if (keflag) {
    int icompute = modify->find_compute(id_temp);
    if (icompute < 0)
      error->all(FLERR,"Could not find compute pressure temperature ID");
    temperature = modify->compute[icompute];
  }

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  delete [] vptrrad;  // add rad
  delete [] vptr_ub;   // add rad
  delete [] vptr_lb;   // add rad
  nvirial = 0;
  nvirialrad = 0; // add rad - counter parallel with nvirial
  vptr = NULL;
  vptrrad = NULL;   // add rad
  vptr_ub = NULL;    // add rad
  vptr_lb = NULL;    // add rad

  if (pairflag && force->pair) nvirial++;
  if (pairflag && force->pair && atom->vrad_flag) nvirialrad++; // add rad

  if (bondflag && atom->molecular && force->bond) nvirial++;
  if (angleflag && atom->molecular && force->angle) nvirial++;
  if (dihedralflag && atom->molecular && force->dihedral) nvirial++;
  if (improperflag && atom->molecular && force->improper) nvirial++;
  if (fixflag)
    for (int i = 0; i < modify->nfix; i++) {
      if (modify->fix[i]->virial_flag) { 
        nvirial++; 
      }
    }

  if (nvirial) {
    vptr = new double*[nvirial];
    vptrrad = new double*[nvirialrad];  // add rad
    vptr_ub = new double*[nvirial_ub];    // add rad
    vptr_lb = new double*[nvirial_lb];    // add rad
    nvirial = 0;
    nvirialrad = 0;   // add rad
    nvirial_ub = 0;    // add rad
    nvirial_lb = 0;    // add rad
    if (pairflag && force->pair) { // add rad only active in pairflag
      vptr[nvirial++] = force->pair->virial;
      // std::cout << "fix here - work! compute_pressure.cpp" << std::endl; 
      if (atom->vrad_flag) vptrrad[nvirialrad++] = force->pair->virialrad;  // nvirialrad++
      if (atom->vrad_flag) vptr_ub[nvirial_ub++]   = force->pair->virial_ub;   // nvirial_ub++
      if (atom->vrad_flag) vptr_lb[nvirial_lb++]   = force->pair->virial_lb;   // nvirial_lb++
    }
    if (bondflag && force->bond) vptr[nvirial++] = force->bond->virial;
    if (angleflag && force->angle) vptr[nvirial++] = force->angle->virial;
    if (dihedralflag && force->dihedral)
      vptr[nvirial++] = force->dihedral->virial;
    if (improperflag && force->improper)
      vptr[nvirial++] = force->improper->virial;
    if (fixflag)
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->virial_flag){
          vptr[nvirial++] = modify->fix[i]->virial;
        }
  }

  // flag Kspace contribution separately, since not summed across procs

  if (kspaceflag && force->kspace) kspace_virial = force->kspace->virial;
  else kspace_virial = NULL;
}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double t;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t = temperature->compute_scalar();
    else t = temperature->scalar;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3,3);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2,2);
    if (keflag)
      scalar = (temperature->dof * boltz * t +
                virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
}
/* ----------------------------------------------------------------------
   compute_scalar for upperbound
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar_ub()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double t_ub;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t_ub = temperature->compute_scalar_ub();
    else t_ub = temperature->scalar_ub;
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3,3);
    if (keflag)
      scalar_ub = (temperature->dof * boltz * t_ub +
                virial_ub[0] + virial_ub[1] + virial_ub[2]) / 3.0 * inv_volume * nktv2p;
    else
      scalar_ub = (virial_ub[0] + virial_ub[1] + virial_ub[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2,2);
    if (keflag)
      scalar_ub = (temperature->dof * boltz * t_ub +
                virial_ub[0] + virial_ub[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar_ub = (virial_ub[0] + virial_ub[1]) / 2.0 * inv_volume * nktv2p;
  }

  // begin debug
  // std::cout << "scalar_ub = " << scalar_ub << std::endl;
  // end debug

  return scalar_ub;
}
/* ----------------------------------------------------------------------
   compute_scalar for lowerbound
------------------------------------------------------------------------- */

double ComputePressure::compute_scalar_lb()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  // invoke temperature if it hasn't been already

  double t_lb;
  if (keflag) {
    if (temperature->invoked_scalar != update->ntimestep)
      t_lb = temperature->compute_scalar_lb();
    else t_lb = temperature->scalar_lb;
  }

  // begin debug
  // std::cout << "temperature->scalar_lb = " << temperature->scalar_lb << std::endl; // is a number
  // for (int j =  0; j < 6; j++) {
  //   std::cout << "virial_lb[" << j << "] = " << virial_lb[j] << std::endl;
  // }
  // end debug

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3,3);
    if (keflag)
      scalar_lb = (temperature->dof * boltz * t_lb +
                virial_lb[0] + virial_lb[1] + virial_lb[2]) / 3.0 * inv_volume * nktv2p;
    else
      scalar_lb = (virial_lb[0] + virial_lb[1] + virial_lb[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2,2);
    if (keflag)
      scalar_lb = (temperature->dof * boltz * t_lb +
                virial_lb[0] + virial_lb[1]) / 2.0 * inv_volume * nktv2p;
    else
      scalar_lb = (virial_lb[0] + virial_lb[1]) / 2.0 * inv_volume * nktv2p;
  }

  // begin debug
  // std::cout << "scalar_lb = " << scalar_lb << std::endl;
  // end debug
  return scalar_lb;
}

/* ----------------------------------------------------------------------
   compute pressure tensor
   assume KE tensor has already been computed
------------------------------------------------------------------------- */

void ComputePressure::compute_vector()
{
  invoked_vector = update->ntimestep;
  if (update->vflag_global != invoked_vector)
    error->all(FLERR,"Virial was not tallied on needed timestep");

  if (force->kspace && kspace_virial && force->kspace->scalar_pressure_flag)
    error->all(FLERR,"Kspace_modify pressure/scalar no required "
               "for components of pressure tensor with kspace_style msm");

  // invoke temperature if it hasn't been already

  double *ke_tensor;
  double *ke_tensor_rad; // add rad
  double *ke_tensor_ub; // add ub
  double *ke_tensor_lb; // add lb
  if (keflag) {
    if (temperature->invoked_vector != update->ntimestep)
      temperature->compute_vector();
      ke_tensor     = temperature->vector;
    // add rad computation - only if 3 flags are on
    if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
      ke_tensor_rad = temperature->vectorrad;    // add rad 
      ke_tensor_ub  = temperature->vector_ub;    // add rad
      ke_tensor_lb  = temperature->vector_lb;    // add rad
      // std::cout << "flag OK, compute ke_tensor_* ..." << std::endl;
    }
  }

  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(6,3);
    if (keflag) {
      for (int i = 0; i < 6; i++) {
        vector[i] = (ke_tensor[i] + virial[i]) * inv_volume * nktv2p;
        // add rad
        if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
          vectorrad[i] = (ke_tensor_rad[i] + virialrad[i]) * inv_volume * nktv2p;   // add rad
          vector_ub[i] = (ke_tensor_ub[i]  + virial_ub[i]) * inv_volume * nktv2p;   // add rad
          vector_lb[i] = (ke_tensor_lb[i]  + virial_lb[i]) * inv_volume * nktv2p;   // add rad
        }

        // begin debug
        // fix here - work!

        // check flag
        // std::cout << "atom->xrad_flag = " << atom->xrad_flag << std::endl;
        // std::cout << "atom->vrad_flag = " << atom->vrad_flag << std::endl;
        // std::cout << "atom->frad_flag = " << atom->frad_flag << std::endl;

        // status: no - vector* not the same 
        // std::cout << "vector_ub[" << i << "] = " << vector_ub[i] << std::endl;
        // std::cout << "vector_lb[" << i << "] = " << vector_lb[i] << std::endl;
        // std::cout << "vector[" << i << "] = " << vector[i] << std::endl;

        // status: no - ke_tensor _ub and _lb all = 0 -- ke_tensor* are not the same
        // std::cout << "ke_tensor_ub[" << i << "] = " << ke_tensor_ub[i] << std::endl;
        // std::cout << "ke_tensor_lb[" << i << "] = " << ke_tensor_lb[i] << std::endl;
        // std::cout << "ke_tensor["    << i << "] = " << ke_tensor[i] << std::endl;

        // std::cout << "virial_ub[" << i << "] = " << virial_ub[i] << std::endl;
        // std::cout << "virial_lb[" << i << "] = " << virial_lb[i] << std::endl;
        // std::cout << "virial[" << i << "] = " << virial[i] << std::endl;
        // end debug
      }
    } else
      for (int i = 0; i < 6; i++) {
        vector[i] = virial[i] * inv_volume * nktv2p;
        // add rad
        if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag)
          vectorrad[i] = virialrad[i] * inv_volume * nktv2p;     // add rad
          vector_ub[i] = virial_ub[i] * inv_volume * nktv2p;    // add rad
          vector_lb[i] = virial_lb[i] * inv_volume * nktv2p;    // add rad
      }
      // finish rad implementation here (only in 3-dimensions)
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(4,2);
    if (keflag) {
      vector[0] = (ke_tensor[0] + virial[0]) * inv_volume * nktv2p;
      vector[1] = (ke_tensor[1] + virial[1]) * inv_volume * nktv2p;
      vector[3] = (ke_tensor[3] + virial[3]) * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    } else {
      vector[0] = virial[0] * inv_volume * nktv2p;
      vector[1] = virial[1] * inv_volume * nktv2p;
      vector[3] = virial[3] * inv_volume * nktv2p;
      vector[2] = vector[4] = vector[5] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputePressure::virial_compute(int n, int ndiag)
{
  int i,j;
  double v[6],*vcomponent;
  double vrad[6],*vcomponentrad;    // add rad
  double v_ub[6],*vcomponent_ub;    // add rad
  double v_lb[6],*vcomponent_lb;    // add rad
  for (i = 0; i < n; i++) v[i] = 0.0;
  for (i = 0; i < n; i++) vrad[i] = 0.0;    // add rad
  for (i = 0; i < n; i++) v_ub[i] = 0.0;    // add rad
  for (i = 0; i < n; i++) v_lb[i] = 0.0;    // add rad

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++){
      v[i] += vcomponent[i];
    }
  }
  /* -------------------------------------------------------------- */
  // add rad - parallel loop for nvirialrad
  for (j = 0; j < nvirialrad; j++) {
    if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
      vcomponentrad = vptrrad[j];   // add rad
    }
    for (i = 0; i < n; i++){
      if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
        vrad[i] += vcomponentrad[i];    // add rad
      }
    }
  }
  /* ---------------------------------------------------------
      add rad - parallel loop for nvirial_ub and nvirial_lb 
  ---------------------------------------------------------- */
  // add _ub - parallel loop for nvirial_ub
  for (j = 0; j < nvirial_ub; j++) {
    if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
      vcomponent_ub = vptr_ub[j];   // add _ub
    }
    for (i = 0; i < n; i++){
      if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
        v_ub[i] += vcomponent_ub[i];    // add _ub
      }
    }
  }

// add _lb - parallel loop for nvirial_lb
  for (j = 0; j < nvirial_lb; j++) {
    if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
      vcomponent_lb = vptr_lb[j];   // add _lb
    }
    for (i = 0; i < n; i++){
      if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
        v_lb[i] += vcomponent_lb[i];    // add _lb
      }
    }
  }
  /* -------------------------------------------------------------- */

  // sum virial across procs

  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);
  // fix here - work!
  if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag) {
    MPI_Allreduce(vrad,virialrad,n,MPI_DOUBLE,MPI_SUM,world);   // add rad
    MPI_Allreduce(v_ub,virial_ub,n,MPI_DOUBLE,MPI_SUM,world);   // add rad
    MPI_Allreduce(v_lb,virial_lb,n,MPI_DOUBLE,MPI_SUM,world);   // add rad
  }
  // KSpace virial contribution is already summed across procs

  if (kspace_virial)
    for (i = 0; i < n; i++) virial[i] += kspace_virial[i];

  // LJ long-range tail correction

  if (force->pair && force->pair->tail_flag) {
    for (i = 0; i < ndiag; i++) {
      virial[i] += force->pair->ptail * inv_volume;
      if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag)
        virialrad[i] += 0; // add rad
    }
  }

  // begin debug
  // status: the same
  // std::cout << "virial_ub = " << *virial_ub << std::endl; 
  // std::cout << "virial_lb = " << *virial_lb << std::endl; 
  // std::cout << "virial = " << *virial << std::endl; 
  // end debug
}

/* ---------------------------------------------------------------------- */

void ComputePressure::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_temp;
  int n = strlen(id_new) + 1;
  id_temp = new char[n];
  strcpy(id_temp,id_new);
}
