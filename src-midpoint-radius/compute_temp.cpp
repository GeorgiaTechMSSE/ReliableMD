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
#include "compute_temp.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "interval.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTemp::ComputeTemp(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute temp command");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[6];
  vectorrad = new double[6]; // add rad
}

/* ---------------------------------------------------------------------- */

ComputeTemp::~ComputeTemp()
{
  delete [] vector;
  delete [] vectorrad; // add rad
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::setup()
{
  fix_dof = -1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::dof_compute()
{
  if (fix_dof) adjust_dof_fix();
  double natoms = group->count(igroup);
  dof = domain->dimension * natoms;
  dof -= extra_dof + fix_dof;
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTemp::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double t = 0.0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * rmass[i];
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) *
          mass[type[i]];
  }

  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  scalar *= tfactor;
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTemp::compute_vector()
{
  int i;
  // std::cout << "used here flag 3" << std::endl;
  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double **vrad = atom->vrad;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone,t[6];
  double trad[6]; // add rad
  for (i = 0; i < 6; i++) {
    t[i] = 0.0;
    trad[i] = 0.0; // init double trad[6] along with double t[6]
  }

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
      /* ---------------------------------------------------------------------------------------------- */
      // add rad
      // NOTE: omit the vrad*vrad part when calculating the 
      // new interval rad because that part (vrad*vrad) belongs
      // to interval mid rather than interval rad
      if (atom->vrad_flag) {
        /* ---------------------------------------------------------------------------------------------- */
        // version 1 - mixed classical/Kaucher arithmetics
        /*trad[0] += massone * ( 2 * v[i][0]*vrad[i][0] + vrad[i][0]*vrad[i][0] );
        trad[1] += massone * ( 2 * v[i][1]*vrad[i][1] + vrad[i][1]*vrad[i][1] );
        trad[2] += massone * ( 2 * v[i][2]*vrad[i][2] + vrad[i][2]*vrad[i][2] );
        trad[3] += massone * ( 2 * v[i][0]*vrad[i][1] + vrad[i][0]*vrad[i][1] );
        trad[4] += massone * ( 2 * v[i][0]*vrad[i][2] + vrad[i][0]*vrad[i][2] );
        trad[5] += massone * ( 2 * v[i][1]*vrad[i][2] + vrad[i][1]*vrad[i][2] );*/
        /* ---------------------------------------------------------------------------------------------- */
        // version 2 - Kaucher arithmetics
        // trad[0] += massone * MathExtra::rad_DoubleMult(v[i][0],vrad[i][0],v[i][0],vrad[i][0]);
        // trad[1] += massone * MathExtra::rad_DoubleMult(v[i][1],vrad[i][1],v[i][1],vrad[i][1]);
        // trad[2] += massone * MathExtra::rad_DoubleMult(v[i][2],vrad[i][2],v[i][2],vrad[i][2]);
        // trad[3] += massone * MathExtra::rad_DoubleMult(v[i][0],vrad[i][0],v[i][1],vrad[i][1]);
        // trad[4] += massone * MathExtra::rad_DoubleMult(v[i][0],vrad[i][0],v[i][2],vrad[i][2]);
        // trad[5] += massone * MathExtra::rad_DoubleMult(v[i][1],vrad[i][1],v[i][2],vrad[i][2]);

        // std::cout << "v[" << i << "][0] = " << v[i][0] << std::endl;
        // std::cout << "v[" << i << "][1] = " << v[i][1] << std::endl;
        // std::cout << "v[" << i << "][2] = " << v[i][2] << std::endl;
        // std::cout << "vrad[" << i << "][0] = " << vrad[i][0] << std::endl;
        // std::cout << "vrad[" << i << "][1] = " << vrad[i][1] << std::endl;
        // std::cout << "vrad[" << i << "][2] = " << vrad[i][2] << std::endl;
        
        // std::cout << "MathExtra::rad_DoubleMult(v["<< i << "][0],vrad["<< i << "][0],v["<< i << "][2],vrad["<< i << "][2]) = " << MathExtra::rad_DoubleMult(v[i][0],vrad[i][0],v[i][2],vrad[i][2]) << std::endl;

        // for (int k = 0; k<6;k++) std::cout << "trad[" << k << "] = " << trad[k] << std::endl;
               
        /* ---------------------------------------------------------------------------------------------- */
        // version 3 - assuming in NVT temperature unchanged
        trad[0] += 0.0;
        trad[1] += 0.0;
        trad[2] += 0.0;
        trad[3] += 0.0;
        trad[4] += 0.0;
        trad[5] += 0.0;
      }
      /* ---------------------------------------------------------------------------------------------- */
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  // add rad
  if (atom->vrad_flag) MPI_Allreduce(trad,vectorrad,6,MPI_DOUBLE,MPI_SUM,world);
  
  for (i = 0; i < 6; i++) {
    vector[i] *= force->mvv2e;
    if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag)
      vectorrad[i] *= force->mvv2e; // add rad
  }
}
