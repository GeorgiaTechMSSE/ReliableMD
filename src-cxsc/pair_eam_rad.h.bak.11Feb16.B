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

#ifdef PAIR_CLASS

PairStyle(eam/rad,PairEAMRad)

#else

#ifndef LMP_PAIR_EAM_RAD_H
#define LMP_PAIR_EAM_RAD_H

#include "stdio.h"
#include "pair.h"
// AT: add CXSC-lib
#include "imatrix.hpp"
#include "ivector.hpp"
#include "interval.hpp"

namespace LAMMPS_NS {

using namespace cxsc; // AT: using namespace cxsc to cut codes short when call interval, ivector, imatrix

class PairEAMRad : public Pair {
 public:
  // add rad flag
  int eamrad_flag;
  int ub_reset_flag; // upperbound reset flag, if 1 then reset
  int lb_reset_flag; // lowerbound reset flag, if 1 then reset

  // public variables so USER-ATC package can access them

  double cutmax;

  // potentials as array data

  int nrho,nr;
  int nfrho,nrhor,nz2r;
  double **frho,**rhor,**z2r;
  double **frho_rad,**rhor_rad,**z2r_rad;
  int *type2frho,**type2rhor,**type2z2r;

  // variables in error generating function
  double zrad_a,zrad_b,rhorad_a,rhorad_b,Frad_a,Frad_b,rho_0;

  // declare ivector/interval variables
  ivector *ivector_rho;
  ivector *ivector_fp;

  // potentials in spline form used for force computation

  double dr,rdr,drho,rdrho,rhomax;
  double ***rhor_spline,***frho_spline,***z2r_spline;

  PairEAMRad(class LAMMPS *);
  virtual ~PairEAMRad();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double single(int, int, int, int, double, double, double, double &);

  virtual int pack_comm(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void swap_eam(double *, double **);

  double type1err_dev(double , double , double );
  double type2err_dev(double , double , double , double );
  void compute_ub(int, int); // compute upper bound of the force
  void compute_lb(int, int); // compute lower bound of the force

  // void compute_kaucher(int, int); // compute force interval based on Kaucher interval arithmetic

 protected:
  int nmax;                   // allocated size of per-atom arrays
  double cutforcesq;

  int signum(double);            // add sign function

  // per-atom arrays

  double *rho,*fp;
  double *rho_ub,*fp_ub; // upperbound
  double *rho_lb,*fp_lb; // lowerbound

  // potentials as file data

  int *map;                   // which element each atom type maps to

  struct Funcfl {
    char *file;
    int nrho,nr;
    double drho,dr,cut,mass;
    double *frho,*rhor,*zr;
    double *frho_rad,*rhor_rad,*zr_rad;
  };
  Funcfl *funcfl;
  int nfuncfl;

  struct Setfl {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,**rhor,***z2r;
    double **frho_rad,**rhor_rad,***z2r_rad;
    double zrad_a,zrad_b,rhorad_a,rhorad_b,Frad_a,Frad_b,rho_0;

  };
  Setfl *setfl;

  struct Fs {
    char **elements;
    int nelements,nrho,nr;
    double drho,dr,cut;
    double *mass;
    double **frho,***rhor,***z2r;
    double **frho_rad,***rhor_rad,***z2r_rad;
  };
  Fs *fs;

  void allocate();
  void array2spline();
  void interpolate(int, double, double *, double **);
  void grab(FILE *, int, double *);

  void read_file(char *);
  void file2array();


};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

*/
