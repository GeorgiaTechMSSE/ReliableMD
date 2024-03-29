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

/* ----------------------------------------------------------------------
   Contributing authors: Stephen Foiles (SNL), Murray Daw (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_eam_rad.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "math_extra.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

PairEAMRad::PairEAMRad(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;

  nmax = 0;
  rho = NULL;
  // rho_ub = NULL; // upperbound
  // rho_lb = NULL; // lowerbound
  fp = NULL;
  // fp_ub = NULL; // upperbound
  // fp_lb = NULL; // lowerbound

  nfuncfl = 0;
  funcfl = NULL;

  setfl = NULL;
  fs = NULL;

  frho = NULL;
  rhor = NULL;
  z2r = NULL;

  frho_rad = NULL;
  rhor_rad = NULL;
  z2r_rad = NULL;

  frho_spline = NULL;
  rhor_spline = NULL;
  z2r_spline = NULL;

  // inherit from pair_eam_alloy.cpp

  one_coeff = 1;
  manybody_flag = 1;

  // set comm size needed by this Pair

  comm_forward = 1;
  comm_reverse = 1;

  // add rad
  atom->frad_flag = 1;
  // ub_reset_flag = 1; // if 1 then reset
  // lb_reset_flag = 1; // if 1 then reset
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairEAMRad::~PairEAMRad()
{
  memory->destroy(rho);
  memory->destroy(fp);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete [] map;
    delete [] type2frho;
    memory->destroy(type2rhor);
    memory->destroy(type2z2r);    
  }

  if (funcfl) {
    for (int i = 0; i < nfuncfl; i++) {
      delete [] funcfl[i].file;
      memory->destroy(funcfl[i].frho);
      memory->destroy(funcfl[i].rhor);
      memory->destroy(funcfl[i].zr);
      if (atom->frad_flag)
      {
        memory->destroy(funcfl[i].frho_rad);
        memory->destroy(funcfl[i].rhor_rad);
        memory->destroy(funcfl[i].zr_rad);
      }
    }
    memory->sfree(funcfl);
  }

  if (setfl) {
    for (int i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    if (atom->frad_flag)
    {
      memory->destroy(setfl->frho_rad);
      memory->destroy(setfl->rhor_rad);
      memory->destroy(setfl->z2r_rad);
    }
    delete setfl;
  }
  
  if (fs) {
    for (int i = 0; i < fs->nelements; i++) delete [] fs->elements[i];
    delete [] fs->elements;
    delete [] fs->mass;
    memory->destroy(fs->frho);
    memory->destroy(fs->rhor);
    memory->destroy(fs->z2r);
    if (atom->frad_flag)
    {
      memory->destroy(fs->frho_rad);
      memory->destroy(fs->rhor_rad);
      memory->destroy(fs->z2r_rad);
    }
    delete fs;
  }

  memory->destroy(frho);
  memory->destroy(rhor);
  memory->destroy(z2r);

  if (atom->frad_flag)
  {
    memory->destroy(frho_rad);
    memory->destroy(rhor_rad);
    memory->destroy(z2r_rad);

  }

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);
}

/* ---------------------------------------------------------------------- */

void PairEAMRad::compute(int eflag, int vflag)
{
  // calculate the upper and lower bounds 
  // ub_reset_flag = 1; 
  // lb_reset_flag = 1; 
  // evaluate force at UPPERBOUND and LOWERBOUND
  compute_ub(eflag,vflag);
  compute_lb(eflag,vflag);  

  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;

  // add rad arguments
  double xtmprad,ytmprad,ztmprad;
  double delxrad,delyrad,delzrad; 
  
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  double **x = atom->x;
  double **xrad = atom->xrad;
  double **f = atom->f;
  double **frad = atom->frad;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    for (i = 0; i < nall; i++){
      rho[i] = 0.0;
      // rho_ub[i] = 0.0; // upperbound
      // rho_lb[i] = 0.0; // lowerbound
    }
  } else for (i = 0; i < nlocal; i++){
      rho[i] = 0.0;
      // rho_ub[i] = 0.0; // upperbound
      // rho_lb[i] = 0.0; // lowerbound
  }

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    // add rad
    xtmprad = xrad[i][0];
    ytmprad = xrad[i][1];
    ztmprad = xrad[i][2];

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      // add  rad
      delxrad = xtmprad - xrad[j][0];
      delyrad = ytmprad - xrad[j][1];
      delzrad = ztmprad - xrad[j][2];

      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        if (newton_pair || j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    if (eflag) {
      phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
      if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      // add rad arguments
      double xtmprad,ytmprad,ztmprad;
      double delxrad,delyrad,delzrad; 
      double fradpair;

      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip;
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip + phip;
        fpair = -psip*recip;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
                
        if (newton_pair || j < nlocal) {
          // action-reaction coupling: apply Newton's law to frad -> create improper force
          f[j][0] -= delx * fpair; 
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        /* ---------- tallying the system ---------- */

        if (eflag) evdwl = phi;
        // if no flag then run this
        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
        
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
  
  // post-force computation
  double **fub = atom->fub;
  double **flb = atom->flb;
  // assign frad based on fub and flb, calculated from ub and lb configuration
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // test
    // std::cout << "f[i][0] = " << f[i][0] << std::endl;
    // std::cout << "f[i][1] = " << f[i][1] << std::endl;
    // std::cout << "f[i][2] = " << f[i][2] << std::endl;
    // std::cout << "fub[i][0] = " << fub[i][0] << std::endl;
    // std::cout << "fub[i][1] = " << fub[i][1] << std::endl;
    // std::cout << "fub[i][2] = " << fub[i][2] << std::endl;
    // std::cout << "flb[i][0] = " << flb[i][0] << std::endl;
    // std::cout << "flb[i][1] = " << flb[i][1] << std::endl;
    // std::cout << "flb[i][2] = " << flb[i][2] << std::endl;

    frad[i][0] = fmin( fub[i][0] - f[i][0], f[i][0] - flb[i][0]);
    frad[i][1] = fmin( fub[i][1] - f[i][1], f[i][1] - flb[i][1]);
    frad[i][2] = fmin( fub[i][2] - f[i][2], f[i][2] - flb[i][2]);

    // std::cout << "frad[i][0] = " << frad[i][0] << std::endl;
    // std::cout << "frad[i][1] = " << frad[i][1] << std::endl;
    // std::cout << "frad[i][2] = " << frad[i][2] << std::endl;

    // for (int j = 0; j < 3; j++) {
    //   std::cout << "flb[" << i << "][" << j << "] = " << flb[i][j] << std::endl;
    //   std::cout << "f[" << i << "][" << j << "] = " << f[i][j] << std::endl;
    //   std::cout << "fub[" << i << "][" << j << "] = " << fub[i][j] << std::endl;
    // }

    // test
    // frad[i][0] = 0.0;
    // frad[i][1] = 0.0;
    // frad[i][2] = 0.0;
  }

  if (atom->xrad_flag && atom->vrad_flag && atom->frad_flag)
  {
    update_virial_rad();
    // std::cout << "update_virial_rad() used" << std::endl;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEAMRad::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  map = new int[n+1];
  for (int i = 1; i <= n; i++) map[i] = -1;

  type2frho = new int[n+1];
  memory->create(type2rhor,n+1,n+1,"pair:type2rhor");
  memory->create(type2z2r,n+1,n+1,"pair:type2z2r");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEAMRad::settings(int narg, char **arg)
{
  if (narg > 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEAMRad::init_style()
{
  // convert read-in file(s) to arrays and spline them

  file2array();
  array2spline();

  neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEAMRad::init_one(int i, int j)
{
  // single global cutoff = max of cut from all files read in
  // for funcfl could be multiple files
  // for setfl or fs, just one file

  if (funcfl) {
    cutmax = 0.0;
    for (int m = 0; m < nfuncfl; m++)
      cutmax = MAX(cutmax,funcfl[m].cut);
  } else if (setfl) cutmax = setfl->cut;
  else if (fs) cutmax = fs->cut;

  cutforcesq = cutmax*cutmax;

  return cutmax;
}

/* ----------------------------------------------------------------------
   read potential values from a DYNAMO single element funcfl file
------------------------------------------------------------------------- */

void PairEAMRad::read_file(char *filename)
{
  // inherits from pair_eam_alloy.cpp

  Setfl *file = setfl;

  // open potential file

  int me = comm->me;
  FILE *fptr;
  char line[MAXLINE];

  if (me == 0) {
    fptr = force->open_potential(filename);
    if (fptr == NULL) {
      char str[128];
      sprintf(str,"Cannot open EAM potential file %s",filename);
      error->one(FLERR,str);
    }
  }

  // read and broadcast header
  // extract element names from nelements line

  int n;
  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    fgets(line,MAXLINE,fptr);
    n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  MPI_Bcast(line,n,MPI_CHAR,0,world);

  sscanf(line,"%d",&file->nelements);
  int nwords = atom->count_words(line);
  if (nwords != file->nelements + 1)
    error->all(FLERR,"Incorrect element names in EAM potential file");

  char **words = new char*[file->nelements+1];
  nwords = 0;
  strtok(line," \t\n\r\f");
  while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

  file->elements = new char*[file->nelements];
  for (int i = 0; i < file->nelements; i++) {
    n = strlen(words[i]) + 1;
    file->elements[i] = new char[n];
    strcpy(file->elements[i],words[i]);
  }
  delete [] words;

  if (me == 0) {
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%d %lg %d %lg %lg",
           &file->nrho,&file->drho,&file->nr,&file->dr,&file->cut);
    std::cout << "file->cut = " << file->cut << std::endl;
  }
 

  MPI_Bcast(&file->nrho,1,MPI_INT,0,world);
  MPI_Bcast(&file->drho,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->nr,1,MPI_INT,0,world);
  MPI_Bcast(&file->dr,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&file->cut,1,MPI_DOUBLE,0,world);

  // read error generating function parameters
  if (me == 0 && atom->frad_flag)
  {
    fgets(line,MAXLINE,fptr);
    sscanf(line,"%lg %lg %lg %lg %lg %lg %lg",
      &file->zrad_a,&file->zrad_b,
      &file->rhorad_a,&file->rhorad_b,
      &file->Frad_a,&file->Frad_b,&file->rho_0);
    std::cout << "file->zrad_a = " << file->zrad_a << std::endl;
    std::cout << "file->zrad_b = " << file->zrad_b << std::endl;
    std::cout << "file->rhorad_a = " << file->rhorad_a << std::endl;
    std::cout << "file->rhorad_b = " << file->rhorad_b << std::endl;
    std::cout << "file->Frad_a = " << file->Frad_a << std::endl;
    std::cout << "file->Frad_b = " << file->Frad_b << std::endl;
    std::cout << "file->rho_0 = " << file->rho_0 << std::endl;
  }

  if (atom->frad_flag)
  {
    MPI_Bcast(&file->zrad_a,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&file->zrad_b,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&file->rhorad_a,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&file->rhorad_b,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&file->Frad_a,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&file->Frad_b,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&file->rho_0,1,MPI_DOUBLE,0,world);
  }    

  file->mass = new double[file->nelements];
  memory->create(file->frho,file->nelements,file->nrho+1,"pair:frho");
  memory->create(file->rhor,file->nelements,file->nr+1,"pair:rhor");
  memory->create(file->z2r,file->nelements,file->nelements,file->nr+1,"pair:z2r");

  if (atom->frad_flag)
  {
    memory->create(file->frho_rad,file->nelements,file->nrho+1,"pair:frho_rad");
    memory->create(file->rhor_rad,file->nelements,file->nr+1,"pair:rhor_rad");
    memory->create(file->z2r_rad,file->nelements,file->nelements,file->nr+1,"pair:z2r_rad");
  }

  int i,j,tmp;
  for (i = 0; i < file->nelements; i++) {
    if (me == 0) {
      fgets(line,MAXLINE,fptr);
      sscanf(line,"%d %lg",&tmp,&file->mass[i]);
    }
    MPI_Bcast(&file->mass[i],1,MPI_DOUBLE,0,world);

    if (me == 0) grab(fptr,file->nrho,&file->frho[i][1]);
    MPI_Bcast(&file->frho[i][1],file->nrho,MPI_DOUBLE,0,world);
    if (me == 0) grab(fptr,file->nr,&file->rhor[i][1]);
    MPI_Bcast(&file->rhor[i][1],file->nr,MPI_DOUBLE,0,world);
  }

  for (i = 0; i < file->nelements; i++)
    for (j = 0; j <= i; j++) {
      if (me == 0) grab(fptr,file->nr,&file->z2r[i][j][1]);
      MPI_Bcast(&file->z2r[i][j][1],file->nr,MPI_DOUBLE,0,world);
    }

  // close the potential file

  if (me == 0) fclose(fptr);
}

/* ----------------------------------------------------------------------
   convert read-in funcfl potential(s) to standard array format
   interpolate all file values to a single grid and cutoff
------------------------------------------------------------------------- */

void PairEAMRad::file2array()
{
  
  // inherits from pair_eam_alloy.cpp
  int i,j,m,n;
  int ntypes = atom->ntypes;

  // set function params directly from setfl file

  nrho = setfl->nrho;
  nr = setfl->nr;
  drho = setfl->drho;
  dr = setfl->dr;
  rhomax = (nrho-1) * drho;

  if (atom->frad_flag)
  {
    zrad_a = setfl->zrad_a;
    zrad_b = setfl->zrad_b;
    rhorad_a = setfl->rhorad_a;
    rhorad_b = setfl->rhorad_b;
    Frad_a = setfl->Frad_a;
    Frad_b = setfl->Frad_b ;
    rho_0 = setfl->rho_0;
  }

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------

  // allocate frho arrays
  // nfrho = # of setfl elements + 1 for zero array

  nfrho = setfl->nelements + 1;
  memory->destroy(frho);
  memory->create(frho,nfrho,nrho+1,"pair:frho");

  // copy each element's frho to global frho

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nrho; m++) {
      frho[i][m] = setfl->frho[i][m];

      // error here
      // if (atom->frad_flag) frho_rad[i][m] = setfl->frho_rad[i][m];
    }

  // add extra frho of zeroes for non-EAM types to point to (pair hybrid)
  // this is necessary b/c fp is still computed for non-EAM atoms

  for (m = 1; m <= nrho; m++) frho[nfrho-1][m] = 0.0;

  // type2frho[i] = which frho array (0 to nfrho-1) each atom type maps to
  // if atom type doesn't point to element (non-EAM atom in pair hybrid)
  // then map it to last frho array of zeroes

  for (i = 1; i <= ntypes; i++)
    if (map[i] >= 0) type2frho[i] = map[i];
    else type2frho[i] = nfrho-1;

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------

  // allocate rhor arrays
  // nrhor = # of setfl elements

  nrhor = setfl->nelements;
  memory->destroy(rhor);
  memory->create(rhor,nrhor,nr+1,"pair:rhor");

  // copy each element's rhor to global rhor

  for (i = 0; i < setfl->nelements; i++)
    for (m = 1; m <= nr; m++) rhor[i][m] = setfl->rhor[i][m];

  // type2rhor[i][j] = which rhor array (0 to nrhor-1) each type pair maps to
  // for setfl files, I,J mapping only depends on I
  // OK if map = -1 (non-EAM atom in pair hybrid) b/c type2rhor not used

  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      type2rhor[i][j] = map[i];

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------

  // allocate z2r arrays
  // nz2r = N*(N+1)/2 where N = # of setfl elements

  nz2r = setfl->nelements * (setfl->nelements+1) / 2;
  memory->destroy(z2r);
  memory->create(z2r,nz2r,nr+1,"pair:z2r");

  // copy each element pair z2r to global z2r, only for I >= J

  n = 0;
  for (i = 0; i < setfl->nelements; i++)
    for (j = 0; j <= i; j++) {
      for (m = 1; m <= nr; m++) z2r[n][m] = setfl->z2r[i][j][m];
      n++;
    }

  // type2z2r[i][j] = which z2r array (0 to nz2r-1) each type pair maps to
  // set of z2r arrays only fill lower triangular Nelement matrix
  // value = n = sum over rows of lower-triangular matrix until reach irow,icol
  // swap indices when irow < icol to stay lower triangular
  // if map = -1 (non-EAM atom in pair hybrid):
  //   type2z2r is not used by non-opt
  //   but set type2z2r to 0 since accessed by opt

  int irow,icol;
  for (i = 1; i <= ntypes; i++) {
    for (j = 1; j <= ntypes; j++) {
      irow = map[i];
      icol = map[j];
      if (irow == -1 || icol == -1) {
        type2z2r[i][j] = 0;
        continue;
      }
      if (irow < icol) {
        irow = map[j];
        icol = map[i];
      }
      n = 0;
      for (m = 0; m < irow; m++) n += m + 1;
      n += icol;
      type2z2r[i][j] = n;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairEAMRad::array2spline()
{
  rdr = 1.0/dr;
  rdrho = 1.0/drho;

  memory->destroy(frho_spline);
  memory->destroy(rhor_spline);
  memory->destroy(z2r_spline);

  memory->create(frho_spline,nfrho,nrho+1,7,"pair:frho");
  memory->create(rhor_spline,nrhor,nr+1,7,"pair:rhor");
  memory->create(z2r_spline,nz2r,nr+1,7,"pair:z2r");

  for (int i = 0; i < nfrho; i++)
    interpolate(nrho,drho,frho[i],frho_spline[i]);

  for (int i = 0; i < nrhor; i++)
    interpolate(nr,dr,rhor[i],rhor_spline[i]);

  for (int i = 0; i < nz2r; i++)
    interpolate(nr,dr,z2r[i],z2r_spline[i]);
}

/* ---------------------------------------------------------------------- */

void PairEAMRad::interpolate(int n, double delta, double *f, double **spline)
{
  for (int m = 1; m <= n; m++) spline[m][6] = f[m];

  spline[1][5] = spline[2][6] - spline[1][6];
  spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
  spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
  spline[n][5] = spline[n][6] - spline[n-1][6];

  for (int m = 3; m <= n-2; m++)
    spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                    8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

  for (int m = 1; m <= n-1; m++) {
    spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
      2.0*spline[m][5] - spline[m+1][5];
    spline[m][3] = spline[m][5] + spline[m+1][5] -
      2.0*(spline[m+1][6]-spline[m][6]);
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;

  for (int m = 1; m <= n; m++) {
    spline[m][2] = spline[m][5]/delta;
    spline[m][1] = 2.0*spline[m][4]/delta;
    spline[m][0] = 3.0*spline[m][3]/delta;
  }
}

/* ----------------------------------------------------------------------
   grab n values from file fp and put them in list
   values can be several to a line
   only called by proc 0
------------------------------------------------------------------------- */

void PairEAMRad::grab(FILE *fptr, int n, double *list)
{
  char *ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n) {
    fgets(line,MAXLINE,fptr);
    ptr = strtok(line," \t\n\r\f");
    list[i++] = atof(ptr);
    while ((ptr = strtok(NULL," \t\n\r\f"))) list[i++] = atof(ptr);
  }
}

/* ---------------------------------------------------------------------- */

double PairEAMRad::single(int i, int j, int itype, int jtype,
                       double rsq, double factor_coul, double factor_lj,
                       double &fforce)
{
  int m;
  double r,p,rhoip,rhojp,z2,z2p,recip,phi,phip,psip;
  double *coeff;

  r = sqrt(rsq);
  p = r*rdr + 1.0;
  m = static_cast<int> (p);
  m = MIN(m,nr-1);
  p -= m;
  p = MIN(p,1.0);

  coeff = rhor_spline[type2rhor[itype][jtype]][m];
  rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = rhor_spline[type2rhor[jtype][itype]][m];
  rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
  coeff = z2r_spline[type2z2r[itype][jtype]][m];
  z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
  z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

  recip = 1.0/r;
  phi = z2*recip;
  phip = z2p*recip - phi*recip;
  psip = fp[i]*rhojp + fp[j]*rhoip + phip;
  fforce = -psip*recip;

  return phi;
}

/* ---------------------------------------------------------------------- */

int PairEAMRad::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = fp[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAMRad::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) fp[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int PairEAMRad::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = rho[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void PairEAMRad::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    rho[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairEAMRad::memory_usage()
{
  double bytes = maxeatom * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);
  bytes += 2 * nmax * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   swap fp array with one passed in by caller
------------------------------------------------------------------------- */

void PairEAMRad::swap_eam(double *fp_caller, double **fp_caller_hold)
{
  double *tmp = fp;
  fp = fp_caller;
  *fp_caller_hold = tmp;
}


/* ----------------------------------------------------------------------
   INHERITS FROM pair_eam_alloy.cpp
------------------------------------------------------------------------- */

void PairEAMRad::coeff(int narg, char **arg)
{
  int i,j;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read EAM setfl file

  if (setfl) {
    for (i = 0; i < setfl->nelements; i++) delete [] setfl->elements[i];
    delete [] setfl->elements;
    delete [] setfl->mass;
    memory->destroy(setfl->frho);
    memory->destroy(setfl->rhor);
    memory->destroy(setfl->z2r);
    if (atom->frad_flag){
      memory->destroy(setfl->frho_rad);
      memory->destroy(setfl->rhor_rad);
      memory->destroy(setfl->z2r_rad);
    }
    delete setfl;
  }
  setfl = new Setfl();
  read_file(arg[2]);

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL

  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < setfl->nelements; j++)
      if (strcmp(arg[i],setfl->elements[j]) == 0) break;
    if (j < setfl->nelements) map[i-2] = j;
    else error->all(FLERR,"No matching element in EAM potential file");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements
  // set mass of atom type if i = j

  int count = 0;
  for (i = 1; i <= n; i++) {
    for (j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j) atom->set_mass(i,setfl->mass[map[i]]);
        count++;
      }
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   type1err_dev calculate derivative of generating error function
   given 2 parameters a & b
------------------------------------------------------------------------- */

double PairEAMRad::type1err_dev(double a1, double b1, double r1)
{
  // test
  // std::cout << "type1err_dev() = " << (-a1 * b1 * exp(-b1*r1) ) << std::endl;
  return (-a1 * b1 * exp(-b1*r1));
  // return (abs(-a1 * b1 * exp(-b1*r1)) );
}

double PairEAMRad::type2err_dev(double a2, double b2, double rho_0, double rho2)
{
  // test
  // if (a2 == 0) std::cout << "type2err_dev() = " << 0 << std::endl;
  // else std::cout << "type2err_dev() = " 
  //               << (a2 * b2 * pow(rho2/rho_0 , b2*rho2 -1) * exp(-b2 * (rho2-rho_0)) * (1-rho2/rho_0) )
  //               << std::endl;

  // begin debug
  // std::cout << "type2err_dev::a2 = " << a2 << std::endl;
  // std::cout << "type2err_dev::b2 = " << b2 << std::endl;
  // std::cout << "type2err_dev::rho_0 = " << rho_0 << std::endl;
  // std::cout << "type2err_dev::rho2 = " << rho2 << std::endl;
  // end debug

  if (rho2 < 0) std::cout << "check type2err_dev() arg rho2" << std::endl;
  if (a2 == 0) return 0;
  // return (a2 * b2 * pow(rho2/rho_0 , b2*rho_0 - 1) * exp(-b2 * (rho2-rho_0)) * abs(1-rho2/rho_0) );
  return (a2 * b2 * pow(rho2/rho_0 , b2*rho_0 -1) * exp(-b2 * (rho2-rho_0)) * (1-rho2/rho_0) );
}


/* ----------------------------------------------------------------------
   ADD UPPERBOUND COMPUTATION
------------------------------------------------------------------------- */

void PairEAMRad::compute_ub(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  i = j = ii = jj = m = inum = jnum = itype = jtype = 0;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup_ub(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  double **x = atom->x;
  double **xub = atom->xub;
  double **xrad = atom->xrad;
  double **fub = atom->fub;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  // std::cout << "nlocal = " << nlocal << std::endl;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    for (i = 0; i < nall; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // xtmp = x[i][0] + xrad[i][0]; // upperbound computation
    // ytmp = x[i][1] + xrad[i][1]; // upperbound computation
    // ztmp = x[i][2] + xrad[i][2]; // upperbound computation
    xtmp = x[i][0]; // upperbound computation
    ytmp = x[i][1]; // upperbound computation
    ztmp = x[i][2]; // upperbound computation

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // delx = xtmp - (x[j][0] + xrad[j][0]); // upperbound computation
      // dely = ytmp - (x[j][1] + xrad[j][1]); // upperbound computation
      // delz = ztmp - (x[j][2] + xrad[j][2]); // upperbound computation

      delx = xtmp - x[j][0]; // upperbound computation
      dely = ytmp - x[j][1]; // upperbound computation
      delz = ztmp - x[j][2]; // upperbound computation

      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        if (newton_pair || j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    // fp[i] += signum(fp[i]) * type2err_dev(Frad_a, Frad_b, rho_0, rho[i]); // upperbound computation 
    fp[i] += type2err_dev(Frad_a, Frad_b, rho_0, rho[i]); // upperbound computation
    // if (eflag) {
    //   phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    //   if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
    //   if (eflag_global) eng_vdwl += phi;
    //   if (eflag_atom) eatom[i] += phi;
    // }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0]; // upperbound computation
    ytmp = x[i][1]; // upperbound computation
    ztmp = x[i][2]; // upperbound computation
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0]; // upperbound computation
      dely = ytmp - x[j][1]; // upperbound computation
      delz = ztmp - x[j][2]; // upperbound computation
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        // rhoip += signum(rhoip) * type1err_dev(rhorad_a, rhorad_b, r); // upperbound computation
        rhoip += type1err_dev(rhorad_a, rhorad_b, r); // upperbound computation
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        // rhojp += signum(rhojp) * type1err_dev(rhorad_a, rhorad_b, r); // upperbound computation
        rhojp += type1err_dev(rhorad_a, rhorad_b, r); // upperbound computation
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        // z2p += signum(z2p) * type1err_dev(zrad_a, zrad_b, r); // upperbound computation
        z2p += type1err_dev(zrad_a, zrad_b, r); // upperbound computation
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip;
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip + phip;
        fpair = -psip*recip;

        fub[i][0] += delx*fpair; // upperbound computation
        fub[i][1] += dely*fpair; // upperbound computation
        fub[i][2] += delz*fpair; // upperbound computation
                
        if (newton_pair || j < nlocal) {
          // action-reaction coupling: apply Newton's law to frad -> create improper force
          fub[j][0] -= delx * fpair; // upperbound computation
          fub[j][1] -= dely * fpair; // upperbound computation
          fub[j][2] -= delz * fpair; // upperbound computation
        }

        if (eflag) evdwl = phi;
        // if no flag then run this
        if (evflag) ev_tally_ub(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  } 
  // int vflag_fdotr_ub = vflag_fdotr;
  if (vflag_fdotr_ub) virial_ub_fdotr_compute();
  ub_reset_flag = 0;
}

/* ----------------------------------------------------------------------
   ADD LOWERBOUND FORCE COMPUTATION
------------------------------------------------------------------------- */

void PairEAMRad::compute_lb(int eflag, int vflag)
{
  int i,j,ii,jj,m,inum,jnum,itype,jtype;
  i = j = ii = jj = m = inum = jnum = itype = jtype = 0;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
  double *coeff;
  
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup_lb(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = eflag_atom = 0;

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    nmax = atom->nmax;
    memory->create(rho,nmax,"pair:rho");
    memory->create(fp,nmax,"pair:fp");
  }

  double **x = atom->x;
  double **xlb = atom->xlb;
  double **xrad = atom->xrad;
  double **flb = atom->flb;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    for (i = 0; i < nall; i++) rho[i] = 0.0;
  } else for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    // xtmp = x[i][0] - xrad[i][0]; // lowerbound computation
    // ytmp = x[i][1] - xrad[i][1]; // lowerbound computation
    // ztmp = x[i][2] - xrad[i][2]; // lowerbound computation
    xtmp = x[i][0]; // lowerbound computation
    ytmp = x[i][1]; // lowerbound computation
    ztmp = x[i][2]; // lowerbound computation

    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // delx = xtmp - (x[j][0] - xrad[j][0]); // lowerbound computation
      // dely = ytmp - (x[j][1] - xrad[j][1]); // lowerbound computation
      // delz = ztmp - (x[j][2] - xrad[j][2]); // lowerbound computation
      delx = xtmp - x[j][0]; // lowerbound computation
      dely = ytmp - x[j][1]; // lowerbound computation
      delz = ztmp - x[j][2]; // lowerbound computation

      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq)*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        if (newton_pair || j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm_pair(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax (e.g. due to close approach of two atoms),
  //   will exceed table, so add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = rho[i]*rdrho + 1.0;
    m = static_cast<int> (p);
    m = MAX(1,MIN(m,nrho-1));
    p -= m;
    p = MIN(p,1.0);
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0]*p + coeff[1])*p + coeff[2];
    // fp[i] -= signum(fp[i]) * type2err_dev(Frad_a, Frad_b, rho_0, rho[i]); // lowerbound computation
    fp[i] -= type2err_dev(Frad_a, Frad_b, rho_0, rho[i]); // lowerbound computation
    // if (eflag) {
    //   phi = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];
    //   if (rho[i] > rhomax) phi += fp[i] * (rho[i]-rhomax);
    //   if (eflag_global) eng_vdwl += phi;
    //   if (eflag_atom) eatom[i] += phi;
    // }
  }

  // communicate derivative of embedding function

  comm->forward_comm_pair(this);

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        r = sqrt(rsq);
        p = r*rdr + 1.0;
        m = static_cast<int> (p);
        m = MIN(m,nr-1);
        p -= m;
        p = MIN(p,1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0]*p + coeff[1])*p + coeff[2];
        // rhoip -= signum(rhoip) * type1err_dev(rhorad_a, rhorad_b, r); // lowerbound computation
        rhoip -= type1err_dev(rhorad_a, rhorad_b, r); // lowerbound computation
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0]*p + coeff[1])*p + coeff[2];
        // rhojp -= signum(rhojp) * type1err_dev(rhorad_a, rhorad_b, r); // lowerbound computation
        rhojp -= type1err_dev(rhorad_a, rhorad_b, r); // lowerbound computation
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0]*p + coeff[1])*p + coeff[2];
        // z2p -= signum(z2p) * type1err_dev(zrad_a, zrad_b, r); // lowerbound computation
        z2p -= type1err_dev(zrad_a, zrad_b, r); // lowerbound computation
        z2 = ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6];

        recip = 1.0/r;
        phi = z2*recip;
        phip = z2p*recip - phi*recip;
        psip = fp[i]*rhojp + fp[j]*rhoip + phip;
        fpair = -psip*recip;

        flb[i][0] += delx*fpair; // lowerbound computation
        flb[i][1] += dely*fpair; // lowerbound computation
        flb[i][2] += delz*fpair; // lowerbound computation
                
        if (newton_pair || j < nlocal) {
          // action-reaction coupling: apply Newton's law to frad -> create improper force
          flb[j][0] -= delx * fpair; // lowerbound computation
          flb[j][1] -= dely * fpair; // lowerbound computation
          flb[j][2] -= delz * fpair; // lowerbound computation
        }
        
        // // begin debug
        // std::cout << "flb[" << i << "][0] = " << flb[i][0] << std::endl;
        // std::cout << "flb[" << i << "][1] = " << flb[i][1] << std::endl;
        // std::cout << "flb[" << i << "][2] = " << flb[i][2] << std::endl;
        // // end debug

        if (eflag) evdwl = phi;
        // if no flag then run this
        if (evflag) ev_tally_lb(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  } 
  
  if (vflag_fdotr_lb) virial_lb_fdotr_compute();
  lb_reset_flag = 0;
}

/* ----------------------------------------------------------------------
   signum function
------------------------------------------------------------------------- */
int PairEAMRad::signum(double value)
{
  int sign;
  if (value>0) sign = 1;
  if (value<0) sign = -1;
  if (value==0) sign = 0;
  return sign;
}