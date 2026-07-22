// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include "compute_temp_eff.h"


#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "universe.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeTempEff::ComputeTempEff(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (!atom->electron_flag)
    error->all(FLERR,"Compute temp/eff requires atom style electron");

  scalar_flag = vector_flag = 1;
  size_vector = 6;
  extscalar = 0;
  extvector = 1;
  tempflag = 1;

  vector = new double[size_vector];
}

/* ---------------------------------------------------------------------- */

ComputeTempEff::~ComputeTempEff()
{
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeTempEff::setup()
{
  dynamic = 0;
  if (dynamic_user || group->dynamic[igroup]) dynamic = 1;
  dof_compute();
}

/* ---------------------------------------------------------------------- */
//计算自由度
void ComputeTempEff::dof_compute()
{
  adjust_dof_fix();
  natoms_temp = group->count(igroup);
  dof = domain->dimension * natoms_temp;
  dof -= extra_dof + fix_dof;
  //fprintf(universe->uscreen,"eff_dof_before:%d %f\n",universe->me, dof);
  //utils::logmesg(lmp,"effdof_before = {:.8}\n",dof);
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (abs(spin[i])==1 || type[i]==1) one++;
    }
  //fprintf(universe->uscreen,"ONE:%d %d\n",universe->me, one);
  int nelectrons;
  MPI_Allreduce(&one,&nelectrons,1,MPI_INT,MPI_SUM,world);
  // Assume 3/2 k T per nucleus

  //dof -= domain->dimension * nelectrons;
  //fprintf(universe->uscreen,"eff_dof_after:%d %f\n",universe->me, dof);
  //utils::logmesg(lmp,"dof_eff, extra_dof, fix_dof = {:d}, {:d}, {:d}\n",dof, extra_dof, fix_dof);
  utils::logmesg(lmp,"dof_eff = {:f}\n",dof);
  if (dof > 0.0) tfactor = force->mvv2e / (dof * force->boltz);
  else tfactor = 0.0;
}

/* ---------------------------------------------------------------------- */

double ComputeTempEff::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;

  double t = 0.0;

  if (mass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
		  t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[type[i]];//计算所有电子的平动动能
        //if (abs(spin[i])==1 || abs(spin[i])==2) t += mefactor*mass[type[i]]*ervel[i]*ervel[i];
      }
    }
  }
//MPI_SUM将所有进程的t相加并返回scalar。
  MPI_Allreduce(&t,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  if (dynamic) dof_compute();
  //fprintf(universe->uscreen,"dof:%d %d\n",universe->me, dof);
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR,"Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  //mv**2/Nk
  //double dof_calc = force->mvv2e / (tfactor * force->boltz);
  //utils::logmesg(lmp,"scalar,dofcalc,dof = {:.8},{:.8},{:.8}\n",scalar,dof_calc,dof);
  //fprintf(universe->uscreen,"eff_scalar:%d %f\n",universe->me, scalar);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeTempEff::compute_vector()
{
  int i;

  invoked_vector = update->ntimestep;

  double **v = atom->v;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;

  double massone,t[6];
  for (i = 0; i < 6; i++) t[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      massone = mass[type[i]];
      t[0] += massone * v[i][0]*v[i][0];
      t[1] += massone * v[i][1]*v[i][1];
      t[2] += massone * v[i][2]*v[i][2];
      t[3] += massone * v[i][0]*v[i][1];
      t[4] += massone * v[i][0]*v[i][2];
      t[5] += massone * v[i][1]*v[i][2];
      if (abs(spin[i])>=1 && abs(spin[i])<=5) {
        t[0] += mefactor*massone*ervel[i]*ervel[i];
        t[1] += mefactor*massone*ervel[i]*ervel[i];
        t[2] += mefactor*massone*ervel[i]*ervel[i];
      }
    }

  MPI_Allreduce(t,vector,6,MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}
