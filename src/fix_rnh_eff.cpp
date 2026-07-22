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


#include "fix_rnh_eff.h"

#include "atom.h"
#include "error.h"
#include "domain.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixRNHEff::FixRNHEff(LAMMPS *lmp, int narg, char **arg) : FixNH(lmp, narg, arg)
{
  if (!atom->electron_flag)
    error->all(FLERR,"Fix nvt/nph/npt/eff requires atom style electron");
}

/* ----------------------------------------------------------------------
   perform half-step update of electron radial velocities
-----------------------------------------------------------------------*/

void FixRNHEff::nve_v()
{
  // standard nve_v velocity update

  //FixNH::nve_v();

  double *erforce = atom->erforce;
  double *ervel = atom->ervel;
  double *mass = atom->mass;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double dtfm;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (abs(spin[i])==1 || (abs(spin[i])>=4 && abs(spin[i])<=7)) {
        dtfm = dtf / mass[type[i]];
		//utils::logmesg(lmp,"mass[type[i]] = {:f}\n",mass[type[i]]);
        ervel[i] += dtfm * erforce[i];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of electron radii
-----------------------------------------------------------------------*/

void FixRNHEff::nve_x()
{
  // standard nve_x position update

  // FixNH::nve_x();

  double *eradius = atom->eradius;
  double *ervel = atom->ervel;
  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (abs(spin[i])>=1 && abs(spin[i])<=7) eradius[i] += dtv * ervel[i];
}

/* ----------------------------------------------------------------------
   perform half-step scaling of electron radial velocities
-----------------------------------------------------------------------*/

void FixRNHEff::nh_v_temp()
{
  // standard nh_v_temp velocity scaling

  //FixNH::nh_v_temp();

  double *ervel = atom->ervel;
  int *spin = atom->spin;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  //utils::logmesg(lmp,"factor_eta_rad = {:f}\n",factor_eta);
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (abs(spin[i])>=1) {
		  ervel[i] *= factor_eta;
		  if(ervel[i] > 0.05)
			  ervel[i] = 0.05;
		  if(ervel[i] < -0.05)
			  ervel[i] = -0.05;
	  }
}
