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
   Contributing author: Andres Jaramillo-Botero
------------------------------------------------------------------------- */

#include "pair_eff_cut.h"
#include "pair_eff_inline.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "min.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "universe.h"
#include "mpi.h"
#include <chrono>

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEffCut::PairEffCut(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  nmax = 0;
  min_eradius = nullptr;
  min_erforce = nullptr;
  nextra = 5;//5个能量 + 4个时间
  pvector = new double[nextra];
  //time1 = time2 = time3 = time4 = 0;
}

/* ---------------------------------------------------------------------- */

PairEffCut::~PairEffCut()
{
  delete [] pvector;
  memory->destroy(min_eradius);
  memory->destroy(min_erforce);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairEffCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,energy;
  double eke,ecoul,epauli,eadd,errestrain,halfcoul,halfpauli;
  double fpair,fx,fy,fz,fadd;
  //用于存储由于环境因子导致的额外受力
  double f1x,f1y,f1z,fh;
  double e1rforce,e2rforce,e1rvirial,e2rvirial;
  double s_fpair, s_e1rforce, s_e2rforce;
  double rsq, rc, r1, r2, cos, p, h;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  int *molecule = atom->molecule;
  //utils::logmesg(lmp,"molecule[0], x[0][0], x[0][1], x[0][2] = {:d}, {:f}, {:f}, {:f}\n",molecule[0], x[0][0], x[0][1], x[0][2]);
  

  energy = eke = epauli = ecoul = errestrain = eadd = 0.0;
  // pvector = [KE, Pauli, ecoul, radial_restraint]
  for (i=0; i<5; i++) pvector[i] = 0.0;

  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *erforce = atom->erforce;
  double *eradius = atom->eradius;
  int *tag = atom->tag;
  int *spin = atom->spin;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  
  
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int timestep = update->ntimestep;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];//包含了所有的nlocal粒子

	//utils::logmesg(lmp,"i = {:d}\n",ilist[ii]);
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
	
    if (abs(spin[i])==1) {
      eke = ecoul = 0.0;
      fpair = e1rforce = e2rforce = 0.0;	
	  if(eradius[i]<0)
		  eradius[i] = -eradius[i];

	  //KinElec(eradius[i],&eke,&e1rforce,2*ke_e1);
	  KinElec(eradius[i],&eke,&e1rforce,2*ke_e1,A_e1,B_e1,C_e1,D_e1);
	  erforce[i] += e1rforce * hhmss2e;//注意单位转换
	  e1rforce = 0.0;
	  Cee_0(eradius[i],&ecoul,&e1rforce,ah_e11,a1_e11,kvv_e11);
      erforce[i] += e1rforce * qqrd2e;
	  // apply unit conversion factors
      eke *= hhmss2e;
      ecoul *= qqrd2e;

      // Sum up contributions
      energy = eke + ecoul;

	  //utils::logmesg(lmp,"evflag,pressure_with_evirials_flag = {:d},{:d}\n",evflag,pressure_with_evirials_flag);
      // Tally energy and compute radial atomic virial contribution
      if (evflag) {
        ev_tally_eff(i,i,nlocal,newton_pair,energy,0.0);
        if (pressure_with_evirials_flag) // iff flexible pressure flag on
          ev_tally_eff(i,i,nlocal,newton_pair,0.0,e1rforce*eradius[i]);
      }
      if (eflag_global) {
        pvector[0] += eke;
        pvector[2] += ecoul;
      }
    }

    for (jj = 0; jj < jnum; jj++) {
	  //utils::logmesg(lmp,"inum, jnum = {:d}, {:d}\n",inum, jnum);
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      rc = sqrt(rsq);

      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        energy = ecoul = epauli = eadd = 0.0;
        fx = fy = fz = fpair = s_fpair = fadd = 0.0;
		f1x = f1y = f1z = fh = 0.0;
        double taper = sqrt(cutsq[itype][jtype]);
        double dist = rc / taper;
        //double spline = cutoff(dist);
        double spline = cutoff(dist);
        double dspline = dcutoff(dist) / taper;

        // nucleus (i) - nucleus (j) Coul interaction
		
        if (q[i] != 0 && q[j] != 0) {
		  e1rforce = e2rforce = s_e1rforce = s_e2rforce = 0.0;
          double qxq = abs(q[i])*abs(q[j]);
		  //utils::logmesg(lmp,"d = {:.8}\n",rc);
          //ElecNucNuc(qxq, rc, &ecoul, &fpair);
		  //ecoul += qxq / rc; 
 		  //fpair += qxq / (rc * rc);
 		  if (force->qqr2e==332.06371)
 		  	{
 		  	if (qxq == 1)
 		  	{
 		  	//ElecNucNuc(qxq, rc, &ecoul, &fpair, khh, thh, phh);
 		  	Unn(rc,ecoul,fpair,kappa_hh,tau_hh,rho_hh,qhh);
			harmonic(rc, 0.4, 200.0, &eadd, &fadd);
			}
		  	else if (q[i] == 8 && q[j] == 1)
		  	{
		  	//ElecNucNuc(qxq, rc, &ecoul, &fpair, kch, tch, pch);
		  	Unn(rc,ecoul,fpair,kappa_oh,tau_oh,rho_oh,qoh);
			harmonic(rc, 0.6, 200.0, &eadd, &fadd);
			}
			else if (q[i] == 1 && q[j] == 8)
		  	{
		  	//ElecNucNuc(qxq, rc, &ecoul, &fpair, kch, tch, pch);
		  	Unn(rc,ecoul,fpair,kappa_oh,tau_oh,rho_oh,qoh);
			harmonic(rc, 0.6, 200.0, &eadd, &fadd);
			}
		  	else if (qxq == 64)
		  	{
		  	//ElecNucNuc(qxq, rc, &ecoul, &fpair, kcc, tcc, pcc);
		  	Uoo(rc,ecoul,fpair,kappa_oo,tau_oo,rho_oo,qoo,rc_wall,A_gauss,B_gauss);
			harmonic(rc, 1.8, 200.0, &eadd, &fadd);
			}
			SmallRForce(delx,dely,delz,rc,fadd,&fx,&fy,&fz);//??
			f[i][0] += fx;
			f[i][1] += fy;
			f[i][2] += fz;
			if (newton_pair || j < nlocal) {
			  f[j][0] -= fx;
			  f[j][1] -= fy;
			  f[j][2] -= fz;
			  }
			pvector[4] += eadd;
 		  	}
		  else if(force->qqr2e==1)
		  	{
		  	ecoul += qxq / rc; 
 		    fpair += qxq / (rc * rc);
		  	}
        }
		
        // nucleus (i) - electron (j) Coul interaction
		// 注意指数项最好比1大一些，否则会导致导数的突变
		//核与外层电子
        else if  (q[i] != 0 && abs(spin[j]) == 1) {
          e1rforce = e2rforce = s_e1rforce = s_e2rforce = 0.0;
		  if(abs(q[i])==1){
			  Ueh(2.0,qhe,rc,eradius[j],&ecoul,&fpair,&e2rforce,ah_e1,a1_e1,a2_e1,a3_e1,a4_e1,a22_e1,a33_e1,a44_e1);
		  }
		  else{
		  	Uec(2.0,qoe,rc,eradius[j],&ecoul,&fpair,&e2rforce,ac_e1,ac1_e1,ac2_e1,ac3_e1,ac4_e1,poe1,ac22_e1,ac33_e1,ac44_e1);
			/*
			harmonic(rc, 0.4, 10000, &eadd, &fadd);//避免外层电子和O核重叠
			SmallRForce(delx,dely,delz,rc,fadd,&fx,&fy,&fz);//??
			f[i][0] += fx;
			f[i][1] += fy;
			f[i][2] += fz;
			if (newton_pair || j < nlocal) {
			  f[j][0] -= fx;
			  f[j][1] -= fy;
			  f[j][2] -= fz;
			  }
			}
			pvector[4] += eadd;
			*/
			}
			//utils::logmesg(lmp,"epauli,s_fpair = {:.12},{:.12}\n",epauli,s_fpair);
			//utils::logmesg(lmp,"e1rforce, s_e1rforce, e2rforce, s_e2rforce = {:.12},{:.12},{:.12},{:.12}\n",e1rforce, s_e1rforce, e2rforce, s_e2rforce);
			e2rforce = spline * qqrd2e * e2rforce;
			erforce[j] += e2rforce;
			
			// Radial electron virial, iff flexible pressure flag set
			if (evflag && pressure_with_evirials_flag) {
			e1rvirial = eradius[i] * e1rforce;
			e2rvirial = eradius[j] * e2rforce;
			ev_tally_eff(i,j,nlocal,newton_pair,0.0,e1rvirial+e2rvirial);
			}
        }

        // electron (i) - nucleus (j) Coul interaction

        else if (abs(spin[i]) == 1 && q[j] != 0) {
          e1rforce = e2rforce = s_e1rforce = s_e2rforce = 0.0;
		  if(abs(q[j])==1){
			  Ueh(2.0,qhe,rc,eradius[i],&ecoul,&fpair,&e1rforce,ah_e1,a1_e1,a2_e1,a3_e1,a4_e1,a22_e1,a33_e1,a44_e1);
		  }
		  else{
		  	Uec(2.0,qoe,rc,eradius[i],&ecoul,&fpair,&e1rforce,ac_e1,ac1_e1,ac2_e1,ac3_e1,ac4_e1,poe1,ac22_e1,ac33_e1,ac44_e1);
			/*
			harmonic(rc, 0.4, 10000, &eadd, &fadd);//避免外层电子和O核重叠
			SmallRForce(delx,dely,delz,rc,fadd,&fx,&fy,&fz);//??
			f[i][0] += fx;
			f[i][1] += fy;
			f[i][2] += fz;
			if (newton_pair || j < nlocal) {
			  f[j][0] -= fx;
			  f[j][1] -= fy;
			  f[j][2] -= fz;
			  }
			}
			pvector[4] += eadd;
			*/
			}
			//utils::logmesg(lmp,"epauli,s_fpair = {:.12},{:.12}\n",epauli,s_fpair);
			e1rforce = spline * qqrd2e * e1rforce;
			erforce[i] += e1rforce;
			//utils::logmesg(lmp,"e1rforce,e2rforce = {:.12},{:.12}\n",e1rforce,e2rforce);
			// Radial electron virial, iff flexible pressure flag set
			if (evflag && pressure_with_evirials_flag) {
			e1rvirial = eradius[i] * e1rforce;
			e2rvirial = eradius[j] * e2rforce;
			ev_tally_eff(i,j,nlocal,newton_pair,0.0,e1rvirial+e2rvirial);
			}
        }

        // electron (i) - electron (j) interactions	
        else if (abs(spin[i]) == 1 && abs(spin[j]) == 1) {
          e1rforce = e2rforce = ecoul = 0.0;
          s_e1rforce = s_e2rforce = 0.0;
          Cee(4.0,rc,eradius[i],eradius[j],&ecoul,&fpair,&e1rforce,&e2rforce,ah_e11,a1_e11,a2_e11,a3_e11,a4_e11,kvv_e11,pvv_e11,a22_e11,a33_e11,a44_e11);
		  pauli1(rc,eradius[i],eradius[j],&epauli,&s_fpair,&s_e1rforce,&s_e2rforce, re_pauli11, rc_pauli11, Kk_11, rk_11, pk_11);
          
		  epauli *= hhmss2e;
          s_fpair *= hhmss2e;
		  //utils::logsmesg(lmp,"rc,eradius[i],eradius[j],ee_e1rforce,e2rforce,ecoul = {:.8},{:.8},{:.8},{:.8},{:.8},{:.8}\n",rc,eradius[i],eradius[j],e1rforce,e2rforce,ecoul);//zzh
          e1rforce = spline * (qqrd2e * e1rforce + hhmss2e * s_e1rforce);
          erforce[i] += e1rforce;
          e2rforce = spline * (qqrd2e * e2rforce + hhmss2e * s_e2rforce);
          erforce[j] += e2rforce;

          // Radial electron virial, iff flexible pressure flag set
          if (evflag && pressure_with_evirials_flag) {
            e1rvirial = eradius[i] * e1rforce;
            e2rvirial = eradius[j] * e2rforce;
            ev_tally_eff(i,j,nlocal,newton_pair,0.0,e1rvirial+e2rvirial);
          }
        }

        // Apply Coulomb conversion factor for all cases
        ecoul *= qqrd2e;
        fpair *= qqrd2e;
		
        // Sum up energy and force contributions
        energy = ecoul + epauli;
        fpair = fpair + s_fpair;

        // Apply cutoff spline
        fpair = fpair * spline - energy * dspline;
        energy = spline * energy;
		
		//eadd, fadd不用添加截断，因为只对近距离的原子对有效
		energy += eadd;

        // Tally cartesian forces
		//utils::logmesg(lmp,"fpair = {:.8}\n",fpair);
        SmallRForce(delx,dely,delz,rc,fpair,&fx,&fy,&fz);//??
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
		//utils::logmesg(lmp,"f[i][2] = {:.8}\n",f[i][2]);
        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
        }

        // Tally energy (in ecoul) and compute normal pressure virials
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,0.0,
                             energy,fx,fy,fz,delx,dely,delz);
        if (eflag_global) {
          if (newton_pair) {
            pvector[1] += spline * epauli;
            pvector[2] += spline * ecoul;
          }
          else {
            halfpauli = 0.5 * spline * epauli;
            halfcoul = 0.5 * spline * ecoul;
            if (i < nlocal) {
              pvector[1] += halfpauli;
              pvector[2] += halfcoul;
            }
            if (j < nlocal) {
              pvector[1] += halfpauli;
              pvector[2] += halfcoul;
            }
          }
        }
      }
	  //utils::logmesg(lmp,"pvector[1] = {:.8}\n",pvector[1]);
    }
	
	// limit electron stifness (size) for periodic systems, to max=half-box-size
	
	if (abs(spin[i]) == 1) {
      double dr, kfactor=hhmss2e*100.0;
      e1rforce = errestrain = 0.0;
	  if (eradius[i] > max_e1radius) {
		  //utils::logmesg(lmp,"eradius[i] = {:.8}\n",eradius[i]);
		  dr = eradius[i]-max_e1radius;
		  errestrain=0.5*kfactor*dr*dr;
		  e1rforce=-kfactor*dr;
		  if (eflag_global) pvector[3] += errestrain;
		  erforce[i] += e1rforce;

		  // Tally radial restrain energy and add radial restrain virial
		  if (evflag) {
			ev_tally_eff(i,i,nlocal,newton_pair,errestrain,0.0);
			if (pressure_with_evirials_flag)  // flexible electron pressure
			  ev_tally_eff(i,i,nlocal,newton_pair,0.0,eradius[i]*e1rforce);
		  }
		}
		else if (eradius[i] < min_e1radius) {
			//utils::logmesg(lmp,"eradius[i] = {:.8}\n",eradius[i]);
			dr = min_e1radius-eradius[i];
			errestrain=0.5*kfactor*dr*dr;
			e1rforce=kfactor*dr;
			if (eflag_global) pvector[3] += errestrain;
			erforce[i] += e1rforce;
			
			if (evflag) {
			ev_tally_eff(i,i,nlocal,newton_pair,errestrain,0.0);
			if (pressure_with_evirials_flag)  // flexible electron pressure
			  ev_tally_eff(i,i,nlocal,newton_pair,0.0,eradius[i]*e1rforce);
		  }
		}	
    }
	
	if (abs(spin[i]) == 3) {
      double dr, kfactor=hhmss2e*100.0;
      e1rforce = errestrain = 0.0;
	  if (eradius[i] > max_e2radius) {
		  dr = eradius[i]-max_e2radius;
		  errestrain=0.5*kfactor*dr*dr;
		  e1rforce=-kfactor*dr;
		  if (eflag_global) pvector[3] += errestrain;
		  erforce[i] += e1rforce;

		  // Tally radial restrain energy and add radial restrain virial
		  if (evflag) {
			ev_tally_eff(i,i,nlocal,newton_pair,errestrain,0.0);
			if (pressure_with_evirials_flag)  // flexible electron pressure
			  ev_tally_eff(i,i,nlocal,newton_pair,0.0,eradius[i]*e1rforce);
		  }
		}
		else if (eradius[i] < min_e2radius) {
			dr = min_e2radius-eradius[i];
			errestrain=0.5*kfactor*dr*dr;
			e1rforce=kfactor*dr;
			if (eflag_global) pvector[3] += errestrain;
			erforce[i] += e1rforce;
			
			if (evflag) {
			ev_tally_eff(i,i,nlocal,newton_pair,errestrain,0.0);
			if (pressure_with_evirials_flag)  // flexible electron pressure
			  ev_tally_eff(i,i,nlocal,newton_pair,0.0,eradius[i]*e1rforce);
		  }
		}	
    }
	
  }

  if (vflag_fdotr) {
    virial_fdotr_compute();
    if (pressure_with_evirials_flag) virial_eff_compute();
  }
  /*
  auto end3 = std::chrono::high_resolution_clock::now();
  auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2);
  time3 += duration3.count();
  if(update->ntimestep % 100 == 0)
	utils::logmesg(lmp,"time2, time3 = {:.6}, {:.6}\n", time2, time3);
*/
}

/* ----------------------------------------------------------------------
   eff-specific contribution to global virial
------------------------------------------------------------------------- */

void PairEffCut::virial_eff_compute()
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  double e_virial;
  int *spin = atom->spin;

  // sum over force on all particles including ghosts

  if (neighbor->includegroup == 0) {
    int nall = atom->nlocal + atom->nghost;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

  // neighbor includegroup flag is set
  // sum over force on initial nfirst particles and ghosts

  } else {
    int nall = atom->nfirst;
    for (int i = 0; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }

    nall = atom->nlocal + atom->nghost;
    for (int i = atom->nlocal; i < nall; i++) {
      if (spin[i]) {
        e_virial = erforce[i]*eradius[i]/3;
        virial[0] += e_virial;
        virial[1] += e_virial;
        virial[2] += e_virial;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into per-atom accumulators
   for virial radial electronic contributions
------------------------------------------------------------------------- */

void PairEffCut::ev_tally_eff(int i, int j, int nlocal, int newton_pair,
                              double energy, double e_virial)
{
  double energyhalf;
  double partial_evirial = e_virial/3.0;
  double half_partial_evirial = partial_evirial/2;

  int *spin = atom->spin;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_pair)
        eng_coul += energy;
      else {
        energyhalf = 0.5*energy;
        if (i < nlocal)
          eng_coul += energyhalf;
        if (j < nlocal)
          eng_coul += energyhalf;
      }
    }
    if (eflag_atom) {
      if (newton_pair || i < nlocal) eatom[i] += 0.5 * energy;
      if (newton_pair || j < nlocal) eatom[j] += 0.5 * energy;
    }
  }

  if (vflag_either) {
    if (vflag_global) {
      if (spin[i] && i < nlocal) {
        virial[0] += half_partial_evirial;
        virial[1] += half_partial_evirial;
        virial[2] += half_partial_evirial;
      }
      if (spin[j] && j < nlocal) {
        virial[0] += half_partial_evirial;
        virial[1] += half_partial_evirial;
        virial[2] += half_partial_evirial;
      }
    }
    if (vflag_atom) {
      if (spin[i]) {
        if (newton_pair || i < nlocal) {
          vatom[i][0] += half_partial_evirial;
          vatom[i][1] += half_partial_evirial;
          vatom[i][2] += half_partial_evirial;
        }
      }
      if (spin[j]) {
        if (newton_pair || j < nlocal) {
          vatom[j][0] += half_partial_evirial;
          vatom[j][1] += half_partial_evirial;
          vatom[j][2] += half_partial_evirial;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairEffCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  /*memory->create(PAULI_RE,"pair:PAULI_RE");
  memory->create(PAULI_RC,"pair:PAULI_RC");
  memory->create(PAULI_RHO,"pair:PAULI_RHO");*/
}

/* ---------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairEffCut::settings(int narg, char **arg)
{
	if (narg < 1)
	  error->all(FLERR,"Illegal pair_style command");
  
	cut_global=utils::numeric(FLERR,arg[0],false,lmp);
	qoo=utils::numeric(FLERR,arg[1],false,lmp);
	qoh=utils::numeric(FLERR,arg[2],false,lmp);
	qoe=utils::numeric(FLERR,arg[3],false,lmp);
	qhh=utils::numeric(FLERR,arg[4],false,lmp);
	qhe=utils::numeric(FLERR,arg[5],false,lmp);
	kappa_hh=utils::numeric(FLERR,arg[6],false,lmp);
	tau_hh=utils::numeric(FLERR,arg[7],false,lmp);
	rho_hh=utils::numeric(FLERR,arg[8],false,lmp);
	kappa_oh=utils::numeric(FLERR,arg[9],false,lmp);
	tau_oh=utils::numeric(FLERR,arg[10],false,lmp);
	rho_oh=utils::numeric(FLERR,arg[11],false,lmp);
	kappa_oo=utils::numeric(FLERR,arg[12],false,lmp);
	tau_oo=utils::numeric(FLERR,arg[13],false,lmp);
	rho_oo=utils::numeric(FLERR,arg[14],false,lmp);
	rc_wall=utils::numeric(FLERR,arg[15],false,lmp);
	A_gauss=utils::numeric(FLERR,arg[16],false,lmp);
	B_gauss=utils::numeric(FLERR,arg[17],false,lmp);
	ah_e1=utils::numeric(FLERR,arg[18],false,lmp);
	a1_e1=utils::numeric(FLERR,arg[19],false,lmp);
	a2_e1=utils::numeric(FLERR,arg[20],false,lmp);
	a3_e1=utils::numeric(FLERR,arg[21],false,lmp);
	a4_e1=utils::numeric(FLERR,arg[22],false,lmp);
	a22_e1=utils::numeric(FLERR,arg[23],false,lmp);
	a33_e1=utils::numeric(FLERR,arg[24],false,lmp);
	a44_e1=utils::numeric(FLERR,arg[25],false,lmp);
	ac_e1=utils::numeric(FLERR,arg[26],false,lmp);
	ac1_e1=utils::numeric(FLERR,arg[27],false,lmp);
	ac2_e1=utils::numeric(FLERR,arg[28],false,lmp);
	ac3_e1=utils::numeric(FLERR,arg[29],false,lmp);
	ac4_e1=utils::numeric(FLERR,arg[30],false,lmp);
	poe1=utils::numeric(FLERR,arg[31],false,lmp);
	ac22_e1=utils::numeric(FLERR,arg[32],false,lmp);
	ac33_e1=utils::numeric(FLERR,arg[33],false,lmp);
	ac44_e1=utils::numeric(FLERR,arg[34],false,lmp);
	ah_e11=utils::numeric(FLERR,arg[35],false,lmp);
	a1_e11=utils::numeric(FLERR,arg[36],false,lmp);
	a2_e11=utils::numeric(FLERR,arg[37],false,lmp);
	a3_e11=utils::numeric(FLERR,arg[38],false,lmp);
	a4_e11=utils::numeric(FLERR,arg[39],false,lmp);
	kvv_e11=utils::numeric(FLERR,arg[40],false,lmp);
	pvv_e11=utils::numeric(FLERR,arg[41],false,lmp);
	a22_e11=utils::numeric(FLERR,arg[42],false,lmp);
	a33_e11=utils::numeric(FLERR,arg[43],false,lmp);
	a44_e11=utils::numeric(FLERR,arg[44],false,lmp);
	re_pauli11=utils::numeric(FLERR,arg[45],false,lmp);
	rc_pauli11=utils::numeric(FLERR,arg[46],false,lmp);
	Kk_11=utils::numeric(FLERR,arg[47],false,lmp);
	rk_11=utils::numeric(FLERR,arg[48],false,lmp);
	pk_11=utils::numeric(FLERR,arg[49],false,lmp);
	ke_e1=utils::numeric(FLERR,arg[50],false,lmp);
	A_e1=utils::numeric(FLERR,arg[51],false,lmp);
	B_e1=utils::numeric(FLERR,arg[52],false,lmp);
	C_e1=utils::numeric(FLERR,arg[53],false,lmp);
	D_e1=utils::numeric(FLERR,arg[54],false,lmp);
	max_e1radius=utils::numeric(FLERR,arg[55],false,lmp);
	min_e1radius=utils::numeric(FLERR,arg[56],false,lmp);
	max_e2radius=utils::numeric(FLERR,arg[57],false,lmp);
	min_e2radius=utils::numeric(FLERR,arg[58],false,lmp);


	pressure_with_evirials_flag = 0;

	// Need to introduce 2 new constants w/out changing update.cpp

	if (force->qqr2e==332.06371) {        // i.e. Real units chosen
	h2e = 627.509;                      // hartree->kcal/mol
	hhmss2e = 175.72044219620075;       // hartree->kcal/mol * (Bohr->Angstrom)^2
	} else if (force->qqr2e==1.0) {        // electron units
	h2e = 1.0;
	hhmss2e = 1.0;
	} else error->all(FLERR,"Check your units");

	// reset cutoffs that have been explicitly set

	if (allocated) {
	int i,j;
	for (i = 1; i <= atom->ntypes; i++)
	  for (j = i; j <= atom->ntypes; j++)
		if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairEffCut::init_style()
{
  // error and warning checks

  if (!atom->q_flag || !atom->spin_flag ||
      !atom->eradius_flag || !atom->erforce_flag)
    error->all(FLERR,"Pair eff/cut requires atom attributes "
               "q, spin, eradius, erforce");

  // add hook to minimizer for eradius and erforce

  if (update->whichflag == 2)
    update->minimize->request(this,1,0.1);//zzh

  // make sure to use the appropriate timestep when using real units

  if (update->whichflag == 1) {
    if (utils::strmatch(update->unit_style,"^real") && update->dt_default)
      error->all(FLERR,"Must lower the default real units timestep for pEFF ");
  }

  // check if any atom's spin = 3 and ECP type was not set

  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  int flag = 0;
  //for (int i = 0; i < nlocal; i++)
    //if (spin[i] == 3) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

  // need a half neigh list

  neighbor->request(this,instance_me);
  //neighbor->add_request(this, NeighConst::REQ_FULL);
    // need a full neighbor list, to store info of Oa and Ob
  //int irequest = neighbor->request(this, instance_me);
  //neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type electron pairs (ECP-only)
------------------------------------------------------------------------- */

void PairEffCut::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if ((strcmp(arg[0],"*") == 0) || (strcmp(arg[1],"*") == 0)) {
    int ilo,ihi,jlo,jhi;
    utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

    double cut_one = cut_global;
    if (narg == 3) cut_one = utils::numeric(FLERR,arg[2],false,lmp);

    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        cut[i][j] = cut_one;
        setflag[i][j] = 1;
         count++;
      }
    }
    if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
  }
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairEffCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEffCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) fwrite(&cut[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEffCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairEffCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairEffCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   returns pointers to the log() of electron radius and corresponding force
   minimizer operates on log(radius) so radius never goes negative
   these arrays are stored locally by pair style
------------------------------------------------------------------------- */

void PairEffCut::min_xf_pointers(int /*ignore*/, double **xextra, double **fextra)
{
  // grow arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(min_eradius);
    memory->destroy(min_erforce);
    nmax = atom->nmax;
    memory->create(min_eradius,nmax,"pair:min_eradius");
    memory->create(min_erforce,nmax,"pair:min_erforce");
  }

  *xextra = min_eradius;
  *fextra = min_erforce;
}

/* ----------------------------------------------------------------------
   minimizer requests the log() of electron radius and corresponding force
   calculate and store in min_eradius and min_erforce
------------------------------------------------------------------------- */

void PairEffCut::min_xf_get(int /*ignore*/)
{
  double *eradius = atom->eradius;
  double *erforce = atom->erforce;
  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (spin[i]) {
      //min_eradius[i] = log(eradius[i]);
      min_eradius[i] = eradius[i];
      min_erforce[i] = erforce[i];
    } else min_eradius[i] = min_erforce[i] = 0.0;
}

/* ----------------------------------------------------------------------
   minimizer has changed the log() of electron radius
   propagate the change back to eradius
------------------------------------------------------------------------- */

void PairEffCut::min_x_set(int /*ignore*/)
{
  double *eradius = atom->eradius;
  int *spin = atom->spin;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    //if (spin[i]) eradius[i] = exp(min_eradius[i]);
    if (spin[i]) eradius[i] = min_eradius[i];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairEffCut::memory_usage()
{
  double bytes = (double)maxeatom * sizeof(double);
  bytes += (double)maxvatom*6 * sizeof(double);
  bytes += (double)2 * nmax * sizeof(double);
  return bytes;
}
