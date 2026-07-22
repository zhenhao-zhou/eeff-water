/* -*- c++ -*- ----------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eff/cut,PairEffCut);
// clang-format on
#else

#ifndef LMP_PAIR_EFF_CUT_H
#define LMP_PAIR_EFF_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairEffCut : public Pair {
 public:
  PairEffCut(class LAMMPS *);
  ~PairEffCut() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void min_pointers(double **, double **);
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  void min_xf_pointers(int, double **, double **) override;
  void min_xf_get(int) override;
  void min_x_set(int) override;
  double memory_usage() override;

 private:
  //double time1, time2, time3, time4;
  // 用来存储每个电子的环境因子与梯度
  double *chi_ke;       // length = atom->nlocal+atom->nghost
  double **dchi_ke;     // [N][3]
  
  // 本进程上所有 proton 和 core 的索引
  std::vector<int> proton_list;
  std::vector<int> core_list;
  std::vector<int> electron_list;
  
  int pressure_with_evirials_flag;
  double cut_global;
  double **cut;

	
	double qoo;
	double qoh;
	double qoe;
	double qhh;
	double qhe;
	double kappa_hh;
	double tau_hh;
	double rho_hh;
	double kappa_oh;
	double tau_oh;
	double rho_oh;
	double kappa_oo;
	double tau_oo;
	double rho_oo;
	double rc_wall;
	double A_gauss;
	double B_gauss;
	double ah_e1;
	double a1_e1;
	double a2_e1;
	double a3_e1;
	double a4_e1;
	double a22_e1;
	double a33_e1;
	double a44_e1;
	double ac_e1;
	double ac1_e1;
	double ac2_e1;
	double ac3_e1;
	double ac4_e1;
	double ac22_e1;
	double ac33_e1;
	double ac44_e1;
	double poe1;
	double ah_e11;
	double a1_e11;
	double a2_e11;
	double a3_e11;
	double a4_e11;
	double kvv_e11;
	double pvv_e11;
	double a22_e11;
	double a33_e11;
	double a44_e11;
	double re_pauli11;
	double rc_pauli11;
	double Kk_11;
	double rk_11;
	double pk_11;
	double ke_e1;
	double A_e1;
	double B_e1;
	double C_e1;
	double D_e1;
	double max_e1radius;
	double min_e1radius;
	double max_e2radius;
	double min_e2radius;



  //ke 
  /*Xq para
  double kq;
  double na;
  double nq;
  */
  double hhmss2e, h2e;
  int nmax;
  double *min_eradius, *min_erforce;

  void allocate();
  void virial_eff_compute();
  void ev_tally_eff(int, int, int, int, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
