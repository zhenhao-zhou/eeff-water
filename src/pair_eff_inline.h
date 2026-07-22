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

/* ----------------------------------------------------------------------
 Contributing authors: Andres Jaramillo-Botero, Hai Xiao, Julius Su (Caltech)
------------------------------------------------------------------------- */
#include "pair_eff_cut.h"

namespace LAMMPS_NS {
/*
#define PAULI_RE 0.9
#define PAULI_RC 1.125
#define PAULI_RHO -0.2
*/
#define ERF_TERMS1 12
#define ERF_TERMS2 7
#define DERF_TERMS 13

// 样条系数
static constexpr double c_proton[6] = {0.5, 0.0, 0.0, 0.625, -0.46875, 0.09375};
static constexpr double c_core1  [6] = {0.0, 0.0, 6.944444, 0.0, -19.290123, 12.860082};
static constexpr double c_core2  [6] = {-8.487654, 49.897119, -100.222908,
									  97.165066, -46.677336,  8.890921};

// smootherstep: S(u) = 6u^5 - 15u^4 + 10u^3
inline double smootherstep(double u) {
    return ((6*u - 15)*u + 10)*u*u*u;
}

// derivative: dS/du
inline double smootherstep_deriv(double u) {
    return (30*u*u*(u*(u - 2) + 1));
}

// zeta_H核（r0 = 0.37，区间 [0, 1.11]）
inline void zeta_proton(double r, double r0, double &z, double &dz) {
    const double r_min = 1e-10;
    const double r_peak = r0;
    const double r_max = 3.0 * r0;

    if (r <= r_min || r >= r_max) {
        z = 0.0;
        dz = 0.0;
    } else if (r <= r_peak) {
        double u = r / r_peak;  // maps [0, r0] → [0, 1]
        z = smootherstep(u);
        dz = smootherstep_deriv(u) / r_peak;
    } else {
        double u = (r_max - r) / (r_max - r_peak);  // maps [r0, 3r0] → [1, 0]
        z = smootherstep(u);
        dz = -smootherstep_deriv(u) / (r_max - r_peak);
    }
}

// zeta_O核（r0 = 0.63，区间 [0, 1.89]）
inline void zeta_core(double r, double r0, double &z, double &dz) {
    const double r_min = 1e-10;
    const double r_peak = r0;
    const double r_max = 3.0 * r0;

    if (r <= r_min || r >= r_max) {
        z = 0.0;
        dz = 0.0;
    } else if (r <= r_peak) {
        double u = r / r_peak;
        z = smootherstep(u);
        dz = smootherstep_deriv(u) / r_peak;
    } else {
        double u = (r_max - r) / (r_max - r_peak);
        z = smootherstep(u);
        dz = -smootherstep_deriv(u) / (r_max - r_peak);
    }
}

//权重取决于共价半径与原子电负性
inline void weight_function(double r, double r0, double chi_A, double &w, double &dw) {
    const double r_max = 3.0 * r0;

    if (r >= r_max) {
        w = 0.0;
        dw = 0.0;
    } else {
        double u = r / r_max;  // u ∈ [0, 1]
        double s = smootherstep(u);
        double dsdu = smootherstep_deriv(u);
        w = chi_A * (1.0 - s);
        dw = - chi_A * dsdu / r_max;  // dw/dr = ds/du * du/dr = ds/du * (1/r_max)
    }
}

inline int samespin(int i, int j) {
    // 判断两个整数是否同号
    if ((i > 0 && j > 0) || (i < 0 && j < 0)) {
        return 1; // 同号
    } else {
        return 0; // 不同号
    }
}

// inline functions for performance

double dist(int i, int j, double **x) {
    double dx = x[i][0] - x[j][0];
    double dy = x[i][1] - x[j][1];
    double dz = x[i][2] - x[j][2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

/* ---------------------------------------------------------------------- */
bool isInArray(int arr[], int size, int value) {  
    for (int i = 0; i < size; i++) {  
        if (arr[i] == value) {  
            return true;  
        }  
    }  
    return false;  
} 


void Unn(double r, double &e, double &f,
         double kappa,
         double tau,
         double rho,
		 double q)
{
    double u = pow(r/rho,3);
    double D = 1.0 + tau*u + u*u;

    // 能量
    e = q/r - kappa / D;

    // 力
    double f1 = q / (r*r);

    double dD = 3.0*(r*r)/(rho*rho*rho) * (tau + 2.0*u);
    double f2 = - kappa * dD / (D*D);

    f = f1 + f2;
}

void Uoo(double r, double &e, double &f,
         double kappa,
         double tau,
         double rho,
         double q,
         double rc_wall, 
         double A_gauss, 
         double B_gauss)
{
    // --- 1. 长程静电 + 屏蔽项 (原始逻辑) ---
    double u = pow(r/rho, 3);
    double D = 1.0 + tau*u + u*u;

    // 长程能量
    e = q/r - kappa / D;

    // 长程力
    double f1 = q / (r*r);
    double dD = 3.0*(r*r)/(rho*rho*rho) * (tau + 2.0*u);
    double f2 = - kappa * dD / (D*D);

    f = f1 + f2;

    // --- 2. 短程指数平滑排斥 (Exponential Tapered Repulsion) ---
    // 形式: E = A * exp(-B*r) * (1 - r/rc)^2
    // 目的: 增强近距离排斥以防止高压下密度过大
    
    if (r < rc_wall) {
        // --- 中间变量计算 ---
        double r_norm = r / rc_wall;          // r / rc
        double taper = 1.0 - r_norm;          // (1 - r/rc)
        double taper2 = taper * taper;        // (1 - r/rc)^2
        
        // 指数衰减项: exp(-B*r)
        double expo = exp( -B_gauss * r ); 

        // --- 能量修正 ---
        double e_rep = A_gauss * expo * taper2;
        e += e_rep;

        // --- 力修正 ---
        // F_rep = -dE/dr
        // 推导: E' = A * [ -B*exp*taper^2 + exp*2*taper*(-1/rc) ]
        //          = -A * exp * taper * [ B*taper + 2/rc ]
        // f_rep = -E' = A * exp * taper * [ B*taper + 2/rc ]
        
        double term_bracket = B_gauss * taper + 2.0 / rc_wall;
        double f_rep = A_gauss * expo * taper * term_bracket;

        f += f_rep;
    }
}

/* ---------------------------------------------------------------------- */

inline void KinElec(
    double r,
    double *eke,
    double *frc,
    // 原 K
    double K,
    // 新增修饰参数
    double A,  // 短程强度
    double B,  // 短程衰减 (r^2 前)
    double C,  // 中程强度
    double D   // 中程衰减 (r 前)
) {
    double r2    = r * r;
    double inv_r2= 1.0 / r2;
    double inv_r3= inv_r2 / r;

    // 修饰因子及其导数
    double e1 = exp(-B * r2);
    double e2 = exp(-D * r);
    double M  = 1.0 + A * e1 + C * e2;
    double dM = -2.0 * A * B * r * e1
                -   C * D      * e2;

    // 势能
    *eke += 0.5 * K * (M * inv_r2);

    // 力
    // F = K*chi_f * [ 2*M/r^3 - (dM)/r^2 ]
    double F = 0.5 * K * (
                   2.0 * M * inv_r3
                 -        dM * inv_r2
               );
    *frc += F;
}

/* ---------------------------------------------------------------------- */

inline void ElecNucNuc(double q, double rc, double *ecoul, double *frc, double knn, double tnn, double pnn)
{
  *ecoul += q / rc - knn/(1+tnn*pow(rc/pnn,3)+pow(rc/pnn,6));
  *frc += knn*(-3*pow(rc, 2)* tnn/pow(pnn, 3) - 6*pow(rc, 5)/pow(pnn, 6))/pow((1 + pow(rc, 3)* tnn/pow(pnn, 3) + 
  	pow(rc, 6)/pow(pnn, 6)),2) + q/pow(rc, 2);
}

/* ---------------------------------------------------------------------- */

inline void Ueh(double multi, double qhe, double rc, double re1, double *ecoul, double *frc, double *fre1,
                double ah, double a1, double a2, double a3, double a4,
                double a22, double a33, double a44)
{
    double ah_plus_re1 = ah + re1;
    /* r0 = rc / (ah + re1) */
    double r0 = rc / ah_plus_re1;

    /* build terms with variable exponents */
    double term1 = pow(a1, 8.0);
    double coeff2 = pow(a2, 8.0 - a22);
    double coeff3 = pow(a3, 8.0 - a33);
    double coeff4 = pow(a4, 8.0 - a44);

    double term2 = coeff2 * pow(r0, a22);
    double term3 = coeff3 * pow(r0, a33);
    double term4 = coeff4 * pow(r0, a44);
    double term5 = pow(r0, 8.0);

    double A = term1 + term2 + term3 + term4 + term5;

    /* A^(-1/8) and A^(-9/8) */
    double A_inv_1_8 = pow(A, -0.125);    /* -1/8 */
    double A_inv_9_8 = pow(A, -1.125);    /* -9/8 */

    /* update energy */
    *ecoul += - multi * qhe * A_inv_1_8 / ah_plus_re1;

    /* derivative dA/dr0 (use general exponents) */
    double dAdr0 = 0.0;
    /* if a22 != 0.0 then term contributes a22 * coeff2 * r0^(a22-1),  
       if a22==0 the coefficient a22 makes contribution zero anyway */
    dAdr0 += coeff2 * a22 * pow(r0, a22 - 1.0);
    dAdr0 += coeff3 * a33 * pow(r0, a33 - 1.0);
    dAdr0 += coeff4 * a44 * pow(r0, a44 - 1.0);
    dAdr0 += 8.0 * pow(r0, 7.0); /* derivative of r0^8 */

    /* dE/drc (note: r0 = rc/(ah+re1) => dr0/drc = 1/(ah+re1))
       Using chain rule and same factoring pattern as original:
       dE_drc (without extra 'multi' factor) = A_inv_9_8/(8 * (ah+re1)^2) * dAdr0
    */
    double denom2 = ah_plus_re1 * ah_plus_re1;
    double denom3 = denom2 * ah_plus_re1;
    double dE_drc = (A_inv_9_8 / (8.0 * denom2)) * dAdr0;

    /* dE/dre1 (two contributions: through A via r0 and through explicit 1/(ah+re1))
       dE_dre1 (without extra 'multi' factor) = - (rc * A_inv_9_8 / (8 * (ah+re1)^3)) * dAdr0
                                                    + A_inv_1_8 / (ah+re1)^2
    */
    double dE_dre1 = - (rc * A_inv_9_8 / (8.0 * denom3)) * dAdr0
                     + A_inv_1_8 / denom2;

    /* Forces are negative gradients; follow original pattern of multiplying by -multi here */
    *frc  += - multi * qhe * dE_drc;
    *fre1 += - multi * qhe * dE_dre1;
}

inline void Uec(double multi, double qoe, double rc, double re1, double *ecoul, double *frc, double *fre1,
                double ah, double a1, double a2, double a3, double a4, double pvo4,
                double a22, double a33, double a44)
{
    double s = ah + re1;                       /* short alias */
    double r0 = rc / s;

    /* protection: avoid s==0 or r0==0 (would cause divisions/pow domain issues) */
    if (s == 0.0 || r0 == 0.0) {
        return;
    }

    /* coefficients with variable exponents */
    double term1 = pow(a1, 8.0);
    double coeff2 = pow(a2, 8.0 - a22);
    double coeff3 = pow(a3, 8.0 - a33);
    double coeff4 = pow(a4, 8.0 - a44);

    /* r0 powers needed */
    double r0_6 = pow(r0, 6.0);
    double r0_7 = pow(r0, 7.0);
    double r0_8 = pow(r0, 8.0);

    /* A = term1 + coeff2 * r0^a22 + coeff3 * r0^a33 + coeff4 * r0^a44 + r0^8 */
    double term2 = coeff2 * pow(r0, a22);
    double term3 = coeff3 * pow(r0, a33);
    double term4 = coeff4 * pow(r0, a44);
    double A = term1 + term2 + term3 + term4 + r0_8;

    /* A^{-1/8} and A^{-9/8} */
    double A_inv_1_8 = pow(A, -0.125);   /* A^{-1/8} */
    double A_inv_9_8 = pow(A, -1.125);   /* A^{-9/8} */

    /* accumulate energy (same structure as before) */
    double pvo4_6 = pow(pvo4, 6.0);
    *ecoul += multi * (-qoe * A_inv_1_8 / s + pvo4_6 / r0_6);

    /* dA/dr0 (general exponents):
       d(term2)/dr0 = coeff2 * a22 * r0^(a22-1), etc.
       plus derivative of r0^8 = 8*r0^7
    */
    double dA_dr0 = 0.0;
    dA_dr0 += coeff2 * a22 * pow(r0, a22 - 1.0);
    dA_dr0 += coeff3 * a33 * pow(r0, a33 - 1.0);
    dA_dr0 += coeff4 * a44 * pow(r0, a44 - 1.0);
    dA_dr0 += 8.0 * r0_7;

    /* --- derivative w.r.t rc ---
       term1 = -qoe * A^{-1/8} / s  -> d/drc = (qoe/8) * A^{-9/8} * (dA/dr0) / s^2
       term2 = pvo4^6 / r0^6       -> d/drc = -6 * pvo4^6 / (r0^7 * s)
    */
    double denom_s_sq = s * s;
    double d_term1_drc = (qoe * 0.125) * A_inv_9_8 * dA_dr0 / denom_s_sq;
    double d_term2_drc = -6.0 * pvo4_6 / (r0_7 * s);

    double dE_drc = multi * (d_term1_drc + d_term2_drc);

    /* *frc is -dE/drc */
    *frc += -dE_drc;

    /* --- derivative w.r.t re1 ---
       For term1:
         d/dre1 = - (qoe/8) * A^{-9/8} * (dA/dr0) * r0 / s^2 + qoe * A^{-1/8} / s^2
         (uses dr0/dre1 = -rc / s^2 = - r0 / s)
       For term2:
         d/dre1 = 6 * pvo4^6 / (r0^6 * s)
         (derivation: dr0/dre1 = -rc/s^2 => overall +6 * pvo4^6 * rc/(r0^7*s^2) = 6*pvo4^6/(r0^6 * s))
    */
    double d_term1_dre1 = - (qoe * 0.125) * A_inv_9_8 * dA_dr0 * r0 / denom_s_sq
                          + qoe * A_inv_1_8 / denom_s_sq;
    double d_term2_dre1 = 6.0 * pvo4_6 / (r0_6 * s);

    double dE_dre1 = multi * (d_term1_dre1 + d_term2_dre1);

    /* *fre1 is -dE/dre1 */
    *fre1 += -dE_dre1;
}


inline void Cee(double multi, double rc, double re1, double re2,
                double *ecoul, double *frc, double *fre1, double *fre2,
                double ah,
                double a1, double a2, double a3, double a4,
                double kvv, double pvv,
                double a22, double a33, double a44)
{
    /* Protective small number to avoid division by zero */
    const double EPS = 1e-16;

    /* compute rlen = sqrt(re1^2 + re2^2) and s = rlen + ah */
    double rlen = sqrt(re1*re1 + re2*re2);
    double s = rlen + ah;

    if (s <= EPS) {
        /* pathological case: return without changing anything */
        return;
    }

    /* r0 = rc / s */
    double r0 = rc / s;

    /* Precompute some basic powers of r0 for later use */
    double r0_5 = pow(r0, 5.0);
    double r0_7 = pow(r0, 7.0);
    double r0_8 = pow(r0, 8.0);

    /* Build A with variable exponents */
    double term1 = pow(a1, 8.0);
    double coeff2 = pow(a2, 8.0 - a22);
    double coeff3 = pow(a3, 8.0 - a33);
    double coeff4 = pow(a4, 8.0 - a44);

    double term2 = coeff2 * pow(r0, a22);
    double term3 = coeff3 * pow(r0, a33);
    double term4 = coeff4 * pow(r0, a44);

    double A = term1 + term2 + term3 + term4 + r0_8;

    /* avoid invalid A (would make pow(A, negative) invalid) */
    if (A <= 0.0) {
        return;
    }

    double A_inv_1_8 = pow(A, -0.125);   /* A^{-1/8} */
    double A_inv_9_8 = pow(A, -1.125);   /* A^{-9/8} */

    /* accumulate energy (same form as before) */
    double B = 1.0 + pow(r0 / pvv, 6.0);    /* assume pvv != 0 */
    *ecoul += multi * ( A_inv_1_8 / s + kvv / B );

    /* dA/dr0 (general exponents):
       d(term2)/dr0 = coeff2 * a22 * r0^(a22-1), etc.
       plus derivative of r0^8 = 8 * r0^7
    */
    double dAdr0 = 0.0;
    dAdr0 += coeff2 * a22 * pow(r0, a22 - 1.0);
    dAdr0 += coeff3 * a33 * pow(r0, a33 - 1.0);
    dAdr0 += coeff4 * a44 * pow(r0, a44 - 1.0);
    dAdr0 += 8.0 * r0_7;

    /* ----- Part 1: derivatives of E1 = multi * A^{-1/8} / s ----- */
    /* ∂E1/∂A = multi * (-1/8) * A^{-9/8} / s */
    /* dr0/drc = 1/s */
    /* => dE1/drc = ∂E1/∂A * dA/dr0 * dr0/drc
                 = -multi*(1/8) * A^{-9/8} * dAdr0 / (s * s) */
    double dE1_drc = - multi * (1.0 / 8.0) * A_inv_9_8 * dAdr0 / (s * s);

    /* dr0/dre1 = - rc * (ds/dre1) / s^2  with ds/dre1 = re1 / rlen (if rlen>0) */
    double dsdre1 = (rlen > EPS) ? (re1 / rlen) : 0.0;
    double dsdre2 = (rlen > EPS) ? (re2 / rlen) : 0.0;
    double dr0_dre1 = - rc * dsdre1 / (s * s);
    double dr0_dre2 = - rc * dsdre2 / (s * s);

    /* ∂E1/∂s = - multi * A^{-1/8} / s^2 */
    double dE1_ds = - multi * A_inv_1_8 / (s * s);

    /* dE1/dre1 = ∂E1/∂A * dA/dr0 * dr0/dre1 + ∂E1/∂s * ds/dre1
       note ∂E1/∂A = multi * (-1/8) * A^{-9/8} / s
    */
    double dE1_dre1 = - multi * (1.0 / 8.0) * A_inv_9_8 * dAdr0 * dr0_dre1 / s + dE1_ds * dsdre1;
    double dE1_dre2 = - multi * (1.0 / 8.0) * A_inv_9_8 * dAdr0 * dr0_dre2 / s + dE1_ds * dsdre2;

    /* ----- Part 2: derivatives of E2 = multi * kvv / B, where B = 1 + (r0/pvv)^6 ----- */
    double pvv6 = pow(pvv, 6.0);  /* assume pvv != 0 */
    double dB_dr0 = 6.0 * r0_5 / pvv6;   /* d/d r0 of (r0/pvv)^6 = 6 * r0^5 / pvv^6 */
    double dE2_dr0 = - multi * kvv * dB_dr0 / (B * B);

    /* chain rules for r0 */
    double dE2_drc = dE2_dr0 * (1.0 / s);
    double dE2_dre1 = dE2_dr0 * dr0_dre1;
    double dE2_dre2 = dE2_dr0 * dr0_dre2;

    /* ----- total derivatives ----- */
    double dE_drc  = dE1_drc  + dE2_drc;
    double dE_dre1 = dE1_dre1 + dE2_dre1;
    double dE_dre2 = dE1_dre2 + dE2_dre2;

    /* forces are negative gradients */
    *frc  += - dE_drc;
    *fre1 += - dE_dre1;
    *fre2 += - dE_dre2;
}



inline void Cee_0(double re1, double *ecoul, double *fre1, double ah, double a1,double kvv)
{
    const double sqrt2 = 1.4142135623730951;

    // a0 = sqrt(2)*re1 + ah
    double a0 = sqrt2 * re1 + ah;

	// 累加能量
	*ecoul += 1 / (a1*a0) + kvv;

	// d(1/a0)/dre1 = -sqrt2 / a0^2
	double dInvA0 = -sqrt2 / (a0 * a0);

	// du/dre1 = 0 （u 为常数）
	// 故 fre1 = prefac * (0 + u * dInvA0)
	*fre1 += - 1 / a1 * dInvA0;
}

inline void pauli1(double rc, double re1, double re2, double* epauli,
	double* frc, double* fre1, double* fre2, double PAULI_RE, double PAULI_RC, double Kk, double rk, double pk)
{
    // safety epsilon to avoid divide-by-zero
    const double EPS = 1e-12;

    // scale parameters as original
    re1 *= PAULI_RE;
    re2 *= PAULI_RE;
    rc  *= PAULI_RC;

    // compute s = re1^2 + re2^2
    double s = re1 * re1 + re2 * re2;
    if (s < EPS) {
        // too small, can't reliably compute reff / y; abort without changing energies/forces
        return;
    }

    double sqrt_s = sqrt(s);

    // reff = rc / sqrt(s)
    double reff = rc / sqrt_s;
	
	if ((reff * reff / rk) > 40.0)
	{
		*epauli += 0;
		*frc += 0;
		*fre1 += 0;
		*fre2 += 0;
		return;
	}

    // To avoid division by zero when computing y = (re1/re2 + re2/re1)/2,
    // use safe copies for direct divides (but keep original re1/re2 values for s)
    double re1s = (fabs(re1) < EPS) ? (re1 >= 0.0 ? EPS : -EPS) : re1;
    double re2s = (fabs(re2) < EPS) ? (re2 >= 0.0 ? EPS : -EPS) : re2;

    // y = (re1/re2 + re2/re1) / 2
    double y = 0.5 * (re1s / re2s + re2s / re1s);

    // precompute common terms
    double y2 = y * y;
    double y3 = y2 * y;
    double reff2 = reff * reff;
    double exp_term = exp(reff2 / rk);

    // M = y^3 * exp(reff^2 / rk) - 1
    double M = y3 * exp_term - 1.0;

    // Denominator D = M * s
    double D = M * s;

    // Protect against D ~ 0 (original formula would divide by D). If D is extremely small,
    // give up computing derivatives to avoid huge / NaN values.
    if (fabs(D) < EPS) {
        return;
    }

    // Numerator N = Kk * ( pk*(y^2 - 1) + reff^2 / rk )
    double N = Kk * ( pk * (y2 - 1.0) + reff2 / rk );

    // Energy add (same as original form; use computed N and D)
    *epauli += 2 * N / D;

    // ---- Now compute partial derivatives ----
    // Need:
    // dE/dx = (dN * D - N * dD) / D^2,  where dD = s * dM + M * ds

    // Derivatives of elementary quantities:

    // dy/dre1, dy/dre2 (using safe re1s,re2s)
    double dy_dre1 = 0.5 * (1.0 / re2s - (re2s / (re1s * re1s)));
    double dy_dre2 = 0.5 * (1.0 / re1s - (re1s / (re2s * re2s)));
    double dy_drc  = 0.0;

    // ds/dre1, ds/dre2, ds/drc
    double ds_dre1 = 2.0 * re1;
    double ds_dre2 = 2.0 * re2;
    double ds_drc  = 0.0;

    // dreff/drc = 1 / sqrt_s
    double dreff_drc = 1.0 / sqrt_s;

    // dreff/dre1 = -rc * re1 / s^(3/2)  (since dreff/dre1 = rc * d(s^-1/2)/dre1)
    double s_times_sqrt_s = s * sqrt_s; // s^(3/2)
    double dreff_dre1 = - (rc * re1) / (s_times_sqrt_s);
    double dreff_dre2 = - (rc * re2) / (s_times_sqrt_s);

    // d(reff^2 / rk)/dx = (2 * reff / rk) * dreff/dx
    double coef_reff2 = 2.0 * reff / rk;
    double d_reff2_rk_drc  = coef_reff2 * dreff_drc;
    double d_reff2_rk_dre1 = coef_reff2 * dreff_dre1;
    double d_reff2_rk_dre2 = coef_reff2 * dreff_dre2;

    // dM/dx = exp_term * ( 3*y^2 * dy/dx + y^3 * d(reff^2/rk)/dx )
    double common_exp = exp_term;
    double dM_drc  = common_exp * ( 3.0 * y2 * dy_drc  + y3 * d_reff2_rk_drc );
    double dM_dre1 = common_exp * ( 3.0 * y2 * dy_dre1 + y3 * d_reff2_rk_dre1 );
    double dM_dre2 = common_exp * ( 3.0 * y2 * dy_dre2 + y3 * d_reff2_rk_dre2 );

    // dN/dx = Kk * ( pk * 2 * y * dy/dx + d(reff^2/rk)/dx )
    double dN_drc  = Kk * ( pk * 2.0 * y * dy_drc  + d_reff2_rk_drc );
    double dN_dre1 = Kk * ( pk * 2.0 * y * dy_dre1 + d_reff2_rk_dre1 );
    double dN_dre2 = Kk * ( pk * 2.0 * y * dy_dre2 + d_reff2_rk_dre2 );

    // dD/dx = s * dM/dx + M * ds/dx
    double dD_drc  = s * dM_drc  + M * ds_drc;
    double dD_dre1 = s * dM_dre1 + M * ds_dre1;
    double dD_dre2 = s * dM_dre2 + M * ds_dre2;

    // dE/dx
    double D2 = D * D;
    double dE_drc  = (dN_drc  * D - N * dD_drc)  / D2;
    double dE_dre1 = (dN_dre1 * D - N * dD_dre1) / D2;
    double dE_dre2 = (dN_dre2 * D - N * dD_dre2) / D2;

    // Forces are negative gradients
    *frc  += -2 * PAULI_RC * dE_drc;
    *fre1 += -2 * PAULI_RE * dE_dre1;
    *fre2 += -2 * PAULI_RE * dE_dre2;

    return;
}

inline void sigmoid_1(int samespin, double rc, double rmin, double *ecoul, double *frc)
{
	if(samespin)
	{
		if(rc < rmin)
		{
			*ecoul += rmin-rc;
			*frc += 1;
		}
		if(rc < (rmin+1.0))
		{
			*ecoul += 0.5 / (1 + exp((rc-rmin)/0.005));
			*frc += 100.0 * exp(200.0 * rc - 200.0 * rmin) / pow(exp(200.0 * rc - 200.0 * rmin) + 1, 2);
		}
		else
		{
			*ecoul += 0;
			*frc += 0;
		}
	}
}

inline void sigmoid_3(int samespin, double rc, double rmin, double *ecoul, double *frc)
{
	if(samespin)
	{
		if(rc < rmin)
		{
			*ecoul += rmin-rc;
			*frc += 1;
		}
		if(rc < (rmin+1.0))
		{
			*ecoul += 0.5 / (1 + exp((rc-rmin)/0.005));
			*frc += 100.0 * exp(200.0 * rc - 200.0 * rmin) / pow(exp(200.0 * rc - 200.0 * rmin) + 1, 2);
		}
		else
		{
			*ecoul += 0;
			*frc += 0;
		}
	}
}

inline void sigmoid_2(double rc, double rmax, double *ecoul, double *frc)
{
	if(rc > rmax + 0.2)
	{
		*ecoul += 0;
		*frc += 0;
	}
	else
	{
		*ecoul += 1.0 / (1 + exp((rmax-rc)/0.005));
		*frc += -200.0 * exp(-200.0 * rc + 200.0 * rmax) / pow(exp(-200.0 * rc + 200.0 * rmax) + 1, 2);
	}
}

inline void sigmoid_add(int samespin, double rc, double rmin, double rnn, double r_ceil, double r_floor, double *ecoul, double *frc)
{
	if(samespin)
	{
		if(rc < rmin)
		{
			*ecoul += rmin-rc;
			*frc += 1;
		}
		if(rc < (rmin+0.6))
		{
			*ecoul += -8 * (rnn - r_ceil) * (rnn - r_floor) / pow(r_ceil - r_floor, 2) / (1 + exp((rc-rmin)/0.05));
			*frc += 40.0 * (4 * r_ceil - 4 * rnn) * (-r_floor + rnn) * exp(20.0 * rc - 20.0 * rmin) / (pow(r_ceil - r_floor, 2) * pow(exp(20.0 * rc - 20.0 * rmin) + 1, 2));
		}
		else
		{
			*ecoul += 0;
			*frc += 0;
		}
	}
}

inline void harmonic(double r, double rmin, double k, double *ecoul, double *frc)
{
	if(r > rmin)
	{
		*ecoul += 0;
		*frc += 0;
	}
	else
	{
		// force & energy
		*frc = -2.0 * k * (r - rmin);
		*ecoul = k * (r - rmin) * (r - rmin);
	}
}
/* ---------------------------------------------------------------------- */

inline void pauli1(int samespin, double rc, double re1, double re2, double chi1, double chi2, double* epauli,
	double* frc, double* fre1, double* fre2, double PAULI_RE, double PAULI_RC, double Kk, double rk, double pk)
{
	re1 *= PAULI_RE;
	re2 *= PAULI_RE;
	rc *= PAULI_RC;
	double PK = pk * (1.4 - 0.4 * (chi1+chi2));
	double reff = rc / pow((re1 * re1 + re2 * re2), 0.5);
	double r2 = pow(re1 / re2 + re2 / re1, 2);
	double r3 = (pow(re1, 2) + pow(re2, 2));
	double y = (re1 / re2 + re2 / re1) / 2;
	if ((reff * reff / rk) > 50.0)
	{
		*epauli += 0;
		*frc += 0;
		*fre1 += 0;
		*fre2 += 0;
	}
	else
	{
		if (samespin)
		{
		if (rc == 0 && y == 1)
			{
				*epauli += 2 * Kk / (re1 * re2);
				*frc += 0;
				*fre1 += 2 * Kk / (re1 * re1 * re1);
				*fre2 += 2 * Kk / (re2 * re2 * re2);
			}
		else{
				*epauli += 2 * Kk * (PK * (y * y - 1) + reff * reff / rk) / (y * y * y * exp(reff * reff / rk) - 1) / (re1 * re1 + re2 * re2);

				*frc += 2 * PAULI_RC *(Kk * rc * (PK * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2) * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - 2 * Kk * rc / (rk * pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

				*fre1 += 2 * PAULI_RE * (2 * Kk * re1 * (PK * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (PK * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re1 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (3.0 / re2 - 3 * re2 / pow(re1, 2)) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (PK * (2.0 / re2 - 2 * re2 / pow(re1, 2)) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re1 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

				*fre2 += 2 * PAULI_RE * (2 * Kk * re2 * (PK * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (PK * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re2 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (-3 * re1 / pow(re2, 2) + 3.0 / re1) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (PK * (-2 * re1 / pow(re2, 2) + 2 / re1) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re2 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));
			}
		}
	}	
}
//dpauli_dchi
inline double dPauli_dchi(int samespin, double rc, double re1, double re2, double PAULI_RE, double PAULI_RC, double Kk, double rk, double pk)
{
	if (samespin){
		re1 *= PAULI_RE;
		re2 *= PAULI_RE;
		rc *= PAULI_RC;
		double reff = rc / pow((re1 * re1 + re2 * re2), 0.5);
		double y = (re1 / re2 + re2 / re1) / 2;
		if ((reff * reff / rk) > 50.0)
		{
			return 0.0;
		}
		else
			return - 2 * Kk * pk * 0.4 * (y * y - 1) / (y * y * y * exp(reff * reff / rk) - 1) / (re1 * re1 + re2 * re2);
	}
	else 
	{
		return 0;
	}
}

inline void pauli1_soft(int samespin, double rc, double re1, double re2, double* epauli,
	double* frc, double* fre1, double* fre2, double PAULI_RE, double PAULI_RC, double Kk, double Ku, double rk, double ru, double pk,double pu)
{
	re1 *= PAULI_RE;
	re2 *= PAULI_RE;
	rc *= PAULI_RC;
	double reff = rc / pow((re1 * re1 + re2 * re2), 0.5);
	double r2 = pow(re1 / re2 + re2 / re1, 2);
	double r3 = (pow(re1, 2) + pow(re2, 2));
	double y = (re1 / re2 + re2 / re1) / 2;
	if ((reff * reff / rk) > 50.0 || (reff * reff / ru) > 50.0)
	{
		*epauli += 0;
		*frc += 0;
		*fre1 += 0;
		*fre2 += 0;
	}
	else
	{
		if (samespin)
		{
		if (rc == 0 && y == 1)
			{
				*epauli += Kk / (re1 * re2);
				*frc += 0;
				*fre1 += Kk / (re1 * re1 * re1);
				*fre2 += Kk / (re2 * re2 * re2);
			}
		else if(rc > 1.0){
				*epauli += Kk * (pk * (y * y - 1) + reff * reff / rk) / (y * y * y * exp(reff * reff / rk) - 1) / (re1 * re1 + re2 * re2);

				*frc += PAULI_RC *(Kk * rc * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2) * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - 2 * Kk * rc / (rk * pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

				*fre1 += PAULI_RE * (2 * Kk * re1 * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re1 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (3.0 / re2 - 3 * re2 / pow(re1, 2)) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (pk * (2.0 / re2 - 2 * re2 / pow(re1, 2)) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re1 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

				*fre2 += PAULI_RE * (2 * Kk * re2 * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re2 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (-3 * re1 / pow(re2, 2) + 3.0 / re1) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (pk * (-2 * re1 / pow(re2, 2) + 2 / re1) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re2 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));
			}
		else
		{
			*epauli += exp(rc-1)*Kk * (pk * (y * y - 1) + reff * reff / rk) / (y * y * y * exp(reff * reff / rk) - 1) / (re1 * re1 + re2 * re2);
			
			*frc +=Kk * rc * (pk * ((re1 / (2 * re2) + re2 / (2 * re1)) * (re1 / re2 + re2 / re1) / 2 - 1) + pow(rc, 2) / (rk * (pow(re1, 2) + pow(re2, 2)))) * pow(re1 / re2 + re2 / re1, 3) * exp(pow(rc, 2) / (rk * (pow(re1, 2) + pow(re2, 2)))) * exp(rc - 1) / (4 * rk * pow(pow(re1, 2) + pow(re2, 2), 2) * pow(pow(re1 / re2 + re2 / re1, 3) * exp(pow(rc, 2) / (rk * (pow(re1, 2) + pow(re2, 2)))) / 8 - 1, 2)) - 2 * Kk * rc * exp(rc - 1) / (rk * pow(pow(re1, 2) + pow(re2, 2), 2) * (pow(re1 / re2 + re2 / re1, 3) * exp(pow(rc, 2) / (rk * (pow(re1, 2) + pow(re2, 2)))) / 8 - 1)) - Kk * (pk * ((re1 / (2 * re2) + re2 / (2 * re1)) * (re1 / re2 + re2 / re1) / 2 - 1) + pow(rc, 2) / (rk * (pow(re1, 2) + pow(re2, 2)))) * exp(rc - 1) / ((pow(re1, 2) + pow(re2, 2)) * (pow(re1 / re2 + re2 / re1, 3) * exp(pow(rc, 2) / (rk * (pow(re1, 2) + pow(re2, 2)))) / 8 - 1));
			
			*fre1 += PAULI_RE * exp(rc-1) * (2 * Kk * re1 * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re1 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (3.0 / re2 - 3 * re2 / pow(re1, 2)) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (pk * (2.0 / re2 - 2 * re2 / pow(re1, 2)) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re1 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

			*fre2 += PAULI_RE * exp(rc-1) * (2 * Kk * re2 * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re2 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (-3 * re1 / pow(re2, 2) + 3.0 / re1) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (pk * (-2 * re1 / pow(re2, 2) + 2 / re1) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re2 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));
			
		}
		}
		else
		{
			*epauli += Ku * (pu * (y * y - 1) + reff * reff / ru) / (y * y * y * exp(reff * reff / ru) + 1) / (re1 * re1 + re2 * re2);
			
			*frc += PAULI_RC *(Ku * rc * (pu * (r2 / 4.0 - 1) + pow(rc, 2) / (ru * r3)) * pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / (4 * ru * pow(r3, 2) * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1, 2)) - 2 * Ku * rc / (ru * pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1)));

			*fre1 += PAULI_RE * (2 * Ku * re1 * (pu * (r2 / 4.0 - 1) + pow(rc, 2) / (ru * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1)) - Ku * (pu * (r2 / 4.0 - 1) + pow(rc, 2) / (ru * r3)) * (pow(rc, 2) * re1 * pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / (4 * ru * pow(r3, 2)) - (3.0 / re2 - 3 * re2 / pow(re1, 2)) * r2 * exp(pow(rc, 2) / (ru * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1, 2)) - Ku * (pu * (2.0 / re2 - 2 * re2 / pow(re1, 2)) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re1 / (ru * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1)));

			*fre2 += PAULI_RE * (2 * Ku * re2 * (pu * (r2 / 4.0 - 1) + pow(rc, 2) / (ru * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1)) - Ku * (pu * (r2 / 4.0 - 1) + pow(rc, 2) / (ru * r3)) * (pow(rc, 2) * re2 * pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / (4 * ru * pow(r3, 2)) - (-3 * re1 / pow(re2, 2) + 3 / re1) * r2 * exp(pow(rc, 2) / (ru * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1, 2)) - Ku * (pu * (-2 * re1 / pow(re2, 2) + 2 / re1) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re2 / (ru * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (ru * r3)) / 8.0 + 1)));
		}
	}	
}
/* ---------------------------------------------------------------------- */
inline void pauli2(double multi, int samespin, double rc, double re1, double re2, double* epauli,
	double* frc, double* fre1, double* fre2, double PAULI_RE, double PAULI_RC, double Kk, double rk, double pk)
{
	re1 *= PAULI_RE;
	re2 *= PAULI_RE;
	rc *= PAULI_RC;
	double reff = rc / pow((re1 * re1 + re2 * re2), 0.5);
	double r2 = pow(re1 / re2 + re2 / re1, 2);
	double r3 = (pow(re1, 2) + pow(re2, 2));
	double y = (re1 / re2 + re2 / re1) / 2;
	if ((reff * reff / rk) > 50.0)
	{
		*epauli += 0;
		*frc += 0;
		*fre1 += 0;
		*fre2 += 0;
	}
	else
	{
		if (samespin)
		{
		if (rc == 0 && y == 1)
			{
				*epauli += multi * Kk / (re1 * re2);
				*frc += 0;
				*fre1 += multi * Kk / (re1 * re1 * re1);
				*fre2 += multi * Kk / (re2 * re2 * re2);
			}
		else{
				*epauli += multi * Kk * (pk * (y * y - 1) + reff * reff / rk) / (y * y * y * exp(reff * reff / rk) - 1) / (re1 * re1 + re2 * re2);

				*frc += multi * PAULI_RC *(Kk * rc * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2) * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - 2 * Kk * rc / (rk * pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

				*fre1 += multi * PAULI_RE * (2 * Kk * re1 * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re1 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (3.0 / re2 - 3 * re2 / pow(re1, 2)) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (pk * (2.0 / re2 - 2 * re2 / pow(re1, 2)) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re1 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));

				*fre2 += multi * PAULI_RE * (2 * Kk * re2 * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) / (pow(r3, 2) * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)) - Kk * (pk * (r2 / 4.0 - 1) + pow(rc, 2) / (rk * r3)) * (pow(rc, 2) * re2 * pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / (4 * rk * pow(r3, 2)) - (-3 * re1 / pow(re2, 2) + 3.0 / re1) * r2 * exp(pow(rc, 2) / (rk * r3)) / 8.0) / (r3 * pow(pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1, 2)) - Kk * (pk * (-2 * re1 / pow(re2, 2) + 2 / re1) * (2 * y) / 4.0 - 2 * pow(rc, 2) * re2 / (rk * pow(r3, 2))) / (r3 * (pow(2 * y, 3) * exp(pow(rc, 2) / (rk * r3)) / 8.0 - 1)));
			}
		}
	}	
}
/*
inline void PauliElecElec(int samespin, double rc, double re1, double re2, double *epauli,
                          double *frc, double *fre1, double *fre2)
{
  double ree, rem;
  double S, t1, t2, tt;
  double dSdr1, dSdr2, dSdr;
  double dTdr1, dTdr2, dTdr;
  double O, dOdS, ratio;

  re1 *= PAULI_RE;
  re2 *= PAULI_RE;
  rc *= PAULI_RC;
  ree = re1 * re1 + re2 * re2;
  rem = re1 * re1 - re2 * re2;

  S = (2.82842712474619 / pow((re2 / re1 + re1 / re2), 1.5)) * exp(-rc * rc / ree);

  t1 = 1.5 * (1 / (re1 * re1) + 1 / (re2 * re2));
  t2 = 2.0 * (3 * ree - 2 * rc * rc) / (ree * ree);
  tt = t1 - t2;

  dSdr1 = (-1.5 / re1) * (rem / ree) + 2 * re1 * rc * rc / (ree * ree);
  dSdr2 = (1.5 / re2) * (rem / ree) + 2 * re2 * rc * rc / (ree * ree);
  dSdr = -2 * rc / ree;
  dTdr1 = -3 / (re1 * re1 * re1) - 12 * re1 / (ree * ree) +
      8 * re1 * (-2 * rc * rc + 3 * ree) / (ree * ree * ree);
  dTdr2 = -3 / (re2 * re2 * re2) - 12 * re2 / (ree * ree) +
      8 * re2 * (-2 * rc * rc + 3 * ree) / (ree * ree * ree);
  dTdr = 8 * rc / (ree * ree);

  if (samespin == 1) {
    O = S * S / (1.0 - S * S) + (1 - PAULI_RHO) * S * S / (1.0 + S * S);
    dOdS = 2 * S / ((1.0 - S * S) * (1.0 - S * S)) +
        (1 - PAULI_RHO) * 2 * S / ((1.0 + S * S) * (1.0 + S * S));
  } else {
    O = -PAULI_RHO * S * S / (1.0 + S * S);
    dOdS = -PAULI_RHO * 2 * S / ((1.0 + S * S) * (1.0 + S * S));
  }

  ratio = tt * dOdS * S;
  *fre1 -= PAULI_RE * (dTdr1 * O + ratio * dSdr1);
  *fre2 -= PAULI_RE * (dTdr2 * O + ratio * dSdr2);
  *frc -= PAULI_RC * (dTdr * O + ratio * dSdr);
  *epauli += tt * O;
}
*/
inline void SmallRForce(double dx, double dy, double dz, double rc, double force, double *fx,
                        double *fy, double *fz)
{
  /* Handles case where rc is small to avoid division by zero */

  if (rc > 1e-10) {
    force /= rc;
    *fx = force * dx;
    *fy = force * dy;
    *fz = force * dz;
  } else {
      *fx = 0.0;
      *fy = 0.0;
      *fz = 0.0;
    //                if (dx < 0) *fx = -*fx;
    //                if (dy < 0) *fy = -*fy;
    //                if (dz < 0) *fz = -*fz;
  }
}

/* ---------------------------------------------------------------------- */

inline double cutoff(double x)
{
  /*  cubic: return x * x * (2.0 * x - 3.0) + 1.0; */
  /*  quintic: return -6 * pow(x, 5) + 15 * pow(x, 4) - 10 * pow(x, 3) + 1; */

  /* Seventh order spline */
  //      return 20 * pow(x, 7) - 70 * pow(x, 6) + 84 * pow(x, 5) - 35 * pow(x, 4) + 1;
  return (((20 * x - 70) * x + 84) * x - 35) * x * x * x * x + 1;
}

/* ---------------------------------------------------------------------- */

inline double dcutoff(double x)
{
  /*  cubic: return (6.0 * x * x - 6.0 * x); */
  /*  quintic: return -30 * pow(x, 4) + 60 * pow(x, 3) - 30 * pow(x, 2); */

  /* Seventh order spline */
  //      return 140 * pow(x, 6) - 420 * pow(x, 5) + 420 * pow(x, 4) - 140 * pow(x, 3);
  return (((140 * x - 420) * x + 420) * x - 140) * x * x * x;
}

}    // namespace LAMMPS_NS
