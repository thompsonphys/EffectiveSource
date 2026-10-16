/*******************************************************************************
 * Analytic m-modes of the equatorial effective source and singular field from
 * the coefficients of effsource_equatorial_ctx_prepare.
 *
 * With dQ = sin(dphib/2), cQ = cos(dphib/2), u = dQ^2, psi = dphib and
 * s^2 = alpha + beta u, every output of effsource_equatorial_ctx_sweep_calc
 * that the m-modes need has the form
 *
 *   f = [P(u) + dQ cQ Q(u)] / s^(2 nu)
 *
 * with nu = 11/2 for src, 7/2 for PhiS and 9/2 for dPhiS/dr and dPhiS/dtheta.
 * The prepare functions of effsource_equatorial_sweep.h store the src
 * coefficients P, Q and alpha as double-double values (src_P + src_P_lo, ...).
 * The PhiS coefficients are the A_R0 and A_R1 fields, and those of the
 * derivatives follow exactly from them and from the dA_dr, dA_dth, dA_dQ and
 * dA_dR fields, all with their low parts. The m-modes
 *
 *   f_m = int_0^{2 pi} f e^{-i m phi} dphi
 *
 * follow from the kernels
 *
 *   K_mu^(nu) = int_0^{2 pi} cos(mu psi) (alpha + beta sin^2(psi/2))^(-nu) dpsi
 *
 * after one of two exact changes of basis of the numerator:
 *
 *   s^2 basis:      P(u) = sum_k e_k s^(2k),  2 u (1 - u) Q(u) = sum_k g_k s^(2k),
 *                   with the kernel orders nu - k,
 *   Fourier basis:  P(u) = sum_l p_l cos(l psi),  dQ cQ Q(u) = sum_l q_l sin(l psi),
 *                   with the kernel order nu only.
 *
 * For src, the contraction of the s^2 basis has a condition number below about
 * 3e3 for C1 = alpha / beta < 0.03 and m <= 40, and that of the Fourier basis
 * grows as C1^-5 (about 1e8 at C1 = 0.04). The kernels, the changes of basis
 * and the contractions use double-double arithmetic, so the result is the
 * exact integral of the coefficients to about 1e-16 relative to each mode.
 * With effsource_equatorial_ctx_prepare_double the low parts are zero,
 * and the rounding of the coefficients, amplified by the condition number of
 * the contraction, limits the small high-m modes (about 1e-10 relative to
 * each mode at C1 ~ 0.05 and m = 20, 1e-15 relative to the largest mode).
 ******************************************************************************/

#ifndef EFFSOURCE_EQUATORIAL_MMODE_H
#define EFFSOURCE_EQUATORIAL_MMODE_H

#include "effsource_equatorial_sweep.h"

/* Number of kernel orders nu = 11/2, 9/2, ..., -5/2 that the s^2 basis uses. */
#define EFFSOURCE_MMODE_NUM_NU 9

/* Largest mu for which the kernels are available. */
#define EFFSOURCE_MMODE_MAX_MU 256

/* Separation C1 = alpha / beta from which effsource_equatorial_ctx_analytic_modes
   uses the Fourier basis instead of the s^2 basis. Between C1 = 3e-3 and this
   value the Fourier basis is also used when the forward recurrence in mu is not
   accurate to m_max, which happens for m_max above about 20. */
#define EFFSOURCE_MMODE_FOURIER_C1 0.03

/* Kernels K[i * (mu_max + 1) + mu] = K_mu^(11/2 - i) for i = 0..num_nu-1 and
   mu = 0..mu_max, computed in double-double arithmetic and rounded to double.
   They come from the forward recurrence in mu if its error at mu_max is below
   about 1e-19, and from the backward (Miller) recurrence otherwise. Returns 0,
   1 if num_nu or mu_max is out of range, or 2 unless alpha > 0 and beta > 0. */
int effsource_equatorial_mmode_kernels(double alpha, double beta, int num_nu, int mu_max, double *K);

/* The m-modes of src, PhiS, dPhiS/dr and dPhiS/dtheta of effsource_equatorial_ctx_sweep_calc
   for m = 0..m_max, where the derivatives are taken at fixed t, theta and phi. Each output that
   is not NULL receives the real and imaginary parts of mode m at indices 2m and 2m + 1, the
   layout of an array of complex double. All outputs share one table of kernels, so each output
   after the first adds only its contraction. Returns 0, 1 if m_max + 7 > EFFSOURCE_MMODE_MAX_MU,
   or 2 if src_alpha <= 0, which happens only on the orbit of the particle. */
int effsource_equatorial_ctx_analytic_modes(const struct effsource_equatorial_ctx *ctx,
  const struct effsource_equatorial_point *point, int m_max, double *src, double *PhiS,
  double *dPhiS_dr, double *dPhiS_dth);

/* The m-modes S_m = re[m] + i im[m] of the src of effsource_equatorial_ctx_sweep_src
   for m = 0..m_max, with the return values of effsource_equatorial_ctx_analytic_modes. */
int effsource_equatorial_ctx_analytic_src_modes(const struct effsource_equatorial_ctx *ctx,
  const struct effsource_equatorial_point *point, int m_max, double *re, double *im);

#endif /* EFFSOURCE_EQUATORIAL_MMODE_H */
