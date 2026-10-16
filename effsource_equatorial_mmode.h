/*******************************************************************************
 * Analytic m-modes of the equatorial effective source from the sweep form.
 *
 * With dQ = sin(dphib/2), cQ = cos(dphib/2), u = dQ^2, psi = dphib and
 * s^2 = alpha + beta u, the prepare functions of effsource_equatorial_sweep.h
 * write the source as
 *
 *   src = [P(u) + dQ cQ Q(u)] / s^11
 *
 * and store P, Q and alpha as double-double values (src_P + src_P_lo, ...).
 * The m-modes
 *
 *   S_m = int_0^{2 pi} src e^{-i m phi} dphi
 *
 * follow from the kernels
 *
 *   K_mu^(nu) = int_0^{2 pi} cos(mu psi) (alpha + beta sin^2(psi/2))^(-nu) dpsi
 *
 * after one of two exact changes of basis of the numerator:
 *
 *   s^2 basis:      P(u) = sum_k e_k s^(2k),  2 u (1 - u) Q(u) = sum_k g_k s^(2k),
 *                   with nu = 11/2 - k,
 *   Fourier basis:  P(u) = sum_l p_l cos(l psi),  dQ cQ Q(u) = sum_l q_l sin(l psi),
 *                   with nu = 11/2 only.
 *
 * The contraction of the s^2 basis has a condition number below about 3e3 for
 * C1 = alpha / beta < 0.03 and m <= 40, and that of the Fourier basis grows as
 * C1^-5 (about 1e8 at C1 = 0.04). The kernels, the changes of basis and the
 * contractions use double-double arithmetic, so the result is the exact
 * integral of the stored coefficients to about 1e-16 relative to each mode.
 ******************************************************************************/

#ifndef EFFSOURCE_EQUATORIAL_MMODE_H
#define EFFSOURCE_EQUATORIAL_MMODE_H

#include "effsource_equatorial_sweep.h"

/* Number of kernel orders nu = 11/2, 9/2, ..., -5/2 that the s^2 basis uses. */
#define EFFSOURCE_MMODE_NUM_NU 9

/* Largest mu for which the kernels are available. */
#define EFFSOURCE_MMODE_MAX_MU 256

/* Separation C1 = alpha / beta from which effsource_equatorial_ctx_sweep_src_modes
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

/* The m-modes S_m = re[m] + i im[m] of the src of effsource_equatorial_ctx_sweep_src
   for m = 0..m_max. Returns 0, 1 if m_max + 7 > EFFSOURCE_MMODE_MAX_MU, or 2 if
   src_alpha <= 0, which happens only on the orbit of the particle. */
int effsource_equatorial_ctx_sweep_src_modes(const struct effsource_equatorial_ctx *ctx,
  const struct effsource_equatorial_sweep *sweep, int m_max, double *re, double *im);

#endif /* EFFSOURCE_EQUATORIAL_MMODE_H */
