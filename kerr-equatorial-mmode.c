/*******************************************************************************
 * Analytic m-modes of the equatorial effective source and singular field. See
 * effsource_equatorial_mmode.h for the formulas.
 *
 * All kernels, basis changes and contractions use double-double arithmetic, so
 * the m-modes are the exact integrals of the numerator coefficients to about
 * 1e-16 relative to each mode. The basis and the recurrence only change the
 * cost and the condition number that multiplies the double-double rounding.
 ******************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "effsource_dd.h"
#include "effsource_equatorial_mmode.h"

/* Largest number of coefficients of P(u), Q(u), G(u) = 2 u (1 - u) Q(u) and of
   the Fourier series of P(u) and sin(psi) Q(u) / 2 over all outputs. */
#define MAX_P 8
#define MAX_Q 7
#define MAX_G (MAX_Q + 2)
#define MAX_F (MAX_Q + 1)

/* Number of outputs: src, PhiS, dPhiS/dr and dPhiS/dtheta. */
#define NUM_OUTPUTS 4

/* Kernel orders nu = 11/2 - i for i = 0..NUM_NU-1 plus one lower order for K_1. */
#define NUM_J (EFFSOURCE_MMODE_NUM_NU + 1)

/* log10 of the relative error of the forward recurrence that is accepted. */
#define FORWARD_LOG10_ERROR (-19.0)

/* log10 of the relative error at which the Miller recurrence starts. */
#define MILLER_LOG10_START (-34.0)

/* Power of two above which the Miller recurrence rescales its values. */
#define MILLER_RESCALE 0x1p600

/* Separation below which the condition number of the Fourier contraction,
   about 1e8 (0.044 / C1)^5, exceeds 1e14 and the s^2 basis is used even where
   its kernels need the Miller recurrence. */
#define FOURIER_MIN_C1 3e-3

static const effsource_dd dd_pi = {3.141592653589793, 1.2246467991473532e-16};

/* a times a power of two. */
static inline effsource_dd dd_scale(effsource_dd a, double power_of_two)
{
  return (effsource_dd){a.hi * power_of_two, a.lo * power_of_two};
}

/* log(rho) with rho = z + sqrt(z^2 - 1) and z = 1 + 2 C1: the rate at which the
   kernels decrease with mu. */
static double log_rho(double C1)
{
  return log1p(2.0 * C1 + 2.0 * sqrt(C1 * (1.0 + C1)));
}

/* The forward recurrence in mu amplifies the double-double rounding by about
   rho^(2 mu) and a polynomial factor in mu. This fit to a 32-digit recurrence for
   nu = 11/2..-5/2 and C1 = 1e-3..0.3 returns whether the error at mu_max stays
   below 10^FORWARD_LOG10_ERROR. */
static int forward_is_accurate(double log_rho_value, int mu_max)
{
  const double polynomial = mu_max > 1 ? fmax(0.0, 7.5 * log10((double)mu_max) - 5.5) : 0.0;
  return -30.0 + 2.0 * mu_max * log_rho_value / M_LN10 + polynomial <= FORWARD_LOG10_ERROR;
}

/* Complete elliptic integrals K(k) and E(k) with the complementary parameter
   m1 = 1 - k^2, by the arithmetic-geometric mean. */
static void elliptic_KE(effsource_dd m1, effsource_dd *K, effsource_dd *E)
{
  effsource_dd a = dd_d(1.0), b = dd_sqrt(m1);
  effsource_dd S = dd_mul_d(dd_sub(dd_d(1.0), m1), 0.5);
  double weight = 1.0;
  for (int n = 0; n < 64; n++) {
    const effsource_dd c = dd_mul_d(dd_sub(a, b), 0.5);
    const effsource_dd a_next = dd_mul_d(dd_add(a, b), 0.5);
    b = dd_sqrt(dd_mul(a, b));
    a = a_next;
    S = dd_add(S, dd_mul_d(dd_mul(c, c), weight));
    weight *= 2.0;
    if (fabs(c.hi) <= 1e-30 * a.hi)
      break;
  }
  *K = dd_div(dd_pi, dd_mul_d(a, 2.0));
  *E = dd_mul(*K, dd_sub(dd_d(1.0), S));
}

/* J[q] = K_0^(q - 7/2) for q = 0..9, from
   K_0^(1/2) = 4 K(k) / sqrt(alpha + beta), K_0^(-1/2) = 4 sqrt(alpha + beta) E(k),
   k^2 = beta / (alpha + beta), and the recurrence in nu
   p D K_0^(p+1) = (2p - 1) A K_0^(p) - (p - 1) K_0^(p-1), D = A^2 - B^2,
   upward to 11/2 and downward to -7/2. Both directions follow the dominant
   solution, so both are stable. */
static void zero_mode_kernels(effsource_dd alpha, double beta, effsource_dd A, effsource_dd *J)
{
  const effsource_dd alpha_beta = dd_add(alpha, dd_d(beta));
  const effsource_dd D = dd_mul(alpha, alpha_beta);
  const effsource_dd root = dd_sqrt(alpha_beta);
  effsource_dd K, E;
  elliptic_KE(dd_div(alpha, alpha_beta), &K, &E);
  J[4] = dd_div(dd_mul_d(K, 4.0), root);
  J[3] = dd_mul_d(dd_mul(root, E), 4.0);
  for (int q = 4; q < NUM_J - 1; q++) {
    const double p = q - 3.5;
    const effsource_dd numerator = dd_sub(dd_mul_d(dd_mul(A, J[q]), 2.0 * p - 1.0), dd_mul_d(J[q - 1], p - 1.0));
    J[q + 1] = dd_div(numerator, dd_mul_d(D, p));
  }
  for (int q = 3; q > 0; q--) {
    const double p = q - 3.5;
    const effsource_dd numerator = dd_sub(dd_mul_d(dd_mul(A, J[q]), 2.0 * p - 1.0), dd_mul_d(dd_mul(D, J[q + 1]), p));
    J[q - 1] = dd_div(numerator, dd_d(p - 1.0));
  }
}

/* The kernels in double-double, same layout as effsource_equatorial_mmode_kernels.
   The recurrence in mu is
   (mu + 1 - nu) K_(mu+1) - 2 mu z K_mu + (mu + nu - 1) K_(mu-1) = 0, z = A / B.
   The kernels are its minimal solution, so the forward direction loses about
   rho^(2 mu) relative to the double-double rounding. If that keeps the error
   small, the forward recurrence starts from K_0 and K_1 = (A K_0^(nu) - K_0^(nu-1)) / B.
   Otherwise the Miller recurrence runs down from zero and one at a guard order
   where the dominant solution is 10^MILLER_LOG10_START smaller, and the result
   is normalised to K_0. */
static void kernels_dd(effsource_dd alpha, double beta, int num_nu, int mu_max, effsource_dd *K)
{
  const double B = 0.5 * beta;
  const effsource_dd A = dd_add(alpha, dd_d(B));
  const effsource_dd z = dd_div_d(A, B);
  effsource_dd J[NUM_J];
  zero_mode_kernels(alpha, beta, A, J);
  const double log_rho_value = log_rho(alpha.hi / beta);
  const int stride = mu_max + 1;
  if (forward_is_accurate(log_rho_value, mu_max)) {
    effsource_dd two_mu_z[EFFSOURCE_MMODE_MAX_MU + 1];
    for (int mu = 1; mu < mu_max; mu++)
      two_mu_z[mu] = dd_mul_d(z, 2.0 * mu);
    for (int i = 0; i < num_nu; i++) {
      const double nu = 5.5 - i;
      effsource_dd *k = K + i * stride;
      k[0] = J[NUM_J - 1 - i];
      if (mu_max < 1)
        continue;
      k[1] = dd_div_d(dd_sub(dd_mul(A, k[0]), J[NUM_J - 2 - i]), B);
      for (int mu = 1; mu < mu_max; mu++) {
        const effsource_dd sum = dd_sub(dd_mul(two_mu_z[mu], k[mu]), dd_mul_d(k[mu - 1], mu + nu - 1.0));
        k[mu + 1] = dd_div_d(sum, mu + 1.0 - nu);
      }
    }
    return;
  }
  const int guard = (int)ceil(-MILLER_LOG10_START * M_LN10 / (2.0 * log_rho_value)) + 20;
  for (int i = 0; i < num_nu; i++) {
    const double nu = 5.5 - i;
    effsource_dd *k = K + i * stride;
    effsource_dd above = dd_d(0.0), current = dd_d(1.0);
    for (int mu = mu_max + guard; mu >= 1; mu--) {
      const effsource_dd sum = dd_sub(dd_mul(dd_mul_d(z, 2.0 * mu), current), dd_mul_d(above, mu + 1.0 - nu));
      above = current;
      current = dd_div_d(sum, mu + nu - 1.0);
      if (mu - 1 <= mu_max)
        k[mu - 1] = current;
      if (fabs(current.hi) > MILLER_RESCALE) {
        above = dd_scale(above, 1.0 / MILLER_RESCALE);
        current = dd_scale(current, 1.0 / MILLER_RESCALE);
        for (int stored = mu - 1; stored <= mu_max; stored++)
          k[stored] = dd_scale(k[stored], 1.0 / MILLER_RESCALE);
      }
    }
    const effsource_dd normalization = dd_div(J[NUM_J - 1 - i], k[0]);
    k[0] = J[NUM_J - 1 - i];
    for (int mu = 1; mu <= mu_max; mu++)
      k[mu] = dd_mul(k[mu], normalization);
  }
}

int effsource_equatorial_mmode_kernels(double alpha, double beta, int num_nu, int mu_max, double *K)
{
  if (num_nu < 1 || num_nu > EFFSOURCE_MMODE_NUM_NU || mu_max < 0 || mu_max > EFFSOURCE_MMODE_MAX_MU)
    return 1;
  if (!(alpha > 0.0 && beta > 0.0))
    return 2;
  effsource_dd kernels[EFFSOURCE_MMODE_NUM_NU * (EFFSOURCE_MMODE_MAX_MU + 1)];
  kernels_dd(dd_d(alpha), beta, num_nu, mu_max, kernels);
  for (int n = 0; n < num_nu * (mu_max + 1); n++)
    K[n] = dd_to_double(kernels[n]);
  return 0;
}

/* The numerator of one output, [P(u) + dQ cQ Q(u)] / s^(2 nu) with nu = 11/2 - top. */
struct numerator {
  int top, num_P, num_Q;
  effsource_dd P[MAX_P], Q[MAX_Q];
};

/* Coefficients of p(w - alpha) in w from those of p(v), in place. */
static void taylor_shift(effsource_dd alpha, int size, effsource_dd *c)
{
  const effsource_dd shift = dd_neg(alpha);
  for (int i = 0; i < size - 1; i++)
    for (int j = size - 2; j >= i; j--)
      c[j] = dd_add(c[j], dd_mul(shift, c[j + 1]));
}

/* Coefficients of sum_j c[j] u^j in powers of s^2 = alpha + beta u, in place. */
static void to_s2_basis(effsource_dd alpha, double beta, int size, effsource_dd *c)
{
  const effsource_dd inverse_beta = dd_div(dd_d(1.0), dd_d(beta));
  effsource_dd scale = dd_d(1.0);
  for (int j = 1; j < size; j++) {
    scale = dd_mul(scale, inverse_beta);
    c[j] = dd_mul(c[j], scale);
  }
  taylor_shift(alpha, size, c);
}

/* Coefficient of cos(l psi) in u^j = sin^2j(psi/2):
   4^-j binomial(2j, j) for l = 0 and 2 (-1)^l 4^-j binomial(2j, j - l) for l > 0. */
static double fourier_coefficient(int j, int l)
{
  double binomial = 1.0;
  for (int i = 1; i <= j - l; i++)
    binomial = binomial * (2 * j - i + 1) / i;
  const double value = ldexp(binomial, -2 * j);
  return l == 0 ? value : (l % 2 ? -2.0 : 2.0) * value;
}

/* even[m] = int cos(m psi) P(u) / s^(2 nu) dpsi and odd[m] = int sin(m psi) sin(psi) Q(u) / (2 s^(2 nu)) dpsi
   in the s^2 basis, with sin(m psi) / sin(psi) = sum_l cos((m - 1 - 2l) psi). K holds the kernels
   K[i * stride + mu] = K_mu^(11/2 - i) for i = 0..n->top + n->num_Q + 1 and mu = 0..m_max. */
static void s2_contraction(const struct numerator *n, effsource_dd alpha, double beta, const effsource_dd *K,
                           int stride, int m_max, double *even, double *odd)
{
  const int num_G = n->num_Q + 2;
  effsource_dd e[MAX_P], g[MAX_G];
  for (int j = 0; j < n->num_P; j++)
    e[j] = n->P[j];
  for (int j = 0; j < num_G; j++) {
    const effsource_dd up = j >= 1 && j - 1 < n->num_Q ? n->Q[j - 1] : dd_d(0.0);
    const effsource_dd down = j >= 2 && j - 2 < n->num_Q ? n->Q[j - 2] : dd_d(0.0);
    g[j] = dd_mul_d(dd_sub(up, down), 2.0);
  }
  to_s2_basis(alpha, beta, n->num_P, e);
  to_s2_basis(alpha, beta, num_G, g);
  const effsource_dd *k_top = K + n->top * stride;
  effsource_dd L_previous[MAX_G], L_current[MAX_G];
  for (int k = 0; k < num_G; k++) {
    L_previous[k] = dd_d(0.0);
    L_current[k] = dd_d(0.0);
  }
  for (int m = 0; m <= m_max; m++) {
    effsource_dd even_sum = dd_d(0.0), odd_sum = dd_d(0.0);
    for (int k = 0; k < n->num_P; k++)
      even_sum = dd_add(even_sum, dd_mul(e[k], k_top[k * stride + m]));
    for (int k = 0; k < num_G; k++) {
      if (m == 1)
        L_current[k] = k_top[k * stride];
      else if (m >= 2) {
        const effsource_dd L_next = dd_add(L_previous[k], dd_mul_d(k_top[k * stride + m - 1], 2.0));
        L_previous[k] = L_current[k];
        L_current[k] = L_next;
      }
      odd_sum = dd_add(odd_sum, dd_mul(g[k], L_current[k]));
    }
    even[m] = dd_to_double(even_sum);
    odd[m] = dd_to_double(odd_sum);
  }
}

/* The same integrals in the Fourier basis P(u) = sum_l p_l cos(l psi) and
   sin(psi) Q(u) / 2 = sum_l q_l sin(l psi), with the kernels k[mu] = K_mu^(nu) of nu = 11/2 - top
   only, for mu = 0..m_max + MAX_F - 1:
   even[m] = sum_l p_l (K_(m+l) + K_|m-l|) / 2, odd[m] = sum_l q_l (K_|m-l| - K_(m+l)) / 2. */
static void fourier_contraction(const struct numerator *n, const effsource_dd *k, int m_max, double *even,
                                double *odd)
{
  effsource_dd p[MAX_F], q[MAX_F], Qc[MAX_Q];
  for (int l = 0; l < MAX_F; l++) {
    p[l] = dd_d(0.0);
    q[l] = dd_d(0.0);
  }
  for (int l = 0; l < n->num_P; l++)
    for (int j = l; j < n->num_P; j++)
      p[l] = dd_add(p[l], dd_mul_d(n->P[j], fourier_coefficient(j, l)));
  for (int l = 0; l < n->num_Q; l++) {
    Qc[l] = dd_d(0.0);
    for (int j = l; j < n->num_Q; j++)
      Qc[l] = dd_add(Qc[l], dd_mul_d(n->Q[j], fourier_coefficient(j, l)));
  }
  if (n->num_Q > 0)
    q[1] = dd_scale(Qc[0], 0.5);
  for (int l = 1; l < n->num_Q; l++) {
    q[l + 1] = dd_add(q[l + 1], dd_scale(Qc[l], 0.25));
    if (l >= 2)
      q[l - 1] = dd_sub(q[l - 1], dd_scale(Qc[l], 0.25));
  }
  const int num_F = n->num_P > n->num_Q + 1 ? n->num_P : n->num_Q + 1;
  for (int l = 1; l < num_F; l++) {
    p[l] = dd_scale(p[l], 0.5);
    q[l] = dd_scale(q[l], 0.5);
  }
  for (int m = 0; m <= m_max; m++) {
    effsource_dd even_sum = dd_mul(p[0], k[m]), odd_sum = dd_d(0.0);
    for (int l = 1; l < num_F; l++) {
      const effsource_dd low = k[abs(m - l)], high = k[m + l];
      even_sum = dd_add(even_sum, dd_mul(p[l], dd_add(low, high)));
      odd_sum = dd_add(odd_sum, dd_mul(q[l], dd_sub(low, high)));
    }
    even[m] = dd_to_double(even_sum);
    odd[m] = dd_to_double(odd_sum);
  }
}

static void clear(struct numerator *n, int top, int num_P, int num_Q)
{
  n->top = top;
  n->num_P = num_P;
  n->num_Q = num_Q;
  for (int j = 0; j < MAX_P; j++)
    n->P[j] = dd_d(0.0);
  for (int j = 0; j < MAX_Q; j++)
    n->Q[j] = dd_d(0.0);
}

/* r(u) += f u^shift a(u) for the polynomial a of size n. */
static void add_scaled(effsource_dd *r, const effsource_dd *a, int n, effsource_dd f, int shift)
{
  for (int j = 0; j < n; j++)
    r[j + shift] = dd_add(r[j + shift], dd_mul(f, a[j]));
}

/* The PhiS-family coefficients of the point structure as double-double polynomials. */
struct field_coefficients {
  effsource_dd A_R0[5], A_R1[5], dA_dr_R0[5], dA_dr_R1[4], dA_dth_R0[4], dA_dth_R1[4];
  effsource_dd dA_dQ_R0[4], dA_dQ_R1[4], dA_dR_R0[5];
};

static void load(const double *hi, const double *lo, int n, effsource_dd *poly)
{
  for (int j = 0; j < n; j++)
    poly[j] = (effsource_dd){hi[j], lo[j]};
}

static void load_field_coefficients(const struct effsource_equatorial_point *point, struct field_coefficients *f)
{
  load(point->A_R0, point->A_R0_lo, 5, f->A_R0);
  load(point->A_R1, point->A_R1_lo, 5, f->A_R1);
  load(point->dA_dr_R0, point->dA_dr_R0_lo, 5, f->dA_dr_R0);
  load(point->dA_dr_R1, point->dA_dr_R1_lo, 4, f->dA_dr_R1);
  load(point->dA_dth_R0, point->dA_dth_R0_lo, 4, f->dA_dth_R0);
  load(point->dA_dth_R1, point->dA_dth_R1_lo, 4, f->dA_dth_R1);
  load(point->dA_dQ_R0, point->dA_dQ_R0_lo, 4, f->dA_dQ_R0);
  load(point->dA_dQ_R1, point->dA_dQ_R1_lo, 4, f->dA_dQ_R1);
  load(point->dA_dR_R0, point->dA_dR_R0_lo, 5, f->dA_dR_R0);
}

/* The src numerator over s^11. */
static void src_numerator(const struct effsource_equatorial_point *point, struct numerator *n)
{
  clear(n, 0, 8, 7);
  load(point->src_P, point->src_P_lo, 8, n->P);
  load(point->src_Q, point->src_Q_lo, 7, n->Q);
}

/* The PhiS numerator over s^7: A = P_A(u) + dR Q_A(u) with dR = 2 dQ cQ. */
static void PhiS_numerator(const struct field_coefficients *f, struct numerator *n)
{
  clear(n, 2, 5, 5);
  add_scaled(n->P, f->A_R0, 5, dd_d(1.0), 0);
  add_scaled(n->Q, f->A_R1, 5, dd_d(2.0), 0);
}

/* The dPhiS/dtheta numerator over s^9: s^2 dA/dtheta - 7 alpha02 dtheta A. */
static void dPhiS_dth_numerator(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                                const struct field_coefficients *f, effsource_dd alpha, struct numerator *n)
{
  clear(n, 1, 5, 5);
  const effsource_dd beta = dd_d(ctx->beta);
  const effsource_dd weight = dd_mul_d(dd_two_sum(point->theta, -ctx->xp.theta), -7.0 * ctx->alpha02);
  add_scaled(n->P, f->dA_dth_R0, 4, alpha, 0);
  add_scaled(n->P, f->dA_dth_R0, 4, beta, 1);
  add_scaled(n->P, f->A_R0, 5, weight, 0);
  add_scaled(n->Q, f->dA_dth_R1, 4, dd_mul_d(alpha, 2.0), 0);
  add_scaled(n->Q, f->dA_dth_R1, 4, dd_mul_d(beta, 2.0), 1);
  add_scaled(n->Q, f->A_R1, 5, dd_mul_d(weight, 2.0), 0);
}

/* The dPhiS/dr numerator over s^9: s^2 DA/Dr - (7/2) A Ds^2/Dr, where the total r
   derivative includes d(dphib)/dr = -c, so that d(dQ)/dr = -c cQ / 2, d(dR)/dr = -c (1 - 2u)
   and Ds^2/Dr = 2 alpha20 dr - c beta dQ cQ. With the partial derivatives in the point structure,
   DA/Dr = D_even(u) + dQ cQ D_odd(u) with
   D_even = dA_dr_R0 - c u (1 - u) dA_dQ_R1 - c (1 - 2u) dA_dR_R0 and D_odd = 2 dA_dr_R1 - c dA_dQ_R0 / 2. */
static void dPhiS_dr_numerator(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                               const struct field_coefficients *f, effsource_dd alpha, struct numerator *n)
{
  clear(n, 1, 7, 5);
  const double c = ctx->c;
  const effsource_dd beta = dd_d(ctx->beta);
  effsource_dd D_even[6], D_odd[4];
  for (int j = 0; j < 6; j++)
    D_even[j] = dd_d(0.0);
  for (int j = 0; j < 4; j++)
    D_odd[j] = dd_d(0.0);
  add_scaled(D_even, f->dA_dr_R0, 5, dd_d(1.0), 0);
  add_scaled(D_even, f->dA_dQ_R1, 4, dd_d(-c), 1);
  add_scaled(D_even, f->dA_dQ_R1, 4, dd_d(c), 2);
  add_scaled(D_even, f->dA_dR_R0, 5, dd_d(-c), 0);
  add_scaled(D_even, f->dA_dR_R0, 5, dd_d(2.0 * c), 1);
  add_scaled(D_odd, f->dA_dr_R1, 4, dd_d(2.0), 0);
  add_scaled(D_odd, f->dA_dQ_R0, 4, dd_d(-0.5 * c), 0);
  const effsource_dd radial = dd_mul_d(dd_two_sum(point->r, -ctx->xp.r), -7.0 * ctx->alpha20);
  const effsource_dd angular = dd_mul_d(beta, 3.5 * c);
  add_scaled(n->P, D_even, 6, alpha, 0);
  add_scaled(n->P, D_even, 6, beta, 1);
  add_scaled(n->P, f->A_R0, 5, radial, 0);
  add_scaled(n->P, f->A_R1, 5, dd_mul_d(angular, 2.0), 1);
  add_scaled(n->P, f->A_R1, 5, dd_mul_d(angular, -2.0), 2);
  add_scaled(n->Q, D_odd, 4, alpha, 0);
  add_scaled(n->Q, D_odd, 4, beta, 1);
  add_scaled(n->Q, f->A_R1, 5, dd_mul_d(radial, 2.0), 0);
  add_scaled(n->Q, f->A_R0, 5, angular, 0);
}

/* The m-modes of the given numerators: even and odd integrals from one kernel table, then the
   phase exp(-i m phi0) with phi0 = phi_p + c dr, written to out[i][2m] and out[i][2m + 1]. */
static void modes(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                  effsource_dd alpha, int m_max, int count, const struct numerator *numerators, double *const *out)
{
  const double beta = ctx->beta;
  const double C1 = point->src_alpha / beta;
  const int s2_basis = C1 < EFFSOURCE_MMODE_FOURIER_C1 && (C1 < FOURIER_MIN_C1 || forward_is_accurate(log_rho(C1), m_max));
  int num_nu = 0;
  for (int i = 0; i < count; i++) {
    const int needed = s2_basis ? numerators[i].top + numerators[i].num_Q + 2 : numerators[i].top + 1;
    num_nu = needed > num_nu ? needed : num_nu;
  }
  const int mu_max = s2_basis ? m_max : m_max + MAX_F - 1;
  effsource_dd K[EFFSOURCE_MMODE_NUM_NU * (EFFSOURCE_MMODE_MAX_MU + 1)];
  kernels_dd(alpha, beta, num_nu, mu_max, K);
  double even[NUM_OUTPUTS][EFFSOURCE_MMODE_MAX_MU + 1], odd[NUM_OUTPUTS][EFFSOURCE_MMODE_MAX_MU + 1];
  for (int i = 0; i < count; i++) {
    if (s2_basis)
      s2_contraction(&numerators[i], alpha, beta, K, mu_max + 1, m_max, even[i], odd[i]);
    else
      fourier_contraction(&numerators[i], K + numerators[i].top * (mu_max + 1), m_max, even[i], odd[i]);
  }
  const double phi0 = ctx->xp.phi + ctx->c * (point->r - ctx->xp.r);
  for (int m = 0; m <= m_max; m++) {
    const effsource_dd angle = dd_two_prod(m, phi0);
    const double c_hi = cos(angle.hi), s_hi = sin(angle.hi);
    const double c = c_hi - angle.lo * s_hi, s = s_hi + angle.lo * c_hi;
    for (int i = 0; i < count; i++) {
      out[i][2 * m] = c * even[i][m] - s * odd[i][m];
      out[i][2 * m + 1] = -(s * even[i][m] + c * odd[i][m]);
    }
  }
}

int effsource_equatorial_ctx_analytic_modes(const struct effsource_equatorial_ctx *ctx,
                                         const struct effsource_equatorial_point *point, int m_max, double *src,
                                         double *PhiS, double *dPhiS_dr, double *dPhiS_dth)
{
  if (m_max < 0 || m_max + MAX_F - 1 > EFFSOURCE_MMODE_MAX_MU)
    return 1;
  if (!(point->src_alpha > 0.0))
    return 2;
  const effsource_dd alpha = {point->src_alpha, point->src_alpha_lo};
  struct numerator numerators[NUM_OUTPUTS];
  double *out[NUM_OUTPUTS];
  int count = 0;
  struct field_coefficients field;
  if (PhiS || dPhiS_dr || dPhiS_dth)
    load_field_coefficients(point, &field);
  if (src) {
    src_numerator(point, &numerators[count]);
    out[count++] = src;
  }
  if (PhiS) {
    PhiS_numerator(&field, &numerators[count]);
    out[count++] = PhiS;
  }
  if (dPhiS_dr) {
    dPhiS_dr_numerator(ctx, point, &field, alpha, &numerators[count]);
    out[count++] = dPhiS_dr;
  }
  if (dPhiS_dth) {
    dPhiS_dth_numerator(ctx, point, &field, alpha, &numerators[count]);
    out[count++] = dPhiS_dth;
  }
  if (count > 0)
    modes(ctx, point, alpha, m_max, count, numerators, out);
  return 0;
}

int effsource_equatorial_ctx_analytic_src_modes(const struct effsource_equatorial_ctx *ctx,
                                             const struct effsource_equatorial_point *point, int m_max,
                                             double *re, double *im)
{
  double src[2 * (EFFSOURCE_MMODE_MAX_MU + 1)];
  const int status = effsource_equatorial_ctx_analytic_modes(ctx, point, m_max, src, NULL, NULL, NULL);
  if (status != 0)
    return status;
  for (int m = 0; m <= m_max; m++) {
    re[m] = src[2 * m];
    im[m] = src[2 * m + 1];
  }
  return 0;
}
