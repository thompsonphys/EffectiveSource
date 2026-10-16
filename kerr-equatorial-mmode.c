/*******************************************************************************
 * Analytic m-modes of the equatorial effective source. See
 * effsource_equatorial_mmode.h for the formulas.
 *
 * All kernels, basis changes and contractions use double-double arithmetic, so
 * the m-modes are the exact integrals of the src coefficients in the sweep
 * structure to about 1e-16 relative to each mode. The basis and the recurrence
 * only change the cost and the condition number that multiplies the
 * double-double rounding.
 ******************************************************************************/

#include <math.h>
#include <stdlib.h>

#include "effsource_dd.h"
#include "effsource_equatorial_mmode.h"

/* Number of coefficients of P(u), Q(u) and G(u) = 2 u (1 - u) Q(u). */
#define NUM_P 8
#define NUM_Q 7
#define NUM_G 9

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

/* even[m] = int cos(m psi) P(u) / s^11 dpsi and odd[m] = int sin(m psi) sin(psi) Q(u) / (2 s^11) dpsi
   in the s^2 basis, with sin(m psi) / sin(psi) = sum_l cos((m - 1 - 2l) psi). */
static void s2_contraction(effsource_dd alpha, double beta, const effsource_dd *P, const effsource_dd *Q,
                           int m_max, double *even, double *odd)
{
  effsource_dd e[NUM_P], g[NUM_G];
  for (int j = 0; j < NUM_P; j++)
    e[j] = P[j];
  for (int j = 0; j < NUM_G; j++) {
    const effsource_dd up = j >= 1 && j - 1 < NUM_Q ? Q[j - 1] : dd_d(0.0);
    const effsource_dd down = j >= 2 && j - 2 < NUM_Q ? Q[j - 2] : dd_d(0.0);
    g[j] = dd_mul_d(dd_sub(up, down), 2.0);
  }
  to_s2_basis(alpha, beta, NUM_P, e);
  to_s2_basis(alpha, beta, NUM_G, g);
  effsource_dd K[EFFSOURCE_MMODE_NUM_NU * (EFFSOURCE_MMODE_MAX_MU + 1)];
  kernels_dd(alpha, beta, NUM_G, m_max, K);
  const int stride = m_max + 1;
  effsource_dd L_previous[NUM_G], L_current[NUM_G];
  for (int k = 0; k < NUM_G; k++) {
    L_previous[k] = dd_d(0.0);
    L_current[k] = dd_d(0.0);
  }
  for (int m = 0; m <= m_max; m++) {
    effsource_dd even_sum = dd_d(0.0), odd_sum = dd_d(0.0);
    for (int k = 0; k < NUM_P; k++)
      even_sum = dd_add(even_sum, dd_mul(e[k], K[k * stride + m]));
    for (int k = 0; k < NUM_G; k++) {
      if (m == 1)
        L_current[k] = K[k * stride];
      else if (m >= 2) {
        const effsource_dd L_next = dd_add(L_previous[k], dd_mul_d(K[k * stride + m - 1], 2.0));
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
   sin(psi) Q(u) / 2 = sum_l q_l sin(l psi), with the kernels for nu = 11/2 only:
   even[m] = sum_l p_l (K_(m+l) + K_|m-l|) / 2, odd[m] = sum_l q_l (K_|m-l| - K_(m+l)) / 2. */
static void fourier_contraction(effsource_dd alpha, double beta, const effsource_dd *P, const effsource_dd *Q,
                                int m_max, double *even, double *odd)
{
  effsource_dd p[NUM_P], q[NUM_P], Qc[NUM_Q];
  for (int l = 0; l < NUM_P; l++) {
    p[l] = dd_d(0.0);
    q[l] = dd_d(0.0);
    for (int j = l; j < NUM_P; j++)
      p[l] = dd_add(p[l], dd_mul_d(P[j], fourier_coefficient(j, l)));
  }
  for (int l = 0; l < NUM_Q; l++) {
    Qc[l] = dd_d(0.0);
    for (int j = l; j < NUM_Q; j++)
      Qc[l] = dd_add(Qc[l], dd_mul_d(Q[j], fourier_coefficient(j, l)));
  }
  q[1] = dd_scale(Qc[0], 0.5);
  for (int l = 1; l < NUM_Q; l++) {
    q[l + 1] = dd_add(q[l + 1], dd_scale(Qc[l], 0.25));
    if (l >= 2)
      q[l - 1] = dd_sub(q[l - 1], dd_scale(Qc[l], 0.25));
  }
  for (int l = 1; l < NUM_P; l++) {
    p[l] = dd_scale(p[l], 0.5);
    q[l] = dd_scale(q[l], 0.5);
  }
  const int mu_max = m_max + NUM_P - 1;
  effsource_dd K[EFFSOURCE_MMODE_MAX_MU + 1];
  kernels_dd(alpha, beta, 1, mu_max, K);
  for (int m = 0; m <= m_max; m++) {
    effsource_dd even_sum = dd_mul(p[0], K[m]), odd_sum = dd_d(0.0);
    for (int l = 1; l < NUM_P; l++) {
      const effsource_dd low = K[abs(m - l)], high = K[m + l];
      even_sum = dd_add(even_sum, dd_mul(p[l], dd_add(low, high)));
      odd_sum = dd_add(odd_sum, dd_mul(q[l], dd_sub(low, high)));
    }
    even[m] = dd_to_double(even_sum);
    odd[m] = dd_to_double(odd_sum);
  }
}

int effsource_equatorial_ctx_sweep_src_modes(const struct effsource_equatorial_ctx *ctx,
                                             const struct effsource_equatorial_sweep *sweep, int m_max,
                                             double *re, double *im)
{
  if (m_max < 0 || m_max + NUM_P - 1 > EFFSOURCE_MMODE_MAX_MU)
    return 1;
  if (!(sweep->src_alpha > 0.0))
    return 2;
  const double beta = ctx->beta;
  const effsource_dd alpha = {sweep->src_alpha, sweep->src_alpha_lo};
  effsource_dd P[NUM_P], Q[NUM_Q];
  for (int j = 0; j < NUM_P; j++)
    P[j] = (effsource_dd){sweep->src_P[j], sweep->src_P_lo[j]};
  for (int j = 0; j < NUM_Q; j++)
    Q[j] = (effsource_dd){sweep->src_Q[j], sweep->src_Q_lo[j]};
  const double C1 = sweep->src_alpha / beta;
  double even[EFFSOURCE_MMODE_MAX_MU + 1], odd[EFFSOURCE_MMODE_MAX_MU + 1];
  if (C1 < EFFSOURCE_MMODE_FOURIER_C1 && (C1 < FOURIER_MIN_C1 || forward_is_accurate(log_rho(C1), m_max)))
    s2_contraction(alpha, beta, P, Q, m_max, even, odd);
  else
    fourier_contraction(alpha, beta, P, Q, m_max, even, odd);
  const double phi0 = ctx->xp.phi + ctx->c * (sweep->r - ctx->xp.r);
  for (int m = 0; m <= m_max; m++) {
    const effsource_dd angle = dd_two_prod(m, phi0);
    const double c_hi = cos(angle.hi), s_hi = sin(angle.hi);
    const double c = c_hi - angle.lo * s_hi, s = s_hi + angle.lo * c_hi;
    re[m] = c * even[m] - s * odd[m];
    im[m] = -(s * even[m] + c * odd[m]);
  }
  return 0;
}
