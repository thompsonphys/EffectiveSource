/*******************************************************************************
 * Check effsource_equatorial_mmode_kernels and effsource_equatorial_ctx_analytic_src_modes.
 *
 *   1. The kernels K_mu^(nu) for nu = 11/2..-5/2 and mu = 0..KERNEL_MU_MAX are
 *      compared with a long double trapezoid rule, which converges
 *      exponentially for the periodic integrand. The error is normalised by
 *      K_0^(nu), the integral of the absolute value of the integrand: the
 *      reference loses relative accuracy where K_mu is far below K_0.
 *   2. The analytic m-modes from effsource_equatorial_ctx_prepare_dd are
 *      compared with the trapezoid rule applied to the src of
 *      effsource_equatorial_ctx_sweep_src with the same coefficients, on a grid
 *      fine enough that aliasing is below 1e-20. The trapezoid rule evaluates
 *      the rational form in long double on nodes of dphib in (-pi, pi], so that
 *      the distance to the peak at dphib = 0 is exact. Both routes integrate the
 *      same function, so the difference measures the kernels and contractions.
 *      The long double evaluation rounds by about 1e-19 I, where
 *      I = int (sum_j |P_j| u^j + |dQ cQ| sum_j |Q_j| u^j) / s^11 dphi, which
 *      exceeds |S_m| by up to 1e7 for the high modes at C1 ~ 0.05. The error of
 *      mode m is normalised by |S_m| + 1e-19 I / MAX_MODE_ERROR, so that it is
 *      below MAX_MODE_ERROR where the difference is below MAX_MODE_ERROR |S_m|
 *      plus the rounding of the reference. With --dump FILE the program also
 *      writes the coefficients and analytic m-modes of every case for
 *      mmode_exact_check.py, which checks each mode against mpmath.
 *   3. The analytic m-modes from the hybrid effsource_equatorial_ctx_prepare
 *      and the numerical m-modes of the hybrid sweep_src route are compared with
 *      the trapezoid rule applied to the long double calc function, where that
 *      grid has at most 2^15 points. The error is normalised by the largest
 *      reference m-mode, as in point-check.
 *   4. The analytic m-modes of PhiS, dPhiS/dr and dPhiS/dtheta from
 *      effsource_equatorial_ctx_analytic_modes are checked twice. First, those from
 *      effsource_equatorial_ctx_prepare_dd are compared mode by mode with
 *      the trapezoid rule applied, on the long double nodes of check 2, to the
 *      chain-rule formulas of effsource_equatorial_ctx_sweep_calc evaluated in
 *      long double from the same coefficients. The error of mode m is
 *      normalised by |f_m| + 1e-19 I / MAX_FIELD_MODE_ERROR with
 *      I = int |f| dphi. This checks the u-polynomial numerators and the
 *      contractions. Second, those from the double-double and the hybrid
 *      prepare are compared with the trapezoid rule applied to the long double
 *      calc function, where that grid has at most 2^17 points, with the error
 *      normalised by the largest reference m-mode of each output. Its nodes
 *      2 pi k / n are rounded to double, which biases the small modes by about
 *      1e-16 |f'|, so this comparison is not made mode by mode.
 *
 * Errors are grouped by decade of C1. The program also reports the time of the
 * analytic m-modes and of the numerical m-modes with the azimuthal grid of the
 * Auto rule of EccentricOrbit. It prints one "check <name> PASS|FAIL <value>"
 * line per check and exits with status 1 if any check fails.
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "effsource_equatorial.h"
#include "effsource_equatorial_mmode.h"
#include "effsource_equatorial_point.h"

#define MMAX 20
#define KERNEL_MU_MAX 40
#define NUM_DECADES 11
#define MAX_KERNEL_ERROR 1e-14
#define MAX_MODE_ERROR 1e-13
#define MAX_SRC_ERROR_RATIO 2.0
#define ERROR_FLOOR 1e-13
#define MAX_FIELD_ERROR 1e-14
#define MAX_FIELD_MODE_ERROR 1e-13
#define FIELD_MAX_POINTS (1 << 17)
#define NUM_FIELDS 3

static const long double pi_long = 3.141592653589793238462643383279502884L;

void effsource_equatorial_ctx_calc_long_double(const struct effsource_equatorial_ctx *ctx,
    const struct coordinate *x, long double *PhiS, long double *dPhiS_dx,
    long double *d2PhiS_dx2, long double *src);

static double now(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

static double log_rho(double C1)
{
  return log1p(2.0 * C1 + 2.0 * sqrt(C1 * (1.0 + C1)));
}

/* Smallest power of two (at least 64) with aliasing below 1e-20 for modes up to
   mu_max, relative to each mode. The kernels decrease as mu^(nu - 1) rho^-mu
   with nu up to 11/2, the Fourier coefficients of the numerator exceed the
   m-modes by up to 1e8, and the highest modes are down to 1e-6 of the largest, so
   n - mu_max >= (46 + 4.5 ln n + 18.4 + 13.8) / ln(rho). */
static int resolved_points(double C1, int mu_max)
{
  const double rate = log_rho(C1);
  int n = 64;
  while (n < 2 * mu_max + 2 || (n - mu_max) * rate < 78.2 + 4.5 * log(n))
    n *= 2;
  return n;
}

/* The Auto rule of EccentricOrbit for mode m. */
static int auto_points(double C1, int m)
{
  const double want = fmax(64.0, fmax(fmax(m, 20) + 16.0 / sqrt(C1), 2.5 * m));
  int n = 64;
  while (n < want)
    n *= 2;
  return n;
}

/* Energy and angular momentum of the prograde circular orbit at radius r (M = 1). */
static void circular_constants(double a, double r, double *E, double *L)
{
  const double v = 1.0 / sqrt(r);
  const double den = sqrt(1.0 - 3.0 * v * v + 2.0 * a * v * v * v);
  *E = (1.0 - 2.0 * v * v + a * v * v * v) / den;
  *L = sqrt(r) * (1.0 - 2.0 * a * v * v * v + a * a * v * v * v * v) / den;
}

static int decade_of(double C1)
{
  return (int)fmin(NUM_DECADES - 1, fmax(0.0, -floor(log10(C1))));
}

/* Largest kernel error for one (alpha, beta). */
static double kernel_error(double alpha, double beta)
{
  enum { NUM_NU = EFFSOURCE_MMODE_NUM_NU, STRIDE = KERNEL_MU_MAX + 1 };
  static double K[NUM_NU * STRIDE];
  static long double reference[NUM_NU * STRIDE];
  if (effsource_equatorial_mmode_kernels(alpha, beta, NUM_NU, KERNEL_MU_MAX, K))
    return INFINITY;
  const int n = resolved_points(alpha / beta, KERNEL_MU_MAX);
  for (int i = 0; i < NUM_NU * STRIDE; i++)
    reference[i] = 0.0L;
  for (int k = 0; k < n; k++) {
    const long double psi = 2.0L * pi_long * (k <= n / 2 ? k : k - n) / n;
    const long double half = sinl(0.5L * psi), s2 = alpha + beta * half * half;
    const long double cos_psi = cosl(psi);
    long double power = 1.0L / (s2 * s2 * s2 * s2 * s2 * sqrtl(s2));
    for (int i = 0; i < NUM_NU; i++) {
      long double c_previous = cos_psi, c = 1.0L;
      for (int mu = 0; mu <= KERNEL_MU_MAX; mu++) {
        reference[i * STRIDE + mu] += power * c;
        const long double c_next = 2.0L * cos_psi * c - c_previous;
        c_previous = c;
        c = c_next;
      }
      power *= s2;
    }
  }
  double worst = 0.0;
  for (int i = 0; i < NUM_NU; i++) {
    const long double scale = fabsl(reference[i * STRIDE]);
    for (int mu = 0; mu <= KERNEL_MU_MAX; mu++)
      worst = fmax(worst, (double)(fabsl(K[i * STRIDE + mu] * n / (2.0L * pi_long) - reference[i * STRIDE + mu]) / scale));
  }
  return worst;
}

/* Numerical m-modes of the rational form of sweep_src with the double-double
   coefficients, evaluated in long double on nodes of psi = dphib in (-pi, pi],
   and the bound I on its rounding. */
static double rational_modes(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                             int n, long double complex *modes)
{
  const long double alpha = (long double)point->src_alpha + point->src_alpha_lo, beta = ctx->beta;
  long double P[8], Q[7];
  for (int j = 0; j < 8; j++)
    P[j] = (long double)point->src_P[j] + point->src_P_lo[j];
  for (int j = 0; j < 7; j++)
    Q[j] = (long double)point->src_Q[j] + point->src_Q_lo[j];
  long double absolute = 0.0L;
  for (int m = 0; m <= MMAX; m++)
    modes[m] = 0.0L;
  for (int k = 0; k < n; k++) {
    const long double psi = 2.0L * pi_long * (k <= n / 2 ? k : k - n) / n;
    const long double half = sinl(0.5L * psi), u = half * half, s2 = alpha + beta * u;
    long double Pu = 0.0L, Qu = 0.0L, Pu_bound = 0.0L, Qu_bound = 0.0L;
    for (int j = 7; j >= 0; j--) {
      Pu = Pu * u + P[j];
      Pu_bound = Pu_bound * u + fabsl(P[j]);
    }
    for (int j = 6; j >= 0; j--) {
      Qu = Qu * u + Q[j];
      Qu_bound = Qu_bound * u + fabsl(Q[j]);
    }
    const long double denominator = s2 * s2 * s2 * s2 * s2 * sqrtl(s2);
    const long double src = (Pu + 0.5L * sinl(psi) * Qu) / denominator;
    const long double complex step = cosl(psi) - I * sinl(psi);
    long double complex phase = 1.0L;
    absolute += (Pu_bound + fabsl(0.5L * sinl(psi)) * Qu_bound) / denominator;
    for (int m = 0; m <= MMAX; m++) {
      modes[m] += src * phase;
      phase *= step;
    }
  }
  const long double phi0 = ctx->xp.phi + ctx->c * (point->r - ctx->xp.r);
  for (int m = 0; m <= MMAX; m++)
    modes[m] *= (cosl(m * phi0) - I * sinl(m * phi0)) * 2.0L * pi_long / n;
  return (double)(absolute * 2.0L * pi_long / n);
}

/* Numerical m-modes of sweep_src on n points, summed in long double, and int |src| dphi. */
static double numerical_modes(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                              int n, long double complex *modes)
{
  long double absolute = 0.0L;
  for (int m = 0; m <= MMAX; m++)
    modes[m] = 0.0L;
  for (int k = 0; k < n; k++) {
    const double phi = 2.0 * M_PI * k / n;
    const long double src = effsource_equatorial_ctx_sweep_src(ctx, point, phi);
    const long double complex step = cexpl(-I * (long double)phi);
    long double complex phase = 1.0L;
    absolute += fabsl(src);
    for (int m = 0; m <= MMAX; m++) {
      modes[m] += src * phase;
      phase *= step;
    }
  }
  for (int m = 0; m <= MMAX; m++)
    modes[m] *= 2.0L * pi_long / n;
  return (double)(absolute * 2.0L * pi_long / n);
}

/* Numerical m-modes of the long double calc function on n points. */
static void reference_modes(const struct effsource_equatorial_ctx *ctx, double r, double theta, int n,
                            long double complex *modes)
{
  for (int m = 0; m <= MMAX; m++)
    modes[m] = 0.0L;
  for (int k = 0; k < n; k++) {
    const double phi = 2.0 * M_PI * k / n;
    struct coordinate x = {.r = r, .theta = theta, .phi = phi, .t = 0.0};
    long double PhiS, d1[4], d2[10], src;
    effsource_equatorial_ctx_calc_long_double(ctx, &x, &PhiS, d1, d2, &src);
    const long double complex step = cexpl(-I * (long double)phi);
    long double complex phase = 1.0L;
    for (int m = 0; m <= MMAX; m++) {
      modes[m] += src * phase;
      phase *= step;
    }
  }
  for (int m = 0; m <= MMAX; m++)
    modes[m] *= 2.0L * pi_long / n;
}

/* Numerical m-modes of PhiS, dPhiS/dr and dPhiS/dtheta of the long double calc function on
   n points, and int |f| dphi of each. */
static void field_reference_modes(const struct effsource_equatorial_ctx *ctx, double r, double theta, int n,
                                  long double complex modes[NUM_FIELDS][MMAX + 1], double *absolute)
{
  long double sum_absolute[NUM_FIELDS] = {0.0L};
  for (int f = 0; f < NUM_FIELDS; f++)
    for (int m = 0; m <= MMAX; m++)
      modes[f][m] = 0.0L;
  for (int k = 0; k < n; k++) {
    const double phi = 2.0 * M_PI * k / n;
    struct coordinate x = {.r = r, .theta = theta, .phi = phi, .t = 0.0};
    long double PhiS, d1[4], d2[10], src;
    effsource_equatorial_ctx_calc_long_double(ctx, &x, &PhiS, d1, d2, &src);
    const long double values[NUM_FIELDS] = {PhiS, d1[1], d1[2]};
    const long double complex step = cexpl(-I * (long double)phi);
    long double complex phase = 1.0L;
    for (int f = 0; f < NUM_FIELDS; f++)
      sum_absolute[f] += fabsl(values[f]);
    for (int m = 0; m <= MMAX; m++) {
      for (int f = 0; f < NUM_FIELDS; f++)
        modes[f][m] += values[f] * phase;
      phase *= step;
    }
  }
  for (int f = 0; f < NUM_FIELDS; f++) {
    absolute[f] = (double)(sum_absolute[f] * 2.0L * pi_long / n);
    for (int m = 0; m <= MMAX; m++)
      modes[f][m] *= 2.0L * pi_long / n;
  }
}

/* Numerical m-modes of PhiS, dPhiS/dr and dPhiS/dtheta from the chain-rule formulas of
   sweep_calc with the double-double coefficients, evaluated in long double on nodes of
   psi = dphib in (-pi, pi], and int |f| dphi of each. */
static void field_rational_modes(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                                 int n, long double complex modes[NUM_FIELDS][MMAX + 1], double *absolute)
{
#define LOAD(name, size) long double name[size]; for (int j = 0; j < size; j++) name[j] = (long double)point->name[j] + point->name##_lo[j]
  LOAD(A_R0, 5); LOAD(A_R1, 5); LOAD(dA_dr_R0, 5); LOAD(dA_dr_R1, 4); LOAD(dA_dth_R0, 4); LOAD(dA_dth_R1, 4);
  LOAD(dA_dQ_R0, 4); LOAD(dA_dQ_R1, 4); LOAD(dA_dR_R0, 5);
#undef LOAD
  const long double alpha = (long double)point->src_alpha + point->src_alpha_lo, beta = ctx->beta, c = ctx->c;
  const long double dr = (long double)point->r - ctx->xp.r, dtheta = (long double)point->theta - ctx->xp.theta;
  long double sum_absolute[NUM_FIELDS] = {0.0L};
  for (int f = 0; f < NUM_FIELDS; f++)
    for (int m = 0; m <= MMAX; m++)
      modes[f][m] = 0.0L;
  for (int k = 0; k < n; k++) {
    const long double psi = 2.0L * pi_long * (k <= n / 2 ? k : k - n) / n;
    const long double dQ = sinl(0.5L * psi), cQ = cosl(0.5L * psi), u = dQ * dQ, dR = sinl(psi), cR = cosl(psi);
#define POLY(name, size) ({ long double value_ = 0.0L; for (int j = size - 1; j >= 0; j--) value_ = value_ * u + name[j]; value_; })
    const long double A = POLY(A_R0, 5) + dR * POLY(A_R1, 5);
    const long double dA_dr = POLY(dA_dr_R0, 5) + dR * POLY(dA_dr_R1, 4);
    const long double dA_dth = POLY(dA_dth_R0, 4) + dR * POLY(dA_dth_R1, 4);
    const long double dA_dQ = dQ * POLY(dA_dQ_R0, 4) + dR * dQ * POLY(dA_dQ_R1, 4);
    const long double dA_dR = POLY(dA_dR_R0, 5);
#undef POLY
    const long double dQ_dr = -c * 0.5L * cQ, dR_dr = -c * cR;
    const long double DA_Dr = dA_dr + dA_dQ * dQ_dr + dA_dR * dR_dr;
    const long double s2 = alpha + beta * u, s = sqrtl(s2), s7 = s2 * s2 * s2 * s, s9 = s7 * s2;
    const long double Ds2_Dr = 2.0L * ctx->alpha20 * dr + 2.0L * beta * dQ * dQ_dr;
    const long double values[NUM_FIELDS] = {A / s7, (-7.0L * Ds2_Dr * A + 2.0L * DA_Dr * s2) / (2.0L * s9),
                                            (-14.0L * ctx->alpha02 * dtheta * A + 2.0L * dA_dth * s2) / (2.0L * s9)};
    const long double complex step = cR - I * dR;
    long double complex phase = 1.0L;
    for (int f = 0; f < NUM_FIELDS; f++)
      sum_absolute[f] += fabsl(values[f]);
    for (int m = 0; m <= MMAX; m++) {
      for (int f = 0; f < NUM_FIELDS; f++)
        modes[f][m] += values[f] * phase;
      phase *= step;
    }
  }
  const long double phi0 = ctx->xp.phi + ctx->c * (point->r - ctx->xp.r);
  for (int f = 0; f < NUM_FIELDS; f++) {
    absolute[f] = (double)(sum_absolute[f] * 2.0L * pi_long / n);
    for (int m = 0; m <= MMAX; m++)
      modes[f][m] *= (cosl(m * phi0) - I * sinl(m * phi0)) * 2.0L * pi_long / n;
  }
}

/* Largest error of the analytic field m-modes relative to the largest reference mode of each
   field, and largest error of a mode relative to itself plus the rounding of the reference. */
static void field_errors(const double analytic[NUM_FIELDS][2 * (MMAX + 1)],
                         long double complex reference[NUM_FIELDS][MMAX + 1], const double *absolute,
                         double *global, double *per_mode)
{
  for (int f = 0; f < NUM_FIELDS; f++) {
    long double scale = 0.0L;
    for (int m = 0; m <= MMAX; m++)
      scale = fmaxl(scale, cabsl(reference[f][m]));
    for (int m = 0; m <= MMAX; m++) {
      const long double difference = cabsl(analytic[f][2 * m] + I * analytic[f][2 * m + 1] - reference[f][m]);
      global[f] = fmax(global[f], (double)(difference / scale));
      per_mode[f] = fmax(per_mode[f],
                         (double)(difference / (cabsl(reference[f][m]) + 1e-19L * absolute[f] / MAX_FIELD_MODE_ERROR)));
    }
  }
}

/* Time per call of the numerical m-modes with the Auto grid, as EccentricOrbit computes them. */
static double numerical_time(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                             int n, double complex *sink)
{
  static double complex table[1 << 22];
  for (int k = 0; k < n; k++)
    table[k] = cexp(-I * (2.0 * M_PI * k / n));
  const int repeats = n > 4096 ? 1 : 4096 / n;
  const double start = now();
  for (int repeat = 0; repeat < repeats; repeat++) {
    for (int k = 0; k < n; k++) {
      const double src = effsource_equatorial_ctx_sweep_src(ctx, point, 2.0 * M_PI * k / n);
      double complex phase = 1.0;
      for (int m = 0; m <= MMAX; m++) {
        sink[m] += src * phase;
        phase *= table[k];
      }
    }
  }
  return (now() - start) / repeats;
}

/* Write the coefficients and analytic m-modes of one case for mmode_exact_check.py. */
static void dump_case(FILE *stream, const char *tag, const struct effsource_equatorial_ctx *ctx,
                      const struct effsource_equatorial_point *point)
{
  double re[MMAX + 1], im[MMAX + 1];
  effsource_equatorial_ctx_analytic_src_modes(ctx, point, MMAX, re, im);
  fprintf(stream, "%s %a %a %a %a\n", tag, point->src_alpha, point->src_alpha_lo, ctx->beta,
          ctx->xp.phi + ctx->c * (point->r - ctx->xp.r));
  fprintf(stream, "P");
  for (int j = 0; j < 8; j++)
    fprintf(stream, " %a %a", point->src_P[j], point->src_P_lo[j]);
  fprintf(stream, "\nQ");
  for (int j = 0; j < 7; j++)
    fprintf(stream, " %a %a", point->src_Q[j], point->src_Q_lo[j]);
  fprintf(stream, "\nS");
  for (int m = 0; m <= MMAX; m++)
    fprintf(stream, " %a %a", re[m], im[m]);
  fprintf(stream, "\n");
}

static int check(const char *name, int pass, double value)
{
  printf("check %s %s %.2e\n", name, pass ? "PASS" : "FAIL", value);
  return pass;
}

int main(int argc, char **argv)
{
  int all_pass = 1;
  FILE *dump = NULL;
  if (argc == 3 && strcmp(argv[1], "--dump") == 0) {
    dump = fopen(argv[2], "w");
    if (dump == NULL) {
      fprintf(stderr, "cannot open %s\n", argv[2]);
      return 1;
    }
  }

  const double kernel_C1[] = {1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 1.0, 10.0, 1e3};
  const double kernel_beta[] = {455.8, 30.0};
  double worst_kernel = 0.0;
  printf("kernels vs long double trapezoid, nu = 11/2..-5/2, mu <= %d\n", KERNEL_MU_MAX);
  for (int ib = 0; ib < 2; ib++)
    for (int ic = 0; ic < (int)(sizeof kernel_C1 / sizeof kernel_C1[0]); ic++) {
      const double error = kernel_error(kernel_C1[ic] * kernel_beta[ib], kernel_beta[ib]);
      printf("  beta %6.1f  C1 %7.0e  worst %.1e\n", kernel_beta[ib], kernel_C1[ic], error);
      worst_kernel = fmax(worst_kernel, error);
    }
  all_pass &= check("kernels", worst_kernel <= MAX_KERNEL_ERROR, worst_kernel);

  const double spins[] = {0.0, 0.9};
  const double particle_r[] = {7.0, 16.0};
  const double particle_ur[] = {0.0, 0.03};
  const double offsets_r[] = {0.0, 0.01, -0.5, 3.0};
  const double offsets_theta[] = {0.5, 0.1, 1e-2, 1e-3, 1e-4};

  double mode_worst[NUM_DECADES] = {0}, src_analytic[NUM_DECADES] = {0}, src_numerical[NUM_DECADES] = {0};
  double time_analytic[NUM_DECADES] = {0}, time_numerical[NUM_DECADES] = {0}, time_prepare[NUM_DECADES] = {0};
  double time_fields[NUM_DECADES] = {0};
  double field_global_dd[NUM_DECADES][NUM_FIELDS] = {{0}}, field_mode_dd[NUM_DECADES][NUM_FIELDS] = {{0}};
  double field_global[NUM_DECADES][NUM_FIELDS] = {{0}}, field_mode[NUM_DECADES][NUM_FIELDS] = {{0}};
  int count[NUM_DECADES] = {0}, count_reference[NUM_DECADES] = {0}, count_fourier[NUM_DECADES] = {0};
  int count_field[NUM_DECADES] = {0};
  double complex sink[MMAX + 1] = {0};

  for (int is = 0; is < 2; is++)
  for (int ir = 0; ir < 2; ir++)
  for (int iu = 0; iu < 2; iu++) {
    double E, L;
    circular_constants(spins[is], particle_r[ir], &E, &L);
    struct effsource_equatorial_ctx *ctx = effsource_equatorial_create(1.0, spins[is]);
    struct coordinate xp = {.r = particle_r[ir], .theta = M_PI_2, .phi = 0.3, .t = 0.0};
    effsource_equatorial_ctx_set_particle(ctx, &xp, E, L, particle_ur[iu]);

    for (int idr = 0; idr < 4; idr++)
    for (int idth = 0; idth < 5; idth++) {
      const double r = particle_r[ir] + offsets_r[idr];
      const double theta = M_PI_2 - offsets_theta[idth];
      struct effsource_equatorial_point point_dd, point;
      effsource_equatorial_ctx_prepare_dd(ctx, r, theta, &point_dd);
      const double C1 = point_dd.src_alpha / ctx->beta;
      const int decade = decade_of(C1);
      count[decade]++;
      count_fourier[decade] += C1 >= EFFSOURCE_MMODE_FOURIER_C1;

      if (dump) {
        struct effsource_equatorial_point point_double;
        effsource_equatorial_ctx_prepare_double(ctx, r, theta, &point_double);
        dump_case(dump, "DD", ctx, &point_dd);
        dump_case(dump, "DBL", ctx, &point_double);
      }
      double re[MMAX + 1], im[MMAX + 1];
      effsource_equatorial_ctx_analytic_src_modes(ctx, &point_dd, MMAX, re, im);
      long double complex numerical[MMAX + 1];
      const int n = resolved_points(C1, MMAX);
      const double absolute = rational_modes(ctx, &point_dd, n, numerical);
      for (int m = 0; m <= MMAX; m++) {
        const long double complex difference = re[m] + I * im[m] - numerical[m];
        const double error = (double)(cabsl(difference) / (cabsl(numerical[m]) + 1e-19L * absolute / MAX_MODE_ERROR));
        mode_worst[decade] = fmax(mode_worst[decade], error);
      }

      double start = now();
      const int repeats = 200;
      for (int repeat = 0; repeat < repeats; repeat++)
        effsource_equatorial_ctx_prepare(ctx, r, theta, &point);
      time_prepare[decade] += (now() - start) / repeats;
      start = now();
      for (int repeat = 0; repeat < repeats; repeat++)
        effsource_equatorial_ctx_analytic_src_modes(ctx, &point, MMAX, re, im);
      time_analytic[decade] += (now() - start) / repeats;
      double fields[NUM_FIELDS][2 * (MMAX + 1)], src_modes[2 * (MMAX + 1)];
      start = now();
      for (int repeat = 0; repeat < repeats; repeat++)
        effsource_equatorial_ctx_analytic_modes(ctx, &point, MMAX, src_modes, fields[0], fields[1], fields[2]);
      time_fields[decade] += (now() - start) / repeats;
      time_numerical[decade] += numerical_time(ctx, &point, auto_points(C1, MMAX), sink);

      if (n <= (1 << 15)) {
        long double complex reference[MMAX + 1], hybrid_numerical[MMAX + 1];
        reference_modes(ctx, r, theta, n, reference);
        numerical_modes(ctx, &point, n, hybrid_numerical);
        long double scale = 0.0L, error_analytic = 0.0L, error_numerical = 0.0L;
        for (int m = 0; m <= MMAX; m++) {
          scale = fmaxl(scale, cabsl(reference[m]));
          error_analytic = fmaxl(error_analytic, cabsl(re[m] + I * im[m] - reference[m]));
          error_numerical = fmaxl(error_numerical, cabsl(hybrid_numerical[m] - reference[m]));
        }
        src_analytic[decade] = fmax(src_analytic[decade], (double)(error_analytic / scale));
        src_numerical[decade] = fmax(src_numerical[decade], (double)(error_numerical / scale));
        count_reference[decade]++;
      }

      double fields_dd[NUM_FIELDS][2 * (MMAX + 1)], unused[NUM_FIELDS];
      effsource_equatorial_ctx_analytic_modes(ctx, &point_dd, MMAX, NULL, fields_dd[0], fields_dd[1], fields_dd[2]);
      {
        long double complex rational[NUM_FIELDS][MMAX + 1];
        double absolute[NUM_FIELDS];
        field_rational_modes(ctx, &point_dd, n, rational, absolute);
        field_errors((const double (*)[2 * (MMAX + 1)])fields_dd, rational, absolute, unused, field_mode_dd[decade]);
      }
      if (n <= FIELD_MAX_POINTS) {
        long double complex reference[NUM_FIELDS][MMAX + 1];
        double absolute[NUM_FIELDS];
        field_reference_modes(ctx, r, theta, n, reference, absolute);
        field_errors((const double (*)[2 * (MMAX + 1)])fields_dd, reference, absolute, field_global_dd[decade],
                     unused);
        field_errors((const double (*)[2 * (MMAX + 1)])fields, reference, absolute, field_global[decade],
                     field_mode[decade]);
        count_field[decade]++;
      }
    }
    effsource_equatorial_free(ctx);
  }

  printf("\nm-modes m <= %d by decade of C1 (cases: all / Fourier basis / with long double reference)\n", MMAX);
  printf("  decade   cases     | analytic vs numerical, same coefficients | hybrid vs long double calc: analytic  numerical"
         " | time per point (us): prepare  analytic  numerical(Auto grid)\n");
  double worst_mode = 0.0, worst_ratio = 0.0;
  for (int d = 0; d < NUM_DECADES; d++) {
    if (count[d] == 0)
      continue;
    printf("  1e-%-2d  %3d/%3d/%3d | %10.1e                               |", d, count[d], count_fourier[d], count_reference[d],
           mode_worst[d]);
    if (count_reference[d])
      printf("                           %9.1e  %9.1e |", src_analytic[d], src_numerical[d]);
    else
      printf("                                   -          - |");
    printf("                     %7.2f  %8.2f  %9.2f\n", 1e6 * time_prepare[d] / count[d], 1e6 * time_analytic[d] / count[d],
           1e6 * time_numerical[d] / count[d]);
    worst_mode = fmax(worst_mode, mode_worst[d]);
    if (count_reference[d] && src_analytic[d] > ERROR_FLOOR)
      worst_ratio = fmax(worst_ratio, src_analytic[d] / fmax(src_numerical[d], ERROR_FLOOR));
  }
  all_pass &= check("modes_vs_numerical", worst_mode <= MAX_MODE_ERROR, worst_mode);
  all_pass &= check("modes_vs_long_double", worst_ratio <= MAX_SRC_ERROR_RATIO, worst_ratio);

  printf("\nPhiS, dPhiS/dr, dPhiS/dtheta m-modes m <= %d vs long double calc, by decade of C1\n", MMAX);
  printf("  decade  cases | vs calc, largest-mode error: prepare_dd | hybrid                             |"
         " vs same coefficients, per mode   | vs calc, per mode, hybrid          | time (us): src  all four\n");
  double worst_field = 0.0, worst_field_mode = 0.0;
  for (int d = 0; d < NUM_DECADES; d++) {
    if (count[d] == 0)
      continue;
    printf("  1e-%-2d  %3d   |", d, count_field[d]);
    for (int f = 0; f < NUM_FIELDS; f++)
      printf(" %9.1e", field_global_dd[d][f]);
    printf("    |");
    for (int f = 0; f < NUM_FIELDS; f++)
      printf(" %9.1e", field_global[d][f]);
    printf("    |");
    for (int f = 0; f < NUM_FIELDS; f++)
      printf(" %9.1e", field_mode_dd[d][f]);
    printf("  |");
    for (int f = 0; f < NUM_FIELDS; f++)
      printf(" %9.1e", field_mode[d][f]);
    printf("  |       %7.2f  %8.2f\n", 1e6 * time_analytic[d] / count[d], 1e6 * time_fields[d] / count[d]);
    for (int f = 0; f < NUM_FIELDS; f++) {
      worst_field = fmax(worst_field, fmax(field_global_dd[d][f], field_global[d][f]));
      worst_field_mode = fmax(worst_field_mode, field_mode_dd[d][f]);
    }
  }
  all_pass &= check("field_modes", worst_field <= MAX_FIELD_ERROR, worst_field);
  all_pass &= check("field_modes_per_mode", worst_field_mode <= MAX_FIELD_MODE_ERROR, worst_field_mode);
  if (creal(sink[0]) == 1.2345)
    printf("\n");
  if (dump)
    fclose(dump);
  return all_pass ? 0 : 1;
}
