/*******************************************************************************
 * Check the numerical m-modes of effsource_equatorial_numerical.h.
 *
 *   1. effsource_equatorial_numerical_num_phi against values of the rule
 *      N = 2^ceil(log2 max(64, max(|m|, 20) + 46 / log(rho), 2.5 |m|)).
 *   2. The source modes from the hybrid effsource_equatorial_ctx_prepare,
 *      with N from the rule, are compared with a trapezoid sum of
 *      effsource_equatorial_ctx_point_src on the coefficients of
 *      effsource_equatorial_ctx_prepare_dd, summed in long double on a grid
 *      fine enough that aliasing is below 1e-20. The error is normalised by
 *      I = int |src| dphi, and the tolerance is 1e-12 I where prepare uses
 *      double-double arithmetic (C1 < EFFSOURCE_PREPARE_DD_C1) and 1e-10 I
 *      elsewhere.
 *   3. The source modes from the coefficients of prepare_dd, with N from the
 *      rule, are compared with the same reference to 1e-13 I. Both sides use
 *      the same coefficients, so this measures the quadrature error of the
 *      rule and the rounding of the sum.
 *   4. The modes of PhiS, dPhiS/dr and dPhiS/dtheta from the hybrid prepare are
 *      compared with a trapezoid sum of the long double calc function where
 *      that grid has at most 2^17 points, and with a trapezoid sum of the
 *      chain-rule formulas of point_calc in long double on the prepare_dd
 *      coefficients elsewhere, to 1e-12 int |f| dphi for each field.
 *   5. A table built for 2^14 nodes and read with stride 16 gives the modes of
 *      a table built for 2^10 nodes.
 *   6. At theta = pi/2 -+ dtheta the modes of src, PhiS and dPhiS/dr are equal
 *      and those of dPhiS/dtheta are opposite, to 1e-13 int |f| dphi. For src
 *      the tolerance grows as 1e-14 / dtheta below dtheta = 0.1, because
 *      M_PI_2 is not exactly pi/2 and the src of calc is asymmetric by up to
 *      about 3e-15 / dtheta.
 *
 * The program also reports the time per field point of prepare plus the source
 * mode at m = 2 against N calls of effsource_equatorial_ctx_calc. It prints one
 * "check <name> PASS|FAIL <value>" line per check, where the value is the
 * largest error in units of the tolerance, and exits with status 1 if any check
 * fails.
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "effsource_equatorial.h"
#include "effsource_equatorial_numerical.h"
#include "effsource_equatorial_point.h"

#define MMAX 20
#define NUM_M 5
#define NUM_FIELDS 3
#define NUM_DECADES 11
#define FIELD_LD_MAX_POINTS (1 << 17)
#define HYBRID_DD_TOLERANCE 1e-12
#define HYBRID_DOUBLE_TOLERANCE 1e-10
#define DD_TOLERANCE 1e-13
#define FIELD_TOLERANCE 1e-12
#define STRIDE_TOLERANCE 1e-15
#define SYMMETRY_TOLERANCE 1e-13

static const int m_list[NUM_M] = {0, 1, 2, 7, 20};
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
   mu_max, relative to each mode (see mmode-check.c). */
static int resolved_points(double C1, int mu_max)
{
  const double rate = log_rho(C1);
  int n = 64;
  while (n < 2 * mu_max + 2 || (n - mu_max) * rate < 78.2 + 4.5 * log(n))
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

/* Numerical m-modes of point_src on n points, summed in long double, and int |src| dphi. */
static double numerical_modes(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                              int n, long double complex *modes)
{
  long double absolute = 0.0L;
  for (int m = 0; m <= MMAX; m++)
    modes[m] = 0.0L;
  for (int k = 0; k < n; k++) {
    const double phi = 2.0 * M_PI * k / n;
    const long double src = effsource_equatorial_ctx_point_src(ctx, point, phi);
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
   point_calc with the double-double coefficients, evaluated in long double on nodes of
   psi = dphib in (-pi, pi], and int |f| dphi of each (see mmode-check.c). */
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

/* Source mode m with a table built for and read at N nodes; NAN on failure. */
static double complex src_mode(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                               int m, size_t N)
{
  struct effsource_equatorial_numerical_table table;
  double re, im;
  if (effsource_equatorial_numerical_table_init(&table, m, N))
    return NAN;
  const int status = effsource_equatorial_ctx_numerical_src_mode(ctx, point, &table, N, &re, &im);
  effsource_equatorial_numerical_table_free(&table);
  return status ? NAN : re + I * im;
}

/* Field modes m with a table built for and read at N nodes; returns the status. */
static int field_modes(const struct effsource_equatorial_ctx *ctx, const struct effsource_equatorial_point *point,
                       int m, size_t N, double complex modes[NUM_FIELDS])
{
  struct effsource_equatorial_numerical_table table;
  double re[NUM_FIELDS], im[NUM_FIELDS];
  if (effsource_equatorial_numerical_table_init(&table, m, N))
    return 1;
  const int status = effsource_equatorial_ctx_numerical_field_modes(ctx, point, &table, N, re, im);
  effsource_equatorial_numerical_table_free(&table);
  for (int f = 0; f < NUM_FIELDS; f++)
    modes[f] = re[f] + I * im[f];
  return status;
}

static int check(const char *name, int pass, double value)
{
  printf("check %s %s %.2e\n", name, pass ? "PASS" : "FAIL", value);
  return pass;
}

int main(void)
{
  int all_pass = 1;

  const struct { double C1; int m; size_t n; } rule[] = {
    {1.0, 2, 64}, {1e-2, 2, 256}, {1e-6, 2, 32768}, {1e-6, -2, 32768},
    {0.25, 40, 128}, {0.0, 2, 0}, {-1.0, 2, 0}, {1e-12, 2, 0}};
  int rule_failures = 0;
  for (size_t i = 0; i < sizeof rule / sizeof rule[0]; i++) {
    const size_t n = effsource_equatorial_numerical_num_phi(rule[i].C1, rule[i].m);
    if (n != rule[i].n) {
      printf("  num_phi(%g, %d) = %zu, expected %zu\n", rule[i].C1, rule[i].m, n, rule[i].n);
      rule_failures++;
    }
  }
  all_pass &= check("num_phi_rule", rule_failures == 0, rule_failures);

  int argument_failures = 0;
  {
    struct effsource_equatorial_numerical_table table;
    argument_failures += effsource_equatorial_numerical_table_init(&table, 2, 100) == 0;
    argument_failures += effsource_equatorial_numerical_table_init(&table, 2, 2) == 0;
  }

  const double spins[] = {0.5, 0.9};
  const double particle_r[] = {7.0, 10.0, 16.0};
  const double particle_ur[] = {0.0, 0.03, -0.03};
  const double offsets_r[] = {0.0, 0.01, -0.5, 3.0};
  const double offsets_theta[] = {0.5, 0.1, 1e-2, 1e-3, 1e-4};

  struct effsource_equatorial_numerical_table stride_big[NUM_M], stride_small[NUM_M];
  for (int im = 0; im < NUM_M; im++)
    if (effsource_equatorial_numerical_table_init(&stride_big[im], m_list[im], 1 << 14)
        || effsource_equatorial_numerical_table_init(&stride_small[im], m_list[im], 1 << 10)) {
      fprintf(stderr, "table allocation failed\n");
      return 1;
    }

  double hybrid_error[NUM_DECADES] = {0}, dd_error[NUM_DECADES] = {0}, field_error[NUM_DECADES][NUM_FIELDS] = {{0}};
  int count[NUM_DECADES] = {0}, count_field_ld[NUM_DECADES] = {0};
  double worst_hybrid = 0.0, worst_dd = 0.0, worst_field = 0.0, worst_stride = 0.0, worst_symmetry = 0.0;
  double min_C1 = INFINITY, time_numerical = 0.0, time_calc = 0.0, sink = 0.0;
  int num_points = 0;

  for (int is = 0; is < 2; is++)
  for (int ir = 0; ir < 3; ir++)
  for (int iu = 0; iu < 3; iu++) {
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
      effsource_equatorial_ctx_prepare(ctx, r, theta, &point);
      const double C1 = effsource_equatorial_point_C1(ctx, &point);
      const int decade = decade_of(C1);
      const double hybrid_tolerance = C1 < EFFSOURCE_PREPARE_DD_C1 ? HYBRID_DD_TOLERANCE : HYBRID_DOUBLE_TOLERANCE;
      min_C1 = fmin(min_C1, C1);
      count[decade]++;

      const int n_reference = resolved_points(C1, MMAX);
      long double complex reference[MMAX + 1];
      const double absolute = numerical_modes(ctx, &point_dd, n_reference, reference);
      long double complex field_reference[NUM_FIELDS][MMAX + 1];
      double field_absolute[NUM_FIELDS];
      if (n_reference <= FIELD_LD_MAX_POINTS) {
        field_reference_modes(ctx, r, theta, n_reference, field_reference, field_absolute);
        count_field_ld[decade]++;
      } else {
        field_rational_modes(ctx, &point_dd, n_reference, field_reference, field_absolute);
      }

      for (int im = 0; im < NUM_M; im++) {
        const int m = m_list[im];
        const size_t N = effsource_equatorial_numerical_num_phi(C1, m);
        const double complex hybrid = src_mode(ctx, &point, m, N);
        const double complex dd = src_mode(ctx, &point_dd, m, N);
        const double e_hybrid = (double)(cabsl(hybrid - reference[m]) / absolute);
        const double e_dd = (double)(cabsl(dd - reference[m]) / absolute);
        hybrid_error[decade] = fmax(hybrid_error[decade], e_hybrid);
        dd_error[decade] = fmax(dd_error[decade], e_dd);
        worst_hybrid = fmax(worst_hybrid, isnan(e_hybrid) ? INFINITY : e_hybrid / hybrid_tolerance);
        worst_dd = fmax(worst_dd, isnan(e_dd) ? INFINITY : e_dd / DD_TOLERANCE);

        double complex fields[NUM_FIELDS];
        const int status = field_modes(ctx, &point, m, N, fields);
        for (int f = 0; f < NUM_FIELDS; f++) {
          const double e = status ? INFINITY : (double)(cabsl(fields[f] - field_reference[f][m]) / field_absolute[f]);
          field_error[decade][f] = fmax(field_error[decade][f], e);
          worst_field = fmax(worst_field, isnan(e) ? INFINITY : e / FIELD_TOLERANCE);
        }

        double re_big, im_big, re_small, im_small;
        if (effsource_equatorial_ctx_numerical_src_mode(ctx, &point, &stride_big[im], 1 << 10, &re_big, &im_big)
            || effsource_equatorial_ctx_numerical_src_mode(ctx, &point, &stride_small[im], 1 << 10, &re_small, &im_small))
          worst_stride = INFINITY;
        else
          worst_stride = fmax(worst_stride, cabs((re_big - re_small) + I * (im_big - im_small)) / absolute / STRIDE_TOLERANCE);

        if (offsets_theta[idth] == 0.1 || offsets_theta[idth] == 1e-3) {
          /* Points exactly symmetric about xp.theta, so that the rounding of
             pi/2 -+ dtheta does not enter. */
          const double dtheta = (xp.theta + offsets_theta[idth]) - xp.theta;
          struct effsource_equatorial_point point_north, point_south;
          effsource_equatorial_ctx_prepare_dd(ctx, r, xp.theta - dtheta, &point_north);
          effsource_equatorial_ctx_prepare_dd(ctx, r, xp.theta + dtheta, &point_south);
          double complex north[NUM_FIELDS], south[NUM_FIELDS];
          const double complex src_north = src_mode(ctx, &point_north, m, N);
          const double complex src_south = src_mode(ctx, &point_south, m, N);
          if (field_modes(ctx, &point_north, m, N, north) || field_modes(ctx, &point_south, m, N, south)) {
            worst_symmetry = INFINITY;
          } else {
            const double signs[NUM_FIELDS] = {1.0, 1.0, -1.0};
            /* The particle sits at xp.theta = M_PI_2, which lies 6e-17 from the
               equator, so the src of calc itself differs between the two points
               by up to about 3e-15 / dtheta of int |src| dphi. */
            double e = cabs(src_north - src_south) / absolute / fmax(1.0, 0.1 / dtheta);
            for (int f = 0; f < NUM_FIELDS; f++)
              e = fmax(e, cabs(north[f] - signs[f] * south[f]) / field_absolute[f]);
            worst_symmetry = fmax(worst_symmetry, isnan(e) ? INFINITY : e / SYMMETRY_TOLERANCE);
          }
        }

        if (m == 2) {
          const int repeats = N > 65536 ? 1 : (int)(65536 / N);
          struct effsource_equatorial_numerical_table table;
          double re, im;
          effsource_equatorial_numerical_table_init(&table, m, N);
          double start = now();
          for (int repeat = 0; repeat < repeats; repeat++) {
            struct effsource_equatorial_point timed;
            effsource_equatorial_ctx_prepare(ctx, r, theta, &timed);
            effsource_equatorial_ctx_numerical_src_mode(ctx, &timed, &table, N, &re, &im);
            sink += re;
          }
          time_numerical += (now() - start) / repeats;
          argument_failures += effsource_equatorial_ctx_numerical_src_mode(ctx, &point, &table, N / 2 * 3, &re, &im) == 0;
          argument_failures += effsource_equatorial_ctx_numerical_src_mode(ctx, &point, &table, 2 * N, &re, &im) == 0;
          effsource_equatorial_numerical_table_free(&table);
          start = now();
          for (int repeat = 0; repeat < repeats; repeat++)
            for (size_t k = 0; k < N; k++) {
              struct coordinate x = {.r = r, .theta = theta, .phi = 2.0 * M_PI * k / N, .t = 0.0};
              double PhiS, d1[4], d2[10], src;
              effsource_equatorial_ctx_calc(ctx, &x, &PhiS, d1, d2, &src);
              sink += src;
            }
          time_calc += (now() - start) / repeats;
          num_points++;
        }
      }
    }
    effsource_equatorial_free(ctx);
  }
  for (int im = 0; im < NUM_M; im++) {
    effsource_equatorial_numerical_table_free(&stride_big[im]);
    effsource_equatorial_numerical_table_free(&stride_small[im]);
  }

  printf("numerical m-modes, m in {0, 1, 2, 7, 20}, N from the rule; errors relative to int |f| dphi (smallest C1 %.1e)\n",
         min_C1);
  printf("  decade  cases (ld calc) | src: hybrid    prepare_dd | PhiS     dPhiS/dr  dPhiS/dtheta (hybrid)\n");
  for (int d = 0; d < NUM_DECADES; d++) {
    if (count[d] == 0)
      continue;
    printf("  1e-%-2d   %3d (%3d)       | %9.1e  %9.1e   | %9.1e %9.1e %9.1e\n", d, count[d], count_field_ld[d],
           hybrid_error[d], dd_error[d], field_error[d][0], field_error[d][1], field_error[d][2]);
  }
  printf("time per field point (us), m = 2: prepare + numerical %.2f, N calls of calc %.2f, ratio %.1f\n",
         1e6 * time_numerical / num_points, 1e6 * time_calc / num_points, time_calc / time_numerical);

  all_pass &= check("bad_arguments", argument_failures == 0, argument_failures);
  all_pass &= check("src_hybrid_vs_reference", worst_hybrid <= 1.0, worst_hybrid);
  all_pass &= check("src_dd_vs_reference", worst_dd <= 1.0, worst_dd);
  all_pass &= check("field_vs_reference", worst_field <= 1.0, worst_field);
  all_pass &= check("stride", worst_stride <= 1.0, worst_stride);
  all_pass &= check("equatorial_symmetry", worst_symmetry <= 1.0, worst_symmetry);
  if (sink == 1.2345)
    printf("\n");
  return all_pass ? 0 : 1;
}
