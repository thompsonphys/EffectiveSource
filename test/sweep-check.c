/*******************************************************************************
 * Check effsource_equatorial_ctx_sweep_calc and effsource_equatorial_ctx_sweep_src
 * against effsource_equatorial_ctx_calc.
 *
 * The calc function and both sweep routes are graded against a long double copy
 * of the calc function at
 * the same double inputs, for a set of particle states and field points with
 * C1 = (alpha20 dr^2 + alpha02 dtheta^2) / beta from about 1 down to 1e-10:
 *   1. every output, sample by sample over the azimuthal grid of the m-mode
 *      projection, normalised by the largest reference value of that output;
 *   2. the m-modes of src for m = 0..MMAX, normalised by the largest
 *      reference m-mode.
 * The sweep_calc route is graded on all outputs, the sweep_src route on src.
 * Errors are grouped by decade of C1. In each decade the worst and the median
 * error of the sweep must be at most MAX_ERROR_RATIO times those of the calc
 * function, or below ERROR_FLOOR. The sweeps use effsource_equatorial_ctx_prepare,
 * which computes the coefficients in double arithmetic where C1 is large. Those
 * coefficients carry errors of 1e-14 to 1e-11, up to three times those of the
 * calc function, so ERROR_FLOOR is set well above them and far below the
 * accuracy the (m,n)-modes need, and the sweep_src route is allowed
 * MAX_SRC_ERROR_RATIO where those coefficients are used. Where C1 is small the
 * coefficients come from double-double arithmetic and the sweep_src errors are
 * orders of magnitude below those of the calc function. Errors are compared by decade rather than
 * case by case, because the error of the calc function in a single case can be
 * far below its typical level by chance: recompiling it with fused
 * multiply-add changes single-case errors by up to a factor of 100 while the
 * decade levels change by less than 1.5. The program also reports the time
 * per evaluation of all three routes. It prints one "check <name> PASS|FAIL <value>"
 * line per check and exits with status 1 if any check fails.
 ******************************************************************************/

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "effsource_equatorial.h"
#include "effsource_equatorial_sweep.h"

#define MMAX 10
#define NUM_OUTPUTS 16
#define NUM_DECADES 11
#define MAX_ERROR_RATIO 2.0
#define MAX_SRC_ERROR_RATIO 4.0
#define MAX_CASES 1024
#define ERROR_FLOOR 1e-12

void effsource_equatorial_ctx_calc_long_double(const struct effsource_equatorial_ctx *ctx,
    const struct coordinate *x, long double *PhiS, long double *dPhiS_dx,
    long double *d2PhiS_dx2, long double *src);

static double now(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}

/* Smallest power of two >= max(64, 20 + 16 / sqrt(C1)), capped at 2^15 for the check. */
static int phi_points(double C1)
{
  const double want = 20.0 + 16.0 / sqrt(C1);
  int n = 64;
  while (n < want && n < (1 << 15))
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

static void unpack(double PhiS, const double *d1, const double *d2, double src, double *out)
{
  out[0] = PhiS;
  for (int i = 0; i < 4; i++)
    out[1 + i] = d1[i];
  for (int i = 0; i < 10; i++)
    out[5 + i] = d2[i];
  out[15] = src;
}

static void unpack_long(long double PhiS, const long double *d1, const long double *d2,
                        long double src, long double *out)
{
  out[0] = PhiS;
  for (int i = 0; i < 4; i++)
    out[1 + i] = d1[i];
  for (int i = 0; i < 10; i++)
    out[5 + i] = d2[i];
  out[15] = src;
}

static double ratio(double sweep_error, double calc_error)
{
  if (sweep_error <= ERROR_FLOOR)
    return 0.0;
  return sweep_error / fmax(calc_error, ERROR_FLOOR);
}

static int compare_doubles(const void *a, const void *b)
{
  const double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

/* Worst and median of the values of `errors` whose decade is `decade`. */
static void decade_stats(const double *errors, const int *decades, int count, int decade,
                         double *worst, double *median)
{
  double selected[MAX_CASES];
  int n = 0;
  for (int i = 0; i < count; i++)
    if (decades[i] == decade)
      selected[n++] = errors[i];
  *worst = 0.0;
  *median = 0.0;
  if (n == 0)
    return;
  qsort(selected, n, sizeof(double), compare_doubles);
  *worst = selected[n - 1];
  *median = n % 2 ? selected[n / 2] : 0.5 * (selected[n / 2 - 1] + selected[n / 2]);
}

int main(void)
{
  const double spins[] = {0.0, 0.5, 0.9};
  const double particle_r[] = {7.0, 10.0, 16.0};
  const double particle_ur[] = {0.0, 0.03, -0.03};
  const double offsets_r[] = {0.0, 0.01, -0.5, 3.0};
  const double offsets_theta[] = {0.5, 0.1, 1e-2, 1e-3, 1e-4};

  static double pointwise_calc[MAX_CASES], pointwise_sweep[MAX_CASES];
  static double mode_calc[MAX_CASES], mode_sweep[MAX_CASES];
  static double src_calc[MAX_CASES], src_rational[MAX_CASES], mode_rational[MAX_CASES];
  static int case_decade[MAX_CASES];
  int nan_mismatch = 0, cases = 0;

  for (int is = 0; is < 3; is++)
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
      const double dr = r - ctx->xp.r, dth = theta - ctx->xp.theta;
      const double C1 = (ctx->alpha20 * dr * dr + ctx->alpha02 * dth * dth) / ctx->beta;
      const int n = phi_points(C1);
      const int decade = (int)fmin(NUM_DECADES - 1, fmax(0.0, -floor(log10(C1))));
      struct effsource_equatorial_point point;
      effsource_equatorial_ctx_prepare(ctx, r, theta, &point);

      double scale[NUM_OUTPUTS] = {0}, error_calc[NUM_OUTPUTS] = {0}, error_sweep[NUM_OUTPUTS] = {0};
      long double complex mode_reference[MMAX + 1] = {0};
      double complex mode_direct[MMAX + 1] = {0}, mode_swept[MMAX + 1] = {0}, mode_ratio[MMAX + 1] = {0};
      double src_scale = 0.0, src_error_calc = 0.0, src_error_rational = 0.0;
      for (int k = 0; k < n; k++) {
        const double phi = 2.0 * M_PI * k / n;
        struct coordinate x = {.r = r, .theta = theta, .phi = phi, .t = 0.0};
        double PhiS, d1[4], d2[10], src, a[NUM_OUTPUTS];
        double sPhiS, sd1[4], sd2[10], ssrc, b[NUM_OUTPUTS];
        long double lPhiS, ld1[4], ld2[10], lsrc, l[NUM_OUTPUTS];
        effsource_equatorial_ctx_calc(ctx, &x, &PhiS, d1, d2, &src);
        effsource_equatorial_ctx_sweep_calc(ctx, &point, phi, &sPhiS, sd1, sd2, &ssrc);
        effsource_equatorial_ctx_calc_long_double(ctx, &x, &lPhiS, ld1, ld2, &lsrc);
        unpack(PhiS, d1, d2, src, a);
        unpack(sPhiS, sd1, sd2, ssrc, b);
        unpack_long(lPhiS, ld1, ld2, lsrc, l);
        const double rsrc = effsource_equatorial_ctx_sweep_src(ctx, &point, phi);
        src_scale = fmax(src_scale, fabsl(lsrc));
        src_error_calc = fmax(src_error_calc, fabsl(src - lsrc));
        src_error_rational = fmax(src_error_rational, fabsl(rsrc - lsrc));
        for (int q = 0; q < NUM_OUTPUTS; q++) {
          if (isnan(a[q]) || isnan(b[q])) {
            nan_mismatch += isnan(a[q]) != isnan(b[q]);
            continue;
          }
          scale[q] = fmax(scale[q], fabsl(l[q]));
          error_calc[q] = fmax(error_calc[q], fabsl(a[q] - l[q]));
          error_sweep[q] = fmax(error_sweep[q], fabsl(b[q] - l[q]));
        }
        for (int m = 0; m <= MMAX; m++) {
          const double complex phase = cexp(-I * m * phi);
          mode_reference[m] += lsrc * phase;
          mode_direct[m] += src * phase;
          mode_swept[m] += ssrc * phase;
          mode_ratio[m] += rsrc * phase;
        }
      }

      double case_calc = 0.0, case_sweep = 0.0;
      for (int q = 0; q < NUM_OUTPUTS; q++) {
        if (scale[q] == 0.0)
          continue;
        case_calc = fmax(case_calc, error_calc[q] / scale[q]);
        case_sweep = fmax(case_sweep, error_sweep[q] / scale[q]);
      }

      double mode_scale = 0.0, mode_error_calc = 0.0, mode_error_sweep = 0.0, mode_error_rational = 0.0;
      for (int m = 0; m <= MMAX; m++) {
        mode_scale = fmax(mode_scale, cabsl(mode_reference[m]));
        mode_error_calc = fmax(mode_error_calc, cabsl(mode_direct[m] - mode_reference[m]));
        mode_error_sweep = fmax(mode_error_sweep, cabsl(mode_swept[m] - mode_reference[m]));
        mode_error_rational = fmax(mode_error_rational, cabsl(mode_ratio[m] - mode_reference[m]));
      }
      case_decade[cases] = decade;
      pointwise_calc[cases] = case_calc;
      pointwise_sweep[cases] = case_sweep;
      mode_calc[cases] = mode_error_calc / mode_scale;
      mode_sweep[cases] = mode_error_sweep / mode_scale;
      mode_rational[cases] = mode_error_rational / mode_scale;
      src_calc[cases] = src_error_calc / src_scale;
      src_rational[cases] = src_error_rational / src_scale;
      cases++;
    }
    effsource_equatorial_free(ctx);
  }

  printf("cases %d\n", cases);
  printf("%-6s | %-23s | %-23s | %-23s | %-23s\n", "C1", "pointwise worst c / s", "pointwise median c / s",
         "m-mode worst c / s", "m-mode median c / s");
  double worst_pointwise_ratio = 0.0, worst_mode_ratio = 0.0;
  for (int decade = 0; decade < NUM_DECADES; decade++) {
    double pcw, pcm, psw, psm, mcw, mcm, msw, msm;
    decade_stats(pointwise_calc, case_decade, cases, decade, &pcw, &pcm);
    decade_stats(pointwise_sweep, case_decade, cases, decade, &psw, &psm);
    decade_stats(mode_calc, case_decade, cases, decade, &mcw, &mcm);
    decade_stats(mode_sweep, case_decade, cases, decade, &msw, &msm);
    if (pcw == 0.0)
      continue;
    printf("1e-%-3d | %.2e / %.2e     | %.2e / %.2e     | %.2e / %.2e     | %.2e / %.2e\n", decade,
           pcw, psw, pcm, psm, mcw, msw, mcm, msm);
    worst_pointwise_ratio = fmax(worst_pointwise_ratio, fmax(ratio(psw, pcw), ratio(psm, pcm)));
    worst_mode_ratio = fmax(worst_mode_ratio, fmax(ratio(msw, mcw), ratio(msm, mcm)));
  }

  printf("%-6s | %-23s | %-23s | %-23s | %-23s\n", "C1", "src worst c / r", "src median c / r",
         "m-mode worst c / r", "m-mode median c / r");
  double worst_src_ratio = 0.0, worst_rational_mode_ratio = 0.0;
  for (int decade = 0; decade < NUM_DECADES; decade++) {
    double pcw, pcm, prw, prm, mcw, mcm, mrw, mrm;
    decade_stats(src_calc, case_decade, cases, decade, &pcw, &pcm);
    decade_stats(src_rational, case_decade, cases, decade, &prw, &prm);
    decade_stats(mode_calc, case_decade, cases, decade, &mcw, &mcm);
    decade_stats(mode_rational, case_decade, cases, decade, &mrw, &mrm);
    if (pcw == 0.0)
      continue;
    printf("1e-%-3d | %.2e / %.2e     | %.2e / %.2e     | %.2e / %.2e     | %.2e / %.2e\n", decade,
           pcw, prw, pcm, prm, mcw, mrw, mcm, mrm);
    worst_src_ratio = fmax(worst_src_ratio, fmax(ratio(prw, pcw), ratio(prm, pcm)));
    worst_rational_mode_ratio = fmax(worst_rational_mode_ratio, fmax(ratio(mrw, mcw), ratio(mrm, mcm)));
  }

  int failed = 0;
  const int nan_ok = nan_mismatch == 0;
  printf("check nan_pattern %s %d\n", nan_ok ? "PASS" : "FAIL", nan_mismatch);
  failed |= !nan_ok;
  const int pointwise_ok = worst_pointwise_ratio <= MAX_ERROR_RATIO;
  printf("check pointwise_error_ratio %s %.2f\n", pointwise_ok ? "PASS" : "FAIL", worst_pointwise_ratio);
  failed |= !pointwise_ok;
  const int mode_ok = worst_mode_ratio <= MAX_ERROR_RATIO;
  printf("check mmode_error_ratio %s %.2f\n", mode_ok ? "PASS" : "FAIL", worst_mode_ratio);
  failed |= !mode_ok;
  const int src_ok = worst_src_ratio <= MAX_SRC_ERROR_RATIO;
  printf("check rational_src_error_ratio %s %.2f\n", src_ok ? "PASS" : "FAIL", worst_src_ratio);
  failed |= !src_ok;
  const int rational_mode_ok = worst_rational_mode_ratio <= MAX_SRC_ERROR_RATIO;
  printf("check rational_mmode_error_ratio %s %.2f\n", rational_mode_ok ? "PASS" : "FAIL",
         worst_rational_mode_ratio);
  failed |= !rational_mode_ok;

  {
    double E, L;
    circular_constants(0.9, 10.0, &E, &L);
    struct effsource_equatorial_ctx *ctx = effsource_equatorial_create(1.0, 0.9);
    struct coordinate xp = {.r = 10.0, .theta = M_PI_2, .phi = 0.0, .t = 0.0};
    effsource_equatorial_ctx_set_particle(ctx, &xp, E, L, 0.02);
    const int points = 200, n = 1024;
    double sum_direct = 0.0, sum_sweep = 0.0;
    double t0 = now();
    for (int j = 0; j < points; j++) {
      struct coordinate x = {.r = 8.0 + 4.0 * j / points, .theta = M_PI_2 - 0.003 - 0.2 * (j % 13) / 13.0,
                             .phi = 0.0, .t = 0.0};
      for (int k = 0; k < n; k++) {
        double PhiS, d1[4], d2[10], src;
        x.phi = 2.0 * M_PI * k / n;
        effsource_equatorial_ctx_calc(ctx, &x, &PhiS, d1, d2, &src);
        sum_direct += src;
      }
    }
    const double t_direct = (now() - t0) / ((double)points * n);
    struct effsource_equatorial_point point;
    double t_prepare = 0.0;
    t0 = now();
    for (int j = 0; j < points; j++) {
      const double tp = now();
      effsource_equatorial_ctx_prepare(ctx, 8.0 + 4.0 * j / points,
                                             M_PI_2 - 0.003 - 0.2 * (j % 13) / 13.0, &point);
      t_prepare += now() - tp;
      for (int k = 0; k < n; k++) {
        double PhiS, d1[4], d2[10], src;
        effsource_equatorial_ctx_sweep_calc(ctx, &point, 2.0 * M_PI * k / n, &PhiS, d1, d2, &src);
        sum_sweep += src;
      }
    }
    const double t_sweep = (now() - t0 - t_prepare) / ((double)points * n);
    double sum_rational = 0.0;
    t0 = now();
    for (int j = 0; j < points; j++) {
      effsource_equatorial_ctx_prepare(ctx, 8.0 + 4.0 * j / points,
                                             M_PI_2 - 0.003 - 0.2 * (j % 13) / 13.0, &point);
      for (int k = 0; k < n; k++)
        sum_rational += effsource_equatorial_ctx_sweep_src(ctx, &point, 2.0 * M_PI * k / n);
    }
    const double t_rational = (now() - t0 - t_prepare) / ((double)points * n);
    printf("time per evaluation: calc %.1f ns, sweep_calc %.1f ns (%.2fx); prepare %.0f ns per field point\n",
           1e9 * t_direct, 1e9 * t_sweep, t_direct / t_sweep, 1e9 * t_prepare / points);
    printf("time per evaluation: sweep_src %.1f ns (%.2fx)\n", 1e9 * t_rational, t_direct / t_rational);
    double t_dd = 0.0, t_double = 0.0;
    t0 = now();
    for (int j = 0; j < points; j++)
      effsource_equatorial_ctx_prepare_dd(ctx, 8.0 + 4.0 * j / points, M_PI_2 - 0.003 - 0.2 * (j % 13) / 13.0, &point);
    t_dd = (now() - t0) / points;
    t0 = now();
    for (int j = 0; j < points; j++)
      effsource_equatorial_ctx_prepare_double(ctx, 8.0 + 4.0 * j / points, M_PI_2 - 0.003 - 0.2 * (j % 13) / 13.0, &point);
    t_double = (now() - t0) / points;
    printf("time per prepare: double-double %.0f ns, double %.0f ns\n", 1e9 * t_dd, 1e9 * t_double);
    printf("sum relative difference: sweep_calc %.2e, sweep_src %.2e\n",
           fabs(sum_direct - sum_sweep) / fabs(sum_direct), fabs(sum_direct - sum_rational) / fabs(sum_direct));
    effsource_equatorial_free(ctx);
  }
  return failed;
}
