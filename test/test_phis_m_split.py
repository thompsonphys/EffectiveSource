"""Validate the kernel split of the m-mode singular field:

    PhiS_m = G0 + GL*ln(alpha),  alpha = alpha20*dr^2 + alpha02*dtheta^2,

with G0, GL analytic across the particle (calc_PhiS_m_split).

1. Recombination identity: G0 + GL*ln(alpha) reproduces calc_PhiS_m_offset
   to near machine precision over a wide range of offsets.
2. Normalization: the log-derivative of PhiS_m along a scaling ray
   s*(dr, dtheta) tends to 2*GL as s -> 0.
3. Analyticity: G0(s), GL(s) along the ray converge exponentially under
   Chebyshev interpolation through s = 0 (including the pure-equatorial ray),
   while PhiS_m itself is log-singular there.

The split is a pure algebraic identity on the ctx coefficients, so it is
exercised both at a physical circular particle (ur = 0) and with ur != 0
(off-shell coefficients are still a valid algebra check).

Run inside the spectre conda env:
    python test/test_phis_m_split.py
"""
import math
import os
import sys

import numpy as np
from numpy.polynomial import chebyshev as C

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import effsource_equatorial as ee

M, a = 1.0, 0.5
rp, thp, php = 10.0, math.pi / 2, math.pi / 3
E = math.sqrt(1 - 2 * M / rp + a**2 / rp**2 - 2 * M * a**2 / rp**3) / math.sqrt(
    1 - 3 * M / rp + 2 * a * math.sqrt(M / rp**3)
)
L = (rp**2 - 2 * M * rp + a * math.sqrt(M * rp)) / (
    rp * math.sqrt(rp - 3 * M + 2 * a * math.sqrt(M / rp))
)

ee.disable_gsl_error_handler()
ctx = ee.EffsourceEquatorialContext(M, a)
xp = ee.make_coordinate(t=1.0, r=rp, theta=thp, phi=php)

failures = 0


def check(name, cond, detail=""):
    global failures
    status = "ok " if cond else "FAIL"
    print(f"[{status}] {name} {detail}")
    if not cond:
        failures += 1


def split_c(m, dr, dth):
    g0, gl = ctx.calc_PhiS_m_split(m, dr, dth)
    return complex(*g0), complex(*gl)


for ur in (0.0, 0.15):
    ctx.set_particle(xp, E, L, ur)
    a20, a02, beta, c = ctx.get_alpha()
    check(f"ur={ur} alpha/beta positive", a20 > 0 and a02 > 0 and beta > 0,
          f"alpha20={a20:.3e} alpha02={a02:.3e} beta={beta:.3e}")

    # --- 1. recombination identity ---------------------------------------
    offsets = [(0.5, 0.2), (0.05, 0.02), (-0.01, 0.0), (1e-4, -3e-5),
               (-1e-6, 1e-6), (1e-8, 0.0), (0.0, 1e-8)]
    for m in (0, 2, 10):
        worst = 0.0
        for dr, dth in offsets:
            alpha = a20 * dr * dr + a02 * dth * dth
            g0, gl = split_c(m, dr, dth)
            rec = g0 + gl * math.log(alpha)
            ref = complex(*ctx.calc_PhiS_m_offset(m, dr, dth))
            rel = abs(rec - ref) / abs(ref)
            worst = max(worst, rel)
        check(f"ur={ur} m={m} recombination", worst < 1e-12,
              f"worst rel={worst:.1e}")

    # --- 2. GL normalization: d PhiS / d ln s -> 2*GL as s -> 0 ----------
    m = 2
    dr0, dth0 = 0.03, 0.012
    s = 1e-5
    p1 = complex(*ctx.calc_PhiS_m_offset(m, s * dr0, s * dth0))
    p2 = complex(*ctx.calc_PhiS_m_offset(m, 2 * s * dr0, 2 * s * dth0))
    _, gl = split_c(m, s * dr0, s * dth0)
    rel = abs((p2 - p1) / (2 * math.log(2.0)) - gl) / abs(gl)
    check(f"ur={ur} GL normalization", rel < 1e-6, f"rel={rel:.1e}")

# --- 3. analyticity of G0, GL through the particle ---------------------------
ctx.set_particle(xp, E, L, 0.15)
H = 0.3
m = 2


def cheb_err(fvals_nodes, fvals_dense, N, dense_x):
    coef = C.chebfit(np.cos(np.pi * (np.arange(N + 1) + 0.5) / (N + 1)),
                     fvals_nodes, N)
    pred = C.chebval(dense_x / H, coef)
    return np.max(np.abs(pred - fvals_dense)) / np.max(np.abs(fvals_dense))


dense = H * np.cos(np.pi * (np.arange(301) + 0.5) / 301)
for ray in ((1.0, 0.4), (1.0, 0.0)):
    dg0 = np.array([split_c(m, s * ray[0], s * ray[1])[0].real for s in dense])
    dgl = np.array([split_c(m, s * ray[0], s * ray[1])[1].real for s in dense])
    dph = np.array([complex(*ctx.calc_PhiS_m_offset(
        m, s * ray[0], s * ray[1])).real for s in dense])
    errs = {}
    for N in (8, 16, 32):
        nodes = H * np.cos(np.pi * (np.arange(N + 1) + 0.5) / (N + 1))
        g0n = np.array([split_c(m, s * ray[0], s * ray[1])[0].real
                        for s in nodes])
        gln = np.array([split_c(m, s * ray[0], s * ray[1])[1].real
                        for s in nodes])
        phn = np.array([complex(*ctx.calc_PhiS_m_offset(
            m, s * ray[0], s * ray[1])).real for s in nodes])
        errs[N] = (cheb_err(g0n, dg0, N, dense),
                   cheb_err(gln, dgl, N, dense),
                   cheb_err(phn, dph, N, dense))
        print(f"    ray={ray} N={N:2d} cheb-err G0={errs[N][0]:.1e} "
              f"GL={errs[N][1]:.1e} PhiS={errs[N][2]:.1e}")
    check(f"ray={ray} G0 analytic", errs[32][0] < 1e-10,
          f"err={errs[32][0]:.1e}")
    check(f"ray={ray} GL analytic", errs[32][1] < 1e-10,
          f"err={errs[32][1]:.1e}")
    check(f"ray={ray} PhiS log-singular (control)", errs[32][2] > 1e-4,
          f"err={errs[32][2]:.1e}")

print(f"\n{failures} failure(s)")
sys.exit(1 if failures else 0)
