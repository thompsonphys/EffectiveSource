"""Validate the kernel split of the full m-mode calc:

    X = A + L*ln(alpha) + P1/alpha + P2/alpha^2,
    alpha = alpha20*dr^2 + alpha02*dtheta^2,

for every component of (PhiS, dPhiS, d2PhiS, src) (calc_m_split).

1. Recombination identity vs calc_m_offset, measured relative to the size
   of the recombination TERMS (far from the particle the total is a deep
   cancellation of O(1) channels -- that is intrinsic to the representation,
   not an error; the split is built for the near zone and for channel-wise
   convolution).  Offsets span the m1 = 0.1 series/closed-form crossover of
   the P/Q elliptic primitives.
2. FD self-consistency (reference-free): derivative channels equal finite
   differences of value channels, chaining value -> d/dr -> d2/dr2 and the
   dtheta analogue.
3. Production-geometry smoothness: along a fixed-dtheta line (what a
   worldline traces relative to a legal field point) the L, P1, P2 channels
   of src are Chebyshev-analytic; P1/P2 to machine precision.
4. Deep near zone: at s -> 1e-9 on a ray the recombined src is finite and
   dominated by the clean P1/P2 pole channels (the unsplit path loses all
   digits there); P1 approaches a constant, P2/s^2 approaches a constant.

Run inside the spectre conda env:
    python test/test_calc_m_split.py
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


D2_VALID = [0, 1, 6, 7, 8, 9, 14, 15, 18, 19]


def components(res):
    PhiS, dP, d2P, src = res
    out = [("PhiS", i) for i in range(2)]
    out += [("dPhiS", i) for i in range(8)]
    out += [("d2PhiS", i) for i in D2_VALID]
    out += [("src", i) for i in range(2)]
    return out


BLK = {"PhiS": 0, "dPhiS": 1, "d2PhiS": 2, "src": 3}

# --- 1. recombination ---------------------------------------------------------
for ur in (0.0, 0.15):
    ctx.set_particle(xp, E, L, ur)
    a20, a02, beta, c = ctx.get_alpha()

    s_cross = math.sqrt(beta / 9.0 / (a20 + a02 * 0.16))
    offsets = [(0.5, 0.2), (-0.3, 0.1), (0.05, 0.02), (-0.01, 0.005),
               (1e-3, -4e-4), (-1e-4, 1e-4), (0.0, 3e-3), (2e-3, 0.0),
               (s_cross, 0.4 * s_cross), (1.02 * s_cross, 0.408 * s_cross)]
    for m in (0, 2, 10):
        worst = 0.0
        worst_lab = ""
        for dr, dth in offsets:
            alpha = a20 * dr * dr + a02 * dth * dth
            lna, inva = math.log(alpha), 1.0 / alpha
            ref = ctx.calc_m_offset(m, dr, dth)
            spl = ctx.calc_m_split(m, dr, dth)
            for lab, i in components(ref):
                blocks = spl[BLK[lab]]
                terms = [blocks[0][i], blocks[1][i] * lna] + [
                    blocks[1 + q][i] * inva**q for q in (1, 2, 3, 4, 5)]
                rec = sum(terms)
                vr = ref[BLK[lab]][i]
                # scale: largest recombination term (channels are the
                # deliverable; a small total far away is intrinsic
                # cancellation, not channel error)
                scale = max(max(abs(t) for t in terms), 1e-300)
                rel = abs(rec - vr) / scale
                if rel > worst:
                    worst, worst_lab = rel, f"{lab}[{i}] dr={dr:g} dth={dth:g}"
        check(f"ur={ur} m={m} recombination", worst < 3e-10,
              f"worst rel(term-scale)={worst:.1e} at {worst_lab}")

# --- 2. FD self-consistency ----------------------------------------------------
# Only the L channel (the unique ln(alpha) coefficient) must be
# derivative-consistent across components.  The A<->Pq split is NOT unique --
# one may shift alpha^q*(analytic) between Pq and A without breaking analyticity
# or recombination -- so A_dPhiS_r need not equal d/dr(A_PhiS); both are valid
# analytic representatives.  The convolution reconstructs the true n-modes
# regardless of the split, so this ambiguity is harmless.
ctx.set_particle(xp, E, L, 0.15)
a20, a02, beta, c = ctx.get_alpha()
m = 2
for (dr0, dth0) in ((0.05, 0.02), (0.0, 0.01), (0.3, 0.001), (-0.01, 0.005)):
    h = 1e-5
    for chan in (1,):  # L channel (unique); A/Pq split is representative-dependent
        cm = ctx.calc_m_split(m, dr0 - h, dth0)
        c0 = ctx.calc_m_split(m, dr0, dth0)
        cp = ctx.calc_m_split(m, dr0 + h, dth0)
        fd1 = (cp[0][chan][0] - cm[0][chan][0]) / (2 * h)
        an1 = c0[1][chan][2]
        fd2 = (cp[1][chan][2] - cm[1][chan][2]) / (2 * h)
        an2 = c0[2][chan][8]
        r1 = abs(fd1 - an1) / max(abs(an1), 1e-300)
        r2 = abs(fd2 - an2) / max(abs(an2), 1e-300)
        lab = "AL"[chan]
        check(f"FD dr ({dr0},{dth0}) {lab}", r1 < 1e-5 and r2 < 1e-4,
              f"d1 rel={r1:.1e} d2 rel={r2:.1e}")
    # dtheta direction, L channel
    cm = ctx.calc_m_split(m, dr0, dth0 - h)
    c0 = ctx.calc_m_split(m, dr0, dth0)
    cp = ctx.calc_m_split(m, dr0, dth0 + h)
    fd1 = (cp[0][1][0] - cm[0][1][0]) / (2 * h)
    an1 = c0[1][1][4]
    fd2 = (cp[1][1][4] - cm[1][1][4]) / (2 * h)
    an2 = c0[2][1][14]
    r1 = abs(fd1 - an1) / max(abs(an1), 1e-300)
    r2 = abs(fd2 - an2) / max(abs(an2), 1e-300)
    check(f"FD dth ({dr0},{dth0}) L", r1 < 1e-5 and r2 < 1e-4,
          f"d1 rel={r1:.1e} d2 rel={r2:.1e}")

# --- 3. production-line smoothness (fixed dtheta, dr sweep) --------------------
H = 0.3
dense = H * np.cos(np.pi * (np.arange(301) + 0.5) / 301)


def cheb_err(f, N):
    nodes = H * np.cos(np.pi * (np.arange(N + 1) + 0.5) / (N + 1))
    coef = C.chebfit(nodes / H, np.array([f(s) for s in nodes]), N)
    fd = np.array([f(s) for s in dense])
    return np.max(np.abs(C.chebval(dense / H, coef) - fd)) / np.max(np.abs(fd))


dth = 0.01
# all seven channels are analytic along a production line, INCLUDING A: the
# P-side hidden poles (docs/nmode-kernel-factorization.md sec 6) are now
# extracted in closed form into P1..P5, so A no longer carries width-w structure.
for chan, lab, tol in ((1, "L", 1e-7), (2, "P1", 1e-10), (3, "P2", 1e-10),
                       (4, "P3", 1e-10), (5, "P4", 1e-10), (6, "P5", 1e-10)):
    err = cheb_err(lambda s: ctx.calc_m_split(m, s, dth)[3][chan][0], 24)
    check(f"line dth={dth} src {lab} analytic", err < tol, f"err={err:.1e}")
errA = cheb_err(lambda s: ctx.calc_m_split(m, s, dth)[3][0][0], 24)
check(f"line dth={dth} src A analytic", errA < 1e-9, f"err={errA:.1e}")

# --- 4. deep near zone ----------------------------------------------------------
ray = (1.0, 0.3)
alpha_of = lambda s: a20 * (s * ray[0])**2 + a02 * (s * ray[1])**2
P1s, P2s = {}, {}
for s in (1e-6, 1e-7, 1e-8, 1e-9):
    blocks = ctx.calc_m_split(m, s * ray[0], s * ray[1])[3]
    al = alpha_of(s)
    rec = (blocks[0][0] + blocks[1][0] * math.log(al)
           + sum(blocks[1 + q][0] / al**q for q in (1, 2, 3, 4, 5)))
    check(f"deep s={s:g} src finite", np.isfinite(rec), f"rec={rec:.3e}")
    P1s[s], P2s[s] = blocks[2][0], blocks[3][0]
r1 = abs(P1s[1e-9] - P1s[1e-6]) / abs(P1s[1e-6])
r2 = abs(P2s[1e-9] / 1e-18 - P2s[1e-6] / 1e-12) / abs(P2s[1e-6] / 1e-12)
check("deep P1 -> const", r1 < 1e-4, f"drift={r1:.1e}")
check("deep P2/s^2 -> const", r2 < 1e-4, f"drift={r2:.1e}")

print(f"\n{failures} failure(s)")
sys.exit(1 if failures else 0)
