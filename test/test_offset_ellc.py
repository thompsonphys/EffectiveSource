"""Validate the offset API + complement-fed elliptic integrals in the
equatorial (eccentric) C code.

1. Consistency: at moderate offsets, the offset entry points agree with the
   coordinate entry points to near machine precision.
2. Accuracy: in the circular limit (ur = 0), the m-mode puncture from the
   offset path matches the BigFloat truth dumped by
   julia/examples/dump_truth.jl even deep in the near zone, while the
   coordinate path degrades. Requires julia/examples/truth_circular.csv.

Run inside the spectre conda env:
    python test/test_offset_ellc.py
"""
import csv
import math
import os
import sys

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
ctx.set_particle(xp, E, L, 0.0)  # circular limit: ur = 0

failures = 0


def check(name, cond, detail=""):
    global failures
    status = "ok " if cond else "FAIL"
    print(f"[{status}] {name} {detail}")
    if not cond:
        failures += 1


# --- 1. offset vs coordinate agreement at moderate offsets ------------------
# Dyadic offsets so rp + dr is exact and both paths see bit-identical
# (dr, dtheta): any disagreement is then a real code difference, not the
# source's legitimate sensitivity to a 1-ulp change of the field point.
for m in (0, 2, 10):
    for (dr, dth) in [(0.125, 0.03125), (-0.0625, 0.015625), (0.25, -0.03125)]:
        x = ee.make_coordinate(t=0.0, r=rp + dr, theta=thp + dth, phi=php)
        pc = ctx.calc_PhiS_m(m, x)
        po = ctx.calc_PhiS_m_offset(m, dr, dth)
        num = abs(complex(*pc) - complex(*po))
        den = abs(complex(*pc))
        rel = num / den if den else num
        check(f"PhiS_m m={m} offset/coord agree", rel < 1e-12, f"rel={rel:.2e}")
        sc = ctx.calc_m(m, x)[3]
        so = ctx.calc_m_offset(m, dr, dth)[3]
        num = abs(complex(*sc) - complex(*so))
        den = abs(complex(*sc))
        rel = num / den if den else num
        check(f"src_m  m={m} offset/coord agree", rel < 1e-10, f"rel={rel:.2e}")

# real-space too
for (dr, dth, dph) in [(0.1, 0.02, 0.03), (-0.05, 0.01, -0.02)]:
    x = ee.make_coordinate(t=0.0, r=rp + dr, theta=thp + dth, phi=php + dph)
    pc = ctx.calc_PhiS(x)
    po = ctx.calc_PhiS_offset(dr, dth, dph)
    rel = abs(pc - po) / abs(pc)
    check("PhiS offset/coord agree", rel < 1e-12, f"rel={rel:.2e}")

# --- 2. near-zone accuracy vs BigFloat truth (circular limit) ---------------
truth_csv = os.path.join(
    os.path.dirname(__file__), "..", "julia", "examples", "truth_circular.csv"
)
if not os.path.exists(truth_csv):
    print(f"truth file missing ({truth_csv}); run julia dump_truth.jl first")
    sys.exit(1)

print("\nnear-zone m-mode puncture: relerr vs BigFloat truth")
print(f"{'m':>3} {'dr':>10} {'dtheta':>10} {'coord':>10} {'offset':>10}")
worst_off = 0.0
with open(truth_csv) as f:
    for row in csv.DictReader(f):
        m = int(row["m"])
        dr, dth = float(row["dr"]), float(row["dtheta"])
        truth_p = complex(float(row["RePhiS"]), float(row["ImPhiS"]))
        x = ee.make_coordinate(t=0.0, r=rp + dr, theta=thp + dth, phi=php)
        pc = complex(*ctx.calc_PhiS_m(m, x))
        po = complex(*ctx.calc_PhiS_m_offset(m, dr, dth))
        ec = abs(pc - truth_p) / abs(truth_p)
        eo = abs(po - truth_p) / abs(truth_p)
        worst_off = max(worst_off, eo)
        s = math.hypot(dr, dth)
        if s < 3e-5:  # print only the interesting near zone
            print(f"{m:>3} {dr:>10.1e} {dth:>10.1e} {ec:>10.1e} {eo:>10.1e}")

check("offset-path puncture matches BigFloat truth", worst_off < 1e-12,
      f"worst rel err = {worst_off:.2e}")

print(f"\n{failures} failure(s)")
sys.exit(1 if failures else 0)
