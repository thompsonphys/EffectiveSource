#!/usr/bin/env python3
"""Drive the compiled `effsource_circular` .so to emit reference values as JSON.

Used by runtests.jl for live parity checking of the Julia port. Prints one JSON
object to stdout: fixed particle params plus a list of field-point records.
"""
import json
import math
import os
import sys

# the compiled module lives at the repo root (two levels up from julia/test/)
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO)

import effsource_circular as es  # noqa: E402

es.disable_gsl_error_handler()

# --- fixed circular orbit (matches test_ctx.py) ----------------------------
M, a = 1.0, 0.5
rp, thp, php, tp = 10.0, math.pi / 2.0, math.pi / 3.0, 1.0
E = math.sqrt(1 - 2 * M / rp + a**2 / rp**2 - 2 * M * a**2 / rp**3) / math.sqrt(
    1 - 3 * M / rp + 2 * a * math.sqrt(M / rp**3)
)
L = (rp**2 - 2 * M * rp + a * math.sqrt(M * rp)) / (
    rp * math.sqrt(rp - 3 * M + 2 * a * math.sqrt(M / rp))
)

ctx = es.EffsourceContext(M, a)
xp = es.make_coordinate(t=tp, r=rp, theta=thp, phi=php)
ctx.set_particle(xp, E, L, 0.0)

MODES = [0, 1, 2, 5, 10, 20]

# grid of field points offset from the particle
records = []
for dr in (0.05, 0.1, -0.1, 0.3):
    for dth in (0.0, 0.01, -0.02):
        for dph in (0.0, 0.02, 0.05):
            x = es.make_coordinate(t=0.0, r=rp + dr, theta=thp + dth, phi=php + dph)
            rec = {
                "x": [0.0, rp + dr, thp + dth, php + dph],
                "PhiS": ctx.calc_PhiS(x),
                "calc": None,
                "phis_m": {},
                "calc_m": {},
            }
            P, d1, d2, s = ctx.calc(x)
            rec["calc"] = {"PhiS": P, "d1": list(d1), "d2": list(d2), "src": s}
            for m in MODES:
                rec["phis_m"][str(m)] = list(ctx.calc_PhiS_m(m, x))
                Pm, d1m, d2m, sm = ctx.calc_m(m, x)
                rec["calc_m"][str(m)] = {
                    "PhiS": list(Pm),
                    "d1": list(d1m),
                    "d2": list(d2m),
                    "src": list(sm),
                }
            records.append(rec)

out = {
    "params": {"M": M, "a": a, "rp": rp, "thp": thp, "php": php, "tp": tp,
               "E": E, "L": L},
    "modes": MODES,
    "records": records,
}


def sanitize(o):
    """Replace non-finite floats with None so the output is strict JSON
    (Julia's JSON parser rejects bare NaN/Infinity tokens)."""
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: sanitize(v) for k, v in o.items()}
    if isinstance(o, list):
        return [sanitize(v) for v in o]
    return o


json.dump(sanitize(out), sys.stdout)
