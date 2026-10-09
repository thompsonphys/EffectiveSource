#!/usr/bin/env python3
"""Dump effsource_equatorial coefficients and kernel outputs as JSON for the Julia parity tests."""
import json
import math
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))

import effsource_equatorial as es

MODES = [0, 1, 2, 5, 10, 20]


def sanitize(o):
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    if isinstance(o, dict):
        return {k: sanitize(v) for k, v in o.items()}
    if isinstance(o, list):
        return [sanitize(v) for v in o]
    return o


def main():
    M, a, rp, E, L, ur = (float(v) for v in sys.argv[1:7])
    ell = sys.argv[7] == "1"
    thp, php, tp = math.pi / 2.0, math.pi / 3.0, 1.0
    es.disable_gsl_error_handler()
    ctx = es.EffsourceEquatorialContext(M, a, ell_from_complement=ell)
    ctx.set_particle(es.make_coordinate(t=tp, r=rp, theta=thp, phi=php), E, L, ur)
    coeffs = {n: getattr(ctx._ctx, n) for n in dir(ctx._ctx)
              if not n.startswith("_") and isinstance(getattr(ctx._ctx, n), float)}
    records = []
    for dr in (0.05, 0.1, -0.1, 0.3):
        for dth in (0.0, 0.01, -0.02):
            for dph in (0.0, 0.02, 0.05):
                x = es.make_coordinate(t=0.0, r=rp + dr, theta=thp + dth, phi=php + dph)
                P, d1, d2, s = ctx.calc(x)
                rec = {"x": [0.0, rp + dr, thp + dth, php + dph],
                       "PhiS": ctx.calc_PhiS(x),
                       "calc": {"PhiS": P, "d1": list(d1), "d2": list(d2), "src": s},
                       "phis_m": {}, "calc_m": {}}
                for m in MODES:
                    rec["phis_m"][str(m)] = list(ctx.calc_PhiS_m(m, x))
                    Pm, d1m, d2m, sm = ctx.calc_m(m, x)
                    rec["calc_m"][str(m)] = {"PhiS": list(Pm), "d1": list(d1m),
                                             "d2": list(d2m), "src": list(sm)}
                records.append(rec)
    near = []
    for d in (1e-5, 1e-6):
        x = es.make_coordinate(t=0.0, r=rp + d, theta=thp + d, phi=php + d)
        near.append({"x": [0.0, rp + d, thp + d, php + d],
                     "phis_m": {str(m): list(ctx.calc_PhiS_m(m, x)) for m in MODES}})
    json.dump(sanitize({"coeffs": coeffs, "modes": MODES, "records": records, "near": near}), sys.stdout)


if __name__ == "__main__":
    main()
