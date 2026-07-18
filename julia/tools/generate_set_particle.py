#!/usr/bin/env python3
"""Generate the Julia coefficient-assignment block for set_particle! from the C
source (lines 841-896: the A#### coefficients plus alpha20/alpha02/beta)."""
from pathlib import Path
from c2jl import convert_expr

HERE = Path(__file__).resolve().parent
SRC_C = (HERE.parent.parent / "kerr-circular-ctx.c").read_text().splitlines()

out = []
for lineno in range(841 - 1, 896):  # 0-based slice covering C lines 841..896
    raw = SRC_C[lineno].strip()
    if not raw.startswith("ctx->"):
        continue
    assert raw.endswith(";"), f"line {lineno+1} not terminated: {raw[:40]}"
    lhs, rhs = raw[:-1].split("=", 1)
    name = lhs.strip()[len("ctx->"):].strip()   # e.g. A006
    jl = convert_expr(rhs.strip())
    # emit into a local variable; set_particle! stores it into the context
    # (rounded to the context type) after computing in the setup precision.
    out.append(f"    {name} = {jl}")

Path(HERE.parent / "src" / "_set_particle_body.jl").write_text("\n".join(out) + "\n")
print(f"Wrote {len(out)} assignments")
