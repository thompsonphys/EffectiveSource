#!/usr/bin/env bash
# Regenerate all machine-translated Julia sources from kerr-circular-ctx.c.
set -euo pipefail
cd "$(dirname "$0")"
python3 c2jl.py                    # transpiler self-test
python3 extract_reei.py            # -> src/reei_data.jl
python3 generate_set_particle.py   # -> src/_set_particle_body.jl
python3 generate_struct.py         # -> src/_ctx_struct.jl, src/_set_particle.jl
python3 generate_functions.py      # -> src/_generated_kernels.jl
echo "Regeneration complete."
