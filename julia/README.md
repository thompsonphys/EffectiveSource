# EffSourceCircular.jl

A Julia port of `kerr-circular-ctx.c` — the scalar singular field and effective
source for a point particle on a **circular equatorial geodesic** in Kerr
spacetime (Barry Wardell's *EffectiveSource*).

The kernels are **generic over the numeric type `T`**, so the same physics runs
in `Float64`, `Float32`, `Float32sr`/`BFloat16sr` (stochastic rounding, via
[StochasticRounding.jl](https://github.com/milankl/StochasticRounding.jl)),
`BigFloat`, etc. This makes it a convenient testbed for numeric-precision
experiments.

## Usage

```julia
using EffSourceCircular

ctx = effsource_create(1.0, 0.5)                      # M, a
xp  = Coordinate{Float64}(r=10.0, theta=π/2, phi=π/3, t=1.0)
E, L = kerr_circular_EL(1.0, 0.5, 10.0)               # circular-geodesic constants
set_particle!(ctx, xp, E, L)

x = Coordinate{Float64}(r=10.1, theta=π/2+0.01, phi=π/3+0.02, t=0.0)

phis(ctx, x)                       # real-space singular field
PhiS, dPhiS, d2PhiS, src = calc(ctx, x)   # + derivatives + effective source
phis_m(ctx, 2, x)                  # m-mode field (Complex)
calc_m(ctx, 2, x)                  # m-mode field, derivatives, source (Complex)
```

Output conventions match the C code:
`dPhiS` is `(∂t, ∂r, ∂θ, ∂φ)`; `d2PhiS` is
`(tt, tr, tθ, tφ, rr, rθ, rφ, θθ, θφ, φφ)`; `src` is the d'Alembertian.

To run at another precision, just build the context and field point in that
type (see `examples/stochastic_rounding.jl`).

**Orbital constants.** `E` and `L` must be those of the circular geodesic at the
particle radius, and `kerr_circular_EL(M, a, r)` returns them. With other values the
effective source is not small, and accuracy statements about it do not carry over
to the physical problem. Rounding the constants to Float64 is itself an error
source near the particle. It leaves an off-shell residual in the real-space
source that grows like d⁻⁴ for an offset d along (1,1,1): relative 1e-12 at
d = 0.1, 1e-4 at 1e-3, above 1 below 1e-4. The m = 2 mode changes by only 3e-10
at d = 1e-6. A BigFloat reference therefore needs the constants in BigFloat,
`kerr_circular_EL(big(M), big(a), big(r))`, and a particle in the equatorial
plane at that precision, `theta=big(pi)/2`. The Float64 value `π/2` converted to
BigFloat lies 6e-17 off the plane, which changes the BigFloat source by a relative
1.7e-6 at d = 1e-3 and 1.7e3 at d = 1e-6. The tests and examples build every
BigFloat context this way, so the Float64 error measured against such a reference
includes the rounding of the constants.

### Reducing numerical cancellation

Two knobs help where this code is ill-conditioned (near the particle, and in the
coefficient setup):

**Offset API — avoid the `r - rp` cancellation.** Near the particle, forming
`dr = r - rp` from absolute coordinates loses digits catastrophically. Pass the
offset directly instead:

```julia
off = Offset{Float64}(dr=1e-6, dtheta=0.0, dphi=0.0)
phis(ctx, off);  calc(ctx, off);  phis_m(ctx, 2, off);  calc_m(ctx, 2, off)
```

Every kernel has both a `Coordinate` and an `Offset` method; they share an
identical core. In the near zone this keeps the puncture and its derivatives at
machine precision (offset 1e-6 along (1,1,1): PhiS 8.6e-11 → 1.1e-16, ∂r PhiS
4.9e-10 → 7.2e-16). It does not rescue the effective source there: with geodesic
constants the source is a near-total cancellation, and its Float64 relative error
at that offset is about 2e11 with either API. `Offset(x, xp)` builds one
from absolute coordinates (no accuracy gain — the subtraction still happens).

**Extended-precision coefficient setup.** The `A####` coefficients are
ill-conditioned alternating sums, computed once per particle. Evaluate them in a
higher precision and store rounded to the working type:

```julia
set_particle!(ctx, xp, E, L; setup=BigFloat)   # or Double64
```

Default `setup` is the context type, reproducing the direct computation exactly.
At the reference orbit, `setup=BigFloat` makes all 53 coefficients
correctly-rounded `Float64` (vs 24/53 with the default) — a last-ULP gain here
that grows in strong-field regimes. Both knobs preserve the working type
(including stochastic-rounding types).

**Compensated summation — `comp=true` (opt-in).** `calc`, `phis_m`, `calc_m`
take a `comp::Bool=false` keyword; `comp=true` sums the m-mode numerator and the
effective-source numerator with compensated (Neumaier) summation:

```julia
calc_m(ctx, 2, off; comp=true)     # compensated source numerator
```

Default is **off**, so results are byte-identical to the naive C reference unless
you ask. The benefit is real but *regime-dependent*: it helps the **source
numerator** where summation error dominates (m = 2 source:
≈40× in Float32 at offset 0.1, ≈6× in Float64 at offset 0.3), gives little for the puncture, and does **not** help
the deep near-zone, where the error is input cancellation rather than summation
(there it can even be neutral or worse — use `Double64`/`BigFloat` or the offset
API instead). It applies only to IEEE floats; stochastic-rounding and `BigFloat`
types fall back to a plain sum (Kahan/Neumaier's error-free transform assumes
round-to-nearest, which SR breaks), so `comp=true` is a safe no-op there.

## How the port is generated

The dense algebraic kernels are **machine-translated** from the C source rather
than hand-transcribed, so thousands of characters of expansion coefficients
carry no transcription risk. The `tools/` scripts convert C → Julia:

| script | output | what it does |
|---|---|---|
| `extract_reei.py` | `src/reei_data.jl` | the `ReEI[21][2][5][27]` elliptic table (5670 coeffs) |
| `generate_set_particle.py` / `generate_struct.py` | `src/_ctx_struct.jl`, `src/_set_particle*.jl` | context struct + coefficient setup |
| `generate_functions.py` | `src/_generated_kernels.jl` | `phis`, `calc`, `phis_m`, `calc_m` |
| `c2jl.py`, `transpile.py` | (library) | the `pow()→^`, float-literal, and body transpiler |

Regenerate everything after touching the C source:

```bash
bash tools/regenerate.sh
```

### Numeric-genericity details

* `pow(x, n)` → `x^n`; **half-integer** powers `pow(x, 3.5)` → `x^3 * sqrt(x)`
  so the result stays in `T` (a `Float64` exponent would promote it).
* `Float64` literals are wrapped as `T(...)`; integer literals are left as-is
  (they promote cleanly). This is what keeps `Float32sr` from decaying to
  `Float64` mid-expression.
* Complete elliptic integrals (used by the `*_m` calls) are computed by the
  **arithmetic–geometric mean** (`_agm_KE` in `src/EffSourceCircular.jl`) rather
  than a special-function library. AGM is pure `+`/`*`/`sqrt`, so it preserves
  `T` exactly — arbitrary precision in `BigFloat`, and stochastic rounding
  propagates through the elliptic path too. It is bit-identical to
  `SpecialFunctions.ellipk/ellipe` in `Float64`. (`SpecialFunctions` is limited
  to `Float32`/`Float64` — no `BigFloat`, and it silently drops SR types — which
  is why it is not used.) The GSL **modulus** convention is kept: the argument is
  `k`, the integral is evaluated at parameter `m = k²`.

## Tests

```bash
julia --project test/runtests.jl
```

`test/ref_dump.py` drives the compiled `effsource_circular` `.so` live and the
Julia `Float64` port is checked against it point-by-point (fields and derivatives to
~1e-12, and the m-mode source, a near-total cancellation for geodesic constants,
to ~1e-4). Further test sets assert: no
`Float64` contamination (a `Float32` input yields `Float32` output throughout);
arbitrary-precision AGM elliptic integrals (Legendre's relation to <1e-70 in
`BigFloat`); the offset API's near-particle accuracy gain; and the extended
`setup=` precision. The `BigFloat` coverage adds a one-time JIT-compilation cost,
so the suite takes a few minutes on a cold run.

Requires a `python3` that can `import effsource_circular` (the CPython-3.12 `.so`
at the repo root). Override the interpreter with `PYTHON=/path/to/python3`.

## Scope

Circular orbits only. The eccentric equatorial code is ported in `EffSourceEquatorial/`, which reuses `Coordinate`, `Offset`, the AGM elliptic integrals and the `ReEI` table from this package.
