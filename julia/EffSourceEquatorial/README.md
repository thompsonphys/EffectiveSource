# EffSourceEquatorial.jl

A Julia port of `kerr-equatorial-ctx.c` and its three coefficient files: the scalar singular field, its derivatives and the effective source for a point particle on an eccentric equatorial geodesic in Kerr spacetime. The kernels and the coefficient set-up are generic over the numeric type `T` (Float64, Float32, Float32sr, Double64, BigFloat).

## Usage

```julia
using EffSourceEquatorial

M, a, p, e = 1.0, 0.5, 10.0, 0.3
E, L = kerr_equatorial_EL(M, a, p, e)
r = 9.0
ur = kerr_equatorial_ur(M, a, E, L, r)

ctx = effsource_equatorial_create(M, a)
set_particle!(ctx, Coordinate{Float64}(r=r, theta=π/2, phi=0.0, t=0.0), E, L, ur)

x = Coordinate{Float64}(r=r + 0.3, theta=π/2 + 0.02, phi=0.1, t=0.0)
phis(ctx, x)
PhiS, dPhiS, d2PhiS, src = calc(ctx, x)
phis_m(ctx, 2, x)
PhiS_m, dPhiS_m, d2PhiS_m, src_m = calc_m(ctx, 2, x)
calc_m(ctx, 2, Offset{Float64}(dr=1e-4, dtheta=0.0))
```

The output layout follows the C code: `dPhiS` is `(∂t, ∂r, ∂θ, ∂φ)` and `d2PhiS` is `(tt, tr, tθ, tφ, rr, rθ, rφ, θθ, θφ, φφ)`. The mixed entries tr, tθ, rθ, rφ and θφ are `NaN`, as in C. `effsource_equatorial_create(M, a; ell_from_complement=true)` evaluates the m-mode elliptic integrals from the complementary parameter. `set_particle!(...; setup=BigFloat)` evaluates the coefficients in BigFloat and stores them rounded to the context type. For Float32 and Float32sr contexts, use `setup=Float64`. In Float32 the set-up overflows: at a = 0.5, p = 10, e = 0.3 and r = 9, 24 of the 280 coefficients are not finite.

## Coefficient set-up runs in the interpreter

The 280 coefficients are generated as one Julia function each, in a submodule compiled with `compile=min infer=false optimize=0`. Compiling the largest coefficient for BigFloat takes 144 s, and the cost grows about quadratically with expression size. Interpreted, one `set_particle!` takes about 0.1 s in Float64 and 0.3 s in BigFloat. A single interpreted function per coefficient file crashes the interpreter, so the per-coefficient layout is required.

## Accuracy notes

- The m-mode kernels use the `ReEI` and `ImEI` tables. Float64, Float32, Float32sr and Double64 use the Float64 values of the C code. BigFloat uses exact rationals derived in Mathematica (`tools/derive_elliptic_tables.wls`), so BigFloat m-modes are exact to working precision: at 256 bits they agree with φ-quadratures of the real-space fields to 1e-68 or better for m ≤ 20. Double64 is limited by the tables. In the circular code at offset (2, 0.3) the m = 20 mode has a table error near 6e-9.
- On the geodesic the puncture equals `notebooks/ScalarSingularField.nb`, and `test/test_independent.jl` checks this to 1e-40. The reference is `TwistedPhiS`, the notebook's re-expansion in the twisted angle Δφ − c·Δr that the code uses. At `ur = 0` it equals `ΦsKerrEquatorial`. It is sampled at nine points on two eccentric orbits (a = 0.5 and 0.9, both radial directions and turning points) and on a circular orbit, and the m-modes are checked against `NIntegrate` of the twisted form.
- The d'Alembertian assumes `M = 1`, as in C.
- The puncture is defined for geodesic `E`, `L` and `ur`. With other values the effective source is not small, and its error statistics do not represent the physical problem.
- Rounding `E`, `L` and `ur` to Float64 is itself an error source near the particle. In the circular code the resulting off-shell residual in the real-space source grows like d⁻⁴ (relative 1e-12 at d = 0.1, 1e-4 at 1e-3, above 1 below 1e-4), while the m = 2 mode changes by only 3e-10 at d = 1e-6. A BigFloat reference therefore needs the constants in BigFloat, for example `kerr_equatorial_EL(big(M), big(a), big(p), big(e))` and `kerr_equatorial_ur` on those values, and the particle at `theta=big(pi)/2`. The Float64 value `π/2` converted to BigFloat lies 6e-17 off the equatorial plane, which in the circular code changes the BigFloat source by a relative 1.7e3 at d = 1e-6. The Float64 error measured against such a reference includes the rounding of the constants.

## Regeneration and tests

```bash
bash ../tools/regenerate.sh
julia --project test/runtests.jl
```

The parity tests run `test/ref_dump.py` against the compiled `effsource_equatorial` module at the repository root. Set `PYTHON` to an interpreter that can import it.
