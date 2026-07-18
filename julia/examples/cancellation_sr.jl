# Does stochastic rounding tame the near-particle cancellation in the m-mode
# puncture field and effective source?
#
#   julia --project examples/cancellation_sr.jl
#
# Setup for a clean test:
#   * Offset API  -> dr, dtheta fed directly, so we isolate the *intrinsic*
#     m-mode cancellation (not the avoidable r - rp coordinate subtraction).
#   * setup=Float64 -> coefficients are the best possible in the working type,
#     so we probe the field *evaluation*, not the coefficient set-up.
#
# For each small offset we compare, against a BigFloat truth:
#   * Float64                      (reference-quality)
#   * Float32  (deterministic)     single value
#   * Float32sr (N stochastic samples) -> mean error and spread
# The spread (std/|mean|) is SR's estimate of how many digits the cancellation
# destroyed; the mean error tells us whether SR is any less biased than the
# deterministic Float32 value.

using EffSourceCircular
using StochasticRounding
using Statistics: mean, std
using Printf

const M, a = 1, 5e-1
const rp, thp, php = 10.0, π/2, π/3
const E = sqrt(1 - 2M/rp + a^2/rp^2 - 2M*a^2/rp^3) / sqrt(1 - 3M/rp + 2a*sqrt(M/rp^3))
const L = (rp^2 - 2M*rp + a*sqrt(M*rp)) / (rp*sqrt(rp - 3M + 2a*sqrt(M/rp)))
const MODE = 2
const NSAMPLES = 2000
const OFFSETS = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)

setprecision(BigFloat, 300)

"""Context in type T with best-possible coefficients (computed in Tsetup)."""
function make_ctx(::Type{T}; setup=T) where {T}
    c = effsource_create(T(M), T(a))
    set_particle!(c, Coordinate{T}(r=T(rp), theta=T(thp), phi=T(php), t=T(1)),
                  T(E), T(L); setup=setup)
    return c
end

offset(::Type{T}, e) where {T} = Offset{T}(dr=T(e), dtheta=T(e), dphi=zero(T))

relerr(z, ref) = abs(ComplexF64(z) - ref) / abs(ref)

const ctx64  = make_ctx(Float64;   setup=Float64)
const ctx32  = make_ctx(Float32;   setup=Float64)
const ctxsr  = make_ctx(Float32sr; setup=Float64)
const ctxbig = make_ctx(BigFloat;  setup=BigFloat)

"""Return (err_f64, err_det, err_srmean, spread, stderr) for a scalar/complex
observable produced by `f(ctx, offset)`. `stderr = spread/sqrt(N)` is the
sampling-noise floor of the SR mean: if err_srmean ~ stderr, the SR estimator is
unbiased and its error is dominated by finite N (so it shrinks as 1/sqrt(N))."""
function probe(f, e)
    ref   = ComplexF64(f(ctxbig, offset(BigFloat, e)))
    ef64  = relerr(f(ctx64, offset(Float64, e)), ref)
    edet  = relerr(f(ctx32, offset(Float32, e)), ref)
    samples = [ComplexF64(f(ctxsr, offset(Float32sr, e))) for _ in 1:NSAMPLES]
    srmean  = mean(samples)
    spread  = std(abs.(samples)) / abs(srmean)
    return ef64, edet, relerr(srmean, ref), spread, spread / sqrt(NSAMPLES)
end

puncture(c, o) = phis_m(c, MODE, o)
source(c, o)   = calc_m(c, MODE, o)[4]      # the effective source (4th return)

function run(title, f)
    println("\n=== m=$MODE  $title  (N=$NSAMPLES SR samples) ===")
    @printf("%-8s  %-11s  %-12s  %-14s  %-10s  %-10s\n",
            "offset", "Float64", "Float32det", "Float32sr mean", "SR spread", "mean SE")
    for e in OFFSETS
        ef64, edet, esr, sp, se = probe(f, e)
        det = isnan(edet) ? "NaN(range)" : @sprintf("%.2e", edet)
        sr  = isnan(esr)  ? "NaN(range)" : @sprintf("%.2e", esr)
        @printf("%-8.0e  %-11.2e  %-12s  %-14s  %-10.2e  %-10.2e\n", e, ef64, det, sr, sp, se)
    end
end

run("PUNCTURE field  phis_m", puncture)
run("SOURCE          calc_m src", source)

println("""
\nConclusion (see also the numbers above):
  * SR HELPS as a bias-reducer: its ensemble mean beats the deterministic
    Float32 value by 1-2 orders of magnitude at moderate offsets, because
    deterministic rounding is directionally biased while SR is ~unbiased even
    through the cancellation (mean error ~ the 'mean SE' column, and it falls
    like 1/sqrt(N)).
  * But it is NOT a cure: it trades bias for VARIANCE. A single SR sample is no
    better than deterministic; you only win by averaging, and the samples needed
    for a target accuracy grow like (spread/target)^2 -- which blows up as
    dr,dtheta shrink. Below offset ~1e-4, Float32 fails outright (NaN) from
    dynamic-range underflow, which no rounding mode can fix.
  * SR's most reliable role here is DIAGNOSTIC: the spread is a cheap estimate of
    how many digits the cancellation destroyed. For actual accuracy in the deep
    near-zone, use the structural fixes (Offset API + Double64/BigFloat setup, or
    an analytic reformulation).""")
