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
const MODE = 2
const NSAMPLES = 2000
const OFFSETS = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)

setprecision(BigFloat, 300)

"""Context in type T with best-possible coefficients (computed in Tsetup)."""
function make_ctx(::Type{T}; setup=T) where {T}
    c = effsource_create(T(M), T(a))
    W = T === BigFloat ? BigFloat : Float64
    th = T === BigFloat ? T(pi)/2 : T(thp)
    set_particle!(c, Coordinate{T}(r=T(rp), theta=th, phi=T(php), t=T(1)),
                  kerr_circular_EL(W(M), W(a), W(rp))...; setup=setup)
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
  * SR HELPS as a bias-reducer at moderate offsets: at offsets 1e-1 and 1e-2 its
    ensemble mean beats the deterministic Float32 value by one to three orders
    of magnitude, for the puncture and for the source.
  * SR is unbiased only at the larger offsets: the mean error is comparable to
    the 'mean SE' column for the puncture down to offset 1e-2 and for the source
    at 1e-1. Closer in the mean error exceeds the mean SE, for the source by
    orders of magnitude at 1e-3, so averaging more samples does not converge to
    the truth there.
  * It is NOT a cure: it trades bias for VARIANCE, and the samples needed for a
    target accuracy grow like (spread/target)^2. The puncture spread grows by
    three orders of magnitude from offset 1e-1 to 1e-3, and the source spread by
    more than one from 1e-1 to 1e-2. From offset 1e-4 down, Float32 fails
    outright (NaN) from dynamic-range underflow, which no rounding mode can fix.
  * As a DIAGNOSTIC the spread bounds the error only where SR is unbiased; for
    the source at offset 1e-3 it understates the error by an order of magnitude.
    For actual accuracy in the deep near-zone, use the structural fixes
    (Offset API + Double64/BigFloat setup, or an analytic reformulation).""")
