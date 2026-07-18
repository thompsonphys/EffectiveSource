# Demonstrate running the effective-source kernel at different numeric
# precisions, including stochastic rounding (Float32sr / BFloat16sr).
#
#   julia --project examples/stochastic_rounding.jl
#
# The point: deterministic low-precision arithmetic carries a fixed rounding
# bias, whereas stochastic rounding is *unbiased* -- averaging many stochastic
# evaluations converges to the high-precision answer. The whole pipeline is
# generic over T -- coefficient setup, field evaluation, and the complete
# elliptic integrals (computed by the arithmetic-geometric mean) -- so the SR
# type propagates through everything, including the m-mode path.

using EffSourceCircular
using StochasticRounding
using Statistics: mean, std
using Printf

# --- fixed circular orbit (same as the test suite) -------------------------
const M, a = 1, 5e-1
const rp, thp, php = 10, π/2, π/3
const Efd = sqrt(1 - 2M/rp + a^2/rp^2 - 2M*a^2/rp^3) / sqrt(1 - 3M/rp + 2a*sqrt(M/rp^3))
const Lfd = (rp^2 - 2M*rp + a*sqrt(M*rp)) / (rp*sqrt(rp - 3M + 2a*sqrt(M/rp)))

"""Build a context + field point in precision T and return phis(ctx, x)."""
function eval_phis(::Type{T}) where {T}
    ctx = effsource_create(T(M), T(a))
    xp = Coordinate{T}(r=T(rp), theta=T(thp), phi=T(php), t=T(1))
    set_particle!(ctx, xp, T(Efd), T(Lfd))
    x = Coordinate{T}(r=T(rp + 0.1), theta=T(thp + 0.01), phi=T(php + 0.02), t=T(0))
    return phis(ctx, x)
end

# high-precision reference
ref = Float64(eval_phis(BigFloat))
@printf("\nReference (BigFloat)      PhiS = %.10f\n\n", ref)

f64 = eval_phis(Float64)
f32 = eval_phis(Float32)
@printf("Float64                   PhiS = %.10f   (rel err %.2e)\n", f64, abs(f64-ref)/abs(ref))
@printf("Float32  (deterministic)  PhiS = %.10f   (rel err %.2e)\n", f32, abs(Float64(f32)-ref)/abs(ref))

# stochastic rounding: sample the same computation many times
function sr_samples(::Type{T}, n) where {T}
    [Float64(eval_phis(T)) for _ in 1:n]
end

for (name, T) in (("Float32sr", Float32sr), ("BFloat16sr", BFloat16sr))
    s = sr_samples(T, 2000)
    m, sd = mean(s), std(s)
    @printf("%-9s (%5d samples) mean = %.10f   (rel err of mean %.2e, spread %.2e)\n",
            name, length(s), m, abs(m-ref)/abs(ref), sd/abs(m))
end

println("""
\nTakeaway: a single deterministic low-precision value hides its rounding error,
while stochastic rounding turns that error into an observable *distribution* --
the reported spread is a direct estimate of how much this kernel's answer can be
trusted at the given precision. (SR is unbiased per operation; for a long
nonlinear pipeline like this one the mean still carries some residual bias, so
the spread is the more honest diagnostic than the mean alone.)""")

# --- m-mode under SR: the AGM elliptic integrals keep the SR type end-to-end
ctx = effsource_create(Float32sr(M), Float32sr(a))
set_particle!(ctx, Coordinate{Float32sr}(r=Float32sr(rp), theta=Float32sr(thp),
                                          phi=Float32sr(php), t=Float32sr(1)),
              Float32sr(Efd), Float32sr(Lfd))
xm = Coordinate{Float32sr}(r=Float32sr(rp+0.1), theta=Float32sr(thp+0.01),
                           phi=Float32sr(php+0.02), t=Float32sr(0))
@printf("\nm=2 mode under Float32sr:  PhiS_m = %s\n", string(phis_m(ctx, 2, xm)))
