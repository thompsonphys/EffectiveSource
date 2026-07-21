"""
    EffSourceCircular

Julia port of Barry Wardell's `kerr-circular-ctx.c` -- the scalar singular field
and effective source for a point particle on a circular equatorial geodesic in
Kerr spacetime.

The kernels are generic over the numeric type `T`, so the same physics can be
run in `Float64`, `Float32`, `Float32sr` (stochastic rounding), `BigFloat`, etc.

Public API:
  * `effsource_create(M, a)`                 -> context
  * `set_particle!(ctx, x_p, E, L, [ur])`    -> fill coefficients
  * `phis(ctx, x)`                           -> real-space singular field
  * `calc(ctx, x)`                           -> (PhiS, dPhiS[4], d2PhiS[10], src)
  * `phis_m(ctx, m, x)`                      -> m-mode field (Complex)
  * `calc_m(ctx, m, x)`                      -> m-mode (PhiS, dPhiS[4], d2PhiS[10], src)

The generated `_*.jl` includes are produced from the C source by the scripts in
`tools/` (see `tools/regenerate.sh`).
"""
module EffSourceCircular

using OffsetArrays
using StaticArrays

export Coordinate, Offset, EffsourceCtx, effsource_create, set_particle!,
       phis, calc, phis_m, calc_m

# --- helpers used by the generated kernels ---------------------------------

"""0-based uninitialised working vector (matches the C fixed-size arrays)."""
@inline _ovec(::Type{T}, n::Integer) where {T} =
    OffsetArray(Vector{T}(undef, n), 0:n-1)

"""0-based vector wrapping an existing element list."""
@inline _ovec0(v::AbstractVector) = OffsetArray(v, 0:length(v)-1)

# ---------------------------------------------------------------------------
# Opt-in compensated (Neumaier) summation for the heavy signed sums (the m-mode
# numerator loops and the effective-source numerator). Enabled per call by the
# `comp=true` keyword on calc/phis_m/calc_m, threaded as a `Val` so the hot
# loops stay type-stable.
#
# Compensation applies ONLY when the mode is Val{true} AND the type is an IEEE
# float: Kahan/Neumaier's error-free transform assumes round-to-nearest, which
# stochastic rounding breaks, and BigFloat does not need it. Every other case
# (mode off, or SR/BigFloat) is an ordinary sum, so the default path and the SR
# experiments are byte-unchanged.

"""Streaming Neumaier step: fold `x` into running sum `s` with compensation `c`."""
@inline function _cadd(::Val{true}, s::T, c::T, x::T) where {T<:Base.IEEEFloat}
    t = s + x
    c += ifelse(abs(s) >= abs(x), (s - t) + x, (x - t) + s)
    return t, c
end
@inline _cadd(::Val, s::T, c::T, x::T) where {T} = (s + x, c)   # off, or SR/BigFloat

"""Compensated sum of a fixed term list (used for the source numerator)."""
@inline function _csum(::Val{true}, xs::Vararg{T,N}) where {T<:Base.IEEEFloat,N}
    s = zero(T); c = zero(T)
    @inbounds for x in xs
        s, c = _cadd(Val(true), s, c, x)
    end
    return s + c
end
@inline _csum(::Val, xs::Vararg{T,N}) where {T,N} = +(xs...)    # off, or SR/BigFloat

# Complete elliptic integrals K(m), E(m) with the GSL *modulus* convention: the
# argument is the modulus k and the integrals are evaluated at parameter m = k^2.
# Computed by the arithmetic-geometric mean -- pure +, *, sqrt, so the working
# type T is preserved exactly (Float64, BigFloat, Float32sr, ...), unlike a
# Float64-only special-function library. Bit-identical to
# SpecialFunctions.ellipk/ellipe in Float64, and unbiased-SR-friendly because
# every operation is a primitive the rounding type can intercept.
#
# Uses the standard AGM recurrence a_{n+1}=(a_n+b_n)/2, b_{n+1}=sqrt(a_n b_n),
# c_{n+1}=(a_n-b_n)/2 with K = pi/(2 a_inf) and
# E = K (1 - sum_{n>=0} 2^{n-1} c_n^2). Converges quadratically.
@inline function _agm_KE(k::T) where {T<:Real}
    m = k*k
    a = one(T)
    b = sqrt(one(T) - m)
    c = k                       # c_0 = sqrt(m) = |k|  (k >= 0 in this code)
    S = c*c/2                   # n = 0 term of the sum
    w = one(T)                  # 2^{n-1}, starting at n = 1
    tol = eps(float(real(one(T))))
    for _ in 1:200
        a, b, c = (a + b)/2, sqrt(a*b), (a - b)/2
        S += w*c*c
        w += w                  # *= 2
        abs(c) <= tol*abs(a) && break
    end
    K = T(pi)/(2a)
    return K, K*(one(T) - S)
end

@inline _ellK(k::Real) = _agm_KE(k)[1]
@inline _ellE(k::Real) = _agm_KE(k)[2]

# Complement-fed variant: K, E of the modulus k with k^2 = 1 - mp, taking the
# complementary parameter mp = 1 - k^2 directly. Near the particle the kernels
# have mp = C1/(1+C1) -> 0 and k -> 1; forming k first rounds away the
# complement (1 - k^2 is only known to ~eps/mp relative), which contaminates
# K, E. Feeding mp keeps b_0 = sqrt(mp) and c_0^2 = 1 - mp fully accurate.
@inline function _agm_KE_c(mp::T) where {T<:Real}
    a = one(T)
    b = sqrt(mp)
    S = (one(T) - mp)/2         # c_0^2/2 without ever forming k
    w = one(T)
    tol = eps(float(real(one(T))))
    for _ in 1:200
        a, b, c = (a + b)/2, sqrt(a*b), (a - b)/2
        S += w*c*c
        w += w
        abs(c) <= tol*abs(a) && break
    end
    K = T(pi)/(2a)
    return K, K*(one(T) - S)
end

# Experiment switch: when true, the m-mode kernels evaluate the elliptic
# integrals from the complementary parameter (cured path); when false they
# reproduce the original k = sqrt(1/(1+C1)) path bit-for-bit.
const ELL_FROM_COMPLEMENT = Ref(false)

@inline function _ellKE_arg(C1::T) where {T<:Real}
    if ELL_FROM_COMPLEMENT[]
        return _agm_KE_c(C1/(one(T) + C1))
    else
        k = sqrt(one(T)/(one(T) + C1))
        return _agm_KE(k)
    end
end

# --- generated code --------------------------------------------------------

include("reei_data.jl")        # REEI_FLAT + reei(T, m,i,j,k)
include("_ctx_struct.jl")      # Coordinate, EffsourceCtx, effsource_create

"""Field-point offset from the particle: `(dr, dtheta, dphi)`.

Passing offsets directly avoids the catastrophic `r - rp` cancellation that
occurs when the field point is close to the particle -- prefer this over
`Coordinate` in the near-zone. Build it from the offsets you already know, e.g.
`Offset{Float64}(dr=1e-6, dtheta=0.0, dphi=0.0)`."""
struct Offset{T}
    dr::T
    dtheta::T
    dphi::T
end
Offset{T}(; dr=zero(T), dtheta=zero(T), dphi=zero(T)) where {T} =
    Offset{T}(T(dr), T(dtheta), T(dphi))
# convenience: derive from absolute coordinates (no accuracy gain -- the
# subtraction still happens here, but handy when you only have x and xp)
Offset(x::Coordinate{T}, xp::Coordinate{T}) where {T} =
    Offset{T}(x.r - xp.r, x.theta - xp.theta, x.phi - xp.phi)

include("_set_particle.jl")    # set_particle!
include("_generated_kernels.jl")  # phis, calc, phis_m, calc_m (Coordinate & Offset)

end # module
