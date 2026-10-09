using Test
using JSON
using StochasticRounding
using DoubleFloats
using LinearAlgebra: norm
using EffSourceCircular
using EffSourceEquatorial

const HERE = @__DIR__
const PY = get(ENV, "PYTHON", "python3")

halfpi(::Type{T}) where {T} = T(pi)/2
thirdpi(::Type{T}) where {T} = T(pi)/3
tq(::Type{T}, x::Rational) where {T} = T(numerator(x))/T(denominator(x))
tq(::Type{T}, x::Integer) where {T} = T(x)

const ORBIT_A = (a=1//2, p=10, e=3//10)
const ORBIT_B = (a=9//10, p=7, e=1//2)

const P2 = merge(ORBIT_A, (r=9, sign=1))
const P3 = merge(ORBIT_A, (r=100//13, sign=0))
const P4 = merge(ORBIT_A, (r=9, sign=-1))
const P5 = merge(ORBIT_B, (r=7, sign=1))

function particle(::Type{T}, s) where {T}
    M, a, p, e, r = one(T), tq(T, s.a), tq(T, s.p), tq(T, s.e), tq(T, s.r)
    E, L = kerr_equatorial_EL(M, a, p, e)
    ur = s.sign == 0 ? zero(T) : kerr_equatorial_ur(M, a, E, L, r; sign=s.sign)
    return (M=M, a=a, r=r, E=E, L=L, ur=ur)
end

const _REF_CACHE = Dict{Any,Any}()

function c_reference(s; ell=false)
    get!(_REF_CACHE, (s, ell)) do
        v = particle(Float64, s)
        args = [string(x) for x in (v.M, v.a, v.r, v.E, v.L, v.ur)]
        JSON.parse(read(`$PY $(joinpath(HERE, "ref_dump.py")) $args $(ell ? 1 : 0)`, String))
    end
end

function ctx_from(::Type{T}, v; setup=T, ell=false) where {T}
    c = effsource_equatorial_create(T(v.M), T(v.a); ell_from_complement=ell)
    xp = Coordinate{T}(r=T(v.r), theta=halfpi(T), phi=thirdpi(T), t=one(T))
    set_particle!(c, xp, T(v.E), T(v.L), T(v.ur); setup)
    return c
end

mkctx(::Type{T}, s; setup=T, ell=false) where {T} = ctx_from(T, particle(T, s); setup, ell)
