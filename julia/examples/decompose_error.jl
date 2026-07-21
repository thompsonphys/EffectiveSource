# Decompose the near-particle Float64 error of the m-mode puncture and
# effective source into its three suspected contributions:
#
#   1. r - rp coordinate cancellation        (Coordinate vs Offset API)
#   2. elliptic-integral modulus rounding     (k = sqrt(1/(1+C1)) vs
#      complement-fed AGM; ELL_FROM_COMPLEMENT switch)
#   3. coefficient set-up cancellation        (setup=Float64 vs setup=BigFloat)
#
#   julia --project examples/decompose_error.jl
#
# Truth is a 300-bit BigFloat evaluation at the *intended* offset. Coordinate
# configs evaluate at Coordinate(rp + dr, ...) like the C code does today, so
# their error includes the inability of absolute coordinates to even represent
# the intended field point -- that is the realistic situation for the eccentric
# code, where dr = r_field - r_p(t) is known far more accurately than the
# rounded difference of the absolute values.

using EffSourceCircular
using Printf
using DelimitedFiles

const M, a = 1.0, 0.5
const rp, thp, php = 10.0, pi/2, pi/3
const E = sqrt(1 - 2M/rp + a^2/rp^2 - 2M*a^2/rp^3) / sqrt(1 - 3M/rp + 2a*sqrt(M/rp^3))
const L = (rp^2 - 2M*rp + a*sqrt(M*rp)) / (rp*sqrt(rp - 3M + 2a*sqrt(M/rp)))

setprecision(BigFloat, 300)

function mkctx(::Type{T}; setup=T) where {T}
    c = effsource_create(T(M), T(a))
    set_particle!(c, Coordinate{T}(r=T(rp), theta=T(thp), phi=T(php), t=T(1)),
                  T(E), T(L); setup=setup)
    return c
end

const CTX64  = mkctx(Float64)                    # default coefficient set-up
const CTXHI  = mkctx(Float64; setup=BigFloat)    # extended-precision coefficients
const CTXBIG = mkctx(BigFloat)

# (label, ctx, elliptic cure on?, API)
const CONFIGS = [
    ("A coord      ", CTX64, false, :coord),   # today's C behaviour
    ("B offset     ", CTX64, false, :off),     # + offset expansion only
    ("C coord+ellc ", CTX64, true,  :coord),   # + elliptic cure only
    ("D offset+ellc", CTX64, true,  :off),     # offset + elliptic cure
    ("E off+ellc+hi", CTXHI, true,  :off),     # + extended-precision coeffs
    ("F offset+hi  ", CTXHI, false, :off),     # coeffs without elliptic cure
]

const MODES = (2, 10)
const DISTS = 10 .^ range(-12, -1; length=45)
const RAYS = [(1.0, 0.0, "radial"), (cos(pi/4), sin(pi/4), "diag")]

relerr(z, truth) = (d = abs(ComplexF64(z) - ComplexF64(truth));
                    d == 0 ? 0.0 : d / abs(ComplexF64(truth)))

function evaluate(ctx, cure, api, m, dr64, dth64)
    EffSourceCircular.ELL_FROM_COMPLEMENT[] = cure
    local p, s
    if api === :off
        o = Offset{Float64}(dr=dr64, dtheta=dth64, dphi=0.0)
        p = phis_m(ctx, m, o)
        s = calc_m(ctx, m, o)[4]
    else
        x = Coordinate{Float64}(r=rp + dr64, theta=thp + dth64, phi=php, t=0.0)
        p = phis_m(ctx, m, x)
        s = calc_m(ctx, m, x)[4]
    end
    EffSourceCircular.ELL_FROM_COMPLEMENT[] = false
    return p, s
end

rows = Any[]
for (cr, ct, rayname) in RAYS, m in MODES, s in DISTS
    dr64  = Float64(s*cr)
    dth64 = Float64(s*ct)
    ob = Offset{BigFloat}(dr=big(dr64), dtheta=big(dth64), dphi=big(0.0))
    pt = phis_m(CTXBIG, m, ob)
    st = calc_m(CTXBIG, m, ob)[4]
    for (label, ctx, cure, api) in CONFIGS
        p, q = evaluate(ctx, cure, api, m, dr64, dth64)
        push!(rows, (rayname, m, s, strip(label), relerr(p, pt), relerr(q, st)))
    end
end

open(joinpath(@__DIR__, "decompose_error.csv"), "w") do io
    println(io, "ray,m,s,config,relerr_phis_m,relerr_src_m")
    for r in rows
        @printf(io, "%s,%d,%.3e,%s,%.3e,%.3e\n", r...)
    end
end

# --- summary table ----------------------------------------------------------
picks = [1e-3, 1e-6, 1e-9, 1e-12]
for (obs, idx) in (("phis_m", 5), ("src_m", 6))
    println("\n== $obs relative error (radial ray) ==")
    for m in MODES
        println("-- m = $m")
        @printf("%-14s", "config")
        for s in picks; @printf("  s=%-8.0e", s); end
        println()
        for (label, _, _, _) in CONFIGS
            @printf("%-14s", strip(label))
            for starget in picks
                _, i = findmin(abs.(log10.(DISTS) .- log10(starget)))
                row = filter(r -> r[1] == "radial" && r[2] == m &&
                                  r[3] == DISTS[i] && r[4] == strip(label), rows)
                @printf("  %-10.1e", row[1][idx])
            end
            println()
        end
    end
end
