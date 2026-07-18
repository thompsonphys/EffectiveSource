# Map the double-precision error of the puncture and effective source along rays
# through the particle in the (r, theta) plane, against a BigFloat truth.
#
#   julia --project examples/ray_error.jl
#
# For each ray angle psi and distance-to-particle s, the field point offset is
#   dr = s*cos(psi),  dtheta = s*sin(psi),  dphi = 0.
# We evaluate the puncture (phis / phis_m) and source (calc src / calc_m src) in
# Float64 two ways -- via the Offset API (intrinsic evaluation error) and via
# absolute Coordinates (adds the r-rp coordinate cancellation) -- and compare
# both to a 256-bit BigFloat truth. Output: ray_error.csv + terminal plots +
# a breakdown-distance summary.

using EffSourceCircular
using UnicodePlots
using Printf

# --- reference orbit -------------------------------------------------------
const M, a = 1, 1/2
const rp, thp, php = 10.0, π/2, π/3
const E = sqrt(1 - 2M/rp + a^2/rp^2 - 2M*a^2/rp^3) / sqrt(1 - 3M/rp + 2a*sqrt(M/rp^3))
const L = (rp^2 - 2M*rp + a*sqrt(M*rp)) / (rp*sqrt(rp - 3M + 2a*sqrt(M/rp)))

const MODES = (2, 10)                       # m-mode values to map
const NANGLES = 12                          # ray angles over [0, 2pi)
const ANGLES = range(0, 2π; length = NANGLES + 1)[1:NANGLES]
# distance-to-particle, log-spaced. Reaches 1e-12 so the real-space fields (which
# are intrinsically robust via the Offset API) still reveal the r-rp coordinate
# cancellation of the absolute-Coordinate path deep in the near zone.
const DISTS = 10 .^ range(-12, 0; length = 49)

setprecision(BigFloat, 256)

# Float64 pipeline (realistic default coefficient set-up) and BigFloat truth.
function make_ctx(::Type{T}) where {T}
    c = effsource_create(T(M), T(a))
    xp = Coordinate{T}(r = T(rp), theta = T(thp), phi = T(php), t = T(1))
    set_particle!(c, xp, T(E), T(L); setup = (T === BigFloat ? BigFloat : T))
    return c
end
const CTX64 = make_ctx(Float64)
const CTXBIG = make_ctx(BigFloat)

# Observables: name => (ctx, offset) -> scalar-or-complex value. calc/calc_m
# return the source in slot 4.
const FIELDS = Any[
    ("phis",   0, (c, o) -> phis(c, o)),
    ("src",    0, (c, o) -> calc(c, o)[4]),
]
for m in MODES
    push!(FIELDS, ("phis_m", m, (c, o) -> phis_m(c, m, o)))
    push!(FIELDS, ("src_m",  m, (c, o) -> calc_m(c, m, o)[4]))
end

offset(::Type{T}, s, psi) where {T} =
    Offset{T}(dr = T(s * cos(psi)), dtheta = T(s * sin(psi)), dphi = zero(T))
coord(s, psi) =
    Coordinate{Float64}(r = rp + s * cos(psi), theta = thp + s * sin(psi), phi = php, t = 0.0)

relerr(z, truth) = (d = abs(ComplexF64(z) - truth); d == 0 ? 0.0 : d / abs(truth))

# --- sweep -----------------------------------------------------------------
# results[(name, m)] = Dict(psi => (relerr_offset::Vector, relerr_abs::Vector, mag::Vector))
results = Dict{Tuple{String,Int},Dict{Float64,NTuple{3,Vector{Float64}}}}()

println("Sweeping $(length(FIELDS)) fields x $NANGLES angles x $(length(DISTS)) distances ",
        "(BigFloat truth, first run compiles ~2 min)...")

open(joinpath(@__DIR__, "ray_error.csv"), "w") do io
    println(io, "field,m,angle,distance,relerr_offset,relerr_abs,mag")
    for (name, m, f) in FIELDS
        byangle = Dict{Float64,NTuple{3,Vector{Float64}}}()
        for psi in ANGLES
            eo, ea, mg = Float64[], Float64[], Float64[]
            for s in DISTS
                truth = ComplexF64(f(CTXBIG, offset(BigFloat, s, psi)))
                vo = f(CTX64, offset(Float64, s, psi))    # intrinsic (offset API)
                va = f(CTX64, coord(s, psi))              # total (absolute coords)
                ro = isfinite(real(ComplexF64(vo))) ? relerr(vo, truth) : Inf
                ra = isfinite(real(ComplexF64(va))) ? relerr(va, truth) : Inf
                push!(eo, ro); push!(ea, ra); push!(mg, abs(truth))
                @printf(io, "%s,%d,%.6f,%.6e,%.6e,%.6e,%.6e\n",
                        name, m, psi, s, ro, ra, abs(truth))
            end
            byangle[psi] = (eo, ea, mg)
        end
        results[(name, m)] = byangle
    end
end
println("Wrote ", joinpath(@__DIR__, "ray_error.csv"), "\n")

# --- terminal plots: relerr (offset API) vs distance, a few angles ---------
fieldlabel(name, m) = m == 0 ? name : "$name(m=$m)"
# indices into ANGLES (step pi/6): outward-r (0), +theta (pi/2), inward-r (pi),
# -theta (3pi/2). Keying by index avoids float-key mismatch.
const SHOWN_IDX = (1, 4, 7, 10)

for (name, m, _) in FIELDS
    by = results[(name, m)]
    plt = nothing
    for (n, i) in enumerate(SHOWN_IDX)
        psi = ANGLES[i]
        y = clamp.(by[psi][1], 1e-18, 1e2)
        lbl = @sprintf("psi=%.2f", psi)
        if n == 1
            plt = lineplot(DISTS, y; xscale = :log10, yscale = :log10, name = lbl,
                           xlabel = "distance to particle", ylabel = "rel err (offset)",
                           title = "Float64 error: $(fieldlabel(name, m))",
                           width = 62, height = 12)
        else
            lineplot!(plt, DISTS, y; name = lbl)
        end
    end
    println(plt)
    println()
end

# --- breakdown-distance summary --------------------------------------------
# Near-particle breakdown distance: the outer edge of the CONTIGUOUS bad region
# that touches s->0. `errs` is indexed with DISTS ascending (errs[1] = closest to
# the particle, where cancellation is worst). Scanning outward, we find the first
# distance at which Float64 recovers the tolerance; everything closer is bad.
# This deliberately ignores isolated relative-error spikes farther out (high-mode
# field zero-crossings, where |value|->0), which are a different phenomenon.
function breakdown(errs, tol)
    j = findfirst(<=(tol), errs)
    j === nothing && return DISTS[end]   # never within tolerance -> bad everywhere
    j == 1 && return 0.0                 # accurate even at the closest sampled point
    return DISTS[j - 1]                  # largest still-bad distance in the inner region
end
# worst-case over ray directions = the direction that fails farthest out
worst(by, col, tol) = maximum(breakdown(by[psi][col], tol) for psi in ANGLES)

println("Breakdown distance -- Float64 error exceeds the tolerance for all field")
println("points CLOSER than this (worst-case over ray directions; larger = worse):\n")
@printf("%-12s %12s %12s   %12s %12s\n", "field",
        "off>1e-6", "off>1e-3", "abs>1e-6", "abs>1e-3")
for (name, m, _) in FIELDS
    by = results[(name, m)]
    @printf("%-12s %12.2e %12.2e   %12.2e %12.2e\n", fieldlabel(name, m),
            worst(by, 1, 1e-6), worst(by, 1, 1e-3),
            worst(by, 2, 1e-6), worst(by, 2, 1e-3))
end

println("""
\nReading it: 'off' isolates the intrinsic evaluation error (Offset API); 'abs'
adds the r-rp coordinate cancellation (absolute Coordinates). abs >= off, and the
gap is the cost of forming dr = r - rp near the particle -- avoidable with the
Offset API. The m-mode source (src_m) breaks down farthest out (most
cancellation-prone); the real-space fields via the Offset API stay accurate
throughout the sampled range (breakdown 0). All are direction dependent -- see
the per-angle curves above and ray_error.csv. (High modes also lose relative
accuracy at field zero-crossings farther out, where |value|->0 -- a separate
effect not counted here; it is visible in the CSV.)""")
