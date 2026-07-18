using Test
using JSON
using StochasticRounding
using EffSourceCircular

# ---------------------------------------------------------------------------
# Live parity check against the compiled Python `.so` reference. We run
# test/ref_dump.py (which drives effsource_circular) and compare the Float64
# Julia port point-by-point. No stored golden file: this exercises the real C
# code each run.
# ---------------------------------------------------------------------------

const HERE = @__DIR__

function load_reference()
    py = get(ENV, "PYTHON", "python3")
    script = joinpath(HERE, "ref_dump.py")
    out = IOBuffer()
    err = IOBuffer()
    try
        run(pipeline(`$py $script`; stdout=out, stderr=err))
    catch e
        @error "ref_dump.py failed" exception=e stderr=String(take!(err))
        rethrow()
    end
    return JSON.parse(String(take!(out)))
end

# non-finite reference values are serialised as JSON null -> `nothing`
refval(x) = x === nothing ? NaN : Float64(x)

"""Relative error with sensible handling of NaN and near-zero magnitudes."""
function relerr(a, b)
    b = refval(b)
    (isnan(a) && isnan(b)) && return 0.0
    (isnan(a) || isnan(b)) && return Inf
    return abs(a - b) / max(abs(a), abs(b), 1e-300)
end

# flatten a complex vector to [re, im, re, im, ...] matching the C output layout
flatci(v) = collect(Iterators.flatten((real(z), imag(z)) for z in v))

ref = load_reference()
p = ref["params"]
M, a = p["M"], p["a"]
E, L = p["E"], p["L"]

ctx = effsource_create(M, a)
xp = Coordinate{Float64}(r=p["rp"], theta=p["thp"], phi=p["php"], t=p["tp"])
set_particle!(ctx, xp, E, L, 0.0)

# accumulate worst-case relative error per channel
worst = Dict{String,Float64}()
bump!(name, e) = (worst[name] = max(get(worst, name, 0.0), e))

@testset "EffSourceCircular parity vs Python .so" begin
    for rec in ref["records"]
        xc = rec["x"]
        x = Coordinate{Float64}(t=xc[1], r=xc[2], theta=xc[3], phi=xc[4])

        # --- real-space PhiS -------------------------------------------------
        bump!("phis", relerr(phis(ctx, x), rec["PhiS"]))

        # --- real-space calc -------------------------------------------------
        P, d1, d2, s = calc(ctx, x)
        c = rec["calc"]
        bump!("calc.PhiS", relerr(P, c["PhiS"]))
        for i in 1:4
            bump!("calc.d1", relerr(d1[i], c["d1"][i]))
        end
        for i in 1:10
            bump!("calc.d2", relerr(d2[i], c["d2"][i]))
        end
        bump!("calc.src", relerr(s, c["src"]))

        # --- m-mode ----------------------------------------------------------
        for m in ref["modes"]
            key = string(m)
            # phis_m
            pm = phis_m(ctx, m, x)
            rpm = rec["phis_m"][key]
            bump!("phis_m", relerr(real(pm), rpm[1]))
            bump!("phis_m", relerr(imag(pm), rpm[2]))
            # calc_m
            Pm, d1m, d2m, sm = calc_m(ctx, m, x)
            cm = rec["calc_m"][key]
            for (jl, py) in ((real(Pm), cm["PhiS"][1]), (imag(Pm), cm["PhiS"][2]))
                bump!("calc_m.PhiS", relerr(jl, py))
            end
            for (jl, py) in zip(flatci(d1m), cm["d1"])
                bump!("calc_m.d1", relerr(jl, py))
            end
            for (jl, py) in zip(flatci(d2m), cm["d2"])
                bump!("calc_m.d2", relerr(jl, py))
            end
            for (jl, py) in ((real(sm), cm["src"][1]), (imag(sm), cm["src"][2]))
                bump!("calc_m.src", relerr(jl, py))
            end
        end
    end

    println("\nWorst-case relative error per channel:")
    for k in sort(collect(keys(worst)))
        println("  ", rpad(k, 14), " ", worst[k])
    end
    println()

    # Fields and derivatives should match to near machine precision. The
    # effective source `src` is a heavily-cancelling d'Alembertian, so its
    # ordering-sensitive error floor is looser (still far below any real bug,
    # which would show up as an O(1) discrepancy).
    @test worst["phis"]        < 1e-11
    @test worst["calc.PhiS"]   < 1e-11
    @test worst["calc.d1"]     < 1e-10
    @test worst["calc.d2"]     < 1e-10
    @test worst["calc.src"]    < 1e-7
    @test worst["phis_m"]      < 1e-10
    @test worst["calc_m.PhiS"] < 1e-10
    @test worst["calc_m.d1"]   < 1e-9
    @test worst["calc_m.d2"]   < 1e-9
    @test worst["calc_m.src"]  < 1e-6
end

@testset "numeric type genericity (no Float64 contamination)" begin
    # If any Float64 literal leaked into the kernels, a Float32 input would be
    # promoted to Float64. Assert the working type is preserved end-to-end.
    ctxf = effsource_create(Float32(M), Float32(a))
    xpf = Coordinate{Float32}(r=Float32(p["rp"]), theta=Float32(p["thp"]),
                              phi=Float32(p["php"]), t=1f0)
    set_particle!(ctxf, xpf, Float32(E), Float32(L))
    xf = Coordinate{Float32}(r=Float32(p["rp"]+0.1), theta=Float32(p["thp"]+0.01),
                             phi=Float32(p["php"]+0.02), t=0f0)

    @test ctxf isa EffsourceCtx{Float32}
    @test phis(ctxf, xf) isa Float32
    Pf, d1f, d2f, sf = calc(ctxf, xf)
    @test Pf isa Float32
    @test eltype(d1f) === Float32
    @test eltype(d2f) === Float32
    @test sf isa Float32
    @test phis_m(ctxf, 2, xf) isa Complex{Float32}
    Pmf, d1mf, d2mf, smf = calc_m(ctxf, 2, xf)
    @test Pmf isa Complex{Float32}
    @test eltype(d1mf) === Complex{Float32}
    @test smf isa Complex{Float32}
end

@testset "AGM elliptic integrals: arbitrary precision + SR propagation" begin
    using EffSourceCircular: _ellK, _ellE

    # Legendre's relation is exact for any modulus k (k' = sqrt(1-k^2)):
    #   E(k) K(k') + E(k') K(k) - K(k) K(k') = pi/2
    # It is independent of the algorithm, so it validates the AGM at any
    # precision. Here _ellK(k)=K(k^2), _ellE(k)=E(k^2) (modulus convention).
    function legendre_residual(::Type{T}, k) where {T}
        k  = T(k)
        kp = sqrt(one(T) - k*k)
        return _ellE(k)*_ellK(kp) + _ellE(kp)*_ellK(k) - _ellK(k)*_ellK(kp) - T(pi)/2
    end

    @test abs(legendre_residual(Float64, 0.6)) < 1e-14

    # arbitrary precision: the residual must shrink to ~BigFloat epsilon, which
    # is far below anything a Float64-limited library could reach.
    setprecision(BigFloat, 256) do
        r = abs(legendre_residual(BigFloat, big"0.6"))
        @test r < 1e-70
    end

    # the m-mode kernel now runs entirely in BigFloat (previously impossible --
    # SpecialFunctions has no BigFloat method) and agrees with Float64.
    Mb, ab = big(M), big(a)
    ctxb = effsource_create(Mb, ab)
    xpb = Coordinate{BigFloat}(r=big(p["rp"]), theta=big(p["thp"]),
                               phi=big(p["php"]), t=big(1))
    set_particle!(ctxb, xpb, big(E), big(L))
    xb = Coordinate{BigFloat}(r=big(p["rp"])+big"0.1", theta=big(p["thp"])+big"0.01",
                              phi=big(p["php"])+big"0.02", t=big(0))
    pm_big = phis_m(ctxb, 2, xb)
    pm_f64 = phis_m(ctx, 2, Coordinate{Float64}(r=p["rp"]+0.1, theta=p["thp"]+0.01,
                                                 phi=p["php"]+0.02, t=0.0))
    @test pm_big isa Complex{BigFloat}
    @test abs(Complex{Float64}(pm_big) - pm_f64) / abs(pm_f64) < 1e-12

    # stochastic rounding now propagates through the elliptic path too:
    # the m-mode result is Float32sr, not a Float64/Float32 fallback.
    T = Float32sr
    ctxs = effsource_create(T(M), T(a))
    set_particle!(ctxs, Coordinate{T}(r=T(p["rp"]), theta=T(p["thp"]),
                                      phi=T(p["php"]), t=T(1)), T(E), T(L))
    xs = Coordinate{T}(r=T(p["rp"]+0.1), theta=T(p["thp"]+0.01),
                       phi=T(p["php"]+0.02), t=T(0))
    @test phis_m(ctxs, 2, xs) isa Complex{Float32sr}
    Pms, d1ms, _, sms = calc_m(ctxs, 2, xs)
    @test Pms isa Complex{Float32sr}
    @test eltype(d1ms) === Complex{Float32sr}
    @test sms isa Complex{Float32sr}
end

@testset "offset API: kills the r-rp cancellation" begin
    setprecision(BigFloat, 300)
    ctxb = effsource_create(big(M), big(a))
    set_particle!(ctxb, Coordinate{BigFloat}(r=big(p["rp"]), theta=big(p["thp"]),
                                             phi=big(p["php"]), t=big(1)), big(E), big(L))

    # near-particle: absolute coords lose digits in r - rp; offsets do not.
    d = 1e-6
    _, _, _, sref = calc(ctxb, Offset{BigFloat}(dr=big(d), dtheta=big(d), dphi=big(d)))
    ref = Float64(sref)
    _, _, _, s_abs = calc(ctx, Coordinate{Float64}(r=p["rp"]+d, theta=p["thp"]+d,
                                                   phi=p["php"]+d, t=0.0))
    _, _, _, s_off = calc(ctx, Offset{Float64}(dr=d, dtheta=d, dphi=d))
    e_abs = abs(s_abs - ref) / abs(ref)
    e_off = abs(s_off - ref) / abs(ref)
    println("\nnear-particle src relerr: absolute=$(e_abs), offset=$(e_off)")
    @test e_off < 1e-12          # offset API stays near machine precision
    @test e_off < e_abs / 50     # and is orders of magnitude better than absolute

    # consistency: at a moderate offset the two APIs agree
    xa = Coordinate{Float64}(r=p["rp"]+0.1, theta=p["thp"]+0.01, phi=p["php"]+0.02, t=0.0)
    oa = Offset{Float64}(dr=0.1, dtheta=0.01, dphi=0.02)
    @test phis(ctx, xa) ≈ phis(ctx, oa) rtol=1e-13
    @test phis_m(ctx, 3, xa) ≈ phis_m(ctx, 3, oa) rtol=1e-12
    # convenience constructor from absolute coordinates
    xp64 = Coordinate{Float64}(r=p["rp"], theta=p["thp"], phi=p["php"], t=1.0)
    o = Offset(xa, xp64)
    @test (o.dr, o.dtheta, o.dphi) == (xa.r - xp64.r, xa.theta - xp64.theta, xa.phi - xp64.phi)
end

@testset "extended-precision coefficient setup" begin
    setprecision(BigFloat, 400)
    cref = effsource_create(big(M), big(a))
    set_particle!(cref, Coordinate{BigFloat}(r=big(p["rp"]), theta=big(p["thp"]),
                                             phi=big(p["php"]), t=big(1)), big(E), big(L))

    xp64 = Coordinate{Float64}(r=p["rp"], theta=p["thp"], phi=p["php"], t=1.0)
    clo = effsource_create(M, a); set_particle!(clo, xp64, E, L)                 # default (Float64)
    chi = effsource_create(M, a); set_particle!(chi, xp64, E, L; setup=BigFloat) # extended

    names = fieldnames(EffsourceCtx)[4:end]   # skip M, a, xp -> the 53 coefficients
    n_hi = count(f -> getfield(chi, f) == Float64(getfield(cref, f)), names)
    n_lo = count(f -> getfield(clo, f) == Float64(getfield(cref, f)), names)
    println("correctly-rounded coefficients: setup=BigFloat $n_hi/$(length(names)), " *
            "default $n_lo/$(length(names))")

    @test n_hi == length(names)      # extended setup -> every coefficient correctly rounded
    @test n_lo < length(names)       # default Float64 setup misses some at the last ULP
    @test eltype(chi) === Float64    # still stored in the context type

    # default setup=Tc is byte-identical to the plain computation (backward compat)
    @test all(getfield(clo, f) === getfield(ctx, f) for f in names)
end

@testset "offset + extended setup preserve the SR type" begin
    T = Float32sr
    c = effsource_create(T(M), T(a))
    set_particle!(c, Coordinate{T}(r=T(p["rp"]), theta=T(p["thp"]),
                                   phi=T(p["php"]), t=T(1)), T(E), T(L); setup=Float64)
    @test eltype(c) === Float32sr            # extended setup keeps the context type
    off = Offset{T}(dr=T(1e-3), dtheta=T(1e-3), dphi=T(1e-3))
    @test phis(c, off) isa Float32sr
    _, _, _, s = calc(c, off)
    @test s isa Float32sr
    @test phis_m(c, 2, off) isa Complex{Float32sr}
end

@testset "opt-in compensated summation (option c)" begin
    using EffSourceCircular: _cadd, _csum

    # default (comp=false) is byte-identical to the plain computation
    x = Coordinate{Float64}(r=p["rp"]+0.1, theta=p["thp"]+0.01, phi=p["php"]+0.02, t=0.0)
    @test calc(ctx, x)[4] === calc(ctx, x; comp=false)[4]
    @test calc_m(ctx, 2, x)[4] === calc_m(ctx, 2, x; comp=false)[4]

    # helpers: Neumaier compensates for IEEE floats, bypasses for SR/BigFloat.
    terms = Float32[(-1f0)^k * (1f0 + (k%7)) * 10f0^(6 - (k%13)) for k in 0:269]
    truth = sum(BigFloat.(terms))
    csum(cm, v) = (s=zero(eltype(v)); c=zero(eltype(v)); for y in v; s,c=_cadd(cm,s,c,y); end; s+c)
    e_off = abs(csum(Val(false), terms) - truth) / abs(truth)
    e_on  = abs(csum(Val(true),  terms) - truth) / abs(truth)
    @test e_on < e_off / 10                    # ~50x on this synthetic sum
    # dispatch proof: the compensation term stays 0 for SR (plain-sum fallback),
    # but captures the lost bits for an IEEE float.
    @test _cadd(Val(true), Float32sr(1f8), zero(Float32sr), Float32sr(1f0))[2] == 0
    @test _cadd(Val(true), 1e16, 0.0, 1.0)[2] != 0

    # measured end-to-end benefit: the source numerator in Float32 at a moderate,
    # summation-limited offset. (The gain is regime-dependent -- negligible for
    # the puncture, and it does not help the input-cancellation-dominated deep
    # near-zone -- so we assert only this robust case and print the numbers.)
    setprecision(BigFloat, 300)
    c32 = effsource_create(Float32(M), Float32(a))
    set_particle!(c32, Coordinate{Float32}(r=Float32(p["rp"]), theta=Float32(p["thp"]),
                  phi=Float32(p["php"]), t=1f0), Float32(E), Float32(L); setup=Float64)
    cbg = effsource_create(big(M), big(a))
    set_particle!(cbg, Coordinate{BigFloat}(r=big(p["rp"]), theta=big(p["thp"]),
                  phi=big(p["php"]), t=big(1)), big(E), big(L))
    o32 = Offset{Float32}(dr=0.1f0, dtheta=0.1f0, dphi=0f0)
    obg = Offset{BigFloat}(dr=big(0.1), dtheta=big(0.1), dphi=zero(BigFloat))
    refs = ComplexF64(calc_m(cbg, 2, obg)[4])
    e_naive = abs(ComplexF64(calc_m(c32, 2, o32)[4]) - refs) / abs(refs)
    e_comp  = abs(ComplexF64(calc_m(c32, 2, o32; comp=true)[4]) - refs) / abs(refs)
    println("\ncompensated source (Float32, offset 0.1): naive=$e_naive comp=$e_comp")
    @test e_comp < e_naive / 2                 # measurable improvement in this regime

    # SR path is untouched: comp=true still returns the SR type and dispatches to
    # the plain-sum fallback (no error-free transform under stochastic rounding).
    T = Float32sr
    csr = effsource_create(T(M), T(a))
    set_particle!(csr, Coordinate{T}(r=T(p["rp"]), theta=T(p["thp"]),
                  phi=T(p["php"]), t=T(1)), T(E), T(L))
    osr = Offset{T}(dr=T(0.05), dtheta=T(0.01), dphi=zero(T))
    @test phis_m(csr, 2, osr; comp=true) isa Complex{Float32sr}
    @test calc_m(csr, 2, osr; comp=true)[4] isa Complex{Float32sr}
end
