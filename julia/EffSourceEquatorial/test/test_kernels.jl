const MODES = (0, 1, 2, 5, 10, 20)

refval(x) = x === nothing ? NaN : Float64(x)

function relerr(a::Real, b)
    b = refval(b)
    (isnan(a) && isnan(b)) && return 0.0
    (isnan(a) || isnan(b)) && return Inf
    return abs(a - b)/max(abs(a), abs(b), 1e-300)
end

function relerr(a::Complex, b::AbstractVector)
    z = complex(refval(b[1]), refval(b[2]))
    (isnan(a) && isnan(z)) && return 0.0
    (isnan(a) || isnan(z)) && return Inf
    return abs(a - z)/max(abs(a), abs(z), 1e-300)
end

function kernel_worst(ctx, ref)
    worst = Dict{String,Float64}()
    bump!(k, e) = (worst[k] = max(get(worst, k, 0.0), e))
    for rec in ref["records"]
        xc = rec["x"]
        x = Coordinate{Float64}(t=xc[1], r=xc[2], theta=xc[3], phi=xc[4])
        bump!("phis", relerr(phis(ctx, x), rec["PhiS"]))
        P, d1, d2, s = calc(ctx, x)
        c = rec["calc"]
        bump!("calc.PhiS", relerr(P, c["PhiS"]))
        foreach(i -> bump!("calc.d1", relerr(d1[i], c["d1"][i])), 1:4)
        foreach(i -> bump!("calc.d2", relerr(d2[i], c["d2"][i])), 1:10)
        bump!("calc.src", relerr(s, c["src"]))
        for m in MODES
            k = string(m)
            bump!("phis_m", relerr(phis_m(ctx, m, x), rec["phis_m"][k]))
            Pm, d1m, d2m, sm = calc_m(ctx, m, x)
            cm = rec["calc_m"][k]
            bump!("calc_m.PhiS", relerr(Pm, cm["PhiS"]))
            foreach(i -> bump!("calc_m.d1", relerr(d1m[i], cm["d1"][2i-1:2i])), 1:4)
            foreach(i -> bump!("calc_m.d2", relerr(d2m[i], cm["d2"][2i-1:2i])), 1:10)
            bump!("calc_m.src", relerr(sm, cm["src"]))
        end
    end
    return worst
end

const TOL = Dict("phis" => 1e-11, "calc.PhiS" => 1e-11, "calc.d1" => 1e-10, "calc.d2" => 1e-10,
                 "phis_m" => 1e-10, "calc_m.PhiS" => 1e-10, "calc_m.d1" => 1e-9, "calc_m.d2" => 1e-9)

function source_errors(ctx, cbig, ref)
    wj = Dict("calc.src" => 0.0, "calc_m.src" => 0.0)
    wc = copy(wj)
    bump!(w, k, e) = (w[k] = max(w[k], Float64(e)))
    for rec in ref["records"]
        xc = rec["x"]
        x = Coordinate{Float64}(t=xc[1], r=xc[2], theta=xc[3], phi=xc[4])
        xb = Coordinate{BigFloat}(t=big(xc[1]), r=big(xc[2]), theta=big(xc[3]), phi=big(xc[4]))
        b = calc(cbig, xb)[4]
        bump!(wj, "calc.src", abs(calc(ctx, x)[4] - b)/abs(b))
        bump!(wc, "calc.src", abs(refval(rec["calc"]["src"]) - b)/abs(b))
        for m in MODES
            bm = calc_m(cbig, m, xb)[4]
            cm = rec["calc_m"][string(m)]["src"]
            bump!(wj, "calc_m.src", abs(calc_m(ctx, m, x)[4] - bm)/abs(bm))
            bump!(wc, "calc_m.src", abs(complex(refval(cm[1]), refval(cm[2])) - bm)/abs(bm))
        end
    end
    return wj, wc
end

@testset "kernels match the C library ($label)" for (label, P, ell) in (
        ("P2", P2, false),
        ("P2 ell_from_complement", P2, true),
        ("P3 turning point", P3, false),
        ("P4 inbound", P4, false),
        ("P5 a=0.9", P5, false))
    ctx = mkctx(Float64, P; ell)
    ref = c_reference(P; ell)
    worst = kernel_worst(ctx, ref)
    for k in sort(collect(keys(TOL)))
        println(rpad("$label $k", 40), worst[k])
        @test worst[k] < TOL[k]
    end
    wj, wc = setprecision(() -> source_errors(ctx, mkctx(BigFloat, P; ell), ref), BigFloat, 256)
    for k in ("calc.src", "calc_m.src")
        println(rpad("$label $k", 40), "vs BigFloat: Julia ", wj[k], "  C ", wc[k])
        @test wj[k] <= 10*wc[k] + 1e-12
    end
end

@testset "ell_from_complement selects the complementary elliptic path near the particle" begin
    ref = c_reference(P2; ell=true)
    cell = mkctx(Float64, P2; ell=true)
    cdef = mkctx(Float64, P2)
    sep = 0.0
    for rec in ref["near"]
        xc = rec["x"]
        x = Coordinate{Float64}(t=xc[1], r=xc[2], theta=xc[3], phi=xc[4])
        for m in MODES
            k = string(m)
            @test relerr(phis_m(cell, m, x), rec["phis_m"][k]) < 1e-12
            @test relerr(calc_m(cell, m, x)[1], rec["phis_m"][k]) < 1e-12
            sep = max(sep, relerr(phis_m(cdef, m, x), rec["phis_m"][k]))
        end
    end
    println("default path vs C ell_from_complement near the particle: ", sep)
    @test sep > 1e-9
end

@testset "mode range and coincident field point" begin
    ctx = mkctx(Float64, P2)
    x = Coordinate{Float64}(r=P2.r + 0.1, theta=pi/2 + 0.01, phi=pi/3, t=0.0)
    @test_throws ErrorException phis_m(ctx, 21, x)
    @test_throws ErrorException phis_m(ctx, -1, x)
    @test_throws ErrorException calc_m(ctx, 21, x)
    @test_throws ErrorException calc_m(ctx, -1, x)
    @test !isfinite(phis(ctx, Offset{Float64}()))
    @test !isfinite(calc(ctx, Offset{Float64}())[4])
    @test !isfinite(calc_m(ctx, 2, Offset{Float64}())[1])
end

@testset "Offset and Coordinate agree away from the particle" begin
    ctx = mkctx(Float64, P2)
    x = Coordinate{Float64}(r=P2.r + 0.3, theta=pi/2 + 0.02, phi=pi/3 + 0.1, t=0.0)
    o = Offset(x, ctx.xp)
    @test phis(ctx, o) ≈ phis(ctx, x) rtol=1e-13
    @test calc(ctx, o)[2] ≈ calc(ctx, x)[2] rtol=1e-12
    @test phis_m(ctx, 3, o) ≈ phis_m(ctx, 3, x) rtol=1e-12
    @test calc_m(ctx, 3, o)[2] ≈ calc_m(ctx, 3, x)[2] rtol=1e-11
end
