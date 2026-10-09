const DT_PAIRS = vcat(
    [(f, Symbol("dAdt", String(f)[2:end]), Symbol("d2Adt2", String(f)[2:end]))
     for f in EffSourceEquatorial.COEFF_FIELDS if occursin(r"^A\d{4}$", String(f))],
    [(:alpha20, :dalphadt20, :d2alphadt220), (:alpha02, :dalphadt02, :d2alphadt202),
     (:beta, :dbetadt, :d2betadt2), (:c, :dcdt, :d2cdt2),
     (:rt, :rtt, nothing), (:urt, :urtt, nothing), (:phit, :phitt, nothing)])

@testset "coefficient fields" begin
    @test length(EffSourceEquatorial.COEFF_FIELDS) == 280
    @test fieldcount(EquatorialCtx{Float64}) == 284
    @test length(DT_PAIRS) == 92
end

@testset "coefficients match the C library ($label)" for (label, P) in
        (("P2", P2), ("P3 turning point", P3), ("P4 inbound", P4), ("P5 a=0.9", P5))
    ref = c_reference(P)["coeffs"]
    c64 = mkctx(Float64, P)
    cbig = mkctx(BigFloat, P)
    bad = String[]
    for f in EffSourceEquatorial.COEFF_FIELDS
        cv = Float64(ref[String(f)])
        for (tag, jv) in (("Float64", getfield(c64, f)), ("BigFloat", Float64(getfield(cbig, f))))
            d = abs(jv - cv)
            (d == 0 || d <= 1e-10*max(abs(cv), abs(jv))) || push!(bad, "$f $tag jl=$jv c=$cv")
        end
    end
    foreach(println, bad)
    @test isempty(bad)
    @test ref["M"] == c64.M && ref["a"] == c64.a
end

@testset "time-derivative coefficients equal finite differences along the orbit (sign $sgn)" for sgn in (1, -1)
    setprecision(BigFloat, 512) do
        s = merge(ORBIT_A, (r=9, sign=sgn))
        P = particle(BigFloat, s)
        c0 = mkctx(BigFloat, s)
        h = big"1e-30"
        function at(t)
            r = P.r + c0.rt*t + c0.rtt*t^2/2
            ur = kerr_equatorial_ur(P.M, P.a, P.E, P.L, r; sign=sgn)
            return EffSourceEquatorial._Coefficients.coefficients(P.M, P.a, r, P.E, P.L, ur)
        end
        vp, v0, vm = at(h), at(zero(h)), at(-h)
        bad = String[]
        for (f, df, d2f) in DT_PAIRS
            d1 = (getfield(vp, f) - getfield(vm, f))/(2h)
            abs(d1 - getfield(c0, df)) <= big"1e-40"*(abs(getfield(c0, df)) + abs(getfield(v0, f))) ||
                push!(bad, "$df")
            d2f === nothing && continue
            d2 = (getfield(vp, f) - 2getfield(v0, f) + getfield(vm, f))/h^2
            abs(d2 - getfield(c0, d2f)) <= big"1e-40"*(abs(getfield(c0, d2f)) + abs(getfield(v0, f))) ||
                push!(bad, "$d2f")
        end
        for (dr0, dth) in ((big"0.3", big"0.1"), (big"-0.2", big"0.05"))
            C1(v, t) = (v.alpha20*(dr0 - c0.rt*t - c0.rtt*t^2/2)^2 + v.alpha02*dth^2)/v.beta
            d1 = (C1(vp, h) - C1(vm, -h))/(2h)
            d2 = (C1(vp, h) - 2C1(v0, zero(h)) + C1(vm, -h))/h^2
            p1 = c0.dC1_dt02*dth^2 + c0.dC1_dt10*dr0 + c0.dC1_dt20*dr0^2
            p2 = c0.d2C1_dt200 + c0.d2C1_dt202*dth^2 + c0.d2C1_dt210*dr0 + c0.d2C1_dt220*dr0^2
            abs(d1 - p1) <= big"1e-40"*abs(d1) || push!(bad, "dC1_dt at $dr0")
            abs(d2 - p2) <= big"1e-40"*abs(d2) || push!(bad, "d2C1_dt2 at $dr0")
        end
        foreach(println, bad)
        @test isempty(bad)
    end
end

@testset "extended-precision setup rounds every coefficient correctly" begin
    hi = mkctx(Float64, P2; setup=BigFloat)
    ref = setprecision(() -> ctx_from(BigFloat, particle(Float64, P2)), BigFloat, 512)
    @test eltype(hi) === Float64
    @test all(getfield(hi, f) == Float64(getfield(ref, f)) for f in EffSourceEquatorial.COEFF_FIELDS)
end

@testset "set_particle! cost" begin
    mkctx(Float64, P2)
    mkctx(BigFloat, P2)
    t64 = @elapsed mkctx(Float64, P2)
    tbig = @elapsed mkctx(BigFloat, P2)
    @info "set_particle! wall time" t64 tbig
    @test t64 < 2.0
    @test tbig < 10.0
end
