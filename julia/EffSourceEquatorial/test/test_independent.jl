function notebook_reference()
    q(s) = (p = split(s, "/"); length(p) == 1 ? BigFloat(parse(BigInt, p[1])) :
            BigFloat(parse(BigInt, p[1]))/BigFloat(parse(BigInt, p[2])))
    parts, real, modes = Dict{String,Any}(), [], []
    for f in split.(readlines(joinpath(HERE, "notebook_reference.txt")))
        if f[1] == "particle"
            parts[String(f[2])] = (a=q(f[3]), p=q(f[4]), e=q(f[5]), r=q(f[6]), sign=parse(Int, f[7]))
        elseif f[1] == "real"
            push!(real, (label=String(f[2]), dr=q(f[3]), dth=q(f[4]), dph=q(f[5]), vals=parse.(BigFloat, f[6:9])))
        else
            push!(modes, (label=String(f[2]), m=parse(Int, f[3]), dr=q(f[4]), dth=q(f[5]),
                          val=complex(parse(BigFloat, f[6]), parse(BigFloat, f[7]))))
        end
    end
    return parts, real, modes
end

function notebook_context(P; ell=false)
    M = one(BigFloat)
    E, L = kerr_equatorial_EL(M, P.a, P.p, P.e)
    ur = P.sign == 0 ? zero(BigFloat) : kerr_equatorial_ur(M, P.a, E, L, P.r; sign=P.sign)
    c = effsource_equatorial_create(M, P.a; ell_from_complement=ell)
    set_particle!(c, Coordinate{BigFloat}(r=P.r, theta=halfpi(BigFloat), phi=zero(BigFloat), t=zero(BigFloat)), E, L, ur)
    return c
end

@testset "puncture and m-modes equal ScalarSingularField.nb along the orbits" begin
    setprecision(BigFloat, 256) do
        parts, real, modes = notebook_reference()
        ctxs = Dict(k => notebook_context(P) for (k, P) in parts)
        cell = Dict(k => notebook_context(P; ell=true) for (k, P) in parts)
        for rec in real
            ctx = ctxs[rec.label]
            off = Offset{BigFloat}(dr=rec.dr, dtheta=rec.dth, dphi=rec.dph)
            d1 = calc(ctx, off)[2]
            for (g, w) in zip((phis(ctx, off), d1[2], d1[3], d1[4]), rec.vals)
                @test abs(g - w) <= big"1e-40"*abs(w)
            end
        end
        for rec in modes, ctx in (ctxs[rec.label], cell[rec.label])
            pm = phis_m(ctx, rec.m, Offset{BigFloat}(dr=rec.dr, dtheta=rec.dth))
            @test abs(pm - rec.val) <= big"1e-40"*abs(rec.val)
        end
    end
end

@testset "orbital kinematics match the geodesic (sign $sgn)" for sgn in (1, -1)
    setprecision(BigFloat, 512) do
        s = merge(ORBIT_A, (r=9, sign=sgn))
        P = particle(BigFloat, s)
        c = mkctx(BigFloat, s)
        M, a, r, E, L = P.M, P.a, P.r, P.E, P.L
        Δ = r^2 - 2M*r + a^2
        ut = ((r^2 + a^2 + 2M*a^2/r)*E - 2M*a*L/r)/Δ
        uphi = ((1 - 2M/r)*L + 2M*a*E/r)/Δ
        V(x) = EffSourceEquatorial._Rpot(M, a, E, L, x)/x^4
        h = big"1e-40"
        @test abs(c.rt - P.ur/ut) < big"1e-60"
        @test abs(c.phit - uphi/ut) < big"1e-60"
        @test abs(c.urt - (V(r + h) - V(r - h))/(4h*ut)) < big"1e-60"
    end
end

@testset "derivatives equal finite differences (BigFloat, 512 bit, sign $sgn)" for sgn in (1, -1)
    setprecision(BigFloat, 512) do
        s = merge(ORBIT_A, (r=9, sign=sgn))
        P = particle(BigFloat, s)
        c0 = mkctx(BigFloat, s)
        h = big"1e-30"
        tol = big"1e-40"
        rel(u, v) = abs(u - v)/abs(v)
        fd1(f) = (f(h) - f(-h))/(2h)
        fd2(f) = (f(h) - 2f(zero(h)) + f(-h))/h^2
        function moved(t)
            r = P.r + c0.rt*t + c0.rtt*t^2/2
            c = effsource_equatorial_create(P.M, P.a)
            xp = Coordinate{BigFloat}(r=r, theta=halfpi(BigFloat),
                                      phi=thirdpi(BigFloat) + c0.phit*t + c0.phitt*t^2/2, t=one(BigFloat) + t)
            set_particle!(c, xp, P.E, P.L, kerr_equatorial_ur(P.M, P.a, P.E, P.L, r; sign=sgn))
            return c
        end
        cp, cm = moved(h), moved(-h)
        for o0 in ((big"0.3", big"0.05", big"0.2"), (big"1e-3", big"1e-3", big"1e-3"), (big"1e-6", big"1e-6", big"1e-6"))
            off(e, n) = Offset{BigFloat}(dr=o0[1] + (n == 1)*e, dtheta=o0[2] + (n == 2)*e, dphi=o0[3] + (n == 3)*e)
            _, d1, d2, _ = calc(c0, off(zero(h), 0))
            for (n, i1, i2) in ((1, 2, 5), (2, 3, 8), (3, 4, 10))
                @test rel(fd1(e -> phis(c0, off(e, n))), d1[i1]) < tol
                @test rel(fd2(e -> phis(c0, off(e, n))), d2[i2]) < tol
            end
            for m in (2, 10)
                _, d1m, d2m, _ = calc_m(c0, m, off(zero(h), 0))
                for (n, i1, i2) in ((1, 2, 5), (2, 3, 8))
                    @test rel(fd1(e -> phis_m(c0, m, off(e, n))), d1m[i1]) < tol
                    @test rel(fd2(e -> phis_m(c0, m, off(e, n))), d2m[i2]) < tol
                end
            end
            x = Coordinate{BigFloat}(r=P.r + o0[1], theta=halfpi(BigFloat) + o0[2],
                                     phi=thirdpi(BigFloat) + o0[3], t=zero(BigFloat))
            _, d1, d2, _ = calc(c0, x)
            @test rel((phis(cp, x) - phis(cm, x))/(2h), d1[1]) < tol
            @test rel((phis(cp, x) - 2phis(c0, x) + phis(cm, x))/h^2, d2[1]) < tol
            @test rel((calc(cp, x)[2][4] - calc(cm, x)[2][4])/(2h), d2[4]) < tol
            for m in (2, 10)
                _, d1m, d2m, _ = calc_m(c0, m, x)
                @test rel((phis_m(cp, m, x) - phis_m(cm, m, x))/(2h), d1m[1]) < tol
                @test rel((phis_m(cp, m, x) - 2phis_m(c0, m, x) + phis_m(cm, m, x))/h^2, d2m[1]) < tol
                @test rel((calc_m(cp, m, x)[2][4] - calc_m(cm, m, x)[2][4])/(2h), d2m[4]) < tol
            end
        end
    end
end

@testset "effective source vanishes linearly at the particle" begin
    setprecision(BigFloat, 512) do
        ctx = mkctx(BigFloat, P2)
        s(d) = abs(calc(ctx, Offset{BigFloat}(dr=d, dtheta=d, dphi=d))[4])
        @test abs(s(big"1e-6")/s(big"1e-3") - big"1e-3") < big"1e-4"
    end
end

@testset "circular limit reproduces EffSourceCircular ($T)" for T in (Float64, BigFloat)
    tolP, tolS = T === Float64 ? (1e-12, 1e-4) : (1e-50, 1e-50)
    M, a, rp = T(1), T(1)/2, T(10)
    E, L = kerr_equatorial_EL(M, a, rp, zero(T))
    xp = Coordinate{T}(r=rp, theta=halfpi(T), phi=thirdpi(T), t=one(T))
    ce = effsource_equatorial_create(M, a)
    set_particle!(ce, xp, E, L, zero(T))
    cc = EffSourceCircular.effsource_create(M, a)
    set_particle!(cc, xp, E, L)
    rel(u, v) = abs(u - v)/max(abs(u), abs(v))
    for (dr, dth, dph) in ((1//10, 1//100, 1//50), (1//2, 1//20, 3//10), (-3//10, 0, 1//10))
        x = Coordinate{T}(r=rp + tq(T, dr), theta=halfpi(T) + tq(T, dth), phi=thirdpi(T) + tq(T, dph), t=zero(T))
        @test rel(phis(ce, x), phis(cc, x)) < tolP
        @test norm(calc(ce, x)[2] - calc(cc, x)[2]) <= 100tolP*norm(calc(cc, x)[2])
        @test rel(calc(ce, x)[4], calc(cc, x)[4]) < tolS
        for m in (0, 2, 10)
            @test rel(phis_m(ce, m, x), phis_m(cc, m, x)) < tolP
            Pe, d1e, _, se = calc_m(ce, m, x)
            Pc, d1c, _, sc = calc_m(cc, m, x)
            @test norm(d1e - d1c) <= 100tolP*norm(d1c)
            @test rel(se, sc) < tolS
        end
    end
end

@testset "m-modes equal phi-quadratures of the real-space fields ($T)" for T in (Float64, BigFloat)
    isbig = T === BigFloat
    tol = isbig ? (P=1e-45, d1=1e-45, tt=1e-45, src=1e-45) : (P=1e-12, d1=1e-11, tt=1e-11, src=1e-6)
    P = particle(T, P2)
    ctx = mkctx(T, P2)
    N = 4096
    geoms = [(T(1)/2, T(1)/20), (T(3)/10, zero(T)), (-T(2)/5, T(1)/10)]
    isbig && push!(geoms, (T(2), T(3)/10))
    modes = isbig ? (0, 2, 10, 20) : (0, 2, 10)
    rel(u, v) = abs(u - v)/abs(v)
    for (dr, dth) in geoms
        φs = [2T(pi)*k/N for k in 0:N-1]
        samples = [calc(ctx, Coordinate{T}(r=P.r + dr, theta=halfpi(T) + dth, phi=φ, t=zero(T))) for φ in φs]
        for m in modes
            w = [cis(-m*φ)*(2T(pi)/N) for φ in φs]
            q(f) = sum(w .* [f(s) for s in samples])
            Pm, d1m, d2m, sm = calc_m(ctx, m, Coordinate{T}(r=P.r + dr, theta=halfpi(T) + dth, phi=zero(T), t=zero(T)))
            @test rel(Pm, q(s -> s[1])) < tol.P
            @test rel(d1m[1], q(s -> s[2][1])) < tol.d1
            @test rel(d1m[2], q(s -> s[2][2])) < tol.d1
            iszero(dth) || @test rel(d1m[3], q(s -> s[2][3])) < tol.d1
            iszero(m) || @test rel(d1m[4], q(s -> s[2][4])) < tol.d1
            @test rel(d2m[1], q(s -> s[3][1])) < tol.tt
            @test rel(sm, q(s -> s[4])) < tol.src
        end
    end
end

@testset "numeric type genericity" begin
    for T in (Float32, Float32sr, Double64, BigFloat)
        ctx = T === Float32sr ? ctx_from(T, particle(Float32, P2)) : mkctx(T, P2)
        x = Coordinate{T}(r=tq(T, 9) + T(1)/10, theta=halfpi(T) + T(1)/100, phi=thirdpi(T) + T(1)/50, t=zero(T))
        @test eltype(ctx) === T
        @test phis(ctx, x) isa T
        Pr, d1, d2, s = calc(ctx, x)
        @test Pr isa T && eltype(d1) === T && eltype(d2) === T && s isa T
        @test phis_m(ctx, 2, x) isa Complex{T}
        Pm, d1m, d2m, sm = calc_m(ctx, 2, x)
        @test Pm isa Complex{T} && eltype(d1m) === Complex{T} && sm isa Complex{T}
    end
    lo = mkctx(Float32, P2; setup=Float64)
    @test all(isfinite(getfield(lo, f)) for f in EffSourceEquatorial.COEFF_FIELDS)
end
