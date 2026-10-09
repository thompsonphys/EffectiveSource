@testset "geodesic constants" begin
    M, a, p, e = 1.0, 0.5, 10.0, 0.3
    E, L = kerr_equatorial_EL(M, a, p, e)
    @test 0 < E < 1
    @test L > 0
    @test abs(EffSourceEquatorial._Rpot(M, a, E, L, p/(1 + e))) < 1e-9
    @test abs(EffSourceEquatorial._Rpot(M, a, E, L, p/(1 - e))) < 1e-9
    @test kerr_equatorial_EL(1, 0.5, 10, 0.3) == (E, L)

    r = 9.0
    ur = kerr_equatorial_ur(M, a, E, L, r)
    Δ = r^2 - 2M*r + a^2
    ut = ((r^2 + a^2 + 2M*a^2/r)*E - 2M*a*L/r)/Δ
    uphi = ((1 - 2M/r)*L + 2M*a*E/r)/Δ
    norm2 = -(1 - 2M/r)*ut^2 - 4M*a/r*ut*uphi + r^2/Δ*ur^2 + (r^2 + a^2 + 2M*a^2/r)*uphi^2
    @test norm2 ≈ -1 rtol=1e-13
    @test kerr_equatorial_ur(M, a, E, L, r; sign=-1) == -ur

    Ec, Lc = kerr_equatorial_EL(M, a, p, 0.0)
    Δp = p^2 - 2M*p + a^2
    utc = ((p^2 + a^2 + 2M*a^2/p)*Ec - 2M*a*Lc/p)/Δp
    uphic = ((1 - 2M/p)*Lc + 2M*a*Ec/p)/Δp
    @test uphic/utc ≈ 1/(p^1.5 + a) rtol=1e-14
    @test abs(EffSourceEquatorial._Rpot(M, a, Ec, Lc, p)) < 1e-9

    setprecision(BigFloat, 256) do
        Eb, Lb = kerr_equatorial_EL(big(M), big(a), big(p), big(e))
        @test abs(EffSourceEquatorial._Rpot(big(M), big(a), Eb, Lb, big(p)/(1 + big(e)))) < big(10)^-60
        @test abs(EffSourceEquatorial._Rpot(big(M), big(a), Eb, Lb, big(p)/(1 - big(e)))) < big(10)^-60
        @test abs(Float64(Eb) - E) < 1e-13
        @test abs(Float64(Lb) - L) < 1e-12
    end
end
