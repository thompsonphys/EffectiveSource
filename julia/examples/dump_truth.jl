# Dump BigFloat-truth values of the m-mode puncture and effective source at a
# set of small offsets from a circular-orbit particle, for validating the
# eccentric (equatorial) C code in its circular limit (ur = 0).
#
#   julia --project examples/dump_truth.jl
#
# Writes examples/truth_circular.csv with columns:
#   m, dr, dtheta, RePhiS, ImPhiS, Resrc, Imsrc

using EffSourceCircular
using Printf

const M, a = 1.0, 0.5
const rp, thp, php = 10.0, pi/2, pi/3
const E = sqrt(1 - 2M/rp + a^2/rp^2 - 2M*a^2/rp^3) / sqrt(1 - 3M/rp + 2a*sqrt(M/rp^3))
const L = (rp^2 - 2M*rp + a*sqrt(M*rp)) / (rp*sqrt(rp - 3M + 2a*sqrt(M/rp)))

setprecision(BigFloat, 300)

ctx = effsource_create(big(M), big(a))
set_particle!(ctx, Coordinate{BigFloat}(r=big(rp), theta=big(thp), phi=big(php), t=big(1)),
              big(E), big(L))

const DISTS = 10 .^ range(-9, -2; length=15)
const RAYS = [(1.0, 0.0), (cos(pi/4), sin(pi/4)), (0.0, 1.0)]
const MODES = (0, 2, 10)

open(joinpath(@__DIR__, "truth_circular.csv"), "w") do io
    println(io, "m,dr,dtheta,RePhiS,ImPhiS,Resrc,Imsrc")
    for m in MODES, (cr, ct) in RAYS, s in DISTS
        dr64  = Float64(s*cr)
        dth64 = Float64(s*ct)
        o = Offset{BigFloat}(dr=big(dr64), dtheta=big(dth64), dphi=big(0.0))
        p = phis_m(ctx, m, o)
        q = calc_m(ctx, m, o)[4]
        @printf(io, "%d,%.17e,%.17e,%.17e,%.17e,%.17e,%.17e\n",
                m, dr64, dth64,
                Float64(real(p)), Float64(imag(p)), Float64(real(q)), Float64(imag(q)))
    end
end
println("wrote truth_circular.csv")
