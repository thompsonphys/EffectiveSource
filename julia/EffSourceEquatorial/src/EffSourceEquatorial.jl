"""
    EffSourceEquatorial

Scalar singular field and effective source for a point particle on an eccentric
equatorial geodesic in Kerr spacetime, generic over the numeric type `T`.
"""
module EffSourceEquatorial

using OffsetArrays
using StaticArrays
using EffSourceCircular: Coordinate, Offset, reei, _agm_KE, _agm_KE_c, _ovec, _ovec0, kerr_circular_EL
import EffSourceCircular: set_particle!, phis, calc, phis_m, calc_m

export Coordinate, Offset, EquatorialCtx, effsource_equatorial_create, set_particle!,
       phis, calc, phis_m, calc_m, kerr_equatorial_EL, kerr_equatorial_ur

include("geodesic.jl")
include("_eq_ctx_struct.jl")

module _Coefficients
Base.Experimental.@compiler_options compile=min infer=false optimize=0
include("_eq_coeffs.jl")
end

include("_eq_set_particle.jl")
include("imei_data.jl")

function _ellKE(ctx::EquatorialCtx{T}, C1::T) where {T}
    ctx.ell_from_complement && return _agm_KE_c(C1/(one(T) + C1))
    return _agm_KE(sqrt(one(T)/(one(T) + C1)))
end

include("_eq_kernels.jl")

end
