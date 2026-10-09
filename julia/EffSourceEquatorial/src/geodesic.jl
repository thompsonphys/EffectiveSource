_Rpot(M, a, E, L, r) = (E*(r^2 + a^2) - a*L)^2 - (r^2 - 2M*r + a^2)*(r^2 + (L - a*E)^2)

function _Rpot_jac(M, a, E, L, r)
    X = E*(r^2 + a^2) - a*L
    Y = L - a*E
    Δ = r^2 - 2M*r + a^2
    return 2X*(r^2 + a^2) + 2a*Δ*Y, -2a*X - 2Δ*Y
end

"""Energy and angular momentum of the prograde equatorial Kerr geodesic with semi-latus rectum `p` and eccentricity `e`."""
function kerr_equatorial_EL(M::T, a::T, p::T, e::T) where {T<:AbstractFloat}
    E, L = kerr_circular_EL(M, a, p)
    iszero(e) && return E, L
    r1, r2 = p/(1 + e), p/(1 - e)
    for _ in 1:100
        F1, F2 = _Rpot(M, a, E, L, r1), _Rpot(M, a, E, L, r2)
        J11, J12 = _Rpot_jac(M, a, E, L, r1)
        J21, J22 = _Rpot_jac(M, a, E, L, r2)
        det = J11*J22 - J12*J21
        dE = (F1*J22 - J12*F2)/det
        dL = (J11*F2 - F1*J21)/det
        E -= dE
        L -= dL
        abs(dE) + abs(dL) <= 4eps(T)*(abs(E) + abs(L)) && break
    end
    return E, L
end

kerr_equatorial_EL(M::Real, a::Real, p::Real, e::Real) =
    kerr_equatorial_EL(promote(float(M), float(a), float(p), float(e))...)

"""Radial four-velocity `u^r` at radius `r` on the equatorial geodesic with constants `E`, `L`; `sign` selects the branch."""
function kerr_equatorial_ur(M::T, a::T, E::T, L::T, r::T; sign::Integer=1) where {T<:AbstractFloat}
    return sign*sqrt(max(_Rpot(M, a, E, L, r), zero(T)))/r^2
end

kerr_equatorial_ur(M::Real, a::Real, E::Real, L::Real, r::Real; sign::Integer=1) =
    kerr_equatorial_ur(promote(float(M), float(a), float(E), float(L), float(r))...; sign)
