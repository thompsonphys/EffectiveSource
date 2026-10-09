"""Energy and angular momentum of the prograde circular equatorial Kerr geodesic of radius `r`."""
function kerr_circular_EL(M::T, a::T, r::T) where {T<:AbstractFloat}
    v = sqrt(M/r)
    den = sqrt(1 - 3v^2 + 2a/M*v^3)
    return (1 - 2v^2 + a/M*v^3)/den, sqrt(M*r)*(1 - 2a/M*v^3 + a^2/M^2*v^4)/den
end

kerr_circular_EL(M::Real, a::Real, r::Real) = kerr_circular_EL(promote(float(M), float(a), float(r))...)
