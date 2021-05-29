export periapsis, apoapsis, semimajoraxis, eccentricity
export meananomaly, trueanomaly, eccentricanomaly
export orbitalperiod, orbitaldistance, orbit

"""
    periapsis(a, e)

Compute the periapsis (closest approach) distance using semi-major axis and eccentricity
"""
periapsis(a, e) = a*(1 - e)

"""
    apoapsis(a, e)

Compute the apoapsis (farthest distance) distance using semi-major axis and eccentricity
"""
apoapsis(a, e) = a*(1 + e)

"""
    semimajoraxis(T, m)

Compute the [semi-major axis](https://en.wikipedia.org/wiki/Semi-major_and_semi-minor_axes) of an orbit
"""
semimajoraxis(T, m) = (𝐆*m*T^2/(4π^2))^(1/3)

"""
    eccentricity(rₚ, rₐ)

Compute [eccentricity](https://en.wikipedia.org/wiki/Orbital_eccentricity)
"""
eccentricity(rₚ, rₐ) = (rₐ - rₚ)/(rₐ + rₚ)

"""
    meananomaly(E, e)

Compute the [mean anomaly](https://en.wikipedia.org/wiki/Mean_anomaly)
"""
meananomaly(E, e) = E - e*sin(E)

"""
    trueanomaly(E, e)

Compute the [true anomaly](https://en.wikipedia.org/wiki/True_anomaly)
"""
function trueanomaly(E, e)
    f = 2*atan(sqrt((1 + e)/(1 - e))*tan(E/2))
    #use [0,2π] instead of [-π,π]
    f < 0 ? f + 2π : f
end

"""
    trueanomaly(t, a, m, e)

Compute the [true anomaly](https://en.wikipedia.org/wiki/True_anomaly)
"""
trueanomaly(t, a, m, e) = trueanomaly(eccentricanomaly(t, a, m, e), e)

"""
    eccentricanomaly(t, a, m, e)

Numerically compute the [eccentric anomaly](https://en.wikipedia.org/wiki/Eccentric_anomaly) using [Kepler's equation](https://en.wikipedia.org/wiki/Kepler%27s_equation)
"""
function eccentricanomaly(t, a, m, e)
    @assert t >= 0 "time must be positive"
    #Kepler's Third Law
    T = orbitalperiod(a, m)
    #definition of mean anomaly
    M = 2π*rem(t, T)/T
    #eccentric anomaly must be found numerically
    E = falseposition(E->meananomaly(E, e) - M, 0, 2π)
    return E
end

"""
    orbitalperiod(a, m)

[Kepler's Third Law](https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Third_law) describing the orbital period of an elliptical orbit
"""
orbitalperiod(a, m) = 2π*√(a^3/(𝐆*m))

"""
    orbitaldistance(a, f, e)

Compute the distance of a planet from its host
"""
orbitaldistance(a, f, e) = a*(1 - e^2)/(1 + e*cos(f))

"""
    orbitaldistance(t, a, m, e)

Compute the distance of a planet from its host
"""
orbitaldistance(t, a, m, e) = orbitaldistance(a, trueanomaly(t, a, m, e), e)

"""
    orbit(a, m, e, N=1000)

Create a distance time-series of `N` points for an elliptical orbit
"""
function orbit(a, m, e, N::Int=1000)
    T = orbitalperiod(a, m)
    t = LinRange(0, T, N+1)[1:end-1]
    f = trueanomaly.(t, a, m, e)
    r = orbitaldistance.(a, f, e)
    return t, r, f
end
