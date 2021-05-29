export substellarlatitude, hourangle
export diurnalfluxfactor, diurnalfluxfactors
export annualfluxfactor, annualfluxfactors

"""
    substellarlatitude(f, γ)

Compute the latitude of the substellar point for a given solar longitude (true anomaly) and obliquity
"""
substellarlatitude(f, γ) = asin(cos(f)*sin(γ))

"""
    hourangle(ϕ, ϕₛ)

Compute the [hour angle](https://en.wikipedia.org/wiki/Hour_angle)
"""
function hourangle(ϕ, ϕₛ)
    x = -sin(ϕ)*sin(ϕₛ)/(cos(ϕ)*cos(ϕₛ))
    if x <= -1
        return π
    elseif x >= 1
        return 0.0
    end
    return acos(x)
end

#ϕ - latitude
#ϕₛ - substellar latitude
"""
    diurnalfluxfactor(ϕ, ϕₛ)

Compute the diurnally averaged fraction of incoming stellar flux received by a point at latitude `ϕ` when the substellar latitude is `ϕₛ`
"""
function diurnalfluxfactor(ϕ, ϕₛ)
    h = hourangle(ϕ, ϕₛ)
    (sin(h)*cos(ϕ)*cos(ϕₛ) + h*sin(ϕ)*sin(ϕₛ))/π
end

#ϕ - latitude
#f - solar longitude
#γ - obliquity
"""
    diurnalfluxfactor(ϕ, f, γ)

Compute the diurnally averaged fraction of incoming stellar flux received by a point at latitude `ϕ` when the planet is at solar longitude (true anomaly) `f`, with obliquity `γ`
"""
diurnalfluxfactor(ϕ, f, γ) = diurnalfluxfactor(ϕ, substellarlatitude(f, γ))

#t - time
#a - semimajor
#m - mass
#e - eccentricity
#ϕ - latitude
#γ - obliquity
#p - precession angle
"""
    diurnalfluxfactor(t, a, m, e, ϕ, γ, p)

Compute the diurnally averaged fraction of incoming stellar flux received by a point at latitude `ϕ` for a general elliptical orbit
"""
function diurnalfluxfactor(t, a, m, e, ϕ, γ, p)
    f = trueanomaly(t, a, m, e)
    r = orbitaldistance(a, f, e)
    diurnalfluxfactor(ϕ, f + p, γ)*(a/r)^2
end

"""
    diurnalfluxfactors(γ; nf=251, nϕ=181)

Compute a grid of diurnally averaged fractions of incoming stellar flux received by a point at latitude `ϕ` for a planet with obliquity `γ` in a circular orbit. Returns  a solar longitude vector (column values), latitude vector (row values), and the grid of flux factors. `nf` indicates the number of points around the orbit and `nϕ` indicates the number of latitudes.
"""
function diurnalfluxfactors(γ; nf::Int=251, nϕ::Int=181)
    ϕ = LinRange(-π/2, π/2, nϕ)
    f = LinRange(0, 2π, nf)
    F, Φ = meshgrid(f, ϕ)
    (f, ϕ, diurnalfluxfactor.(Φ, F, γ))
end

"""
    diurnalfluxfactors(a, m, e, γ, p; nt=251, nϕ=181)

Compute a grid of diurnally averaged fractions of incoming stellar flux for a planet in a general elliptical orbit. Returns a time vector (column values), latitude vector (row values), and the grid of flux factors. `nt` indicates the number of time samples around the orbit and `nϕ` indicates the number of latitudes.
"""
function diurnalfluxfactors(a, m, e, γ, p; nt::Int=251, nϕ::Int=181)
    t = LinRange(0, orbitalperiod(a, m), nt)
    ϕ = LinRange(-π/2, π/2, nϕ)
    T, Φ = meshgrid(t, ϕ)
    (t, ϕ, diurnalfluxfactor.(T, a, m, e, Φ, γ, p))
end

"""
    annualfluxfactor(a, m, e, ϕ, γ, p)

Compute the annually averaged flux factor for a latitude `ϕ` on a planet in a general elliptical orbit
"""
function annualfluxfactor(a, m, e, ϕ, γ, p)
    T = orbitalperiod(a, m)
    F, err = hquadrature(t->diurnalfluxfactor(t, a, m, e, ϕ, γ, p), 0, T)
    return F/T
end

"""
    annualfluxfactor(a, m, e, ϕ, γ, p; nϕ=181)

Compute a range of annually averaged flux factors for a planet in a general elliptical orbit. Returns a latitude vector (row values) and a vector of flux factors. `nϕ` indicates the number of latitude samples.
"""
function annualfluxfactors(a, m, e, γ, p; nϕ::Int=181)
    ϕ = LinRange(-π/2, π/2, nϕ)
    F = annualfluxfactor.(a, m, e, ϕ, γ, p)
    return ϕ, F
end
