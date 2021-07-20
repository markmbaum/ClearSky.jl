export substellarlatitude, hourangle
export diurnalfluxfactor, diurnalfluxfactors
export annualfluxfactor, annualfluxfactors

"""
    substellarlatitude(f, γ)

Compute the latitude of the substellar point for a given solar longitude `f` (true anomaly) and obliquity `γ`
"""
substellarlatitude(f, γ) = asin(cos(f)*sin(γ))

"""
    hourangle(θ, θₛ)

Compute the [hour angle](https://en.wikipedia.org/wiki/Hour_angle)
"""
function hourangle(θ, θₛ)
    x = -sin(θ)*sin(θₛ)/(cos(θ)*cos(θₛ))
    if x <= -1
        return π
    elseif x >= 1
        return 0.0
    end
    return acos(x)
end

#θ - latitude
#θₛ - substellar latitude
"""
    diurnalfluxfactor(θ, θₛ)

Compute the diurnally averaged fraction of incoming stellar flux received by a point at latitude `θ` when the substellar latitude is `θₛ`
"""
function diurnalfluxfactor(θ, θₛ)
    h = hourangle(θ, θₛ)
    return (sin(h)*cos(θ)*cos(θₛ) + h*sin(θ)*sin(θₛ))/π
end

#θ - latitude
#f - solar longitude
#γ - obliquity
"""
    diurnalfluxfactor(θ, f, γ)

Compute the diurnally averaged fraction of incoming stellar flux received by a point at latitude `θ` when the planet is at solar longitude (true anomaly) `f`, with obliquity `γ`
"""
diurnalfluxfactor(θ, f, γ) = diurnalfluxfactor(θ, substellarlatitude(f, γ))

"""
    diurnalfluxfactor(t, a, m, e, θ, γ, p)

Compute the diurnally averaged fraction of incoming stellar flux received by a point at latitude `θ` for a general elliptical orbit
"""
function diurnalfluxfactor(t, a, m, e, θ, γ, p)
    f = trueanomaly(t, a, m, e)
    r = orbitaldistance(a, f, e)
    return diurnalfluxfactor(θ, f - p, γ)*(a/r)^2
end

"""
    diurnalfluxfactors(γ; nf=251, nθ=181)

Compute a grid of diurnally averaged fractions of incoming stellar flux received by a point at latitude `θ` for a planet with obliquity `γ` in a circular orbit. Returns  a solar longitude vector (column values), latitude vector (row values), and the grid of flux factors. `nf` indicates the number of points around the orbit and `nθ` indicates the number of latitudes.
"""
function diurnalfluxfactors(γ; nf::Int=251, nθ::Int=181)
    θ = LinRange(-π/2, π/2, nθ)
    f = LinRange(0, 2π, nf)
    F, Θ = meshgrid(f, θ)
    return (f, θ, diurnalfluxfactor.(Θ, F, γ))
end

"""
    diurnalfluxfactors(a, m, e, γ, p; nt=251, nθ=181)

Compute a grid of diurnally averaged fractions of incoming stellar flux for a planet in a general elliptical orbit. Returns a time vector (column values) over one orbital period, latitude vector (row values), and the grid of flux factors. `nt` indicates the number of time samples around the orbit and `nθ` indicates the number of latitudes.
"""
function diurnalfluxfactors(a, m, e, γ, p; nt::Int=251, nθ::Int=181)
    t = LinRange(0, orbitalperiod(a, m), nt)
    θ = LinRange(-π/2, π/2, nθ)
    T, Θ = meshgrid(t, θ)
    return (t, θ, diurnalfluxfactor.(T, a, m, e, Θ, γ, p))
end

"""
    annualfluxfactor(e, θ, γ, p)

Compute the annually averaged flux factor for a latitude `θ` on a planet in a general elliptical orbit.
"""
function annualfluxfactor(e, θ, γ, p; tol::Float64=1e-4)
    T = orbitalperiod(1.0, 1.0)
    f(t) = diurnalfluxfactor(t, 1.0, 1.0, e, θ, γ, p)
    F, _ = hquadrature(f, 0, T, reltol=tol, abstol=tol)
    return F/T
end

"""
    annualfluxfactors(e, γ, p; nθ=181)

Compute a range of annually averaged flux factors for a planet in a general elliptical orbit. Returns a latitude vector (row values) and a vector of flux factors. `nθ` indicates the number of latitude samples.
"""
function annualfluxfactors(e, γ, p; nθ::Int=181)
    θ = LinRange(-π/2, π/2, nθ)
    F = annualfluxfactor.(e, θ, γ, p)
    return θ, F
end
