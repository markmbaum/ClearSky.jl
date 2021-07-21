#-------------------------------------------------------------------------------
export opticaldepth

"""
    opticaldepth(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute monochromatic [optical depths](https://en.wikipedia.org/wiki/Optical_depth#Atmospheric_sciences) (τ) between two pressure levels

# Arguments

* `P₁`: first pressure level [Pa]
* `P₂`: second pressure level [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `θ`: angle [radians] of path, must be ∈ [0,π/2), where 0 is straight up/down
* `absorbers`: at least one gas object and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

Returns a vector of optical depths across all wavenumbers stored in gas objects. The `tol` keyword argument adjusts integrator error tolerance.
"""
function opticaldepth(P₁::Real,
                      P₂::Real,
                      g::Real,
                      fT::Q,
                      fμ::R,
                      θ::Real,
                      absorbers...;
                      tol::Float64=1e-5
                      )::Vector{Float64} where {Q,R}
    #initialization
    P₁, P₂ = max(P₁, P₂), min(P₁, P₂)
    A = setupintegration(P₁, P₂, absorbers)
    checkazimuth(θ)
    #integrate wavenumbers in parallel
    τ = zeros(Float64, A.nν)
    m = 1/cos(θ)
    @threads for i ∈ eachindex(A.ν)
        #integrate
        τ[i] = depth(dτdω, P2ω(P₁), P2ω(P₂), A, i, g, m, fT, fμ, tol)
    end
    return τ
end

#-------------------------------------------------------------------------------
#the basic transmittance method exp(-τ) is already exported

"""
    transmittance(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute monochromatic [transmittances](https://en.wikipedia.org/wiki/Transmittance.) between two pressure levels

Accepts the same arguments as [`opticaldepth`](@ref) and returns a vector of transmittances across all wavenumbers stored in gas objects.
"""
transmittance(X...; kwargs...) = transmittance.(opticaldepth(X...; kwargs...))

#-------------------------------------------------------------------------------
export outgoing

"""
    outgoing(Pₛ, Pₜ, g, fT, fμ, absorbers; nstream=5, tol=1e-5)

Compute outgoing monochromatic radiative fluxes [W/m``^2``/cm``^{-1}``], line-by-line. Integrates the [`schwarzschild`](@ref) equation from `Pₛ` to `Pₜ` at each wavenumber in the provided gas object(s) using any number of streams/angles. Total flux [W/m``^2``] can be evaluated with the [`trapz`](@ref) function.

# Arguments
* `Pₛ`: surface pressure [Pa]
* `Pₜ`: top of atmopshere pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `absorbers`: at least one [gas object](gas_objects.md) and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

The keyword argument `nstream` specifies how many independent streams, or beam angles through the atmosphere, to integrate. The keyword argument `tol` is a numerical error tolerance passed to the [`radau`](https://github.com/wordsworthgroup/ScalarRadau.jl) integrator.
"""
function outgoing(Pₛ::Real,
                  Pₜ::Real,
                  g::Real,
                  fT::Q,
                  fμ::R,
                  absorbers...;
                  nstream::Int64=5, #number of streams to use
                  tol::Float64=1e-5 #integrator tolerance
                  )::Vector{Float64} where {Q,R}
    #initialization
    ca = setupintegration(Pₛ, Pₜ, absorbers)
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    #surface temperature
    T₀ = fT(Pₛ)
    #integrate wavenumbers in parallel
    F = zeros(Float64, ca.nν)
    @threads for i ∈ eachindex(ca.ν)
        I₀ = planck(ca.ν[i], T₀)
        F[i] = streams(dIdω, I₀, ω₁, ω₂, ca, i, g, nstream, fT, fμ, tol)
    end
    return F
end

#-------------------------------------------------------------------------------
export topfluxes, topimbalance

"""
    topfluxes(Pₛ, Pₜ, g, fT, fμ, absorber::UnifiedAbsorber; nstream=5, θ=0.952, tol=1e-5)

Compute the upward and downward monochromatic fluxes at the top of the atmopshere.

# Arguments
* `Pₛ`: surface pressure [Pa]
* `Pₜ`: top of atmopshere pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `absorber`: a `UnifiedAbsorber`, which is either a [`GroupedAbsorber`](@ref) or an [`AcceleratedAbsorber`](@ref)
"""
function topfluxes(Pₛ::Real,
                   Pₜ::Real,
                   g::Real,
                   fT::Q, # fT(P)
                   fμ::R, # fμ(T,P)
                   fstellar::S, # fstellar(ν) [W/m^2]
                   albedo::Real,
                   absorber::UnifiedAbsorber;
                   nstream::Int64=5, #number of streams to use
                   θ::Float64=0.952, #corresponds to cos(θ) = 0.66
                   tol::Float64=1e-5)::NTuple{2,Vector{Float64}} where {Q,R,S}
    #setup
    checkpressures(absorber, Pₛ, Pₜ)
    checkazimuth(θ)
    ν, nν = absorber.ν, absorber.nν
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    ι₁, ι₂ = P2ι(Pₜ, Pₛ)

    #incoming stellar flux at TOA
    Fₜ⁻ = fstellar.(ν)*cos(θ)
    #downward stellar flux at the ground/surface
    Fₛ⁻ = zeros(nν)
    @threads for i ∈ eachindex(ν)
        τ = depth(dτdι, ι₁, ι₂, absorber, i, g, 1/cos(θ), fT, fμ, tol)
        Fₛ⁻[i] = Fₜ⁻[i]*exp(-τ)
    end
    #some of it gets reflected back upward at the surface
    Iₛ⁺ = Fₛ⁻*albedo/π #Lambertian reflection
    #surface temperature
    Tₛ = fT(Pₛ)
    #radiation streams up out of the atmosphere
    Fₜ⁺ = zeros(nν)
    @threads for i ∈ eachindex(ν)
        I₀ = Iₛ⁺[i] + planck(ν[i], Tₛ)
        Fₜ⁺[i] = streams(dIdω, I₀, ω₁, ω₂, absorber, i, g, nstream, fT, fμ, tol)
    end

    return Fₜ⁻, Fₜ⁺
end

"""
    topimbalance(Pₛ, Pₜ, g, fT, fμ, absorber::UnifiedAbsorber; nstream=5, θ=0.952, tol=1e-5)

Compute the difference between the total upward and downward radiative flux at the top of the atmopshere.

# Arguments
* `Pₛ`: surface pressure [Pa]
* `Pₜ`: top of atmopshere pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `absorber`: a `UnifiedAbsorber`, which is either a [`GroupedAbsorber`](@ref) or an [`AcceleratedAbsorber`](@ref)
"""
function topimbalance(X...; kwargs...)::Float64
    #compute monochromatic fluxes at top of atmosphere
    Fₜ⁻, Fₜ⁺ = topfluxes(X...; kwargs...)
    #get the wavenumber vector, knowing that last arg is UnifiedAbsorber
    ν = X[end].ν
    #integrate across wavenumbers
    trapz(ν, Fₜ⁻) - trapz(ν, Fₜ⁺)
end

#-------------------------------------------------------------------------------

