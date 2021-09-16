#-------------------------------------------------------------------------------
export opticaldepth

"""
    opticaldepth(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute monochromatic [optical depths](https://en.wikipedia.org/wiki/Optical_depth#Atmospheric_sciences) (``\\tau``) between two pressure levels

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
    A = unifyabsorbers(absorbers)
    checkpressures(A, P₁, P₂)
    ν, _ = A.ν, A.nν
    ω₁, ω₂ = P2ω(P₁), P2ω(P₂)
    checkazimuth(θ)
    m = 1/cos(θ)
    #spawn integrations in parallel, dynamic schedule
    tasks = [@spawn depth(dτdω, ω₁, ω₂, A, i, g, m, fT, fμ, tol) for i ∈ eachindex(ν)]
    #fetch the results
    [fetch(task) for task ∈ tasks]
end

#-------------------------------------------------------------------------------
#the basic transmittance method exp(-τ) is already exported

"""
    transmittance(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute monochromatic [transmittances](https://en.wikipedia.org/wiki/Transmittance). between two pressure levels

Accepts the same arguments as [`opticaldepth`](@ref) and returns a vector of transmittances across all wavenumbers stored in gas objects.
"""
transmittance(args...; kwargs...) = transmittance.(opticaldepth(args...; kwargs...))

#-------------------------------------------------------------------------------
export outgoing

"""
    outgoing(Pₛ, g, fT, fμ, absorbers; Pₜ=1.0, nstream=5, tol=1e-5)

Compute outgoing monochromatic radiative fluxes [W/m``^2``/cm``^{-1}``], line-by-line. Integrates the [`schwarzschild`](@ref) equation from `Pₛ` to `Pₜ` at each wavenumber in the provided gas object(s) using any number of streams/angles. Total flux [W/m``^2``] can be evaluated with the [`trapz`](@ref) function. This function does not include reflected stellar radiation, which is accounted for in [`topfluxes`](@ref) and [`topimbalance`](@ref)

# Arguments
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `absorbers`: at least one [gas object](gas_objects.md) and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

# Keywords
* `Pₜ`: top of atmopshere pressure [Pa]
* `nstream`: the number of atmospheric streams (radiation angles) to calculate and integrate
* `tol`: integrator error tolerance

The keyword argument `nstream` specifies how many independent streams, or beam angles through the atmosphere, to integrate. The keyword argument `tol` is a numerical error tolerance passed to the [`radau`](https://github.com/markmbaum/ScalarRadau.jl) integrator.
"""
function outgoing(Pₛ::Real,
                  g::Real,
                  fT::Q,
                  fμ::R,
                  absorbers...;
                  Pₜ::Real=1.0, #top of atmosphere pressure, 1 Pa by default
                  nstream::Int=5, #number of streams to use
                  tol::Real=1e-5 #integrator tolerance
                  ) where {Q,R}
    #initialization
    A = unifyabsorbers(absorbers)
    ν, nν = A.ν, A.nν
    checkpressures(A, Pₛ, Pₜ)
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    #surface temperature
    T₀ = fT(Pₛ)
    #integrate in parallel, dynamic schedule
    tasks = Vector{Task}(undef, nν)
    for (i,νᵢ) ∈ enumerate(ν)
        I₀ = planck(νᵢ, T₀)
        tasks[i] = @spawn streams(dIdω, I₀, ω₁, ω₂, A, i, g, fT, fμ, nstream, tol)
    end
    [fetch(task) for task ∈ tasks]
end

#-------------------------------------------------------------------------------
export netflux, heating

function netflux(P::AbstractVector{<:Real},
                 g::Real,
                 fT::Q,
                 fμ::R,
                 fS::S,
                 fα::U,
                 absorbers...;
                 nstream::Int=5,
                 θₛ::Float64=0.841,
                 tol::Float64=1e-5
                 )::Vector{Float64} where {Q,R,S,U}
    #setup
    A = unifyabsorbers(absorbers)
    #use ascending pressure coordinates
    idx = sortperm(P)
    P = P[idx]
    checkpressures(A, P[1], P[end])
    np = length(P)
    checkazimuth(θₛ)
    ν, nν = A.ν, A.nν
 
    #big blocks of flux
    M⁺ = zeros(np, nν) #monochromatic upward fluxes
    M⁻ = zeros(np, nν) #monochromatic downward fluxes
    #asynchronous, parallel integrations
    tasks = Vector{Task}(undef, A.nν)
    for i ∈ eachindex(A.ν)
        Mᵥ⁻ = view(M⁻,:,i)
        Mᵥ⁺ = view(M⁺,:,i)
        tasks[i] = @spawn monochromaticfluxes!(Mᵥ⁻, Mᵥ⁺, P, A.β[i], g, fT, fμ, fS, fα, nstream, θₛ, tol)
    end
    [fetch(task) for task ∈ tasks]

    Fₙ = zeros(np)
    @threads for i ∈ eachindex(P)
        #the upward flux order must be reversed to match
        Mᵥ⁺ = @view M⁺[np-i+1,:]
        #downward fluxes are already in the same order as P
        Mᵥ⁻ = @view M⁻[i,:]
        #take the difference of total flux at each level
        Fₙ[i] = trapz(ν, Mᵥ⁺) - trapz(ν, Mᵥ⁻)
    end

    return Fₙ[sortperm(idx)]
end

function heating(P::AbstractVector{<:Real},
                 g::Real,
                 fT::Q,
                 fμ::R,
                 fS::S,
                 fα::U,
                 fcₚ::V, #heat capacity cₚ(T,P) [J/kg/K]
                 absorbers...;
                 kwargs...
                 )::Vector{Float64} where {Q,R,S,U,V}
    #get pressures sorted in ascending order
    idx = sortperm(P)
    P = P[idx]
    #insert points for high-accuracy finite differences
    Φ, δ = insertdiff(P)
    #compute the net fluxes
    Fₙ = netflux(Φ, g, fT, fμ, fS, fα, absorbers...; kwargs...)
    #take dFₙ/dP
    ∂ = evaldiff(Fₙ, δ)
    #apply thermodynamical properties
    H = zeros(Float64, length(∂))
    for i ∈ eachindex(∂)
        #atmospheric temperature
        T = fT(P[i])
        #heat capacity
        cₚ = fcₚ(T,P[i])
        #heating rate [Kelvin/s]
        H[i] = ∂[i]*(g/cₚ)
    end
    return H[sortperm(idx)]
end