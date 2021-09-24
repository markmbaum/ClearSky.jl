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
                      ) where {Q,R}
    #initialization
    P₁, P₂ = max(P₁, P₂), min(P₁, P₂)
    𝔸, ν, nν = unifyabsorbers(absorbers)
    checkpressures(𝔸, P₁, P₂)
    ω₁, ω₂ = P2ω(P₁), P2ω(P₂)
    checkazimuth(θ)
    𝓂 = 1/cos(θ)
    #spawn integrations in parallel, dynamic schedule
    tasks = [@spawn depth(dτdω, ω₁, ω₂, 𝔸, i, g, 𝓂, fT, fμ, tol) for i ∈ 1:nν]
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
    𝔸, ν, nν = unifyabsorbers(absorbers)
    checkpressures(𝔸, Pₛ, Pₜ)
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    #surface temperature
    T₀ = fT(Pₛ)
    #integrate in parallel, dynamic schedule
    tasks = Vector{Task}(undef, nν)
    for (i,νᵢ) ∈ enumerate(ν)
        I₀ = planck(νᵢ, T₀)
        tasks[i] = @spawn streams(dIdω, I₀, ω₁, ω₂, 𝔸, i, g, fT, fμ, nstream, tol)
    end
    [fetch(task) for task ∈ tasks]
end

#-------------------------------------------------------------------------------
export monochromaticfluxes, fluxes, netfluxes, netfluxderivs, heating

function monochromaticfluxes(P::AbstractVector{<:Real},
                             g::Real,
                             fT::Q,
                             fμ::R,
                             fS::S,
                             fα::U,
                             absorbers...;
                             nstream::Int=5,
                             θₛ::Float64=0.841,
                             tol::Float64=1e-4) where {Q,R,S,U}
    #setup
    𝔸, ν, nν = unifyabsorbers(absorbers)
    #use ascending pressure coordinates
    idx = sortperm(P)
    P = P[idx]
    checkpressures(𝔸, P[end], P[1])
    np = length(P)
    checkazimuth(θₛ)
    #transformed coordinates
    ω, ι = P2ω.(P), P2ι.(P)
    reverse!(ω)

    #big blocks of flux
    M⁺ = zeros(eltype(P), np, nν) #monochromatic upward fluxes
    M⁻ = zeros(eltype(P), np, nν) #monochromatic downward fluxes
    #asynchronous, parallel integrations
    tasks = Vector{Task}(undef, nν)
    for i ∈ eachindex(ν)
        Mᵢ⁻ = @view M⁻[:,i]
        Mᵢ⁺ = @view M⁺[:,i]
        tasks[i] = @spawn monoflux!(Mᵢ⁻, Mᵢ⁺, P, ω, ι, 𝔸, i, g, fT, fμ, fS, fα, nstream, θₛ, tol)
    end
    [fetch(task) for task ∈ tasks]

    #make ordering consistent with input pressure ordering
    idx = sortperm(idx)
    M⁻ = M⁻[idx,:]
    M⁺ = M⁺[idx,:]
    #return the big blocks of monochromatic fluxes
    return M⁻, M⁺
end

function fluxes(P::AbstractVector{<:Real},
                g::Real,
                fT::Q,
                fμ::R,
                fS::S,
                fα::U,
                absorbers...;
                kwargs...) where {Q,R,S,U}
    #setup
    𝔸, ν, nν = unifyabsorbers(absorbers)
    #get monochromatic fluxes
    M⁻, M⁺ = monochromaticfluxes(P, g, fT, fμ, fS, fα, 𝔸; kwargs...)
    #integrate over wavenumber
    F⁻ = similar(M⁻, size(M⁻, 1))
    F⁺ = similar(M⁺, size(M⁺, 1))
    @threads for i ∈ eachindex(P)
        F⁻[i] = trapz(ν, view(M⁻,i,:))
        F⁺[i] = trapz(ν, view(M⁺,i,:))
    end
    return F⁻, F⁺
end

function netfluxes(P::AbstractVector{<:Real},
                   g::Real,
                   fT::Q,
                   fμ::R,
                   fS::S,
                   fα::U,
                   absorbers...;
                   kwargs...) where {Q,R,S,U}
    #wavenumber integrated fluxes [W/m^2] at each pressure level
    F⁻, F⁺ = fluxes(P, g, fT, fμ, fS, fα, absorbers...; kwargs...)
    #net flux
    [(F⁺[i] - F⁻[i]) for i ∈ eachindex(P)]
end

function netfluxderivs(P::AbstractVector{<:Real},
                       g::Real,
                       fT::Q,
                       fμ::R,
                       fS::S,
                       fα::U,
                       absorbers...;
                       kwargs...) where {Q,R,S,U}
    #ensure pressures are sorted in ascending order
    idx = sortperm(P)
    P = P[idx]
    #insert points for 2nd order finite differencing in ln(P) space
    Φ, δ = insertdiff(P)
    #compute the net fluxes
    Fₙ = netfluxes(Φ, g, fT, fμ, fS, fα, absorbers...; kwargs...)
    #evaluate ∂Fₙ/∂P
    ∂ = evaldiff(Fₙ, δ)
    #match input ordering and return
    ∂[sortperm(idx)]
end

function heating(P::AbstractVector{<:Real},
                 g::Real,
                 fT::Q,
                 fμ::R,
                 fS::S,
                 fα::U,
                 fcₚ::V, #heat capacity cₚ(T,P) [J/kg/K]
                 absorbers...;
                 kwargs...) where {Q,R,S,U,V}
    #evaluate derivative of net flux w/r/t pressure
    ∂ = netfluxderivs(P, g, fT, fμ, fS, fα, absorbers...; kwargs...)
    #compute heating rates
    H = similar(∂)
    for i ∈ eachindex(∂)
        #atmospheric temperature
        T = fT(P[i])
        #heat capacity
        cₚ = fcₚ(T,P[i])
        #heating rate [Kelvin/s]
        H[i] = ∂[i]*(g/cₚ)
    end
    return H
end