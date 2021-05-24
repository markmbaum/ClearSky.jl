#for organizing a group of absorbers
struct ParsedAbsorbers
    G::Tuple #all gases
    C::Tuple #CIA objects
    F::Tuple #functions σ(ν, T, P)
end

#splits a group of gas, cia, & functions objects into their own tuples
function ParsedAbsorbers(absorbers::Tuple)
    #can't be empty
    @assert length(absorbers) > 0 "no absorbers... use stefanboltzmann?"
    #check for dups
    @assert length(absorbers) == length(unique(absorbers)) "duplicate absorbers"
    #types of absorbers
    T = map(typeof, absorbers)
    #check for unexpected types
    for t ∈ T
        if !((t <: AbstractGas) | (t == CIATables) | (t <: Function))
            throw("absorbers must only be gasses (<: AbstractGas), CIA objects, or functions in the form σ(ν, T, P)")
        end
    end
    #all gases
    G = tuple(absorbers[findall(t->t<:AbstractGas, T)]...)
    #cia tables, pairing with the correct gases in the process
    C = tuple([CIA(x, G...) for x ∈ absorbers[findall(t->t==CIATables, T)]]...)
    #functions in the form σ(ν, T, P)
    F = tuple(absorbers[findall(t->!(t<:AbstractGas) & !(t==CIATables), T)]...)
    #construct the object
    ParsedAbsorbers(G, C, F)
end

#very important part of dealing with Vararg absorber inputs efficiently
# https://discourse.julialang.org/t/tuple-indexing-taking-time/58309/3
function σgas(gases::Tuple, i::Int, T, P)::Float64
    isempty(gases) && return 0.0
    first(gases)(i, T, P) + σgas(tail(gases), i, T, P)
end

#same as above but with a wavenumber argument instead of an index
function σgen(absorbers::Tuple, ν::Float64, T, P)::Float64
    isempty(absorbers) && return 0.0
    first(absorbers)(ν, T, P) + σgen(tail(absorbers), ν, T, P)
end

function (A::ParsedAbsorbers)(i, ν, T, P)::Float64
    σgas(A.G, i, T, P) + σgen(A.C, ν, T, P) + σgen(A.F, ν, T, P)
end

#-------------------------------------------------------------------------------
#more internal functions

#sets up gaussian quadrature for averaging multiple streams
function streamnodes(nstream::Int)
    #gauss-legendre quadrature points and weights in [-1,1]
    x, w = gauss(nstream)
    #map to cos(θ) ∈ [0,1]
    x = (x .+ 1)./2
    w ./= 2
    #get θ in case it's useful
    θ = acos.(x)
    #precompute 1/cos(θ), use m because μ is for gas molar mass
    m = 1 ./ x
    return m, w, θ
end

function checkpressures(G::Tuple, pressures...)
    #check pressure within domains of gases
    Pmin = minimum(map(g->g.Ω.Pmin, G))
    Pmax = maximum(map(g->g.Ω.Pmax, G))
    for P ∈ pressures
        @assert P >= Pmin "Pressure $P Pa too low, domain minimum is $Pmin"
        @assert P <= Pmax "Pressure $P Pa too low, domain minimum is $Pmax"
    end
end

function checkazimuth(θ)
    @assert 0 <= θ < π/2 "azimuth angle θ must be ∈ [0,π/2)"
end

#checks for identical wavenumber sampling across different gases
function setupwavenumbers(G::AbstractGas...)::Tuple{Vector{Float64},Int64}
    ν₁ = G[1].ν
    for g ∈ G
        @assert ν₁ == g.ν "gases must have identical wavenumber vectors"
    end
    return ν₁, length(ν₁)
end

#operations common to setting up high-level radiative functions
function setupintegration(Pₛ, Pₜ, absorbers)
    #split gas objects from cia objects
    A = ParsedAbsorbers(absorbers)
    #check pressures in order
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    #check pressures against AtmosphericDomains
    checkpressures(A.G, Pₛ, Pₜ)
    #insist on identical wavenumber vectors in each gas
    ν, nν = setupwavenumbers(A.G...)
    return A, ν, nν
end

#-------------------------------------------------------------------------------
#internal for integrating τ with proper parameters and coordinates

function dτdP(P::Float64, τ::Float64, param::Tuple)::Float64
    #unpack parameters
    A, i, ν, g, m, fT, fμ = param
    #temperature from given profile
    T = fT(P)
    #mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = A(i, ν, T, P)
    #compute dτ/dlnP, scaled by the angle m = 1/cos(θ)
    m*dτdP(σ, g, μ)
end

#-------------------------------------------------------------------------------
#internals for integrating Schwarzschild with proper parameters and coordinates

function dIdP(P::Float64, I::Float64, param::Tuple)::Float64
    #unpack parameters
    A, i, ν, g, m, fT, fμ = param
    #compute temperature from given profile
    T = fT(P)
    #compute mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = A(i, ν, T, P)
    #compute dI/dlnP, scaled by the angle m = 1/cos(θ)
    m*schwarzschild(I, ν, σ, g, μ, T)
end

function dIdω(ω::Float64, I::Float64, param::Tuple)::Float64
    P = ω2P(ω)
    P*dIdP(P, I, param)
end

function dIdι(ι::Float64, I::Float64, param::Tuple)::Float64
    P = ι2P(ι)
    P*dIdP(P, I, param)
end

#-------------------------------------------------------------------------------
#internal to integrate up the atmosphere for an outgoing stream

function outgoingstream(ω₁, #transformed pressure at start of integration
                        ω₂, #end of integration
                        A::ParsedAbsorbers,
                        i, #index for gas wavenumber and opaacity table
                        ν, #wavenumber [cm^-1]
                        g, #gravity [m/s^2]
                        m, #1/cos(θ), where θ is the stream angle
                        fT::Q, #temperature profile fT(P)
                        fμ::R, #mean molar mass as a function of T,P
                        tol #integrator tolerance
                        )::Float64 where {Q,R}
    #initial intensity
    I₀ = planck(ν, fT(ω2P(ω₁)))
    #pack parameters
    param = (A, i, ν, g, m, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coordinates and return
    radau(dIdω, I₀, ω₁, ω₂, param, atol=tol, rtol=tol)
end

function outgoingstream!(I, #vector output irradiance
                         ω, #vector output coordinates
                         ω₁, #transformed pressure at start of integration
                         ω₂, #end of integration
                         A::ParsedAbsorbers,
                         i, #index for gas wavenumber and opaacity table
                         ν, #wavenumber [cm^-1]
                         g, #gravity [m/s^2]
                         m, #1/cos(θ), where θ is the stream angle
                         fT::Q, #temperature profile fT(P)
                         fμ::R, #mean molar mass as a function of T,P
                         tol #integrator tolerance
                         )::Float64 where {Q,R}
    @assert length(I) == length(ω)
    #initial intensity
    I₀ = planck(ν, fT(ω2P(ω₁)))
    #pack parameters
    param = (A, i, ν, g, m, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coordinates and return
    radau!(I, ω, dIdω, I₀, ω₁, ω₂, param, atol=tol, rtol=tol)
end

#-------------------------------------------------------------------------------
#internal to integrate down the atmosphere for an incoming stream

function incomingstream(ι₁, #transformed pressure at start of integration
                        ι₂, #end of integration
                        A::ParsedAbsorbers,
                        i, #index for gas wavenumber and opaacity table
                        ν, #wavenumber [cm^-1]
                        g, #gravity [m/s^2]
                        m, #1/cos(θ), where θ is the stream angle
                        fT::Q, #temperature profile fT(P)
                        fμ::R, #mean molar mass as a function of T,P
                        fIₜ::S, #initial irradiance as a function of wavenumber fIₜ(ν)
                        tol #integrator tolerance
                        )::Float64 where {Q,R,S}
    #initial intensity
    I₀ = fIₜ(ν)
    #pack parameters
    param = (A, i, ν, g, m, fT, fμ)
    #integrate
    radau(dIdι, I₀, ι₁, ι₂, param, atol=tol, rtol=tol)
end

function incomingstream!(I, #vector output irradiance
                         ι, #output coordinates
                         ι₁, #transformed pressure at start of integration
                         ι₂, #end of integration
                         A::ParsedAbsorbers,
                         i, #index for gas wavenumber and opaacity table
                         ν, #wavenumber [cm^-1]
                         g, #gravity [m/s^2]
                         m, #1/cos(θ), where θ is the stream angle
                         fT::Q, #temperature profile fT(P)
                         fμ::R, #mean molar mass as a function of T,P
                         fIₜ::S, #initial irradiance as a function of wavenumber fIₜ(ν)
                         tol #integrator tolerance
                         )::Float64 where {Q,R,S}
    #initial intensity
    I₀ = fIₜ(ν)
    #pack parameters
    param = (A, i, ν, g, m, fT, fμ)
    #integrate
    radau!(I, ι, dIdι, I₀, ι₁, ι₂, param, atol=tol, rtol=tol)
end

#-------------------------------------------------------------------------------
export opticaldepth

"""
    opticaldepth(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute the optical depth (τ) between two pressure levels

# Arguments

* `P₁`: first pressure level [Pa]
* `P₂`: second pressure level [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `θ`: angle [rad] of path, must be ∈ [0,π/2), where 0 is straight up
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
                      tol::Float64=1e-5)::Vector{Float64} where {Q,R}
    #initialization
    if P₁ > P₂
        Pₛ, Pₜ = P₁, P₂
    else
        Pₜ, Pₛ = P₁, P₂
    end
    A, ν, nν = setupintegration(Pₛ, Pₜ, absorbers)
    checkazimuth(θ)
    m = 1/cos(θ)
    #integrate wavenumbers in parallel
    τ = zeros(Float64, nν)
    @threads for i = 1:nν
        #pack parameters
        param = (A, i, ν[i], g, m, fT, fμ)
        #integrate
        τ[i] = radau(dτdP, 0.0, Pₜ, Pₛ, param, atol=tol, rtol=tol)
    end
    return τ
end

#-------------------------------------------------------------------------------
#the basic transmittance method exp(-τ) is already exported

"""
    transmittance(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute the transmittance between two pressure levels

# Arguments

* `P₁`: first pressure level [Pa]
* `P₂`: second pressure level [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `θ`: angle [rad] of path, must be ∈ `[0,π/2)`, where 0 is straight up
* `absorbers`: at least one gas object and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

Returns a vector of transmittances across all wavenumbers stored in gas objects. The `tol` keyword argument adjusts integrator error tolerance.
"""
transmittance(X...; kwargs...) = transmittance.(opticaldepth(X...; kwargs...))

#-------------------------------------------------------------------------------
export outgoing

"""
    outgoing(Pₛ, Pₜ, g, fT, fμ, absorbers; nstream=3, tol=1e-5)

Compute outgoing radiative flux [W/m``^2``/cm``^{-1}``] at each wavenumber defined in the provided gas object(s). Total flux [W/m``^2``] can be evaluated with the [`trapz`](@ref) function.

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
                  fT::R,
                  fμ::Q,
                  absorbers...;
                  nstream::Int=3, #number of streams to use
                  tol::Float64=1e-5 #integrator tolerance
                  )::Vector{Float64} where {Q,R}
    #initialization
    A, ν, nν = setupintegration(Pₛ, Pₜ, absorbers)
    m, w, θ = streamnodes(nstream)
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    #integrate wavenumbers in parallel
    F = zeros(Float64, nν)
    @threads for i = 1:nν
        for j = 1:nstream
            #outgoing irradiance at Pₜ
            Iⱼ = outgoingstream(ω₁, ω₂, A, i, ν[i], g, m[j], fT, fμ, tol)
            #contribute to gauss quad
            #F[i] += w[j]*(2*π*sin(θ[j])*cos(θ[j])*Iⱼ) #why not this???
            F[i] += w[j]*π*Iⱼ
        end
    end
    #for hemisphere flux, ∫∫ I cos(θ) sin(θ) dθ dϕ = πI, assuming isotropy
    #F = π.*I
    return F
end

#-------------------------------------------------------------------------------
export incoming

"""
    incoming(Pₜ, Pₛ, g, fT, fμ, absorbers; nstream=3, tol=1e-5)

Compute incoming radiative flux [W/m``^2``/cm``^{-1}``] at each wavenumber defined in the provided gas object(s). Total flux [W/m``^2``] can be evaluated with the [`trapz`](@ref) function.

# Arguments
* `Pₜ`: top of atmopshere pressure [Pa]
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `absorbers`: at least one [gas object](gas_objects.md) and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

The keyword argument `nstream` specifies how many independent streams, or beam angles through the atmosphere, to integrate. The keyword argument `tol` is a numerical error tolerance passed to the [`radau`](https://github.com/wordsworthgroup/ScalarRadau.jl) integrator.
"""
function incoming(Pₜ::Real,
                  Pₛ::Real,
                  g::Real,
                  fT::Q,
                  fμ::R,
                  fIₜ::S,
                  θ::Real,
                  absorbers...;
                  tol::Float64=1e-5
                  )::Vector{Float64} where {Q,R,S}
    #initialization
    A, ν, nν = setupintegration(Pₛ, Pₜ, absorbers)
    ι₁, ι₂ = P2ι(Pₜ, Pₛ)
    checkazimuth(θ)
    m = 1/cos(θ)
    #integrate wavenumbers in parallel
    I = zeros(Float64, nν)
    @threads for i = 1:nν
        #incoming stream at only one angle
        I[i] = incomingstream(ι₁, ι₂, A, i, ν[i], g, m, fT, fμ, fIₜ, tol)
    end
    #single stream, flux is simply scaled by the angle
    F = cos(θ).*I
    return F
end

#-------------------------------------------------------------------------------
export heating

function heating(Pₛ,
                 Pₜ,
                 g,
                 fT::Q,
                 fμ::R,
                 fIₜ::S,
                 θ::Real,
                 outgoingabsorbers::Tuple,
                 incomingabsorbers::Tuple;
                 tol::Float64=1e-5) where {Q,R,S}
    #initialize outgoing stuff
    Aₒ, νₒ, nvₒ = setupintegration(Pₛ, Pₜ, outgoingabsorbers)
    mₒ, w, _ = streamnodes(nstream)
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    #initialize incoming stuff
    Aᵢ, νᵢ, nνᵢ = setupintegration(Pₛ, Pₜ, incomingabsorbers)
    ι₁, ι₂ = P2ι(Pₜ, Pₛ)
    checkazimuth(θ)
    mᵢ = 1/cos(θ)
    #
    #Iᵢ
    #Iₒ
    #@threads for i = 1:nν
end
