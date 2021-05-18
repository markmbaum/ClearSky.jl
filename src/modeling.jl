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
    σ::Float64 = 0.0
    σ += σgas(A.G, i, T, P)
    σ += σgen(A.C, ν, T, P)
    σ += σgen(A.F, ν, T, P)
    return σ
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
#internals to integrate up the atmosphere for an outgoing stream

function outgoingstream(ω₁, #transformed pressure at start of integration
                        ω₂, #end of integration
                        A, #ParsedAbsorbers
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

#-------------------------------------------------------------------------------
#internals to integrate down the atmosphere for an incoming stream

function incomingstream(ι₁, #transformed pressure at start of integration
                        ι₂, #end of integration
                        A, #ParsedAbsorbers
                        i, #index for gas wavenumber and opaacity table
                        ν, #wavenumber [cm^-1]
                        g, #gravity [m/s^2]
                        m, #1/cos(θ), where θ is the stream angle
                        fT::Q, #temperature profile fT(P)
                        fμ::R, #mean molar mass as a function of T,P
                        fI₀::S, #initial irradiance as a function of wavenumber fI₀(ν)
                        tol #integrator tolerance
                        )::Float64 where {Q,R,S}
    #initial intensity
    I₀ = fI₀(ν)
    #pack parameters
    param = (A, i, ν, g, m, fT, fμ)
    #integrate
    radau(dIdι, I₀, ι₁, ι₂, param, atol=tol, rtol=tol)
end

#-------------------------------------------------------------------------------
export opticaldepth

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

function transmittance(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)
    transmittance.(opticaldepth(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=tol))
end

#-------------------------------------------------------------------------------
export outgoing

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
    I = zeros(Float64, nν)
    @threads for i = 1:nν
        for j = 1:nstream
            Iⱼ = outgoingstream(ω₁, ω₂, A, i, ν[i], g, m[j], fT, fμ, tol)
            I[i] += π*w[j]*Iⱼ #contribute to gauss quad
        end
    end
    return I
end

#-------------------------------------------------------------------------------
export incoming

function incoming(Pₜ::Real,
                  Pₛ::Real,
                  g::Real,
                  fT::Q,
                  fμ::R,
                  fI₀::S,
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
        I[i] = incomingstream(ι₁, ι₂, A, i, ν[i], g, m, fT, fμ, fI₀, tol)
    end
    return I
end
