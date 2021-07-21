#-------------------------------------------------------------------------------
# general gets and checks

function pressurelimits(gases)::NTuple{2,Float64}
    #largest minimum pressure in gas atmospheric domains
    Pmin = maximum(map(g->g.Ω.Pmin, gases))
    #smallest maximum pressure in gas atmospheric domains
    Pmax = minimum(map(g->g.Ω.Pmax, gases))
    return Pmin, Pmax
end

function checkpressures(gases, pressures...)::Nothing
    #pressure bounds
    Pmin, Pmax = pressurelimits(gases)
    #demand all pressures within the range
    for P ∈ pressures
        @assert P >= Pmin "Pressure $P Pa too low, domain minimum is $Pmin"
        @assert P <= Pmax "Pressure $P Pa too low, domain minimum is $Pmax"
    end
    nothing
end

function checkazimuth(θ)::Nothing
    @assert 0 <= θ < π/2 "azimuth angle θ must be ∈ [0,π/2)"
    nothing
end

function getwavenumbers(absorbers::Tuple)::Vector{Float64}
    G = absorbers[findall(a -> typeof(a) <: AbstractGas, absorbers)]
    @assert length(G) > 0 "no gas objects found"
    getwavenumbers(G...)
end

#checks for identical wavenumber sampling across different gases
function getwavenumbers(G::AbstractGas...)::Vector{Float64}
    ν₁ = G[1].ν
    for g ∈ G
        @assert ν₁ == g.ν "gases must have identical wavenumber vectors"
    end
    return ν₁
end

#-------------------------------------------------------------------------------
# super for both consolidated absorber types

abstract type UnifiedAbsorber end

export UnifiedAbsorber, noabsorption, getσ

#-------------------------------------------------------------------------------
# specialized container for absorbing objects and functions

export GroupedAbsorber

"""
    GroupedAbsorber(absorbers...)

A struct for consolidating absorbers. Construct with any number of [gas objects](gas_objects.md), functions in the form `σ(ν, T, P)` and [`CIATables`](@ref).
"""
struct GroupedAbsorber{T,U,V} <: UnifiedAbsorber
    #tuple of gas objects
    gas::T
    #tuple of CIA objects
    cia::U
    #tuple of functions in the for f(ν, T, P)
    fun::V
    #wavenumber vector [cm^-1], must be identical for all gases
    ν::Vector{Float64}
    #length of wavenumber vector
    nν::Int64
    #flag indicating whether cia and fun are both empty
    gasonly::Bool
    #flags indicating where all gas objects are empty
    gasempty::Vector{Bool}
end

GroupedAbsorber(absorbers...) = GroupedAbsorber(absorbers)

#splits a group of gas, cia, & functions objects into their own tuples
function GroupedAbsorber(absorbers::Tuple)
    #can't be empty
    @assert length(absorbers) > 0 "no absorbers... no need for modeling?"
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
    G = absorbers[findall(t->t<:AbstractGas, T)]
    #cia tables, pairing with the correct gases in the process
    C = tuple([CIA(x, G...) for x ∈ absorbers[findall(t->t==CIATables, T)]]...)
    #functions in the form σ(ν, T, P)
    F = absorbers[findall(t->!(t<:AbstractGas) & !(t==CIATables), T)]
    #wavenumber vector, must be identical for all gases
    ν = getwavenumbers(G...)
    nν = length(ν)
    #flag indicating whether there are only gases present, no cias orfunctions
    gasonly = isempty(C) & isempty(F)
    #flags indicating whether all gases are empty at each wavenumber
    gasempty = zeros(Bool, nν)
    for i ∈ eachindex(ν)
        gasempty[i] = all(ntuple(j->G[j].Π[i].empty, length(G)))
    end
    return GroupedAbsorber(G, C, F, ν, nν, gasonly, gasempty)
end

#https://discourse.julialang.org/t/tuple-indexing-taking-time/58309/18?u=markmbaum
function σrecur(A::Q, x, T, P)::Float64 where {Q}
    isempty(A) && return 0.0
    first(A)(x, T, P) + σrecur(tail(A), x, T, P)
end

function getσ(ga::GroupedAbsorber, i::Int, T, P)::Float64
    ν = ga.ν[i]
    σ = σrecur(ga.gas, i, T, P) + σrecur(ga.cia, ν, T, P) + σrecur(ga.fun, ν, T, P)
    return σ
end

(ga::GroupedAbsorber)(i::Int, T, P)::Float64 = getσ(ga, i, T, P)

function (ga::GroupedAbsorber)(T::Real, P::Real)::Vector{Float64}
    [ga(i, T, P) for i ∈ eachindex(ga.ν)]
end

#check whether integration is pointless because there's no absorption
noabsorption(ga::GroupedAbsorber, i::Int)::Bool = ga.gasonly && ga.gasempty[i]

function noabsorption(ga::GroupedAbsorber)::Vector{Bool}
    [noabsorption(ga, i) for i ∈ eachindex(ga.ν)]
end

function checkpressures(ga::GroupedAbsorber, pressures...)
    checkpressures(ga.gas, pressures...)
end

#-------------------------------------------------------------------------------
# accelerated interpolation of cross-sections

export AcceleratedAbsorber, update!

"""
    AcceleratedAbsorber(ga, P, T)

An accelerated struct for getting cross-sections from groups of absorbers. Pressure and temperature coordinates must be provided. 
"""
struct AcceleratedAbsorber <: UnifiedAbsorber
    #cross-section interpolators
    ϕ::Vector{LinearInterpolator{Float64,WeakBoundaries}}
    #wavenumber vector [cm^-1], must be identical for all gases
    ν::Vector{Float64}
    #length of wavenumber vector
    nν::Int64
    #flag indicating whether there is no absorption
    empty::Vector{Bool}
    #original pressures
    P::Vector{Float64}
    #reference to GroupedAbsorber
    ga::GroupedAbsorber
end

function AcceleratedAbsorber(ga::GroupedAbsorber, P::Vector{Float64}, T::Vector{Float64})
    ν, nν = ga.ν, ga.nν
    logP = log.(P)
    ϕ = Vector{LinearInterpolator{Float64,WeakBoundaries}}(undef, nν)
    empty = zeros(Bool, nν)
    @threads for i ∈ eachindex(ν)
        empty[i] = noabsorption(ga, i)
        if !empty[i]
            σ = ga.(i, T, P)
            σ[σ .< 1e-200] .= 1e-200
            ϕ[i] = LinearInterpolator(logP, log.(σ), WeakBoundaries())
        end
    end
    AcceleratedAbsorber(ϕ, ν, nν, empty, P, ga)
end

#also a method specifically for interpolators, P vs log(σ)
getσ(aa::AcceleratedAbsorber, i, _, P)::Float64 = exp(aa.ϕ[i](log(P)))

(aa::AcceleratedAbsorber)(i::Int, P)::Float64 = getσ(aa, i, nothing, P)

function (aa::AcceleratedAbsorber)(P::Real)::Vector{Float64}
    [aa(i, P) for i ∈ eachindex(aa.ν)]
end

noabsorption(aa::AcceleratedAbsorber, i::Int)::Bool = aa.empty[i]

function noabsorption(ga::AcceleratedAbsorber)::Vector{Bool}
    [noabsorption(aa, i) for i ∈ eachindex(ga.ν)]
end

function checkpressures(aa::AcceleratedAbsorber, pressures...)
    checkpressures(aa.ga.gas, pressures...)
end

function update!(aa::AcceleratedAbsorber, T::Vector{Float64})
    @assert length(T) == length(aa.P)
    @threads for i ∈ eachindex(aa.ϕ)
        if !aa.empty[i]
            for j ∈ eachindex(aa.P)
                aa.ϕ[i].r.y[j] = max(log(aa.ga(i, T[j], aa.P[j])), 1e-200)
            end
        end
    end
end

#-------------------------------------------------------------------------------
#function and cache for gaussian quadrature of multiple streams over the azimuth

const NODECACHE = Dict{Int64,NTuple{2,Vector{Float64}}}()

function streamnodes(n::Int)::NTuple{2,Vector{Float64}}
    #pedantic with the key type
    n = convert(Int64, n)
    #getting the gauss nodes is fast but not trivial
    if !haskey(NODECACHE, n)
        if n < 4
            @warn "careful! using nstream < 4 is likely to be inaccurate!" maxlog=1
        end
        #gauss-legendre quadrature points and weights in [-1,1]
        x, w = gauss(n)
        #map angles and weights to θ ∈ [0,π/2]
        θ = @. (π/2)*(x + 1)/2
        w .*= (π/2)/2
        #sines and cosines of θ
        c = cos.(θ)
        s = sin.(θ)
        #precompute 2π*cos(θ)*sin(θ)*wzxc
        W = @. 2π*w*c*s
        #precompute 1/cos(θ) using "m" because μ is for gas molar masses
        m = 1 ./ c
        #store these values
        NODECACHE[n] = (m, W)
    end
    return NODECACHE[n]
end

#-------------------------------------------------------------------------------
# functions for setting up absorbers, coordinates, etc. for integration

#operations common to setting up high-level radiative functions
function setupintegration(Pₛ, Pₜ, absorbers)
    #split gas objects from cia objects
    ga = GroupedAbsorber(absorbers)
    #check pressures in order
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    #check pressures against AtmosphericDomains
    checkpressures(ga.gas, Pₛ, Pₜ)
    return ga
end

#-------------------------------------------------------------------------------
# core differential equations with Tuples of parameters

function dτdP(P::Float64, τ::Float64, param::Tuple)::Float64
    #unpack parameters
    A, i, g, m, fT, fμ = param
    #temperature from given profile
    T = fT(P)
    #mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = getσ(A, i, T, P)
    #compute dτ/dlnP, scaled by the angle m = 1/cos(θ)
    m*dτdP(σ, g, μ) #no Planck emission
end

function dIdP(P::Float64, I::Float64, param::Tuple)::Float64
    #unpack parameters
    A, i, g, m, fT, fμ = param
    #compute temperature from given profile
    T = fT(P)
    #compute mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = getσ(A, i, T, P)
    #compute dI/dlnP, scaled by the angle m = 1/cos(θ)
    m*schwarzschild(I, A.ν[i], σ, g, μ, T)
end

#-------------------------------------------------------------------------------
# wrappers for log pressure coordinates

function dτdι(ι::Float64, τ::Float64, param::Tuple)::Float64
    P = ι2P(ι)
    P*dτdP(P, τ, param)
end

function dτdω(ω::Float64, τ::Float64, param::Tuple)::Float64
    P = ω2P(ω)
    P*dτdP(P, τ, param)
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
# functions for optical depth paths

function depth(dτdx::Q,
               x₁::Real,
               x₂::Real,
               A::UnifiedAbsorber,
               idx::Int,
               g::Real,
               m::Real, # 1/cos(θ)
               fT::R,
               fμ::S,
               tol::Float64
               )::Float64 where {Q,R,S}
    #if zero absorption, don't integrate
    noabsorption(A, idx) && return 0.0
    #pack parameters
    param = (A, idx, g, m, fT, fμ)
    #integrate with the ODE solver (appears to be faster than quadrature)
    radau(dτdx, 0.0, x₁, x₂, param, atol=tol, rtol=tol)
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere, no storage

function stream(dIdx::Q, #version of schwarzschild equation
                I₀::Real, #initial irradiance
                x₁::Real, #initial pressure coordinate
                x₂::Real, #final pressure coordinate
                A::UnifiedAbsorber,
                idx::Int, #index for wavenumber and opacity table
                g::Real, #gravity [m/s^2]
                m::Real, #1/cos(θ), where θ is the stream angle
                fT::R, #temperature profile fT(P)
                fμ::S, #mean molar mass μ(T,P)
                tol::Float64 #integrator error tolerance
                )::Float64 where {Q,R,S}
    #if zero absorption, don't integrate
    noabsorption(A, idx) && return I₀
    #pack parameters
    param = (A, idx, g, m, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords and return
    radau(dIdx, I₀, x₁, x₂, param, atol=tol, rtol=tol)
end

function streams(dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 x₁::Real, #initial pressure coordinate
                 x₂::Real, #final pressure coordinate
                 A::UnifiedAbsorber,
                 idx::Int, #index for wavenumber and opacity table
                 g::Real, #gravity [m/s^2]
                 nstream::Int,
                 fT::R, #temperature profile fT(P)
                 fμ::S, #mean molar mass μ(T,P)
                 tol::Float64 #integrator error tolerance
                 )::Float64 where {Q,R,S}
    #setup gaussian quadrature nodes
    m, W = streamnodes(nstream)
    #solve schwarzschild w multiple streams, integrating over hemisphere
    F = 0.0
    for i ∈ 1:nstream
        I = stream(dIdx, I₀, x₁, x₂, A, idx, g, m[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        F += W[i]*I #W = 2π*w*cos(θ)*sin(θ), precomputed
    end
    return F
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere with in-place storage

function stream!(dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 I::AbstractVector{Float64}, #output/solution vector
                 x::Vector{Float64}, #output/solution coordinates
                 x₁::Real, #initial pressure coordinate
                 x₂::Real, #final pressure coordinate
                 A::UnifiedAbsorber,
                 idx::Int, #index for wavenumber and interpolator
                 g::Real, #gravity [m/s^2]
                 m::Real, #1/cos(θ), where θ is the stream angle
                 fT::R, #temperature profile fT(P)
                 fμ::S, #mean molar mass μ(T,P)
                 tol::Float64 #integrator error tolerance
                 )::Nothing where {Q,R,S}
    #if zero absorption, don't integrate
    if noabsorption(A, idx)
        I .= I₀
        return nothing
    end
    #pack parameters
    param = (A, idx, g, m, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords, in-place
    radau!(I, x, dIdx, I₀, x₁, x₂, param, atol=tol, rtol=tol)
    return nothing
end

function streams!(dIdx::Q, #version of schwarzschild equation
                  I₀::Real, #initial irradiance
                  I::AbstractVector{Float64}, #temporary irradiance vector
                  F::AbstractVector{Float64}, #output/solution vector
                  x::Vector{Float64}, #output/solution coordinates
                  x₁::Real, #initial pressure coordinate
                  x₂::Real, #final pressure coordinate
                  A::UnifiedAbsorber,
                  idx::Int, #index for wavenumber and opacity table
                  g::Real, #gravity [m/s^2]
                  nstream::Int,
                  fT::R, #temperature profile fT(P)
                  fμ::S, #mean molar mass μ(T,P)
                  tol::Float64 #integrator error tolerance
                  )::Nothing where {Q,R,S}
    @assert length(I) == length(F) == length(x)
    #setup gaussian quadrature nodes
    m, W = streamnodes(nstream)
    #wipe any pre-existing values
    F .= 0.0
    #solve schwarzschild w multiple streams, integrating over hemisphere
    for i ∈ 1:nstream
        stream!(dIdx, I₀, I, x, x₁, x₂, A, idx, g, m[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        for j ∈ eachindex(F)
            F[j] += W[i]*I[j] #W = 2π*w*cos(θ)*sin(θ), precomputed
        end
    end
    return nothing
end


