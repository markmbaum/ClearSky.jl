#-------------------------------------------------------------------------------
# general gets and checks

function pressurelimits(gases)::NTuple{2,Float64}
    #largest minimum pressure in gas atmospheric domains
    Pmin = maximum(map(g->g.Ω.Pmin, gases))
    #smallest maximum pressure in gas atmospheric domains
    Pmax = minimum(map(g->g.Ω.Pmax, gases))
    return Pmin, Pmax
end

function checkpressures(gases, Pₛ, Pₜ)::Nothing
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    #pressure bounds
    Pmin, Pmax = pressurelimits(gases)
    #demand all pressures within the range
    for P ∈ (Pₛ, Pₜ)
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
    gas = absorbers[findall(t->t<:AbstractGas, T)]
    #cia tables, pairing with the correct gases in the process
    cia = tuple([CIA(x, gas) for x ∈ absorbers[findall(t->t==CIATables, T)]]...)
    #functions in the form σ(ν, T, P)
    fun = absorbers[findall(t->!(t<:AbstractGas) & !(t==CIATables), T)]
    #wavenumber vector, must be identical for all gases
    ν = getwavenumbers(gas...)
    nν = length(ν)
    #flag indicating whether there are only gases present, no cias orfunctions
    gasonly = isempty(cia) & isempty(fun)
    #flags indicating whether all gases are empty at each wavenumber
    gasempty = zeros(Bool, nν)
    for i ∈ eachindex(ν)
        gasempty[i] = all(ntuple(j->gas[j].Π[i].empty, length(gas)))
    end
    return GroupedAbsorber(gas, cia, fun, ν, nν, gasonly, gasempty)
end

#https://discourse.julialang.org/t/tuple-indexing-taking-time/58309/18?u=markmbaum
function σrecur(A::Q, x, T, P)::Float64 where {Q}
    isempty(A) && return 0.0
    first(A)(x, T, P) + σrecur(tail(A), x, T, P)
end

function getσ(G::GroupedAbsorber, i::Int, T, P)::Float64
    ν = G.ν[i]
    σ = σrecur(G.gas, i, T, P) + σrecur(G.cia, ν, T, P) + σrecur(G.fun, ν, T, P)
    return σ
end

(G::GroupedAbsorber)(i::Int, T, P)::Float64 = getσ(G, i, T, P)

function (G::GroupedAbsorber)(T::Real, P::Real)::Vector{Float64}
    [G(i, T, P) for i ∈ eachindex(G.ν)]
end

#check whether integration is pointless because there's no absorption
noabsorption(G::GroupedAbsorber, i::Int)::Bool = G.gasonly && G.gasempty[i]

function noabsorption(G::GroupedAbsorber)::Vector{Bool}
    [noabsorption(G, i) for i ∈ eachindex(G.ν)]
end

function checkpressures(G::GroupedAbsorber, Pₛ, Pₜ)
    checkpressures(G.gas, Pₛ, Pₜ)
end

#-------------------------------------------------------------------------------
# accelerated interpolation of cross-sections

export AcceleratedAbsorber, update!

"""
    AcceleratedAbsorber(P, T, G::GroupedAbsorber)
    AcceleratedAbsorber(P, T, absorbers...)

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
    G::GroupedAbsorber
end

function AcceleratedAbsorber(P::Vector{Float64},T::Vector{Float64}, G::GroupedAbsorber)
    #pull out wavenumber info
    ν, nν = G.ν, G.nν
    #flip pressures if they're not ascending (for interpolators)
    if P[1] > P[end]
        P = reverse(P)
    end
    #log pressure coordinates as usual
    logP = log.(P)
    ϕ = Vector{LinearInterpolator{Float64,WeakBoundaries}}(undef, nν)
    empty = zeros(Bool, nν)
    @threads for i ∈ eachindex(ν)
        empty[i] = noabsorption(G, i)
        if !empty[i]
            σ = G.(i, T, P)
            #σ[σ .< 1e-200] .= 1e-200
            ϕ[i] = LinearInterpolator(logP, log.(σ), WeakBoundaries())
        end
    end
    AcceleratedAbsorber(ϕ, ν, nν, empty, P, G)
end

function AcceleratedAbsorber(P::Vector{Float64}, T::Vector{Float64}, absorbers...)
    AcceleratedAbsorber(P, T, GroupedAbsorber(absorbers))
end

#also a method specifically for interpolators, P vs log(σ)
function getσ(A::AcceleratedAbsorber, i, _, P)::Float64
    A.empty[i] && return 0.0
    exp(A.ϕ[i](log(P)))
end

(A::AcceleratedAbsorber)(i::Int, P)::Float64 = getσ(A, i, nothing, P)

function (A::AcceleratedAbsorber)(P::Real)::Vector{Float64}
    [A(i, P) for i ∈ eachindex(A.ν)]
end

noabsorption(A::AcceleratedAbsorber, i::Int)::Bool = A.empty[i]

function noabsorption(G::AcceleratedAbsorber)::Vector{Bool}
    [noabsorption(A, i) for i ∈ eachindex(G.ν)]
end

function checkpressures(A::AcceleratedAbsorber, Pₛ, Pₜ)
    checkpressures(A.G.gas, Pₛ, Pₜ)
end

"""
    update!(A::AcceleratedAbsorber, T::Vector{Float64})

Update the cross-section interpolators underlyig an `AcceleratedAbsorber` with a new set of temperatures. The new temperatures should correspond to the pressure levels used when originally constructing the `AcceleratedAbsorber`.  
"""
function update!(A::AcceleratedAbsorber, T::Vector{Float64})
    @assert length(T) == length(A.P)
    @threads for i ∈ eachindex(A.ϕ)
        if !A.empty[i]
            for j ∈ eachindex(A.P)
                A.ϕ[i].r.y[j] = max(log(A.G(i, T[j], A.P[j])), 1e-200)
            end
        end
    end
end

#-------------------------------------------------------------------------------
#making sense of variable absorber inputs

function unifyabsorbers(absorbers::Tuple)::UnifiedAbsorber
    length(absorbers) == 0 && error("no absorbers")
    if (length(absorbers) == 1) && (typeof(absorbers[1]) <: UnifiedAbsorber)
        return absorbers[1]
    end
    GroupedAbsorber(absorbers)
end

#-------------------------------------------------------------------------------
#function and cache for gaussian quadrature of multiple streams over the azimuth

const NODECACHE = Dict{Int64,NTuple{2,Tuple}}()

function _streamnodes(n::Int64)::NTuple{2,Tuple}
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
    return tuple(m...), tuple(W...)
end

function streamnodes(n::Int64)::NTuple{2,Tuple}
    #too few streams is likely problematic
    if n < 4
        @warn "careful! using nstream < 4 is likely to be inaccurate!" maxlog=1
    end
    #getting the gauss nodes is pretty fast but not trivial
    if !haskey(NODECACHE, n)
        #store these values
        NODECACHE[n] = _streamnodes(n)
    end
    return NODECACHE[n]
end

#load some stream nodes in serial
for n ∈ 1:16
    NODECACHE[n] = _streamnodes(n)
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
    #F .= 0.0
    #solve schwarzschild w multiple streams, integrating over hemisphere
    for i ∈ 1:nstream
        stream!(dIdx, I₀, I, x, x₁, x₂, A, idx, g, m[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        for j ∈ eachindex(F)
            @inbounds F[j] += W[i]*I[j] #W = 2π*w*cos(θ)*sin(θ), precomputed
        end
    end
    return nothing
end
