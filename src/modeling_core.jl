#-------------------------------------------------------------------------------
# general gets and checks

function pressurelimits(gases::Tuple)::NTuple{2,Float64}
    #largest minimum pressure in gas atmospheric domains
    Pmin = maximum(map(g->g.Ω.Pmin, gases))
    #smallest maximum pressure in gas atmospheric domains
    Pmax = minimum(map(g->g.Ω.Pmax, gases))
    return Pmin, Pmax
end

function checkpressures(gases::Tuple, Pₛ, Pₜ)::Nothing
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
    G = absorbers[findall(a -> typeof(a) <: Gas, absorbers)]
    @assert length(G) > 0 "no gas objects found"
    getwavenumbers(G...)
end

#checks for identical wavenumber sampling across different gases
function getwavenumbers(G::Gas...)::Vector{Float64}
    ν₁ = G[1].ν
    for g ∈ G
        @assert ν₁ == g.ν "gases must have identical wavenumber vectors"
    end
    return ν₁
end

#-------------------------------------------------------------------------------
# super for both consolidated absorber types

abstract type AbstractAbsorber end

export AbstractAbsorber, getσ

#-------------------------------------------------------------------------------
# specialized container for absorbing objects and functions

export UnifiedAbsorber

"""
    UnifiedAbsorber(absorbers...)

A struct for consolidating absorbers. Construct with any number of [gas objects](gas_objects.md), functions in the form `σ(ν, T, P)`, and [`CIATables`](@ref).
"""
struct UnifiedAbsorber{T,U,V} <: AbstractAbsorber
    #tuple of Gas objects
    gas::T
    #tuple of Gas objects
    cia::U
    #tuple of Gas objects
    fun::V
    #wavenumber vector [cm^-1], must be identical for all gases
    ν::Vector{Float64}
    #length of wavenumber vector
    nν::Int64
end

UnifiedAbsorber(absorbers...) = UnifiedAbsorber(absorbers)

#splits a group of gas, cia, & functions objects into their own tuples
function UnifiedAbsorber(absorbers::Tuple)
    #can't be empty
    @assert length(absorbers) > 0 "no absorbers... nothing to group"
    #check for dups
    @assert length(absorbers) == length(unique(absorbers)) "duplicate absorbers"
    #types of absorbers
    T = map(typeof, absorbers)
    #check for unexpected types
    for t ∈ T
        if !((t <: Gas) | (t == CIATables) | (t <: Function))
            throw("absorbers must only be gases (<: Gas), CIA objects, or functions in the form σ(ν, T, P)")
        end
    end
    #all gases
    gas = absorbers[findall(t -> t <: Gas, T)]
    #cia tables, pairing with the correct gases in the process
    cia = tuple([CIA(x, gas) for x ∈ absorbers[findall(t -> t == CIATables, T)]]...)
    #functions in the form σ(ν, T, P)
    fun = absorbers[findall(t -> !(t <: Gas) & !(t == CIATables), T)]
    #wavenumber vector, must be identical for all gases
    ν = getwavenumbers(gas...)
    nν = length(ν)
    #construct the UnifiedAbsorber
    UnifiedAbsorber(gas, cia, fun, ν, nν)
end

#https://discourse.julialang.org/t/tuple-indexing-taking-time/58309/18?u=markmbaum
#also see "applychain" here for a similar example: https://github.com/FluxML/Flux.jl/blob/dbb9f82ef8d4e196259ff1af56aeddc626159bf3/src/layers/basic.jl#L46
σchain(::Tuple{}, x, T, P) = zero(T)

σchain(A::Tuple, x, T, P) = first(A)(x, T, P) + σchain(tail(A), x, T, P)

function σchain(U::UnifiedAbsorber, i::Int, ν, T, P)
    σchain(U.gas, i, T, P) + σchain(U.cia, ν, T, P) + σchain(U.fun, ν, T, P)
end

(U::UnifiedAbsorber)(i::Int, T, P) = σchain(U, i, U.ν[i], T, P)

(U::UnifiedAbsorber)(T, P) = [U(i, T, P) for i ∈ eachindex(U.ν)]

getσ(U::UnifiedAbsorber, i::Int, T, P) = @inbounds σchain(U, i, U.ν[i], T, P)

checkpressures(U::UnifiedAbsorber, P...) = checkpressures(U.gas, P...)

#-------------------------------------------------------------------------------
# accelerated interpolation of cross-sections

export AcceleratedAbsorber, update!

"""
    AcceleratedAbsorber(P, T, G::UnifiedAbsorber)
    AcceleratedAbsorber(P, T, absorbers...)

An accelerated struct for getting cross-sections from groups of absorbers. Pressure and temperature coordinates must be provided. 
"""
struct AcceleratedAbsorber <: AbstractAbsorber
    #cross-section interpolators
    ϕ::Vector{LinearInterpolator{Float64, NoBoundaries}}
    #wavenumber vector [cm^-1], must be identical for all gases
    ν::Vector{Float64}
    #length of wavenumber vector
    nν::Int64
    #original pressures
    P::Vector{Float64}
    #reference to UnifiedAbsorber
    U::UnifiedAbsorber
end

function AcceleratedAbsorber(T, P, U::UnifiedAbsorber)
    #pull out wavenumber info
    ν, nν = U.ν, U.nν
    #flip vectors if pressure is not ascending (interpolators require this)
    idx = sortperm(P)
    P = P[idx]
    T = T[idx]
    #log pressure coordinates as usual
    logP = log.(P)
    #prepare interpolators
    ϕ = Vector{LinearInterpolator{Float64, NoBoundaries}}(undef, nν)
    for i ∈ eachindex(ν)
        #cross-sections from the UnifiedAbsorber
        σ = U.(i, T, P)
        #fix values so small they could wreak havoc (specifically, zeros)
        #set all the tinies to the super low threshold value
        @. σ[σ < TINY] = TINY
        #store an interpolating object
        ϕ[i] = LinearInterpolator(logP, log.(σ), NoBoundaries())
    end
    AcceleratedAbsorber(ϕ, ν, nν, P, U)
end

function AcceleratedAbsorber(T, P, absorbers...)
    AcceleratedAbsorber(T, P, UnifiedAbsorber(absorbers))
end

(A::AcceleratedAbsorber)(i::Int, P) = exp(A.ϕ[i](log(P)))

(A::AcceleratedAbsorber)(P) = [A(i, P) for i ∈ eachindex(A.ν)]

getσ(A::AcceleratedAbsorber, i, T, P) = @inbounds exp(A.ϕ[i](log(P)))

checkpressures(A::AcceleratedAbsorber, P...) = checkpressures(A.U, P...)

"""
    update!(A::AcceleratedAbsorber, T::Vector{Float64})

Update the cross-section interpolators underlying an `AcceleratedAbsorber` with a new set of temperatures. The new temperatures should correspond to the pressure levels used when originally constructing the `AcceleratedAbsorber`.  
"""
function update!(A::AcceleratedAbsorber, T::AbstractVector{<:Real})
    @assert length(T) == length(A.P)
    minlnσ = log(floatmin(Float64))
    for i ∈ eachindex(A.β)
        σ = values(A.β[i].ϕ)
        for j ∈ eachindex(σ)
            @inbounds σ[j] = max(log(A.U(i, T[j], A.P[j])), minlnσ)
        end
    end
end

#-------------------------------------------------------------------------------
#making sense of variable absorber inputs

function unifyabsorbers(absorbers::Tuple)::AbstractAbsorber
    length(absorbers) == 0 && error("no absorbers")
    if (length(absorbers) == 1) & (typeof(absorbers[1]) <: AbstractAbsorber)
        return absorbers[1]
    end
    UnifiedAbsorber(absorbers)
end

#-------------------------------------------------------------------------------
#function and cache for gaussian quadrature of multiple streams over the azimuth

const NODECACHE = Dict{Int64,NTuple{2,Vector{Float64}}}()

function newstreamnodes(n::Int64)::NTuple{2,Vector{Float64}}
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
    return m, W
end

function streamnodes(n::Int64)::NTuple{2,Vector{Float64}}
    #too few streams is likely problematic
    n < 4 && @warn "careful! using nstream < 4 is likely to be inaccurate!" maxlog=1
    #getting the gauss nodes is pretty fast but not trivial
    if !haskey(NODECACHE, n)
        #store these values
        NODECACHE[n] = newstreamnodes(n)
    end
    return NODECACHE[n]
end

#-------------------------------------------------------------------------------
# core differential equations with Tuples of parameters

function dτdP(P, τ, param::Tuple)
    #unpack parameters
    A, idx, g, m, fT, fμ = param
    #temperature from given profile
    T = fT(P)
    #mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = getσ(A, idx, T, P)
    #compute dτ/dlnP, scaled by the angle m = 1/cos(θ)
    m*dτdP(σ, g, μ) #no Planck emission
end

function dIdP(P, I, param::Tuple)
    #unpack parameters
    A, idx, g, m, fT, fμ = param
    #compute temperature from given profile
    T = fT(P)
    #compute mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = getσ(A, idx, T, P)
    #pull out wavenumber
    ν = @inbounds A.ν[idx]
    #compute dI/dlnP, scaled by the angle m = 1/cos(θ)
    m*schwarzschild(I, ν, σ, g, μ, T)
end

#-------------------------------------------------------------------------------
# wrappers for log pressure coordinates

function dτdι(ι, τ, param::Tuple)
    P = ι2P(ι)
    P*dτdP(P, τ, param)
end

function dτdω(ω, τ, param::Tuple)
    P = ω2P(ω)
    P*dτdP(P, τ, param)
end

function dIdω(ω, I, param::Tuple)
    P = ω2P(ω)
    P*dIdP(P, I, param)
end

function dIdι(ι, I, param::Tuple)
    P = ι2P(ι)
    P*dIdP(P, I, param)
end

#-------------------------------------------------------------------------------
# functions for optical depth paths

function depth(dτdx::Q,
               x₁::Real,
               x₂::Real,
               A::R,
               idx::Int,
               g::Real,
               m::Real, # 1/cos(θ)
               fT::S,
               fμ::U,
               tol::Float64
               ) where {Q,R<:AbstractAbsorber,S,U}
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
                A::R,
                idx::Int,
                g::Real, #gravity [m/s^2]
                m::Real, #1/cos(θ), where θ is the stream angle
                fT::S, #temperature profile fT(P)
                fμ::U, #mean molar mass μ(T,P)
                tol::Real #integrator error tolerance
                ) where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (A, idx, g, m, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords and return
    radau(dIdx, I₀, x₁, x₂, param, atol=tol, rtol=tol)
end

function streams(dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 x₁::Real, #initial pressure coordinate
                 x₂::Real, #final pressure coordinate
                 A::R,
                 idx::Int,
                 g::Real, #gravity [m/s^2]
                 fT::S, #temperature profile fT(P)
                 fμ::U, #mean molar mass μ(T,P)
                 nstream::Int,
                 tol::Real #integrator error tolerance
                 ) where {Q,R<:AbstractAbsorber,S,U}
    #setup gaussian quadrature nodes
    m, W = streamnodes(nstream)
    #solve schwarzschild w multiple streams, integrating over hemisphere
    M = zero(I₀)
    for i ∈ 1:nstream
        I = stream(dIdx, I₀, x₁, x₂, A, idx, g, m[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        M += W[i]*I #W = 2π*w*cos(θ)*sin(θ), precomputed
    end
    return M
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere with in-place storage

function stream!(I, #output/solution vector
                 x, #output/solution coordinates
                 dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 A::R,
                 idx::Int,
                 g::Real, #gravity [m/s^2]
                 m::Real, #1/cos(θ), where θ is the stream angle
                 fT::S, #temperature profile fT(P)
                 fμ::U, #mean molar mass μ(T,P)
                 tol::Real #integrator error tolerance
                 )::Nothing where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (A, idx, g, m, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords, in-place
    radau!(I, x, dIdx, I₀, x[1], x[end], param, atol=tol, rtol=tol)
    return nothing
end

function streams!(M, #output/solution vector
                  x, #output/solution coordinates
                  dIdx::Q, #version of schwarzschild equation
                  I₀::R, #initial irradiance
                  A::S,
                  idx::Int,
                  g::Real, #gravity [m/s^2]
                  fT::U, #temperature profile fT(P)
                  fμ::V, #mean molar mass μ(T,P)
                  nstream::Int,
                  tol::Real #integrator error tolerance
                  )::Nothing where {Q,R<:Real,S<:AbstractAbsorber,U,V}
    @assert length(M) == length(x)
    L = length(M)
    #setup gaussian quadrature nodes
    m, W = streamnodes(nstream)
    #temporary irradiance vector
    I = Vector{R}(undef,L) #allocating shouldn't usually be that costly, overall
    #solve schwarzschild w multiple streams, integrating over hemisphere
    for i ∈ 1:nstream
        stream!(I, x, dIdx, I₀, A, idx, g, m[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        for j ∈ eachindex(M)
            @inbounds M[j] += W[i]*I[j] #W = 2π*w*cos(θ)*sin(θ), precomputed
        end
    end
    return nothing
end

#-------------------------------------------------------------------------------
# core function for whole atmosphere upward and downward monochromatic fluxes

export monochromaticfluxes!

function monochromaticfluxes!(M⁻, #downward monochromatic fluxes [W/m^2/cm^-1]
                              M⁺, 
                              P, #pressure coordinates of output
                              A::Q,
                              idx::Int,
                              g::Real, #gravity [m/s^2]
                              fT::R, #temperature profile fT(P)
                              fμ::S, #mean molar mass μ(T,P)
                              fS::U, #incoming stellar radiation fS(ν) [W/m^2]
                              fα::V, #surface albedo fα(ν)
                              nstream::Int=5, #number of streams to integrate in both directions
                              θₛ::Real=0.841, #stellar radiation angle, corresponds to cos(θ) = 2/3
                              tol::Real=1e-4
                              ) where {Q<:AbstractAbsorber,R,S,U,V}
    @assert length(M⁻) == length(M⁺) == length(P)
    #surface pressure assuming ascending pressures
    Pₛ = P[end]
    #surface temperature
    Tₛ = fT(Pₛ)
    #wavenumber
    ν = A.ν[idx]
    #angle factor for incoming stellar radiation
    m = 1/cos(θₛ)
    #downward stellar irradiance at ν
    Iₜ⁻ = fS(ν)
    #transformed pressure coordinates
    ω = P2ω.(P)
    reverse!(ω)
    ι = P2ι.(P)
    #downward stellar irradiance throughout atmosphere
    stream!(M⁻, ι, dIdι, Iₜ⁻, A, idx, g, m, fT, fμ, tol)
    #convert to flux
    M⁻ .*= cos(θₛ)
    #add the atmospheric contribution to downward flux
    streams!(M⁻, ι, dIdι, zero(Iₜ⁻), A, idx, g, fT, fμ, nstream, tol)
    #some of the downward stellar flux is reflected
    Iₛ⁺ = M⁻[end]*fα(ν)/π #Lambertian
    #and the surface emits some radiation
    Iₛ⁺ += planck(ν, Tₛ)
    #upward radiation streams
    streams!(M⁺, ω, dIdω, Iₛ⁺, A, idx, g, fT, fμ, nstream, tol)
    #reverse the upward flux to match the coordinate ordering of P and ι
    reverse!(M⁺)

    nothing
end

export monochromaticfluxes

function monochromaticfluxes(P, #pressure coordinates of output
                             A::Q, #AbstractAbsorber
                             g::Real, #gravity [m/s^2]
                             fT::R, #temperature profile fT(P)
                             fμ::S, #mean molar mass μ(T,P)
                             fS::U, #incoming stellar radiation fS(ν) [W/m^2]
                             fα::V, #surface albedo fα(ν)
                             nstream::Int=5, #number of streams to integrate in both directions
                             θₛ::Real=0.841, #stellar radiation angle, corresponds to cos(θ) = 2/3
                             tol::Real=1e-4
                             ) where {Q<:AbstractAbsorber,R,S,U,V}
    #allocate space for results
    M⁻ = zeros(length(P), length(A.ν))
    M⁺ = zeros(length(P), length(A.ν))
    #asynchronous execution 'cause spectrally variable optical depth
    tasks = Vector{Task}(undef, A.nν)
    for i ∈ eachindex(A.ν)
        Mᵥ⁻ = view(M⁻,:,i)
        Mᵥ⁺ = view(M⁺,:,i)
        tasks[i] = @spawn monochromaticfluxes!(Mᵥ⁻, Mᵥ⁺, P, A, i, g, fT, fμ, fS, fα, nstream, θₛ, tol)
    end
    [fetch(task) for task ∈ tasks]

    return M⁻, M⁺
end