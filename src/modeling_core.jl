#-------------------------------------------------------------------------------
# general gets and checks

function pressurelimits(gases::Tuple)::NTuple{2,Float64}
    g = gases[findall(g -> typeof(g) <: Gas, gases)]
    isempty(g) && return(0.0, Inf)
    #largest minimum pressure in gas atmospheric domains
    Pmin = maximum(map(g->g.Ω.Pmin, g))
    #smallest maximum pressure in gas atmospheric domains
    Pmax = minimum(map(g->g.Ω.Pmax, g))
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
    G = absorbers[findall(a -> typeof(a) <: AbstractGas, absorbers)]
    @assert length(G) > 0 "no gas objects found"
    getwavenumbers(G...)
end

#checks for identical wavenumber sampling across different gases
function getwavenumbers(G::AbstractGas...)::Vector{Float64}
    @assert all(g -> g.ν == G[1].ν, G) "gases must have identical wavenumber vectors"
    return G[1].ν
end

#-------------------------------------------------------------------------------
# super for both consolidated absorber types

abstract type AbstractAbsorber end

export AbstractAbsorber

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
        if !((t <: AbstractGas) | (t == CIATables) | (t <: Function))
            throw("absorbers must only be gases (<: Gas), CIA objects, or functions in the form σ(ν, T, P)")
        end
    end
    #all gases
    gas = absorbers[findall(t -> t <: AbstractGas, T)]
    isempty(gas) && error("must have at least one Gas object, which specifies wavenumber samples")
    #cia tables, pairing with the correct gases in the process
    cia = tuple([CIA(x, gas) for x ∈ absorbers[findall(t -> t == CIATables, T)]]...)
    #functions in the form σ(ν, T, P)
    fun = absorbers[findall(t -> !(t <: AbstractGas) & !(t == CIATables), T)]
    #wavenumber vector, must be identical for all gases
    ν = getwavenumbers(gas...)
    nν = length(ν)
    #construct the UnifiedAbsorber
    UnifiedAbsorber(gas, cia, fun, ν, nν)
end

#https://discourse.julialang.org/t/tuple-indexing-taking-time/58309/18?u=markmbaum
#also see "applychain" here for a similar example: https://github.com/FluxML/Flux.jl/blob/dbb9f82ef8d4e196259ff1af56aeddc626159bf3/src/layers/basic.jl#L46
σchain(::Tuple{}, x, T, P) = 0

σchain(A::Tuple, x, T, P) = first(A)(x, T, P) + σchain(tail(A), x, T, P)

function σchain(U::UnifiedAbsorber, i::Int, ν, T, P)
    ( σchain(U.gas, i, T, P)
    + σchain(U.cia, ν, T, P)
    + σchain(U.fun, ν, T, P))
end

#internal
𝝈(U::UnifiedAbsorber, i::Int, T, P) = @inbounds σchain(U, i, U.ν[i], T, P)

(U::UnifiedAbsorber)(i::Int, T, P) = σchain(U, i, U.ν[i], T, P)

(U::UnifiedAbsorber)(T, P) = [U(i, T, P) for i ∈ eachindex(U.ν)]

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
    lnP = log.(P)
    #initalize an AcceleratedAbsorber with empty interpolators
    ϕ = Vector{LinearInterpolator{Float64, NoBoundaries}}(undef, nν)
    for i ∈ eachindex(ν)
        ϕ[i] = LinearInterpolator(lnP, similar(lnP), NoBoundaries())
    end
    A = AcceleratedAbsorber(ϕ, ν, nν, P, U)
    #then update the cross-sections in-place
    update!(A, T)
    #and return the updated AcceleratedAbsorber
    return A
end

"""
    update!(A::AcceleratedAbsorber, T)

Update the cross-section interpolators underlying an `AcceleratedAbsorber` with a new set of temperatures. The new temperatures should correspond to the pressure levels used when originally constructing the `AcceleratedAbsorber`.  
"""
function update!(A::AcceleratedAbsorber, T)::Nothing
    #check lengths
    @assert length(T) == length(A.P)
    #log of smallest float
    lntiny = log(TINY)
    #update each interpolators
    for (i,ϕ) ∈ enumerate(A.ϕ)
        #update each value
        for j ∈ eachindex(T)
            #retrieve cross-section from UnifiedAbsorber
            lnσⱼ = log(𝝈(A.U, i, T[j], A.P[j]))
            #set the new value
            ϕ[j] = (lnσⱼ < lntiny) ? lntiny : lnσⱼ
        end
    end
    nothing
end

function AcceleratedAbsorber(T, P, absorbers...)
    AcceleratedAbsorber(T, P, UnifiedAbsorber(absorbers))
end

#internal
𝝈(A::AcceleratedAbsorber, i::Int, ::Any, P) = @inbounds exp(A.ϕ[i](log(P)))

(A::AcceleratedAbsorber)(i::Int, P) = exp(A.ϕ[i](log(P)))

(A::AcceleratedAbsorber)(P) = [𝛔(A, i, nothing, P) for i ∈ eachindex(A.ν)]

checkpressures(A::AcceleratedAbsorber, P...) = checkpressures(A.U, P...)

#-------------------------------------------------------------------------------
#making sense of variable absorber inputs

function unifyabsorbers(absorbers::Tuple)
    L = length(absorbers)
    L == 0 && error("no absorbers")
    T = typeof(absorbers[1])
    if (L == 1) & (T <: AbstractAbsorber)
        A = absorbers[1]
    else
        A = UnifiedAbsorber(absorbers)
    end
    return A, A.ν, A.nν
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
    #pre-scaled weights and flipped cosines
    𝒲 = @. 2π*w*cos(θ)*sin(θ)
    𝓂 = @. 1/cos(θ)
    return 𝓂, 𝒲
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

function dτdP(P, ::Any, param::Tuple)
    #unpack parameters
    A, idx, g, 𝓂, fT, fμ = param
    #temperature from given profile
    T = fT(P)
    #mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = 𝝈(A, idx, T, P)
    #compute dτ/dlnP, scaled by the angle m = 1/cos(θ)
    𝓂*dτdP(σ, g, μ) #no Planck emission
end

function dIdP(P, I, param::Tuple)
    #unpack parameters
    A, idx, g, 𝓂, fT, fμ = param
    #compute temperature from given profile
    T = fT(P)
    #compute mean molar mass
    μ = fμ(T, P)
    #sum of all cross-sections
    σ = 𝝈(A, idx, T, P)
    #pull out wavenumber
    ν = @inbounds A.ν[idx]
    #compute dI/dlnP, scaled by the angle m = 1/cos(θ)
    𝓂*schwarzschild(I, ν, σ, g, μ, T)
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
               𝔸::R,
               idx::Int,
               g::Real,
               𝓂::Real, # 1/cos(θ)
               fT::S,
               fμ::U,
               tol::Float64
               ) where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (𝔸, idx, g, 𝓂, fT, fμ)
    #integrate with the ODE solver (appears to be faster than quadrature)
    radau(dτdx, 0.0, x₁, x₂, param, atol=tol, rtol=tol)
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere, no storage

function stream(dIdx::Q, #version of schwarzschild equation
                I₀::Real, #initial irradiance
                x₁::Real, #initial pressure coordinate
                x₂::Real, #final pressure coordinate
                𝔸::R,
                idx::Int,
                g::Real, #gravity [m/s^2]
                𝓂::Real, #1/cos(θ), where θ is the stream angle
                fT::S, #temperature profile fT(P)
                fμ::U, #mean molar mass μ(T,P)
                tol::Real #integrator error tolerance
                ) where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (𝔸, idx, g, 𝓂, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords and return
    radau(dIdx, I₀, x₁, x₂, param, atol=tol, rtol=tol)
end

function streams(dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 x₁::Real, #initial pressure coordinate
                 x₂::Real, #final pressure coordinate
                 𝔸::R,
                 idx::Int, #index of wavenumber
                 g::Real, #gravity [m/s^2]
                 fT::S, #temperature profile fT(P)
                 fμ::U, #mean molar mass μ(T,P)
                 nstream::Int,
                 tol::Real #integrator error tolerance
                 ) where {Q,R<:AbstractAbsorber,S,U}
    #setup gaussian quadrature nodes
    𝓂, 𝒲 = streamnodes(nstream)
    #solve schwarzschild w multiple streams, integrating over hemisphere
    M = zero(I₀)
    for i ∈ 1:nstream
        I = stream(dIdx, I₀, x₁, x₂, 𝔸, idx, g, 𝓂[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        M += 𝒲[i]*I #W = 2π*w*cos(θ)*sin(θ), precomputed
    end
    return M
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere with in-place storage

function stream!(I, #output/solution vector
                 x, #output/solution coordinates
                 dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 𝔸::R,
                 idx::Int,
                 g::Real, #gravity [m/s^2]
                 𝓂::Real, #1/cos(θ), where θ is the stream angle
                 fT::S, #temperature profile fT(P)
                 fμ::U, #mean molar mass μ(T,P)
                 tol::Real #integrator error tolerance
                 )::Nothing where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (𝔸, idx, g, 𝓂, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords, in-place
    radau!(I, x, dIdx, I₀, x[1], x[end], param, atol=tol, rtol=tol)
    return nothing
end

function streams!(M, #output/solution vector
                  x, #output/solution coordinates
                  dIdx::Q, #version of schwarzschild equation
                  I₀::R, #initial irradiance
                  𝔸::S,
                  idx::Int,
                  g::Real, #gravity [m/s^2]
                  fT::U, #temperature profile fT(P)
                  fμ::V, #mean molar mass μ(T,P)
                  nstream::Int,
                  tol::Real #integrator error tolerance
                  )::Nothing where {Q,R<:Real,S<:AbstractAbsorber,U,V}
    @assert length(M) == length(x)
    #setup gaussian quadrature nodes
    𝓂, 𝒲 = streamnodes(nstream)
    #solve schwarzschild w multiple streams, integrating over hemisphere
    for i ∈ 1:nstream-1
        stream!(M, x, dIdx, I₀, 𝔸, idx, g, 𝓂[i], fT, fμ, tol)
        M .*= 𝒲[i]/𝒲[i+1]
    end
    stream!(M, x, dIdx, I₀, 𝔸, idx, g, 𝓂[end], fT, fμ, tol)
    M .*= 𝒲[end]
    return nothing
end

#-------------------------------------------------------------------------------
# core function for whole atmosphere upward and downward monochromatic fluxes

function monoflux!(M⁻, #downward monochromatic fluxes [W/m^2/cm^-1]
                   M⁺, #upward monochromatic fluxes [W/m^2/cm^-1]
                   P, #pressure coordinates of output
                   ω, #transformed pressure coords
                   ι, #transformed pressure coords
                   𝔸::Q,
                   idx::Int,
                   g::Real, #gravity [m/s^2]
                   fT::R, #temperature profile fT(P)
                   fμ::S, #mean molar mass μ(T,P)
                   fS::U, #incoming stellar radiation fS(ν) [W/m^2]
                   fα::V, #surface albedo fα(ν)
                   nstream::Int, #number of streams to integrate in both directions
                   θₛ::Real, #stellar radiation angle, corresponds to cos(θ) = 2/3
                   tol::Real) where {Q<:AbstractAbsorber,R,S,U,V}
    #setup
    @assert length(M⁻) == length(M⁺) == length(P)
    #surface pressure assuming ascending pressures
    Pₛ = P[end]
    #surface temperature
    Tₛ = fT(Pₛ)
    #wavenumber
    ν = 𝔸.ν[idx]
    #angle factor for incoming stellar radiation
    𝓂 = 1/cos(θₛ)
    #downward stellar irradiance at ν
    Iₜ⁻ = fS(ν)

    #===================================
    downgoing flux throughout atmosphere
    ===================================#
    #cosine of the stellar zenith angle
    c = cos(θₛ)
    #atmospheric contribution to downward flux
    streams!(M⁻, ι, dIdι, zero(Iₜ⁻), 𝔸, idx, g, fT, fμ, nstream, tol)
    #divide by c before adding the stellar irradiance in-place
    M⁻ ./= c
    #add downward stellar irradiance
    stream!(M⁻, ι, dIdι, Iₜ⁻, 𝔸, idx, g, 𝓂, fT, fμ, tol)
    #multiply everything by c to get the true flux
    M⁻ .*= c

    #===================================
    upgoing flux throughout atmosphere
    ===================================#
    #some of the downward stellar flux is reflected
    Iₛ⁺ = M⁻[end]*fα(ν)/π #Lambertian
    #and the surface emits some radiation
    Iₛ⁺ += planck(ν, Tₛ)
    #upward radiation streams
    streams!(M⁺, ω, dIdω, Iₛ⁺, 𝔸, idx, g, fT, fμ, nstream, tol)
    #reverse the upward flux to match the coordinate ordering of P and ι
    reverse!(M⁺)

    nothing
end