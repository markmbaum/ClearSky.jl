#-------------------------------------------------------------------------------
# super for both consolidated absorber types

export AbstractAbsorber

abstract type AbstractAbsorber end

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
    #tuple of CIA objects
    cia::U
    #tuple of functions
    fun::V
    #wavenumber vector [cm^-1], must be identical for all gases
    ν::Vector{Float64}
    #length of wavenumber vector
    nν::Int64
end

function Base.show(io::IO, U::UnifiedAbsorber)
    print(io, "UnifiedAbsorber")
    if !isempty(U.gas)
        s = join(map(g->g.formula, U.gas), ", ", " and ")
        print(io, "\n  $(length(U.gas)) gas(es): $s")
    end
    if !isempty(U.cia)
        s = join(map(c->c.name, U.cia), ", ", " and ")
        print(io, "\n  $(length(U.cia)) CIA(s): $s")
    end
    if !isempty(U.fun)
        print(io, "\n  $(length(U.fun)) other function(s)")
    end

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
    #real gases, ignoring Gray
    realgas = gas[findall(g -> isa(g, Gas), gas)]
    #cia tables, pairing with the correct gases in the process
    cia = tuple([CIA(x, realgas) for x ∈ absorbers[findall(t -> t == CIATables, T)]]...)
    #functions in the form σ(ν, T, P)
    fun = absorbers[findall(t -> !(t <: AbstractGas) & !(t == CIATables), T)]
    #wavenumber vector, must be identical for all gases
    ν = getwavenumbers(gas...)
    nν = length(ν)
    #construct the UnifiedAbsorber
    UnifiedAbsorber(gas, cia, fun, ν, nν)
end

#blank update function for similarity with AcceleratedAbsorber
function update!(A::UnifiedAbsorber, T)::Nothing end

#https://discourse.julialang.org/t/tuple-indexing-taking-time/58309/18?u=markmbaum
#also see "applychain" here for a similar example: https://github.com/FluxML/Flux.jl/blob/dbb9f82ef8d4e196259ff1af56aeddc626159bf3/src/layers/basic.jl#L46
σchain(::Tuple{}, x, T, P) = 0

σchain(a::Tuple, x, T, P) = first(a)(x, T, P) + σchain(tail(a), x, T, P)

function σchain(U::UnifiedAbsorber, i::Int, ν, T, P)
    ( σchain(U.gas, i, T, P)
    + σchain(U.cia, ν, T, P)
    + σchain(U.fun, ν, T, P))
end

#internal
Σ(U::UnifiedAbsorber, i::Int, T, P) = @inbounds σchain(U, i, U.ν[i], T, P)

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
struct AcceleratedAbsorber{V,Q<:UnifiedAbsorber} <: AbstractAbsorber
    #cross-section interpolators
    ϕ::Vector{LinearInterpolator{V, NoBoundaries}}
    #wavenumber vector [cm^-1], must be identical for all gases
    ν::Vector{Float64}
    #length of wavenumber vector
    nν::Int64
    #original temperatures
    T::Vector{Float64}
    #original pressures
    P::Vector{Float64}
    #stored reference to UnifiedAbsorber
    U::Q
end

function Base.show(io::IO, A::AcceleratedAbsorber)
    print(io, "AcceleratedAbsorber ($(length(A.P)) pressure samples)\n")
    print("generated from ")
    print(A.U)
end

function AcceleratedAbsorber(T::AbstractArray{V},
                             P::AbstractArray,
                             U::Q) where {V<:Real,Q<:UnifiedAbsorber}
    #pull out wavenumber info
    ν, nν = U.ν, U.nν
    #flip vectors if pressure is not ascending (interpolators require this)
    idx = sortperm(P)
    P = P[idx]
    T = T[idx]
    np = length(P)
    #log pressure coordinates as usual
    lnP = log.(P)
    #initalize an AcceleratedAbsorber with empty interpolators
    ϕ = Vector{LinearInterpolator{V, NoBoundaries}}(undef,nν)
    for i ∈ eachindex(ν)
        ϕ[i] = LinearInterpolator(lnP, Vector{V}(undef,np), NoBoundaries())
    end
    A = AcceleratedAbsorber{V,Q}(ϕ, ν, nν, T, P, deepcopy(U))
    #then update the cross-sections in-place
    update!(A, T)
    #and return the updated AcceleratedAbsorber
    return A
end

function AcceleratedAbsorber(T, P, absorbers...)
    AcceleratedAbsorber(T, P, UnifiedAbsorber(absorbers))
end

function AcceleratedAbsorber(::Any, P, A::AcceleratedAbsorber)
    @assert all(P .== A.P) "cannot change AcceleratedAbsorber's pressure coordinates after construction"
    return A
end

"""
    update!(A::AcceleratedAbsorber, T)

Update the cross-section interpolators underlying an `AcceleratedAbsorber` with a new set of temperatures. The new temperatures should correspond to the pressure levels used when originally constructing the `AcceleratedAbsorber`.  
"""
function update!(A::AcceleratedAbsorber{V,Q}, T::AbstractVector)::Nothing where {V,Q}
    #check lengths
    @assert length(T) == length(A.P)
    #not parallel b/c BichebyshevInterpolator not thread safe
    for i ∈ eachindex(T)
        update!(A, T[i], i)
    end
    return nothing
end

function update!(A::AcceleratedAbsorber{V,Q}, T::Real, idx::Int)::Nothing where {V,Q}
    #uniform type (or error)
    T = convert(V, T)
    #smallest allowable number
    logtiny = log(floatmin(T))
    #pressure level
    P = A.P[idx]
    #update each interpolator at index idx
    for (i,ϕ) ∈ enumerate(A.ϕ)
        #log of cross-section
        lnσ = log(Σ(A.U, i, T, P))
        #set value inside interpolator
        ϕ[idx] = (lnσ < logtiny) ? logtiny : lnσ
    end
    #remember the temperature
    A.T[idx] = T
    return nothing
end

#internal
Σ(A::AcceleratedAbsorber, i::Int, ::Any, P) = @inbounds exp(A.ϕ[i](log(P)))

(A::AcceleratedAbsorber)(i::Int, P) = exp(A.ϕ[i](log(P)))

(A::AcceleratedAbsorber)(P) = [A(i, P) for i ∈ eachindex(A.ν)]

checkpressures(A::AcceleratedAbsorber, P...) = checkpressures(A.U, P...)

#-------------------------------------------------------------------------------
#making sense of absorber inputs to modeling functions

unifyabsorbers(x::Tuple{UnifiedAbsorber}) = x[1], x[1].ν, x[1].nν

unifyabsorbers(x::Tuple{AcceleratedAbsorber}) = x[1], x[1].ν, x[1].nν

function unifyabsorbers(x::Tuple)
    U = UnifiedAbsorber(x)
    return U, U.ν, U.nν
end

unifyabsorbers(::Tuple{}) = error("no absorbers")

#checks for identical wavenumber sampling across different gases
function getwavenumbers(G::AbstractGas...)::Vector{Float64}
    @assert all(g -> g.ν == G[1].ν, G) "gases must have identical wavenumber vectors"
    return G[1].ν
end

function getwavenumbers(absorbers::Tuple)::Vector{Float64}
    G = absorbers[findall(a -> typeof(a) <: AbstractGas, absorbers)]
    @assert length(G) > 0 "no gas objects found"
    getwavenumbers(G...)
end

function checkpressures(gases::Tuple, Pₛ, Pₜ)
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    #pressure bounds
    Pmin, Pmax = pressurelimits(gases)
    #demand all pressures within the range
    for P ∈ (Pₛ, Pₜ)
        @assert P >= Pmin "Pressure $P Pa too low, domain minimum is $Pmin"
        @assert P <= Pmax "Pressure $P Pa too low, domain minimum is $Pmax"
    end
end

function pressurelimits(gases::Tuple)::NTuple{2,Float64}
    g = gases[findall(g -> typeof(g) <: Gas, gases)]
    isempty(g) && return(0.0, Inf)
    #largest minimum pressure in gas atmospheric domains
    Pmin = maximum(map(g->g.Ω.Pmin, g))
    #smallest maximum pressure in gas atmospheric domains
    Pmax = minimum(map(g->g.Ω.Pmax, g))
    return Pmin, Pmax
end

function temperaturelimits(gases::Tuple)::NTuple{2,Float64}
    g = gases[findall(g -> typeof(g) <: Gas, gases)]
    isempty(g) && return(0.0, Inf)
    #largest minimum pressure in gas atmospheric domains
    Tmin = maximum(map(g->g.Ω.Tmin, g))
    #smallest maximum pressure in gas atmospheric domains
    Tmax = minimum(map(g->g.Ω.Tmax, g))
    return Tmin, Tmax
end

temperaturelimits(U::UnifiedAbsorber) = temperaturelimits(U.gas)

temperaturelimits(A::AcceleratedAbsorber) = temperaturelimits(A.U)