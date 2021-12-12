#-------------------------------------------------------------------------------
# defining the radiative domain and associated functions

export AtmosphericDomain

"""
Structure defining the temperature and pressure ranges over which absorption cross-sections are generated when constructing gas objects. `AtmosphericDomain` objects store the temperature and pressure coordinates of cross-section interpolation grids. More points lead to higher accuracy interpolation. Generally, about 12 temperature points and 24 pressure points results in maximum error of ~1 % and much smaller average error.

| Field | Type | Description |
| ----- | :--- | :---------- |
| `T` | `Vector{Float64}` | temperature coordinates of grid [K] |
| `Tmin` | `Float64` | lowest temperature value |
| `Tmax` | `Float64` | highest temperature value |
| `nT` | `Int64` | number of temperature coordinates |
| `P` | `Vector{Float64}` | pressure coordinates of grid [Pa] |
| `Pmin` | `Float64` | lowest pressure value |
| `Pmax` | `Float64` | highest pressure value |
| `nP` | `Int64` | number of pressure coordinates |

# Constructor

    AtmosphericDomain(Trange, nT, Prange, nP)

Creates a domain with the given temperature/pressure ranges and numbers of points. `Trange` and `Prange` should be tuples of two values. `nT` and `nP` indicate the number of points to use.
"""
struct AtmosphericDomain
    #temperature samples for interpolator [K]
    T::Vector{Float64}
    Tmin::Float64
    Tmax::Float64
    nT::Int64
    #pressure samples for interpolator [atm]
    P::Vector{Float64}
    Pmin::Float64
    Pmax::Float64
    nP::Int64
end

function Base.show(io::IO, Ω::AtmosphericDomain)
    print(io, "AtmosphericDomain:\n")
    print(io, "  $(Ω.nT) temperature nodes ∈ [$(Ω.Tmin),$(Ω.Tmax)] K\n")
    print(io, "  $(Ω.nP) pressure nodes ∈ [$(Ω.Pmin),$(Ω.Pmax)] Pa")
end

function AtmosphericDomain(Trange::NTuple{2,Real}, nT::Int,
                           Prange::NTuple{2,Real}, nP::Int)
    #check for negatives
    @assert all(Trange .> 0) "temperature range must be positive"
    @assert all(Prange .> 0) "pressure range must be positive"
    #check the Qref/Q range
    @assert all(Trange .>= TMIN) "minimum temperature with Qref/Q accuracy is $TMIN K"
    @assert all(Trange .<= TMAX) "maximum temperature with Qref/Q accuracy is $TMAX K"
    #order
    @assert Trange[1] < Trange[2] "Trange[1] ($(Trange[1])) can't be greater than Trange[2] ($(Trange[2]))"
    @assert Prange[1] < Prange[2] "Prange[1] ($(Prange[1])) can't be greater than Prange[2] ($(Prange[2]))"
    #generate grid points
    T = chebygrid(Trange[1], Trange[2], nT)
    P = exp.(chebygrid(log(Prange[1]), log(Prange[2]), nP))
    #assemble!
    AtmosphericDomain(T, Trange[1], Trange[2], nT, P, Prange[1], Prange[2], nP)
end

#-------------------------------------------------------------------------------
#wrapper type for BichebyshevInterpolators used for cross-sections

export OpacityTable

struct OpacityTable{M,N}
    Φ::BichebyshevInterpolator{M,N,Float64}
end

# T: temperature grid coordinates [K]
# P: pressure grid coordinates [Pa]
# σ: cross-section grid [cm^2/molecule]
function OpacityTable(T, P, σ)
    tiny = floatmin(eltype(T))
    ζ = all(σ .<= tiny) #whether all values are effectively zero
    lnP = log.(P)
    lnσ = ζ ? fill(log(tiny), size(σ)) : log.(σ)
    Φ = BichebyshevInterpolator(T, lnP, lnσ)
    OpacityTable(Φ)
end

#gets cross-section out of interpolators, un-logged [cm^2/molecule]
(Π::OpacityTable)(T, P) = exp(Π.Φ(T, log(P)))

#-------------------------------------------------------------------------------
#function for building gas opacity tables

function checkν(ν::AbstractVector)::Nothing
    #check wavenumbers for problems
    @assert all(diff(ν) .> 0) "wavenumbers must be unique and in ascending order"
    @assert all(ν .>= 0) "wavenumbers must be positive"
    nothing
end

function bake(sl::SpectralLines,
              fC::F,
              shape!::G,
              Δνcut::Real,
              ν::Vector{Float64},
              Ω::AtmosphericDomain,
              progress::Bool=true
              ) where {F, G<:Function}
    checkν(ν)
    #number of wavenumbers
    nν = length(ν)
    #create a single block of cross-sections
    σ = zeros(nν, Ω.nT, Ω.nP)
    #initialize progress meter if desired
    if progress
        prg = Progress(Ω.nT*Ω.nP, 0.1, "evaluating $(sl.formula) cross-sections ")
    end
    #fill σ by evaluating in batches of wavenumbers (slow part)
    @threads for i ∈ eachindex(Ω.T)
        for j ∈ eachindex(Ω.P)
            #get a view into the big σ array
            σᵢⱼ = @view σ[:,i,j]
            #get temperature, pressure, concentration
            T = Ω.T[i]
            P = Ω.P[j]
            C = fC(T, P)
            #make sure concentration isn't wacky
            @assert 0 <= C <= 1 "gas molar concentrations must be in [0,1], not $C (encountered @ $T K, $P Pa)"
            #evaluate line shapes (slow part)
            shape!(σᵢⱼ, ν, sl, T, P, C*P, Δνcut)
            #update progress meter
            progress && next!(prg)
        end
    end
    #check for weirdness
    z = zeros(Bool, nν)
    for i ∈ eachindex(ν)
        σᵥ = @view σ[i,:,:]
        if (minimum(σᵥ) == 0) & (maximum(σᵥ) > 0)
            z[i] = true
        end
    end
    if any(z)
        @info "Zero cross-section values are mixed with non-zero values for the following wavenumbers for $(sl.name):\n\n$(ν[z])\n\n Likely, absorption is extremely weak in these regions, causing underflow. Absorption is being set to zero for all temperatures and pressures at those wavenumbers to avoid non-smooth and inaccurate interpolation tables."
        σ[z,:,:] .= 0.0
    end
    #split the block and create opacity tables for each ν
    [OpacityTable(Ω.T, Ω.P, σ[i,:,:]) for i ∈ eachindex(ν)]
end

#-------------------------------------------------------------------------------
#for testing opacity table errors

export opacityerror

function opacityerror(Π::OpacityTable,
                      Ω::AtmosphericDomain,
                      sl::SpectralLines,
                      ν::Real,
                      C::F, #C(T,P)
                      shape::G=voigt,
                      N::Int=50) where {F,G}
    #create T and P grids from the domain
    T = LinRange(Ω.Tmin, Ω.Tmax, N)
    P = 10 .^ LinRange(log10(Ω.Pmin), log10(Ω.Pmax), N)
    #compute exact and approximate cross-sections
    σop = zeros(N,N)
    σex = zeros(N,N)
    @threads for i = 1:N
        for j = 1:N
            σop[i,j] = Π(T[i], P[j])
            σex[i,j] = shape(ν, sl, T[i], P[j], C(T[i],P[j])*P[j])
        end
    end
    #compute error and relative error
    aerr = σop .- σex
    rerr = aerr./σex
    return T, P, aerr, rerr
end

#-------------------------------------------------------------------------------
#defining Gas struct and access to cross-sections

abstract type AbstractGas end

export Gas
export rawσ, concentration, reconcentrate

"""
Gas type for real-gas, radiatively active, atmospheric constituents. Must be constructed from a `.par` file or a [`SpectralLines`](@ref) object.

# Constructors

    Gas(sl::SpectralLines, fC, ν, Ω, shape!=voigt!, Δνcut=25; progress=true)

* `sl`: a [`SpectralLines`](@ref) object
* `fC`: molar concentration of the constituent [mole/mole] as a function of temperature and pressure `fC(T,P)`
* `ν`: vector of wavenumber samples [cm``^{-1}``]
* `Ω`: [`AtmosphericDomain`](@ref)
* `shape!`: line shape to use, must be the in-place version ([`voigt!`](@ref), [`lorentz!`](@ref), etc.)
* `Δνcut`: profile truncation distance [cm``^{-1}``]
* `progress`: whether to display the progress meter


    Gas(par::String, fC, ν, Ω, shape!=voigt!, Δνcut=25; progress=true)

Same arguments as the first constructor, but reads a `par` file directly into the gas object. Keyword arguments are passed through to [`readpar`](@ref).
"""
struct Gas{T,F} <: AbstractGas
    name::String
    formula::String
    μ::Float64 #mean molar mass
    ν::Vector{Float64}
    Ω::AtmosphericDomain
    Π::Vector{T} #cross-section interpolators
    fC::F #molar concentration [mole/mole] function fC(T,P)
end

function Base.show(io::IO, g::Gas)
    print(io, "$(g.name) ($(g.formula))\n")
    print(io, "  μ = $(round(g.μ, sigdigits=8)) kg/mole\n")
    print(io, "  ν: $(minimum(g.ν)) - $(maximum(g.ν)) cm^-1 ($(length(g.ν)) samples)\n")
    Tmin, Tmax = g.Ω.Tmin, g.Ω.Tmax
    Pmin, Pmax = g.Ω.Pmin, g.Ω.Pmax
    print(io, "  T: $Tmin - $Tmax K ($(g.Ω.nT) samples)\n")
    print(io, "  P: $Pmin - $Pmax Pa ($(g.Ω.nP) samples)\n")
end

function Gas(sl::SpectralLines,
             fC::F,
             ν::AbstractVector{<:Real},
             Ω::AtmosphericDomain,
             shape!::Function=voigt!,
             Δνcut::Real=25;
             progress::Bool=true) where {F}
    @assert length(ν) > 0
    μ = sum(sl.A .* sl.μ)/sum(sl.A) # mean molar mass [kg/mole]
    ν = collect(Float64, ν)
    Π = bake(sl, fC, shape!, Δνcut, ν, Ω, progress)
    T = eltype(Π)
    Gas{T,F}(sl.name, sl.formula, μ, ν, Ω, Π, fC)
end

function Gas(filename::String,
             fC::F,
             ν::AbstractVector{<:Real},
             Ω::AtmosphericDomain,
             shape!::Function=voigt!,
             Δνcut=25;
             kwargs...) where {F}
    sl = SpectralLines(filename; kwargs...)
    Gas(sl, fC, ν, Ω, shape!, Δνcut; kwargs...)
end

"""
    rawσ(g::Gas, i, T, P)

Retrieve the absorption cross-section at wavenumber index `i`, temperature `T` [K], and pressure `P` [Pa], **without** the concentration scaling.
"""
rawσ(g::Gas, i::Int, T, P) = g.Π[i](T, P)

"""
    rawσ(g::Gas, T, P)

Retrieve the absorption cross-section for all wavenumbers at temperature `T` [K], and pressure `P` [Pa], **without** the concentration scaling.
"""
rawσ(g::Gas, T, P) = [Π(T, P) for Π ∈ g.Π]

"""
    concentration(g::Gas, T, P)

Furnishes the molar concentration [mole/mole] of a [`Gas`](@ref) object at a particular temperature and pressure. Identical to `g.C(T,P)`.
"""
concentration(g::Gas, T, P) = g.fC(T,P)

function concentration(g::Gas, T::AbstractVector, P::AbstractVector)
    @assert length(T) == length(P)
    [concentration(g, T[i], P[i]) for i ∈ eachindex(P)]
end

#single concentration-scaled cross-section
(g::Gas)(i::Int, T, P) = concentration(g, T, P)*rawσ(g, i, T, P)

#for full vectors of concentration-scaled cross-sections with whatever gas
(g::Gas)(T, P) = concentration(g, T, P)*rawσ(g, T, P)

"""
    reconcentrate(g::Gas, fC)

Create a copy of a [`Gas`](@ref) object with a new molar concentration function, `fC(T,P)` [mole/mole].

!!! warning

    The self-broadening component of the line shape is not recomputed when using the `reconcentrate` function. This component is generally small when partial pressure is low, but may be appreciable if the concentration changes significantly.
"""
function reconcentrate(g::Gas, fC)::Gas
    #check for invalid concentrations or unexpected input
    for P ∈ g.Ω.P, T ∈ g.Ω.T
        C = try
            fC(T,P)
        catch
            error("concentration function must have the form fC(T,P), where T is temperature [K] and P is pressure [Pa], but something went wrong when the function was evaluated")
        end
        @assert 0 <= C <= 1.0 "gas molar concentrations must be in [0,1], not $C, which was encountered at T=$T P=$P"
    end
    #construct a new gas WITHOUT COPYING
    Gas(g.name, g.formula, g.μ, g.ν, g.Ω, g.Π, fC)
end

"""
    reconcentrate(g::Gas, C)

Create a copy of a [`Gas`](@ref) object with a new, **uniform** molar concentration [mole/mole].

!!! warning

    The self-broadening component of the line shape is not recomputed when using the `reconcentrate` function. This component is generally small when partial pressure is low, but may be appreciable if the concentration changes significantly.
"""
function reconcentrate(g::Gas, C::Real)::Gas where {U,V}
    #check for invalid concentrations
    @assert 0 <= C <= 1.0 "gas molar concentrations must be in [0,1], not $C"
    #construct a new gas WITHOUT COPYING
    Gas(g.name, g.formula, g.μ, g.ν, g.Ω, g.Π, (T,P)->float(C))
end

function Base.getindex(g::Gas, x)
    Gas(
        g.name[:],
        g.formula[:],
        g.μ,
        g.ν[x],
        g.Ω,
        g.Π[x],
        g.fC
    )
end

Base.firstindex(g::Gas) = 1

Base.lastindex(g::Gas) = length(g.ν)

#---------------------------------------

export GrayGas

struct GrayGas{T} <: AbstractGas
    name::String
    formula::String
    μ::Float64
    ν::Vector{Float64}
    σ::T
end

function Base.show(io::IO, g::GrayGas)
    print(io, "Gray Gas\n")
    print(io, "  ν: $(minimum(g.ν)) - $(maximum(g.ν)) cm^-1 ($(length(g.ν)) samples)\n")
    print(io, "  σ = $(g.σ) cm^2/molecule")
end

function GrayGas(σ::Real, ν::AbstractVector{<:Real})
    GrayGas("Gray", "Gray", NaN, collect(Float64, ν), σ)
end

(g::GrayGas)(x...) = g.σ

#---------------------------------------

export SemiGrayGas

struct SemiGrayGas{T} <: AbstractGas
    name::String
    formula::String
    μ::Float64
    ν::Vector{Float64}
    νcut::Float64
    σ::T
end

function Base.show(io::IO, g::SemiGrayGas)
    println(io, "Semi-Gray Gas")
    println(io, "  ν: $(minimum(g.ν)) - $(maximum(g.ν)) cm^-1 ($(length(g.ν)) samples)")
    println(io, "  σ = $(g.σ) cm^2/molecule")
    print(io, "  absorption cutoff = $(g.νcut) cm^-1")
end

function SemiGrayGas(σ::Real, ν::AbstractVector{<:Real}, νcut::Real)
    SemiGrayGas("SemiGray", "SemiGray", NaN, collect(Float64, ν), Float64(νcut), σ)
end

(g::SemiGrayGas)(i::Int, _...) = (g.ν[i] <= g.νcut) ? g.σ : zero(g.σ)
