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

# Constructors

    AtmosphericDomain(Trange, nT, Prange, nP)

Creates a domain with the given temperature/pressure ranges and numbers of points. `Trange` and `Prange` should be tuples of two values. `nT` and `nP` indicate the number of points to use.

    AtmosphericDomain()

For convenience, creates a domain with 12 temperature points in `[25, 550]` K and 24 pressure points in `[1,1e6]` Pa.
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

function AtmosphericDomain(Trange::Tuple{Real,Real}, nT::Int,
                           Prange::Tuple{Real,Real}, nP::Int)
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

function AtmosphericDomain()
    Trange = (25, 550)
    Prange = (1, 1e6)
    nT = 12
    nP = 24
    AtmosphericDomain(Trange, nT, Prange, nP)
end

#-------------------------------------------------------------------------------
#wrapper type for the BichebyshevInterpolators used for cross-sections

export OpacityTable

"""
An `OpacityTable` is a simple object wrapping a [BichebyshevInterpolator](https://wordsworthgroup.github.io/BasicInterpolators.jl/stable/2d/#BasicInterpolators.BichebyshevInterpolator). Inside, the interpolator stores a grid of `log` cross-section values along `log` pressure coordinates and temperature coordinates. An `OpacityTable` behaves like a function, recieving a temperature and pressure. When called, it retrieves a cross-section from the interpolator, undoes the `log`, and returns it. When constructing a gas object, each wavenumber is allocated a unique `OpacityTable` for fast and accurate cross-section evaluation at any temperature and pressure inside the `AtmosphericDomain`. Generally, `OpacityTable` objects should be used indirectly through gas objects.
"""
struct OpacityTable
    itp::BichebyshevInterpolator
    emp::Bool
end

function OpacityTable(T::AbstractVector{<:Real},
                      P::AbstractVector{<:Real},
                      σ::AbstractArray{<:Real,2})
    if all(σ .== 0)
        #avoid evaluating log(0) and passing -Infs to the interp constructor
        itp = BichebyshevInterpolator(T, log.(P), fill(0.0, size(σ)))
        emp = true
    else
        itp = BichebyshevInterpolator(T, log.(P), log.(σ))
        emp = false
    end
    OpacityTable(itp, emp)
end

#gets cross sections out of interpolators, un-logged, cm^2/molecule
#also explicitly handles empty tables
(Π::OpacityTable)(T, P)::Float64 = Π.emp ? 0.0 : exp(Π.itp(T, log(P)))

#-------------------------------------------------------------------------------
#function for building gas opacity tables

function bake(sl::SpectralLines,
              Cfun::F,
              shape!::G,
              ν::AbstractVector{<:Real},
              Ω::AtmosphericDomain
              )::Vector{OpacityTable} where {F, G<:Function}
    #check wavenumbers for problems
    @assert all(diff(ν) .> 0) "wavenumbers must be unique and in ascending order"
    @assert all(ν .>= 0) "wavenumbers must be positive"
    #number of wavenumbers
    nν = length(ν)
    #create a single block of cross-sections
    σ = zeros(nν, Ω.nT, Ω.nP)
    #fill it by evaluating in batches of wavenumbers (slow part)
    @threads for i = 1:Ω.nT # @distributed??
        for j = 1:Ω.nP
            #get a view into the big σ array
            σᵢⱼ = view(σ,:,i,j)
            #get temperature, pressure, concentration
            T = Ω.T[i]
            P = Ω.P[j]
            C = Cfun(T, P)
            #make sure concentration isn't wacky
            @assert 0 <= C <= 1 "gas molar concentrations must be in [0,1], not $C (encountered @ $T K, $P Pa)"
            #evaluate line shapes (slow part)
            shape!(σᵢⱼ, ν, sl, T, P, C*P)
        end
    end
    #check for weirdness
    z = zeros(Bool, nν)
    for i = 1:nν
        σᵥ = view(σ,i,:,:)
        if (minimum(σᵥ) == 0) & (maximum(σᵥ) > 0)
            z[i] = true
        end
    end
    if any(z)
        @info "Zero cross-section values are mixed with non-zero values for the following wavenumbers for $(sl.name):\n\n$(ν[z])\n\n Likely, absorption is extremely weak in these regions. Absorption is being set to zero for all T and P at those wavenumbers to avoid non-smooth and inaccurate interpolation tables."
        σ[z,:,:] .= 0.0
    end
    #split the block and create interpolators for each ν
    Π = Vector{OpacityTable}(undef, nν)
    @threads for i = 1:nν
        Π[i] = OpacityTable(Ω.T, Ω.P, σ[i,:,:])
    end
    #fresh out the oven
    return Π
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
#defining absorbers and access to cross-sections

abstract type AbstractGas end

export AbstractGas
export WellMixedGas, VariableGas
export concentration, reconcentrate

#abundance weighted average molar mass
meanmolarmass(sl::SpectralLines) = sum(sl.A .* sl.μ)/sum(sl.A)

#-------------------------------

struct WellMixedGas <: AbstractGas
    name::String
    formula::String
    μ::Float64 #mean molar mass [kg/mole]
    C::Float64 #constant fraction [mole/mole] of dry gas constituents
    ν::Vector{Float64}
    Ω::AtmosphericDomain
    Π::Vector{OpacityTable} #cross-section interpolators
end

function WellMixedGas(sl::SpectralLines,
                      C::Real,
                      ν::AbstractVector{<:Real},
                      Ω::AtmosphericDomain,
                      shape!::Function=voigt!)
    μ = meanmolarmass(sl)
    ν = collect(Float64, ν)
    Π = bake(sl, (T,P)->C, shape!, ν, Ω)
    WellMixedGas(sl.name, sl.formula, μ, C, ν, Ω, Π)
end

function WellMixedGas(par, C, ν, Ω, shape!::Function=voigt!; kwargs...)
    sl = SpectralLines(par; kwargs...)
    WellMixedGas(sl, C, ν, Ω, shape!)
end

#-------------------------------

struct VariableGas{F} <: AbstractGas
    name::String
    formula::String
    C::F #concentration [mole/mole] from temperature and pressure, C(T,P)
    ν::Vector{Float64}
    Ω::AtmosphericDomain
    Π::Vector{OpacityTable} #cross-section interpolators
end

function VariableGas(sl::SpectralLines,
                     C::Q,
                     ν::AbstractVector{<:Real},
                     Ω::AtmosphericDomain,
                     shape!::Function=voigt!) where {Q}
    μ = meanmolarmass(sl)
    ν = collect(Float64, ν)
    Π = bake(sl, C, shape!, ν, Ω)
    VariableGas(sl.name, sl.formula, μ, C, ν, Ω, Π)
end

function VariableGas(par, C, ν, Ω, shape!::Function=voigt!; kwargs...)
    sl = SpectralLines(par; kwargs...)
    VariableGas(sl, C, ν, Ω, shape!)
end

#-------------------------------

concentration(g::WellMixedGas) = g.C

concentration(g::WellMixedGas, X...)::Float64 = g.C

concentration(g::VariableGas, T, P)::Float64 = g.C(T,P)

function reconcentrate(g::WellMixedGas, C::Real)::WellMixedGas
    @assert 0 <= C <= 1 "gas molar concentrations must be in [0,1], not $C"
    Ω = deepcopy(g.Ω)
    Π = deepcopy(g.Π)
    WellMixedGas(g.name[:], g.formula[:], g.μ, C, g.ν, Ω, Π)
end

#-------------------------------------------------------------------------------

#single cross-section
(g::AbstractGas)(i::Int, T, P)::Float64 = concentration(g, T, P)*g.Π[i](T,P)

#for full vectors of cross-sections with whatever gas
function (g::AbstractGas)(T::Real, P::Real)::Vector{Float64}
    Float64[g(i, T, P) for i ∈ 1:length(g.ν)]
end
