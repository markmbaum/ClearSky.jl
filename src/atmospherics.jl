#minimum pressure in temperature profiles and floor for hydrostatic profile
const PMIN = 1e-9

#-------------------------------------------------------------------------------
#constructing generalized hydrostatic pressure profiles and inverting for z

export scaleheight, hydrostatic, altitude
export Hydrostatic

"""
    scaleheight(g, μ, T)

Evaluate the [atmospheric scale height](https://en.wikipedia.org/wiki/Scale_height),

``\\frac{RT}{μg}``

where ``R`` is the [universial gas constant](https://en.wikipedia.org/wiki/Gas_constant).

# Arguments
* `g`: gravitational acceleration [m/s``^s``]
* `μ`: mean molar mass [kg/mole]
* `T`: temperature [K]
"""
scaleheight(g, μ, T) = 𝐑*T/(μ*g)

#parameterized hydrostatic relation in log coordinates
function dlnPdz(z, lnP, param::Tuple)
    #unpack parameters
    Pₛ, g, fT, fμ = param
    #evaluate temperature and mean molar mass [kg/mole]
    P = exp(lnP)
    P < PMIN && return zero(lnP)
    P = min(P, Pₛ) #don't allow tiny numerical dips below Pₛ
    T = fT(P)
    μ = fμ(T, P)
    #evaluate derivative
    -μ*g/(𝐑*T)
end

"""
    hydrostatic(z, Pₛ, g, fT, fμ)

Compute the hydrostatic pressure [Pa] at a specific altitude using arbitrary atmospheric profiles of temperature and mean molar mass. This function integrates the hydrostatic relation,

``\\frac{dP}{dz} = \\frac{\\mu g}{R T}

from the surface to a height of ``z``, where ``R`` is the [universial gas constant](https://en.wikipedia.org/wiki/Gas_constant).

# Arguments

* `z`: altitude [m] to compute pressure at
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure, `fT(P)`
* `fμ`: mean molar mass [kg/mole] as a function of pressure and temperature `fμ(T,P)`
"""
function hydrostatic(z, Pₛ, g, fT::T, fμ::U) where {T,U}
    @assert z >= 0 "cannot compute pressure at negative altitude $z m"
    @assert Pₛ > PMIN "pressure cannot be less than $PMIN Pa"
    #parameters
    param = (Pₛ, g, fT, fμ)
    #integrate in log coordinates
    lnP = radau(dlnPdz, log(Pₛ), 0, z, param)
    #convert
    exp(lnP)
end

"""
    altitude(P, Pₛ, g, fT, fμ)

Compute the altitude [m] at which a specific hydrostatic pressure occurs using arbitrary atmospheric profiles of temperature and mean molar mass. This function applies a root finder to the [`hydrostatic`](@ref) function.

# Arguments

* `P`: pressure [Pa] to compute altitude at
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure, `fT(P)`
* `fμ`: mean molar mass [kg/mole] as a function of pressure and temperature `fμ(T,P)`
"""
function altitude(P, Pₛ, g, fT::T, fμ::U) where {T,U}
    @assert P < Pₛ "surface pressure must be greater than pressure aloft"
    #pressure decreases monotonically, find altitudes bracketing Pₜ
    z₁ = 0.0
    z₂ = 1e2
    P₁ = Pₛ
    P₂ = hydrostatic(z₂, Pₛ, g, fT, fμ)
    while P₂ > P
        z₁ = z₂
        z₂ *= 2
        P₁ = P₂
        P₂ = hydrostatic(z₂, Pₛ, g, fT, fμ)
    end
    #find precise altitude where P = Pₜ
    regulafalsi((z,p) -> log(hydrostatic(z, Pₛ, g, fT, fμ)) - log(P), z₁, z₂)
end

"""
[Function-like type](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) for initializing and evaluating a hydrostatic pressure profile with arbitrary temperature and mean molar mass profiles. A `Hydrostatic` object maps altitude to pressure. Internally, a pressure vs altitude profile is generated and used for interpolation.

# Constructor

    Hydrostatic(Pₛ, Pₜ, g, fT, fμ, N=250)

* `Pₛ`: surface pressure [Pa]
* `Pₜ`: top of profile pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of presssure, `fT(P)`
* `fμ`: mean molar mass [kg/mole] as a function of temperature and pressure, `fμ(T,P)`
* `N`: optional, number of interpolation nodes

For a constant molar mass or temperature, you can use [anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions) directly. For example, to construct a hydrostatic pressure profile for a crude Earth-like atmosphere:

```@example
#moist adiabatic temperature profile
M = MoistAdiabat(288, 1e5, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Ptropo=1e4);
#hydrostatic pressure profile with constant mean molar mass
H = Hydrostatic(1e5, 1, 9.8, M, (T,P)->0.029);
#evaluate pressures at a few different altitudes
H.([0, 1e3, 1e4])
```
"""
struct Hydrostatic
    ϕ::LinearInterpolator{Float64,WeakBoundaries}
    zₜ::Float64
end

function Hydrostatic(Pₛ, Pₜ, g, fT::T, fμ::U, N::Int=250) where {T,U}
    #find the altitude corresponding to Pₜ
    zₜ = altitude(Pₜ, Pₛ, g, fT, fμ)
    #interpolation knots and output array
    z = logrange(0, zₜ, N)
    lnP = zeros(Float64, N)
    #integrate to get a full pressure profile
    radau!(lnP, z, dlnPdz, log(Pₛ), 0.0, zₜ, (Pₛ, g, fT, fμ))
    #construct and return
    Hydrostatic(LinearInterpolator(z, lnP, WeakBoundaries()), zₜ)
end

(H::Hydrostatic)(z)::Float64 = exp(H.ϕ(z))

"""
    altitude(H::Hydrostatic, P)

Compute the altitude at which a specific pressure occurs in a [`Hydrostatic`](@ref) pressure profile. A root finder is applied to the object.
"""
function altitude(H::Hydrostatic, P::Real)::Float64
    regulafalsi((z,p) -> log(H(z)) - log(P), 0.0, H.zₜ)
end

#-------------------------------------------------------------------------------

#general function for adiabat with one condensible in bulk non-condensible
function dTdω(ω, T, cₚₙ, cₚᵥ, Rₙ, Rᵥ, L, psat::F) where {F}
    #molar mixing ratio of condensible
    α = psat(T)/ω2P(ω)
    #whole expression at once
    -T*(Rₙ/cₚₙ)*(1 + α*L/(Rₙ*T))/(1 + α*(cₚᵥ/cₚₙ + (L/(T*Rᵥ) - 1)*L/(cₚₙ*T)))
end

dTdω(ω, T, param::Tuple) = dTdω(ω, T, param...)

#------------------------------------

abstract type AbstractAdiabat end

export MoistAdiabat, DryAdiabat
export tropopause

function checkadiabat(Tₛ, Pₛ, Pₜ, Tstrat, Ptropo)
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    @assert Pₜ > 0 "Pₜ must be greater than 0"
    if Tstrat > 0
        @assert Tstrat < Tₛ "Tstrat cannot be greater than Tₛ"
    end
    if (Tstrat != 0) & (Ptropo != 0)
        throw("Cannot have nonzero Tstrat and Ptropo, must use one or the other")
    end
end

#------------------------------------

"""
[Function-like type](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) for initializing and evaluating a dry adiabatic temperature profile. Optionally, a uniform upper atmospheric temperature can be set below a specified temperature or pressure.

# Constructor

    DryAdiabat(Tₛ, Pₛ, cₚ, μ; Tstrat=0.0, Ptropo=0.0, Pₜ=$PMIN)

* `Tₛ`: surface temperature [K]
* `Pₛ`: surface pressure [K]
* `cₚ`: specific heat of the atmosphere [J/kg/K]
* `μ`: molar mass of the atmosphere [kg/mole]

If `Tstrat` is greater than zero, the temperature profile will never drop below that temperature. If `Ptropo` is greater than zero, the temperature profile at pressures lower than `Ptropo` will be equal to the temperature at exactly `Ptropo`. `Tstrat` and `Ptropo` cannot both greater than zero.

`Pₜ` defines the highest pressure [Pa] in the temperature profile. If should generally be small but cannot be zero.

# Example

Once constructed, use a `DryAdiabat` like a function to compute temperature at a given pressure.

```@example
Tₛ = 288; #surface temperature [K]
Pₛ = 1e5; #surface pressure [Pa]
cₚ = 1040; #specific heat of air [J/kg/K]
μ = 0.029; #molar mass of air [kg/mole]

#construct the dry adiabat with an upper atmosphere temperature of 190 K
D = DryAdiabat(Tₛ, Pₛ, cₚ, μ, Tstrat=190);

#temperatures at 40-10 kPa
D.([4e4, 3e4, 2e4, 1e4])
```
"""
struct DryAdiabat <: AbstractAdiabat
    Tₛ::Float64
    Pₛ::Float64
    Pₜ::Float64
    cₚ::Float64
    μ::Float64
    Tstrat::Float64
    Ptropo::Float64
end

function DryAdiabat(Tₛ, Pₛ, cₚ, μ; Tstrat=0, Ptropo=0, Pₜ=PMIN)
    checkadiabat(Tₛ, Pₛ, Pₜ, Tstrat, Ptropo)
    DryAdiabat(Tₛ, Pₛ, Pₜ, cₚ, μ, Tstrat, Ptropo)
end

#------------------------------------

"""
[Function-like type](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) for initializing and evaluating a moist adiabatic temperature profile. Optional uniform upper atmospheric temperature below a specified temperature or pressure.

# Constructor

    MoistAdiabat(Tₛ, Pₛ, cₚₙ, cₚᵥ, μₙ, μᵥ, L, psat; Tstrat=0, Ptropo=0, N=1000, Pₜ=$PMIN)

* `Tₛ`: surface temperature [K]
* `Pₛ`: surface pressure [K]
* `cₚₙ`: specific heat of the non-condensible atmospheric component (air) [J/kg/K]
* `cₚᵥ`: specific heat of the condensible atmospheric component [J/kg/K]
* `μₙ`: molar mass of the non-condensible atmospheric component (air) [kg/mole]
* `μᵥ`: molar mass of the condensible atmospheric component [kg/mole]
* `L`: condsible component's latent heat of vaporization [J/kg]
* `psat`: function defining the saturation vapor pressure for a given temperature, `psat(T)`

If `Tstrat` is greater than zero, the temperature profile will never drop below that temperature. If `Ptropo` is greater than zero, the temperature profile at pressures lower than `Ptropo` will be equal to the temperature at exactly `Ptropo`. `Tstrat` and `Ptropo` cannot both greater than zero.

`Pₜ` defines the highest pressure [Pa] in the temperature profile. If should generally be small but cannot be zero.

The profile is evaluated along a number of pressure values in the atmosphere set by `N`. Those points are then used to construct a cubic spline interpolator for efficient and accurate temperature calculation. Experience indicates that 1000 points is very accurate and also fast.

# Example

Once constructed, use a `MoistAdiabat` like a function to compute temperature at a given pressure.

```@example
Tₛ = 288; #surface temperature [K]
Pₛ = 1e5; #surface pressure [Pa]
cₚₙ = 1040; #specific heat of air [J/kg/K]
cₚᵥ = 1996; #specific heat of H2O [J/kg/K]
μₙ = 0.029; #molar mass of air [kg/mole]
μᵥ = 0.018; #molar mass of H2O [kg/mole]
L = 2.3e6; #H2O latent heat of vaporization [J/kg]

#a saturation vapor pressure function for H2O is built in
psat = psatH2O;

#construct the moist adiabat with a tropopause pressure of 1e4 Pa
M = MoistAdiabat(Tₛ, Pₛ, cₚₙ, cₚᵥ, μₙ, μᵥ, L, psat, Ptropo=1e4);

#temperatures at 30-5 kPa
M.([3e4, 2e4, 1e4, 5e3])
```
"""
struct MoistAdiabat <: AbstractAdiabat
    ϕ::LinearInterpolator{Float64,WeakBoundaries}
    Pₛ::Float64
    Pₜ::Float64
    Tstrat::Float64
    Ptropo::Float64
end

function MoistAdiabat(Tₛ, Pₛ, cₚₙ, cₚᵥ, μₙ, μᵥ, L, psat::F;
                      Tstrat=0.0,
                      Ptropo=0.0,
                      N::Int=100,
                      Pₜ=PMIN) where {F}
    checkadiabat(Tₛ, Pₛ, Pₜ, Tstrat, Ptropo)
    #interpolation knots and output vector
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    ω = logrange(ω₁, ω₂, N)
    T = zeros(Float64, N)
    #pack the parameters
    param = (cₚₙ, cₚᵥ, 𝐑/μₙ, 𝐑/μᵥ, L, psat)
    #integrate with in-place dense output
    radau!(T, ω, dTdω, Float64(Tₛ), ω[1], ω[end], param)
    #natural spline in log pressure coordinates
    ϕ = LinearInterpolator(ω, T, WeakBoundaries())
    MoistAdiabat(ϕ, Pₛ, Pₜ, Tstrat, Ptropo)
end

#------------------------------------
#general operations with an adiabat

#direct calculation without temperature/pressure floors
temperature(Γ::DryAdiabat, P) = Γ.Tₛ*(P/Γ.Pₛ)^(𝐑/(Γ.μ*Γ.cₚ))

#direct calculation without temperature/pressure floors
temperature(Γ::MoistAdiabat, P) = Γ.ϕ(P2ω(P))

#find the pressure corresponding to a temperature, ignoring Tstrat & Ptropo
function pressure(Γ::AbstractAdiabat, T)
    Tₛ = temperature(Γ, Γ.Pₛ)
    Tₜ = temperature(Γ, Γ.Pₜ)
    @assert Tₛ >= T >= Tₜ "temperature $T K out of adiabat range [$(Tₛ),$(Tₜ)] K"
    regulafalsi((P,p) -> temperature(Γ, P) - T, Γ.Pₛ, Γ.Pₜ)
end

function (Γ::AbstractAdiabat)(P)
    #check if pressure is below tropopause
    P < Γ.Ptropo && return temperature(Γ, Γ.Ptropo)
    #what the temperature would be without any floor
    T = temperature(Γ, P)
    #apply the floor, if desired
    T < Γ.Tstrat && return Γ.Tstrat
    #ensure positive
    @assert T > 0 "non-positive temperature ($T K) encountered in adiabat at $P Pa"
    return T
end

"""
    tropopause(Γ::AbstractAdiabat)

Compute the temperature [K] and pressure [Pa] at which the tropopause occurs in an adiabatic temperature profile. This function can be called on a `DryAdiabat` or a `MoistAdiabat` if it was constructed with nonzero `Tstrat` or `Ptropo`. Returns the tuple `(T,P)`.
"""
function tropopause(Γ::AbstractAdiabat)
    Γ.Ptropo != 0 && return temperature(Γ, Γ.Ptropo), Γ.Ptropo
    Γ.Tstrat != 0 && return Γ.Tstrat, pressure(Γ, Γ.Tstrat)
    error("no stratosphere temperature or pressure has been defined (Tstrat/Ptropo)")
end

#-------------------------------------------------------------------------------
export psatH2O, tsatCO2, ozonelayer

"""
    psatH2O(T)

Compute the saturation partial pressure of water vapor at a given temperature using expressions from

* [Murphy, D. M. & Koop, T. Review of the vapour pressures of ice and supercooled water for atmospheric applications. Q. J. R. Meteorol. Soc. 131, 1539–1565 (2005).](https://rmets.onlinelibrary.wiley.com/doi/10.1256/qj.04.94)

This function uses equation 10 in the paper above when ``T >= 273.15`` K and equation 7 otherwise.
"""
function psatH2O(T)
    a = log(T)
    b = 1/T
    if T >= 273.15
        #equation 10
        c = 53.878 - 1331.22*b - 9.44523*a + 0.014025*T
        d = c*tanh(0.0415*(T - 218.8))
        P = exp(54.842763 - 6763.22*b - 4.21*a + 3.67e-4*T + d)
    else
        #equation 7
        P = exp(9.550426 - 5723.265*b + 3.53068*a - 0.00728332*T)
    end
    return P
end

"""
    tsatCO2(P)

Compute the saturation pressure of carbon dioxide at a certain pressure using equation 19 from

* [Fanale, F. P., Salvail, J. R., Bruce Banerdt, W. & Steven Saunders, R. Mars: The regolith-atmosphere-cap system and climate change. Icarus 50, 381–407 (1982)](https://doi.org/10.1016/0019-1035(82)90131-2)

The equation is inverted to express temperature as a function of pressure.
"""
function tsatCO2(P)
    @assert P <= 518000.0 "Pressure cannot be above 518000 Pa for CO2 saturation temperature"
    A = 1.2264e12 #[Pa]
    B = -3167.8 #[K]
    B/log(P/A)
end

"""
    ozonelayer(P, Cmax=8e-6)

Approximate the molar concentration of ozone in Earth's ozone layer using an 8 ppm peak at 1600 Pa which falls to zero at 100 Pa and 25500 Pa. Peak concentration is defined by `Cmax`. This approximation is discussed in

* Jacob, D. Introduction to Atmospheric Chemistry. (Princeton University Press, 1999).

"""
function ozonelayer(P, Cmax=8e-6)
    P = log(P)
    P₁ = 10.146433731146518 #ln(25500) 
    P₂ = 7.3777589082278725 #ln(1600)   
    P₃ = 4.605170185988092  #ln(100)   
    if P₂ <= P <= P₁
        return Cmax*(P₁ - P)/(P₁ - P₂)
    elseif P₃ <= P <= P₂
        return Cmax*(P - P₃)/(P₂ - P₃)
    end
    return 0.0*P
end

#-------------------------------------------------------------------------------
#miscellaneous stuff

export condensibleprofile, haircut!, rayleighCO2

"""
    condensibleprofile(Γ::AbstractAdiabat, fPₛ)

Create a function defining concentration vs pressure for a condensible with uniform upper-atmosphere (stratosphere) concentration. The new concentration profile is created with reference to an existing adiabatic profile ([`DryAdiabat`](@ref) or [`MoistAdiabat`](@ref)), which must have `Ptropo != 0` or `Tstrato != 0`. Lower atmospheric concentration is determined by the temperature dependent partial pressure function `fPₛ(T)`. The concentration is `P/(fPₛ + P)`, where `P` is the dry/non-condensible pressure.

"""
function condensibleprofile(Γ::AbstractAdiabat, fPₛ::F)::Function where {F}
    #insist on an isothermal stratosphere
    @assert ((Γ.Ptropo != 0) | (Γ.Tstrat != 0)) "adiabat must have isothermal stratosphere"
    #compute tropopause pressure and temperature
    Tₜ, Pₜ = tropopause(Γ)
    #compute saturation partial pressure at tropopause
    Pₛₜ = fPₛ(Tₜ)
    #create concentration function and return it
    let (Pₜ, Pₛₜ) = (Pₜ, Pₛₜ)
        function(T, P)
            if P >= Pₜ
                Pₛ = fPₛ(T)
                C = Pₛ/(Pₛ + P)
            else
                C = Pₛₜ/(Pₜ + Pₛₜ)
            end
            return C
        end
    end
end

#this is already exported in gases.jl
"""
    reconcentrate(G::Gas, Γ::AbstractAdiabat, fPₛ)

Create a new concentration function for a [`Gas`](@ref) and use it to [`reconcentrate`](@the gas). This does not automatically copy gas data.
"""
function reconcentrate(G::Gas, Γ::AbstractAdiabat, fPₛ)::Gas
    #construct the concentration function
    C = condensibleprofile(Γ, fPₛ)
    #assign the gas a new concentration profile
    reconcentrate(G, C)
end

export haircut!
"""
    haircut!(T, P, fTₛ)

Put a temperature floor on a temperature profile using the saturation temperature function `fTₛ(P)`.
"""
function haircut!(T, P, fTₛ::F)::Nothing where {F}
    @assert length(T) == length(P)
    for i ∈ eachindex(T)
        Tₛ = fTₛ(P[i])
        if Tₛ > T[i]
            T[i] = Tₛ
        end
    end
    nothing
end

function rayleighCO2(ν, Pₛ, g, θ)
    #convert to wavelength in micrometers
    λ = ν2λ(ν)*1e6
    #formula 2.32 from
    #  Hansen and Travis, Light Scattering in Planetary Atmospheres (1974)
    τ₀ = 1.527*(1/λ^4)*(1 + 0.013*(1/λ^2))
    #conversion from Venus conditions
    τ₀ *= 8.7/(93*𝐀)
    #final optical depth
    τ = (Pₛ/g)*τ₀

    γ = 0.75
    μ = cos(θ)
    β = 1 - exp(-τ/μ)
    f = γ*τ
    R⁻ = ((0.5 - γ*μ)*β + f)/(1 + f)
    R⁺ = f/(1 + f)
    R = 1 - (1 - R⁺)*(1 - R⁻)/((1 - R⁻))
end