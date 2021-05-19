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
scaleheight(g, μ, T)::Float64 = 𝐑*T/(μ*g)

#parameterized hydrostatic relation in log coordinates
function dlnPdz(z, lnP, param::Tuple)::Float64
    #unpack parameters
    Pₛ, g, fT, fμ = param
    #evaluate temperature and mean molar mass [kg/mole]
    P = exp(lnP)
    if P < PMIN
        return 0.0
    end
    P = min(P, Pₛ) #don't allow tiny numerical dips below Pₛ
    T = fT(P)
    μ = fμ(T, P)
    #evaluate derivative
    -μ*g/(𝐑*T)
end

"""
    hydrostatic(z, Pₛ, g, fT, fμ)

Compute the hydrostatic pressure [Pa] at a specific altitude using arbitrary atmospheric profiles of temperature of mean molar mass

# Arguments

* `z`: altitude [m] to compute pressure at
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature as a function of pressure, `fT(P)`
* `fμ`: mean molar mass as a function of pressure and temperature `fμ(T,P)`
"""
function hydrostatic(z, Pₛ, g, fT::T, fμ::U)::Float64 where {T,U}
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

Compute the altitude [m] at which a specific hydrostatic pressure occurs using arbitrary atmospheric profiles of temperature of mean molar mass

# Arguments

* `P`: pressure [Pa] to compute altitude at
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature as a function of pressure, `fT(P)`
* `fμ`: mean molar mass as a function of pressure and temperature `fμ(T,P)`
"""
function altitude(P, Pₛ, g, fT::T, fμ::U)::Float64 where {T,U}
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
    falseposition(z -> log(hydrostatic(z, Pₛ, g, fT, fμ)) - log(P), z₁, z₂)
end

"""
[Function-like type](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) for initializing and evaluating a hydrostatic pressure profile with arbitrary temperature and mean molar mass profiles. A `Hydrostatic` object maps altitude to pressure.

# Constructor

    Hydrostatic(Pₛ, Pₜ, g, fT, fμ)

* `Pₛ`: surface pressure [Pa]
* `Pₜ`: top of profile pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature as a function of presssure, `fT(P)`
* `fμ`: mean molar mass as a function of temperature and pressure, `fμ(T,P)`

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
    ϕ::LinearInterpolator
    zₜ::Float64
end

function Hydrostatic(Pₛ, Pₜ, g, fT::T, fμ::U, N::Int=1000) where {T,U}
    #find the altitude corresponding to Pₜ
    zₜ = altitude(Pₜ, Pₛ, g, fT, fμ)
    #interpolation knots and output array
    z = logrange(0, zₜ, N)
    lnP = zeros(Float64, N)
    #integrate to get a full pressure profile
    radau!(lnP, z, dlnPdz, log(Pₛ), 0, zₜ, (Pₛ, g, fT, fμ))
    #construct and return
    Hydrostatic(LinearInterpolator(z, lnP), zₜ)
end

function (H::Hydrostatic)(z)::Float64
    @assert 0 <= z <= H.zₜ "altitude $z out of range [0,$(H.zₜ)]"
    exp(H.ϕ(z, false))
end

"""
    altitude(H::Hydrostatic, P)

Compute the altitude at which a specific pressure occurs in a hydrostatic pressure profile.
"""
altitude(H::Hydrostatic, P)::Float64 = falseposition(z -> log(H(z)) - log(P), 0.0, H.zₜ)

#-------------------------------------------------------------------------------

#general function for adiabat with one condensable in bulk non-condensable
function dTdω(ω, T, cₚₙ, cₚᵥ, Rₙ, Rᵥ, L, psat::F)::Float64 where {F}
    #molar mixing ratio of condensible
    α = psat(T)/ω2P(ω)
    #whole expression at once
    -T*(Rₙ/cₚₙ)*(1 + α*L/(Rₙ*T))/(1 + α*(cₚᵥ/cₚₙ + (L/(T*Rᵥ) - 1)*L/(cₚₙ*T)))
end

dTdω(ω, T, param::Tuple)::Float64 = dTdω(ω, T, param...)

#------------------------------------

abstract type AbstractAdiabat end

export MoistAdiabat, DryAdiabat
export tropopause

#------------------------------------

"""
[Function-like type](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) for initializing and evaluating a dry adiabatic temperature profile. Optional uniform upper atmospheric temperature below a specified temperature or pressure.

# Constructor

    DryAdiabat(Tₛ, Pₛ, cₚ, μ, Pₜ=$PMIN; Tstrat=0.0, Ptropo=0.0)

* `Tₛ`: surface temperature [K]
* `Pₛ`: surface pressure [K]
* `Pₜ`: highest allowable pressure (can be small but not zero) [Pa]
* `cₚ`: specific heat of the atmosphere [J/kg/K]
* `μ`: molar mass of the atmosphere [kg/mole]
* `Pₜ`: highest pressure in the temperature profile (should generally be small to avoid evaluating out of range) [Pa]

If `Tstrat` is greater than zero, the temperature profile will not drop below that temperature. If `Ptropo` is greater than zero, the temperature profile at pressures lower than `Ptropo` will be equal to the temperature at exactly `Ptropo`. `Tstrat` and `Ptropo` cannot be greater than zero simultaneously.

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

function DryAdiabat(Tₛ, Pₛ, cₚ, μ, Pₜ=PMIN; Tstrat=0, Ptropo=0)
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    @assert Pₜ > 0 "Pₜ must be greater than 0"
    DryAdiabat(Tₛ, Pₛ, Pₜ, cₚ, μ, Tstrat, Ptropo)
end

#------------------------------------

"""
[Function-like type](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) for initializing and evaluating a moist adiabatic temperature profile. Optional uniform upper atmospheric temperature below a specified temperature or pressure.

# Constructor

    MoistAdiabat(Tₛ, Pₛ, cₚₙ, cₚᵥ, μₙ, μᵥ, L, psat, Pₜ=$PMIN; Tstrat=0, Ptropo=0, N=1000)

* `Tₛ`: surface temperature [K]
* `Pₛ`: surface pressure [K]
* `cₚₙ`: specific heat of the non-condensible atmospheric component (air) [J/kg/K]
* `cₚᵥ`: specific heat of the condensible atmospheric component [J/kg/K]
* `μₙ`: molar mass of the non-condensible atmospheric component (air) [kg/mole]
* `μᵥ`: molar mass of the condensible atmospheric component [kg/mole]
* `L`: condsible component's latent heat of vaporization [J/kg]
* `psat`: function defining the saturation vapor pressure for a given temperature, `psat(T)`
* `Pₜ`: highest pressure in the temperature profile (should generally be small to avoid evaluating pressures out of range) [Pa]

If `Tstrat` is greater than zero, the temperature profile will not drop below that temperature. If `Ptropo` is greater than zero, the temperature profile at pressures lower than `Ptropo` will be equal to the temperature at exactly `Ptropo`. `Tstrat` and `Ptropo` cannot be greater than zero simultaneously.

The profile is evaluated along a number of pressure values in the atmosphere set by `N`. REVISE. Those points are then used to construct a cubic spline interpolator for efficient and accurate temperature calculation. Experience indicates that 1000 points is very accurate and also fast.

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
    ϕ::LinearInterpolator
    Pₛ::Float64
    Pₜ::Float64
    Tstrat::Float64
    Ptropo::Float64
end

function MoistAdiabat(Tₛ, Pₛ, cₚₙ, cₚᵥ, μₙ, μᵥ, L, psat::F,
                      Pₜ=PMIN;
                      Tstrat=0.0,
                      Ptropo=0.0,
                      N::Int=1000) where {F<:Function}
    if (Tstrat != 0) & (Ptropo != 0)
        throw("Cannot have nonzero Tstrat and Ptropo, must use one or the other")
    end
    @assert Pₛ > Pₜ "Pₛ must be greater than Pₜ"
    @assert Pₜ > 0 "Pₜ must be greater than 0"
    #interpolation knots and output vector
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    ω = logrange(ω₁, ω₂, N)
    T = zeros(Float64, N)
    #pack the parameters
    param = (cₚₙ, cₚᵥ, 𝐑/μₙ, 𝐑/μᵥ, L, psat)
    #integrate with in-place dense output
    radau!(T, ω, dTdω, Tₛ, ω[1], ω[end], param)
    #natural spline in log pressure coordinates
    MoistAdiabat(LinearInterpolator(ω, T), Pₛ, Pₜ, Tstrat, Ptropo)
end

#------------------------------------
#general operations with an adiabat

#direct calculation
temperature(Γ::DryAdiabat, P)::Float64 = Γ.Tₛ*(P/Γ.Pₛ)^(𝐑/(Γ.μ*Γ.cₚ))

#coordinate conversion and interpolation
temperature(Γ::MoistAdiabat, P)::Float64 = Γ.ϕ(P2ω(P), false)

#find the pressure corresponding to a temperature (ignores Tstrat, Ptropo)
function pressure(Γ::AbstractAdiabat, T)::Float64
    Tₛ = Γ(Γ.Pₛ)
    Tₜ = Γ(Γ.Pₜ)
    @assert Tₛ >= T >= Tₜ "temperature $T K out of adiabat range [$(Tₛ),$(Tₜ)] K"
    falseposition(P -> Γ(P) - T, Γ.Pₛ, Γ.Pₜ)
end

function (Γ::AbstractAdiabat)(P)::Float64
    #check bounds
    if (P < Γ.Pₜ) && !(P ≈ Γ.Pₜ)
        throw("Adiabat defined within $(Γ.Pₜ) and $(Γ.Pₛ) Pa, $P Pa is too low.")
    end
    if (P > Γ.Pₛ) && !(P ≈ Γ.Pₛ)
        throw("Adiabat defined within $(Γ.Pₜ) and $(Γ.Pₛ) Pa, $P Pa is too high.")
    end
    #check if pressure is below tropopause
    if P < Γ.Ptropo
        return temperature(Γ, Γ.Ptropo)
    end
    #what the temperature would be without any floor
    T = temperature(Γ, P)
    #apply the floor, if desired
    if T < Γ.Tstrat
        return Γ.Tstrat
    end
    #ensure positive
    @assert T > 0 "non-positive temperature ($T K) encountered in adiabat at $P Pa"
    return T
end

"""
    tropopause(Γ::AbstractAdiabat)

Compute the temperature and pressure at which the tropopause occurs in an adiabatic temperature profile. This function can be called on a `DryAdiabat` or `MoistAdiabat` if it was constructed with nonzero `Tstrat` or `Ptropo`. Returns a tuple, `(T,P)`.
"""
function tropopause(Γ::AbstractAdiabat)::Tuple{Float64,Float64}
    if Γ.Ptropo != 0
        return temperature(Γ, Γ.Ptropo), Γ.Ptropo
    end
    if Γ.Tstrat != 0
        return Γ.Tstrat, pressure(Γ, Γ.Tstrat)
    end
    throw("no stratosphere temperature or pressure has been defined (Tstrat/Ptropo)")
end

#-------------------------------------------------------------------------------

export psatH2O, tsatCO2, ozonelayer

"""
    psatH2O(T)

Compute the saturation partial pressure of water vapor at a certain temperature using expressions from

* Murphy, D. M. & Koop, T. Review of the vapour pressures of ice and supercooled water for atmospheric applications. Q. J. R. Meteorol. Soc. 131, 1539–1565 (2005).
"""
function psatH2O(T)::Float64
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

Compute the saturation pressure of carbon dioxide at a certain pressure using expression from Fanale et al. (1982)
"""
function tsatCO2(P)::Float64
    @assert P <= 518000.0 "Pressure cannot be above 518000 Pa for CO2 saturation temperature"
    -3167.8/(log(0.01*P) - 23.23)
end

"""
    ozonelayer(P)

Approximate the molar concentration of ozone in Earth's ozone layer using an 8 ppm peak at 1600 Pa
"""
function ozonelayer(P)::Float64 #yields molar concentration with peak at 8 ppm
    P = log(P)
    P₁ = 10.146433731146518 #ln(25500)
    P₂ = 7.3777589082278725 #ln(1600)
    P₃ = 4.605170185988092  #ln(100)
    C₂ = 8e-6
    if P₂ <= P <= P₁
        return C₂*(P₁ - P)/(P₁ - P₂)
    elseif P₃ <= P <= P₂
        return C₂*(P - P₃)/(P₂ - P₃)
    end
    return 0
end
