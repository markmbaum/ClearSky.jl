#-------------------------------------------------------------------------------
#conversions between spectral units

export ν2f, f2ν, ν2λ, λ2ν, λ2f, f2λ

"""
Convert wavenumber [cm``^{-1}``] to frequency [1/s]
"""
ν2f(ν) = 100.0*𝐜*ν

"""
Convert frequency [1/s] to wavenumber [cm``^{-1}``]
"""
f2ν(f) = f/(100.0*𝐜)

"""
Convert wavenumber [cm``^{-1}``] to wavelength [m]
"""
ν2λ(ν) = 0.01/ν

"""
Convert wavelength [m] to wavenumber [cm``^{-1}``]
"""
λ2ν(λ) = 0.01/λ

"""
Convert wavelength [m] to frequency [1/s]
"""
λ2f(λ) = 𝐜/λ

"""
Convert frequency [1/s] to wavelength [m]
"""
f2λ(f) = f/𝐜

#-------------------------------------------------------------------------------
export planck, normplanck, dplanck, stefanboltzmann, equilibriumtemperature

"""
    planck(ν, T)

Compute black body intensity [W/m``^2``/cm``^{-1}``/sr] using [Planck's law](https://en.wikipedia.org/wiki/Planck%27s_law)

# Arguments
* `ν`: wavenumger [cm``^{-1}``]
* `T`: temperature [Kelvin]
"""
function planck(ν, T)
    νₘ = 100.0*ν #convert from cm^-1 to m^-1
    x = 𝐡*𝐜*νₘ/(𝐤*T) #exponent
    p = 2*𝐡*𝐜^2*νₘ^3 #prefactor
    #result, converting back to cm^-1 with the factor of 100
    100.0*p/(exp(x) - 1.0)
end

"""
    normplanck(ν, T)

Compute black body intensity [W/m``^2``/cm``^{-1}``/sr] using [Planck's law](https://en.wikipedia.org/wiki/Planck%27s_law), normalized by the power emitted per unit area at the given temperature ([`stefanboltzmann`](@ref)),

```
B(ν,T)/σT^4
```

yielding units of 1/cm``^{-1}``/sr.

# Arguments
* `ν`: wavenumger [cm``^{-1}``]
* `T`: temperature [Kelvin]
"""
normplanck(ν, T) = planck(ν, T)/stefanboltzmann(T)

"""
    dplanck(ν, T)

Evaluates ``\\frac{\\partial B}{\\partial T}``
"""
function dplanck(ν, T)
    νₘ = 100.0*ν #convert from cm^-1 to m^-1
    x = 𝐡*𝐜*νₘ/(𝐤*T) #exponent
    y = exp(x)
    if isinf(y)
        return zero(y) #fail accurately, Inf/(Inf - 1)^2 → 0
    end
    p = 2*𝐡^2*𝐜^3*νₘ^4/(𝐤*T^2) #prefactor
    #result, converting back to cm^-1 with the factor of 100
    100.0*p*y/(y - 1)^2
end

"""
    stefanboltzmann(T)

Compute black body radiation power using the [Stefan-Boltzmann](https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law) law, ``σT^4`` [W/m``^2``].
"""
stefanboltzmann(T) = 𝛔*(T^4)

"""
    equilibriumtemperature(F, A)

Compute the [planetary equilibrium temperature](https://en.wikipedia.org/wiki/Planetary_equilibrium_temperature), or equivalent blackbody temperature of a planet.

``(\\frac{(1 - A)F}{4\\sigma})^{1/4}``

# Arguments
* `F`: stellar flux [W/m``^2``]
* `A`: albedo
"""
equilibriumtemperature(F, A) = ((1 - A)*F/(4*𝛔))^(1/4)

"""
    equilibriumtemperature(L, A, R)

Compute the [planetary equilibrium temperature](https://en.wikipedia.org/wiki/Planetary_equilibrium_temperature), or equivalent blackbody temperature of a planet.

``(\\frac{(1 - A)L}{16 \\sigma \\pi R^2})^{1/4}``

# Arguments
* `L`: stellar luminosity [W]
* `A`: albedo
* `R`: orbital distance [m]
"""
equilibriumtemperature(L, A, R) = (L*(1 - A)/(16*𝛔*π*R^2))^(1/4)

#-------------------------------------------------------------------------------
export dτdP, transmittance, schwarzschild, absorption, emission

"""
    dτdP(σ, g, μ)

Evaluate the differential increase in optical depth in pressure coordinates, equivalent to the [`schwarzschild`](@ref) equation without Planck emission.

``\\frac{dτ}{dP} = σ\\frac{\\textrm{N}_A}{g μ}``

where ``N_A`` is Avogadro's number.

# Arguments
* `σ`: absorption cross-section [cm``^2``/molecule]
* `g`: gravitational acceleration [m/s``^2``]
* `μ`: mean molar mass [kg/mole]
"""
dτdP(σ, g, μ) = 1e-4*σ*𝐍𝐚/(μ*g)

"""
    transmittance(τ)

Evaluate transmittance from optical depth, ``t = e^{-τ}``
"""
transmittance(τ) = exp(-τ)

"""
    schwarzschild(I, ν, σ, T, P)

Evaluate the [Schwarzschild differential equation](https://en.wikipedia.org/wiki/Schwarzschild%27s_equation_for_radiative_transfer) for radiative transfer with units of length/height [m] and assuming the ideal gas law.

``\\frac{dI}{dz} = σ\\frac{P}{k_B T}[B_ν(T) - I]``

where ``B_ν`` is [`planck`](@ref)'s law.

# Arguments
* `I`: radiative intensity [W/m``^2``/cm``^{-1}``/sr]
* `ν`: radiation wavenumber [cm``^{-1}``]
* `σ`: absorption cross-section [cm``^2``/molecule]
* `T`: temperature [K]
* `P`: pressure [Pa]
"""
schwarzschild(I, ν, σ, T, P) = 1e-4*σ*(P/(𝐤*T))*(planck(ν,T) - I)

"""
    schwarzschild(I, ν, σ, g, μ, T)

Evaluate the [Schwarzschild differential equation](https://en.wikipedia.org/wiki/Schwarzschild%27s_equation_for_radiative_transfer) for radiative transfer with pressure units [Pa] and assuming the ideal gas law.

``\\frac{dI}{dP} = σ\\frac{\\textrm{N}_A}{g μ}[B_ν(T) - I]``

where ``B_ν`` is [`planck`](@ref)'s law and ``N_A`` is Avogadro's number.

# Arguments
* `I`: radiative intensity [W/m``^2``/cm``^{-1}``/sr]
* `ν`: radiation wavenumber [cm``^{-1}``]
* `σ`: absorption cross-section [cm``^2``/molecule]
* `g`: gravitational acceleration [m/s``^2``]
* `μ`: mean molar mass [kg/mole]
* `T`: temperature [K]
"""
schwarzschild(I, ν, σ, g, μ, T) = 1e-4*σ*(𝐍𝐚/(μ*g))*(planck(ν,T) - I)

"""
    absorption(I, σ, g, μ)

Evaluate the [Schwarzschild differential equation](https://en.wikipedia.org/wiki/Schwarzschild%27s_equation_for_radiative_transfer) for radiative transfer with pressure units [Pa] and assuming the ideal gas law, but **without Planck emission**. This can be useful if the wavenumber is far away from the emission region of the temperature and the Planck emission is known to be negligible.

``\\frac{dI}{dP} = -σ\\frac{\\textrm{N}_A}{g μ}I``

where ``N_A`` is Avogadro's number.

# Arguments
* `I`: radiative intensity [W/m``^2``/cm``^{-1}``/sr]
* `σ`: absorption cross-section [cm``^2``/molecule]
* `g`: gravitational acceleration [m/s``^2``]
* `μ`: mean molar mass [kg/mole]
"""
absorption(I, σ, g, μ) = -1e-4*σ*(𝐍𝐚/(μ*g))*I

"""
    emission(I, ν, σ, g, μ, T)

Evaluate the [Schwarzschild differential equation](https://en.wikipedia.org/wiki/Schwarzschild%27s_equation_for_radiative_transfer) for radiative transfer with pressure units [Pa] and assuming the ideal gas law, but **without Planck emission**. This can be useful if the wavenumber is far away from the emission region of the temperature and the Planck emission is known to be negligible.

``\\frac{dI}{dP} = -σ\\frac{\\textrm{N}_A}{g μ}B_ν(T)``

where ``B_ν`` is [`planck`](@ref)'s law

# Arguments
* `I`: radiative intensity [W/m``^2``/cm``^{-1}``/sr]
* `σ`: absorption cross-section [cm``^2``/molecule]
* `g`: gravitational acceleration [m/s``^2``]
* `μ`: mean molar mass [kg/mole]
"""
emission(ν, σ, g, μ, T) = 1e-4*σ*(𝐍𝐚/(μ*g))*planck(ν,T)

