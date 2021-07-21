# Line Shapes

Line shapes are computed following the definitions and equations in [HITRAN](https://hitran.org/docs/definitions-and-units/) (and elsewhere).

In `ClearSky`, line shapes can be computed from [`SpectralLines`](@ref) objects and built in shape functions. The following shape functions are implemented with the necessary supporting functions:

* [`voigt`](@ref)
* [`lorentz`](@ref), pressure-broadening
* [`doppler`](@ref)
* [`PHCO2`](@ref), the Perrin & Hartman sublorentzian shape for carbon dioxide

Each of these functions has three methods for computing cross-sections:

1. at a single wavenumber, temperature, and pressure
2. over a sorted vector of wavenumbers, a single temperature, and a single pressure
3. the same as item 2, but computing cross-sections in-place

Using multiple dispatch, the arguments supplied to the functions determine the behavior. For example, calling the [`voigt`](@ref) function with a single number for the `ν` argument executes the function for a single cross-section. Calling the same function with a vector of wavenumbers in the `ν` argument executes the version of the function optimized for that scenario.

Another way to get cross-sections is through [gas objects](gas_objects.md), which are used for higher-level modeling.

##### A Note on Handling TIPS

Evaluating line shapes requires evaluating the [temperature dependence of line intensities](https://hitran.org/docs/definitions-and-units/#mjx-eqn-eqn-intensity-temperature-dependence). To compute this scaling, the ratio of total internal partition functions (TIPS),

``Q(T_{ref})/Q(T)``

must be evaluated. The necessary information is provided by HITRAN [for every isotopologue](https://hitran.org/docs/iso-meta/) and computing the ratio requires interpolating a range of ``Q(T)`` values for the appropriate temperature.

`ClearSky` evaluates the ratio accurately and automatically inside the [`scaleintensity`](@ref) function.

To facilitate this, the [`molparam.py`](https://github.com/markmbaum/ClearSky.jl/blob/main/scripts/molparam.py) script was used to download Q data for each isotopologue, generate high-accuracy interpolating Chebyshev polynomials for each one, and write the information to a Julia source file called [`molparam.jl`](https://github.com/markmbaum/ClearSky.jl/blob/main/src/molparam.jl). The pre-computed interpolating coefficients are defined directly in source code, allowing rapid and accurate evaluation of the TIPS ratio. The interpolating functions are guaranteed to reproduce the provided data with less than 1 % error between 25 and 1000 K.

-----

## Voigt Profile

```@docs
fvoigt
voigt
voigt!
```

-----

## Lorentz Profile

```@docs
γlorentz
florentz
lorentz
lorentz!
```

-----

## Doppler Profile

```@docs
αdoppler
fdoppler
doppler
doppler!
```

-----

## Perrin & Hartman Sublorentzian CO2 Profile

```@docs
ΧPHCO2
PHCO2
PHCO2!
```

-----

## Other

```@docs
scaleintensity
```
