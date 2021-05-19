# Computing Line Shapes

Line shapes are computed following the definitions and equations in [HITRAN](https://hitran.org/docs/definitions-and-units/) (and elsewhere).

In `ClearSky`, line shapes can be computed from [`SpectralLines`](@ref) objects and built in shape functions. The following shape functions are implemented:

* [`voigt`](@ref)
* [`lorentz`](@ref), pressure-broadening
* [`doppler`](@ref)
* [`PHCO2`](@ref), the Perrin & Hartman sublorentzian shape for carbon dioxide

Each of these functions has three versions/methods for computing cross-sections:

1. at a single wavenumber
2. with vectors of wavenumbers
3. in-place with vectors of wavenumbers

The methods for vectors of wavenumbers are optimized for that scenario and will generally be much faster than simply broadcasting the single-wavenumber method.

##### A Note on Handling TIPS

Evaluating line shapes requires evaluating the [temperature dependence of line intensities](https://hitran.org/docs/definitions-and-units/#mjx-eqn-eqn-intensity-temperature-dependence). To compute this scaling, the ratio of total internal partition functions (TIPS),

``Q(T_{ref})/Q(T)``

must be evaluated. The necessary information is provided by HITRAN [for each isotopologue](https://hitran.org/docs/iso-meta/) and computing the ratio requires interpolating a range of ``Q(T)`` values for the appropriate temperature.

`ClearSky` does this automatically. The [`molparam.py`](https://github.com/wordsworthgroup/ClearSky.jl/blob/main/scripts/molparam.py) script downloads Q data for each isotopologue, generates high-accuracy interpolating Chebyshev polynomials for each one, and writes an array with the information to a file in Julia syntax called [`molparam.jl`](https://github.com/wordsworthgroup/ClearSky.jl/blob/main/src/molparam.jl). These pre-computed interpolating coefficients are used directly as source code in the model, allowing rapid and accurate evaluation of the TIPS ratio. The interpolating functions are guaranteed to reproduce the provided data with less than 1 % error between 25 and 1000 K.

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
