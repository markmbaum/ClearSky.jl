# Atmospheric Profiles

`ClearSky` contains functions for common atmospheric profiles and quantities.

## Temperature Profiles

Many radiative transfer calculations require an atmospheric temperature profile. `ClearSky` is designed to work with any arbitrary temperature profile if it can be defined as a function of pressure, `T(P)`.

For convenience, dry and moist adiabatic profiles are available through the [`DryAdiabat`](@ref) and [`MoistAdiabat`](@ref) types, which are [function-like types](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects). There is also a [`tropopause`](@ref) function.

```@docs
DryAdiabat
MoistAdiabat
tropopause
```

## Pressure Profiles

In case a pressure profile with constant scale height isn't sufficient, hydrostatic profiles with arbitrary temperature and mean molar mass functions are available through the [`Hydrostatic`](@ref) type and related functions.

```@docs
Hydrostatic
hydrostatic
scaleheight
altitude
```

## Other Functions

```@docs
psatH2O
tsatCO2
ozonelayer
condensibleprofile
haircut!
rayleighCO2
```
