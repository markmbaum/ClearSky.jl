# Gas Objects

Gas objects are high-level representations of greenhouse gases that allow fast retrieval of absorption cross-sections over continuous ranges of temperature and pressure.

## Creating Gases

Before creating a gas object, you should start your Julia session with all available threads on your system. For example, if your computer has 8 threads available, use
```shell
julia --threads 8
```

then you can define the

1. vector of wavenumbers
2. temperature and pressure ranges ([`AtmosphericDomain`](@ref))

over which its absorption cross-sections will be defined. For example,

```julia
using ClearSky
ν = LinRange(1, 2500, 2500);
Ω = AtmosphericDomain((100,350), 12, (1,1e5), 24);
```

defines 2500 evenly spaced wavenumber samples over a typical window of longwave wavenumbers and an atmospheric domain between 100-350 K and 1-1e5 Pa. The numbers 12 and 24 define the number of interpolation nodes along the temperature and pressure axes, respectively.

Now you can create a gas object directly from a `par` file containing the spectral line data from HITRAN. For example, to load carbon dioxide from the file `"CO2.par"` and assign a well-mixed concentration of 400 ppm,
```julia
co2 = WellMixedGas("CO2.par", 400e-6, ν, Ω)
```
In the background, `ClearSky` does the following
1. reads the line data
2. computes absorption cross-sections for each wavenumber, temperature, and pressure point defined by `ν` and `Ω` (using the [`voigt!`](@ref) profile by default)
3. generates high-accuracy interpolation functions for the temperature-pressure grid at each wavenumber
4. stores concentration information

Consequently, loading gases may take some time. It will be faster with more threads and with fewer wavenumber, temperature, and pressure points.

## Retrieving Cross-Sections

Gases are [function-like objects](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects). They can be used like functions to retrieve **concentration-scaled** cross-sections at any temperature and pressure within the atmospheric domain. Computing cross-sections at a specific temperature and pressure, 250 K and 10000 Pa for example, is as simple as

```julia
co2(250, 1e4)
```

This returns a vector of cross-section values [cm``^2``/molecule] at each wavenumber point the gas was created with. The cross sections are scaled by the gas molar concentration that was used when constructing the gas.

If you only need a cross-section for one of the specific wavenumber points in the gas, you must pass the index of that wavenumber before the temperature and pressure. For example, to get the cross-section corresponding to `ν[600]`,

```julia
co2(600, 250, 1e4)
```

To retrieve cross-sections that are not concentration-scaled, use the [`rawσ`](@ref) function.

## Storing Gases

Creating gas objects may take some time if you have few threads, a huge number of wavenumbers, and/or a dense temperature-pressure grid in your [`AtmosphericDomain`](@ref). To avoid loading the same gas twice, you can use Julia's built-in [serialization functions](https://docs.julialang.org/en/v1/stdlib/Serialization/) to save gases to files and quickly reload them. For example, assuming you have a gas object named `co2`, the following code will write the gas to file, then reload it.
```julia
using Serialization
#write the gas to a file called "co2"
serialize("co2", co2);
#reload the same gas from that file
co2 = deserialize("co2");
```

-----

```@docs
AtmosphericDomain
WellMixedGas
VariableGas
rawσ
concentration
reconcentrate
```
