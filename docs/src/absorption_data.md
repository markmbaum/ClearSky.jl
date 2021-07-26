# Absorption Data

## Spectral Lines

The model is designed to use the [HITRAN](https://hitran.org/) database of spectral line data for various atmospheric gases. To download line data, you must register for a free account then [search](https://hitran.org/lbl/) for lines. It's generally best to download all lines for a single gas into a single file, which will have a `.par` extension. These are text files in fixed-width format and they can be read by functions in `ClearSky`.

To work with spectral line data directly, use the [`readpar`](@ref) function to load the data. This function simply parses the fixed-width file into a dictionary of vectors with the appropriate data types. For more information, see the [`readpar`](@ref) documentation.

If you want to compute line shapes, read par files into [`SpectralLines`](@ref) objects. The constructor reads files using [`readpar`](@ref) then does some rearranging for fast line shape calculations. Unnecessary information is dropped and the molecule name, formula, and molar masses are assigned. To compute line shapes, see [Line Shapes](line_shapes.md).

For only high-level calculations, `par` files can also be loaded directly into [gas objects](gas_objects.md).

```@docs
readpar
SpectralLines
```

## Collision Induced Absorption (CIA)

The model also makes it easy to include CIA data from HITRAN. These files can be  [downloaded directly](https://hitran.org/cia/) or all at once using the [`download_cia.py`](https://github.com/markmbaum/ClearSky.jl/blob/main/scripts/download_cia.py) script. Each file contains potentially many tables of absorption data at different wavenumbers and temperatures.

Like the line data, there is a function for reading these CIA files without doing anything else. The [`readcia`](@ref) function reads a `cia` file into a vector of dictionaries. Each dictionary represents a table of absorption data. This is the raw data, but it is relatively hard to work with.

A [`CIATables`](@ref) object arranges each table of absorption data into an interpolator and makes it easy to compute the CIA absorption coefficient at any wavenumber and temperature.

A [`CIA`](@ref) struct links a [`CIATables`](@ref) object with the gasses that are inducing absorption. Because concentrations are defined in the gas objects, [`CIA`](@ref) objects can be used to compute induced absorption cross-sections using only the wavenumber, temperature, and pressure values.

The [`cia`](@ref) function also provides a number of methods for computing induced cross-sections.

-----

```@docs
readcia
CIATables
CIA
cia
cia!
```
