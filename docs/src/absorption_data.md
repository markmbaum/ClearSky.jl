# Absorption Data

## Spectral Lines

The model is designed to use the [HITRAN](https://hitran.org/) database of spectral line data for various atmospheric gases. To download line data, you must register for a free account then [search](https://hitran.org/lbl/) for lines. It is generally best to download all lines for a single gas into a single file, which will have a `.par` extension. These are text files in fixed-width format and they can be read by functions in `ClearSky`.

To work with spectral line data directly, use the [`readpar`](@ref) function to load the data. This function simply parses the fixed-width file into a dictionary of vectors with the appropriate data types. For more information, see the [`readpar`](@ref) documentation.

If you plan to compute line shapes directly, read par files into [`SpectralLines`](@ref) objects. The constructor reads files using [`readpar`](@ref) then rearranges it for line shape calculations. Unnecessary information is dropped and the molecule name, formula, and molar masses are assigned. To compute line shapes, see [Computing Line Shapes](computing_line_shapes.md).

For only high-level calculations, `par` files can be loaded directly into gas objects, as described [[]]

```@docs
readpar
SpectralLines
```

## Collision Induced Absorption (CIA)

The model also makes it easy to include CIA data from HITRAN. These files can be  [downloaded directly](https://hitran.org/cia/) or all at once using the [`download_cia.py`](https://github.com/wordsworthgroup/ClearSky.jl/blob/main/scripts/download_cia.py) script. Each file contains potentially many tables of absorption data at different wavenumbers and temperatures.

Like the line data, there is a function for reading these CIA files without doing anything else. The [`readcia`](@ref) function reads a `cia` file into a vector of dictionaries. Each dictionary represents a table of absorption data. This is the raw data, but it is relatively hard to work with.

A [`CIATables`](@ref) object arranges each table of absorption data into an interpolator and makes it easy to compute the CIA absorption coefficient at any wavenumber and temperature. Also, in combination with the [`cia`](@ref) function, a [`CIATables`](@ref) can be used to compute absorption cross-sections from provided wavenumber, temperature, and partial pressures.

-----

```@docs
readcia
CIATables
```
