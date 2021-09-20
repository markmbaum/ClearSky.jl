#-------------------------------------------------------------------------------
# reading and constructing CIA

export readcia, CIATables

function nextspace(s::String, i::Int, L::Int)
    while (i <= L) && (s[i] != ' ')
        i += 1
    end
    return i
end

function nextnonspace(s::String, i::Int, L::Int)
    while (i <= L) && (s[i] == ' ')
        i += 1
    end
    return i
end

"""
    readcia(filename)

Read a collision induced absorption data file. These files are available from [HITRAN](https://hitran.org/cia/) and desribed by [this readme](https://hitran.org/data/CIA/CIA_Readme.pdf). A vector of dictionaries is returned. Each dictionary represents absorption coefficients for a single temperature over a range of wavenumbers. Each dictionary contains:

| Key | Type | Description |
| --- | :--- | :---------- |
| `symbol` | `String` | chemical symbol |
| `νmin` | `Float64` | minimum wavenumber of absorption range [cm``^{-1}``] |
| `νmax` | `Float64` | maximum wavenumber of absorption range [cm``^{-1}``] |
| `npts` | `Int64` | number of points
| `T` | `Float64` | temperature for absorption data [K] |
| `ν` | `Vector{Float64}` | wavenumber samples [cm``^{-1}``]
| `k` | `Vector{Float64}` | absorption coefficients [cm``^5``/molecule``^2``]
| `maxcia` | `Float64` | maximum absorption coefficient [cm``^5``/molecule``^2``]
| `res` | `Float64` | ? |
| `comments` | `String` | miscelleneous comments |
| `reference` | `Int64` | indices of data references |
"""
function readcia(filename::String)
    @assert filename[end-3:end] == ".cia" "expected file with .cia extension downloaded from https://hitran.org/cia/"
    lines = readlines(filename)
    #lengths of all lines
    L = length.(lines)
    @assert maximum(L) == 100 "unexpected maximum line length in cia file, expected 100 but got $(maximum(L))"
    #find locations of header lines
    hidx = findall(len->len == 100, L)
    #allocate a big array of dictionaries for each range of data
    data = Vector{Dict{String,Any}}(undef, length(hidx))
    #parse the tables
    push!(hidx, length(lines) + 1)
    for i = 1:length(hidx) - 1
        #start the dictionary for this table
        data[i] = Dict{String,Any}()
        #line indices for the ith table
        ia = hidx[i]
        ib = hidx[i+1]
        #parse the header values
        line = lines[ia]
        data[i]["symbol"]    = strip(line[1:20])
        data[i]["νmin"]      = parse(Float64, line[21:30])
        data[i]["νmax"]      = parse(Float64, line[31:40])
        data[i]["npts"]      = parse(Int64, line[41:47])
        data[i]["T"]         = parse(Float64, line[48:54])
        data[i]["maxcia"]    = parse(Float64, line[55:64])
        data[i]["res"]       = parse(Float64, line[65:70])
        data[i]["comments"]  = strip(line[71:97])
        data[i]["reference"] = parse(Int64, line[98:100])
        #read the data columns
        table = lines[ia+1:ib-1]
        L = length(table)
        ν = zeros(Float64, L)
        k = zeros(Float64, L)
        for j = 1:L
            #string representing row of table
            line = table[j]
            N = length(line)
            #index of first non-space character
            na = nextnonspace(line, 1, N)
            #next space
            nb = nextspace(line, na, N)
            #and the next non-space
            nc = nextnonspace(line, nb, N)
            #finally, the next space or line end
            nd = nextspace(line, nc, N)
            #parse the values into arrays
            ν[j] = parse(Float64, line[na:nb-1])
            k[j] = parse(Float64, line[nc:nd-1])
        end
        #add to the dictionary
        data[i]["ν"] = ν
        data[i]["k"] = k
    end
    return data
end

"""
Organizing type for collision induced absorption data, with data tables loaded into [interpolators](https://markmbaum.github.io/BasicInterpolators.jl/dev/).

| Field | Type | Description |
| ----- | :--- | :---------- |
| `name` | `String` | molecular symbol, i.e. `"CO2-H2"` |
| `formulae` | `Tuple{String,String}` | split molecular formulae, i.e `("CO2", "H2")` |
| `Φ` | `Vector{BilinearInterpolator}` | interpolators for each grid of absorption coefficients |
| `ϕ` | `Vector{LinearInterpolator}` | interpolators for isolated ranges of absorption coefficients |
| `T` | `Vector{Float64}` | temperatures [K] for single ranges in `ϕ` |
| `extrapolate` | `Bool` | whether to extrapolate using flat boundaries from the coefficient grids in `Φ` |
| `singles` | `Bool` | whether to use the single ranges in `ϕ` at all |

The interpolator objects are described in the [`BasicInterpolators.jl`](https://markmbaum.github.io/BasicInterpolators.jl/dev/) documentation.

# Constructors

    CIATables(cia::Vector{Dict}; extrapolate=false, singles=false, verbose=true)

Construct a `CIATables` object from a dictionary of coefficient data, which can be read from files using [`readcia`](@ref). Keywords `extrapolate` and `singles` are used to set those fields of the returned object.

    CIATables(filename; extrapolate=false, singles=false, verbose=true)

Construct a `CIATables` object directly from file, using [`readcia`](@ref) along the way.

# Examples
A `CIATables` object is [function-like](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects). To retrieve the absorption coefficient at a given wavenumber [cm``^{-1}``] and temperature [K], use the object like a function.

```julia
co2co2 = CIATables("data/cia/CO2-CO2_2018.cia"); #read data
ν = 100; #wavenumber [cm^-1]
T = 288; #temperature [K]
k = co2co2(ν, T) #absorption coefficient [cm^5/molecule^2]
```

The object interpolates and sums all data tables that contain `ν` and `T`. If `extrapolate` is `true`, boundary values are included whenever the temperature is out of range. If `singles` is `true`, data ranges for a single temperature are included whenever they contain `ν`.

A `CIATables` can be passed to the [`cia`](@ref) function to compute an absorption cross-section with different temperatures and pressures.

```julia
co2ch4 = CIATables("data/cia/CO2-CH4_2018.cia"); #read data
ν = 250; #wavenumber [cm^-1]
T = 310 #temperature [K]
Pa = 1e5; #air pressure [Pa]
Pco2 = 40; #CO2 partial pressure [Pa]
Pch4 = 0.1; #CH4 partial pressure [Pa]
σ = cia(ν, co2ch4, T, Pa, Pco2, Pch4) #absorption cross-section [cm^2/molecule]
```
"""
struct CIATables
    #molecular symbol (CO2-CO2, H2-H2, etc.)
    name::String
    #a Set object of the two gases involved
    formulae::Tuple{String,String}
    #interpolation structs
    Φ::Vector{BilinearInterpolator{Float64,NoBoundaries}}
    ϕ::Vector{LinearInterpolator{Float64,NoBoundaries}}
    #temperatures for singles in ϕ
    T::Vector{Float64}
    #whether to extrapolate (!)
    extrapolate::Bool
    #whether to use single-temperature CIATables ranges
    singles::Bool
end

function CIATables(data::Vector{Dict{String,Any}};
                   extrapolate::Bool=false,
                   singles::Bool=false,
                   verbose::Bool=true)
    verbose && println("creating CIATables")
    #pull out wavenumber ranges and temperatures for each grid
    νmin = map(x->x["νmin"], data)
    νmax = map(x->x["νmax"], data)
    T = map(x->x["T"], data)
    #select unique wavenumber ranges and sort 'em
    νranges = sort(unique(zip(νmin, νmax)), by=x->x[1])
    n = length(νranges)
    #now get the interpolation tables and temperature ranges for each νrange
    Φ = Vector{BilinearInterpolator{Float64,NoBoundaries}}()
    ϕ = Vector{LinearInterpolator{Float64,NoBoundaries}}()
    τ = Vector{Float64}()
    for i = 1:n
        νmin[i], νmax[i] = νranges[i]
        #pull out the data for this wavenumber range
        idx = findall(x->(x["νmin"] ≈ νmin[i]) & (x["νmax"] ≈ νmax[i]), data)
        T = map(x->x["T"], data[idx])
        ν = map(x->x["ν"], data[idx])
        k = map(x->x["k"], data[idx])
        #if there is only one range, make a linear interpolation
        if length(T) == 1
            ν, k = ν[1], k[1]
            k[k .<= 0.0] .= 0.0
            push!(ϕ, LinearInterpolator(ν, log.(k), NoBoundaries()))
            push!(τ, T[1])
            if verbose
                println("  ? single temperature CIA range found at $(T[1]) K, $(minimum(ν)) - $(maximum(ν)) cm^-1")
            end
        else
            #assert that all ν vectors are identical, then take one of them
            for j = 2:length(ν)
                @assert sum(ν[1] .- ν[j]) ≈ 0.0 "wavenumber sample within a wavenumber range appear to be different"
            end
            ν = ν[1]
            #ensure sorted by temperature
            idx = sortperm(T)
            T, k = T[idx], k[idx]
            #put the k values into a 2d array
            k = convert.(Float64, hcat(k...))
            #replace bizarre negative values with tiny values
            k[k .<= 0.0] .= floatmin(Float64)
            #construct an interpolator WITH THE LOG OF K for accuracy
            push!(Φ, BilinearInterpolator(ν, T, log.(k), NoBoundaries()))
        end
    end
    #make sure symbols are all the same and get the individual gas strings
    symbols = unique(map(x->x["symbol"], data))
    @assert length(symbols) == 1
    symbol = symbols[1]
    formulae = Tuple(map(String, split(symbol, '-')))
    #construct
    tables = CIATables(symbol, formulae, Φ, ϕ, τ, extrapolate, singles)
    #print some info if desired
    if verbose
        println("  formulae: $(tables.formulae[1]) & $(tables.formulae[2])")
        Φ = tables.Φ
        n = length(Φ)
        ϕ = tables.ϕ
        m = length(ϕ)
        println("  $(n+m) absorption region(s)")
        for i = 1:n
            println("    $i) ν = $(Φ[i].G.xa) - $(Φ[i].G.xb) cm^-1")
            println("       T = $(Φ[i].G.ya) - $(Φ[i].G.yb) K")
        end
        for i = 1:m
            println("    $(i+n)) ν = $(ϕ[i].r.xa) - $(ϕ[i].r.xb) cm^-1")
            println("       T = $(tables.T[i]) (single)")
        end
    end
    return tables
end

function CIATables(fn::String;
                   extrapolate::Bool=false,
                   singles::Bool=false,
                   verbose::Bool=true)
    CIATables(readcia(fn), extrapolate=extrapolate, singles=singles, verbose=verbose)
end

#-------------------------------------------------------------------------------
# interpolating k, the raw CIA values

function (tables::CIATables)(ν, T)
    k = zero(T)
    #look at each grid
    for Φ ∈ tables.Φ
        if Φ.G.xa <= ν <= Φ.G.xb
            #inside wavenumber range
            if Φ.G.ya <= T <= Φ.G.yb
                #interpolate inside the grid of data
                k += exp(Φ(ν, T))
            elseif tables.extrapolate
                #otherwise extrapolate using flat boundary values, if desired
                k += exp(Φ(ν, T > Φ.G.yb ? Φ.G.yb : Φ.G.ya))
            end
        end
    end
    #optionally include the weird ranges at a single temperature
    if tables.singles
        for ϕ ∈ tables.ϕ
            #wavenumber range
            if ϕ.r.xa <= ν <= ϕ.r.xb
                k += exp(ϕ(ν))
            end
        end
    end
    return k
end

#-------------------------------------------------------------------------------
# computing collision induced absorption cross-sections

export cia, cia!

#the Loschmidt number, but in molecules/cm^3 then squared [molecules^2/cm^6]
const Locmsq = 7.21879268e38

"""
    cia(k, T, Pₐ, P₁, P₂)

Compute a collision induced absorption cross-section

# Arguments
* `k`: absorption coefficient [cm``^5``/molecule``^2``]
* `T`: temperature [K]
* `Pₐ`: total air pressure [Pa]
* `P₁`: partial pressure of first gas [Pa]
* `P₂`: partial pressure of second gas [Pa]
"""
function cia(k, T, Pₐ, P₁, P₂)
    #number densities of gases, in amagats
    ρ₁ = (P₁/𝐀)*(𝐓₀/T)
    ρ₂ = (P₂/𝐀)*(𝐓₀/T)
    #number density of air, in molecules/cm^3
    ρₐ = 1e-6*Pₐ/(𝐤*T)
    #σ in cm^2/molecule, converting k from cm^5/molecule^2 to cm^-1/amagat^2
    (k*Locmsq)*ρ₁*ρ₂/ρₐ
end

"""
    cia(ν, x::CIATables, T, Pₐ, P₁, P₂)

Compute a collision induced absorption cross-section after retrieving the total absorption coefficient from a [`CIATables`](@ref) object

# Arguments
* `ν`: wavenumber [cm``^{-1}``]
* `x`: [`CIATables`](@ref) object
* `T`: temperature [K]
* `Pₐ`: total air pressure [Pa]
* `P₁`: partial pressure of first gas [Pa]
* `P₂`: partial pressure of second gas [Pa]
"""
function cia(ν::Real, x::CIATables, T, Pₐ, P₁, P₂)
    #first retrieve the absorption coefficient from the interpolator
    k = x(ν, T) #cm^5/molecule^2
    #then compute the cross-section
    cia(k, T, Pₐ, P₁, P₂) #cm^2/molecule
end

"""
    cia!(σ::AbstractVector, ν::AbstractVector, x::CIATables, T, Pₐ, P₁, P₂)

Compute a vector of collision induced absorption cross-sections in-place, retrieving the total absorption coefficient from a [`CIATables`](@ref) object.

# Arguments
* `σ`: vector to store computed cross-sections
* `ν`: vector of wavenumbers [cm``^{-1}``]
* `x`: [`CIATables`](@ref) object
* `T`: temperature [K]
* `Pₐ`: total air pressure [Pa]
* `P₁`: partial pressure of first gas [Pa]
* `P₂`: partial pressure of second gas [Pa]
"""
function cia!(σ::AbstractVector, ν::AbstractVector, x::CIATables, T, Pₐ, P₁, P₂)
    @assert length(σ) == length(ν)
    for i = 1:length(σ)
        @inbounds σ[i] += cia(ν[i], x, T, P₁, P₂, Pₐ)
    end
end

"""
    cia(ν::AbstractVector, x::CIATables, T, Pₐ, P₁, P₂)

Compute a vector of collision induced absorption cross-sections, retrieving the total absorption coefficient from a [`CIATables`](@ref) object.

# Arguments
* `ν`: vector of wavenumbers [cm``^{-1}``]
* `x`: [`CIATables`](@ref) object
* `T`: temperature [K]
* `Pₐ`: total air pressure [Pa]
* `P₁`: partial pressure of first gas [Pa]
* `P₂`: partial pressure of second gas [Pa]
"""
function cia(ν::AbstractVector, x::CIATables, T, Pₐ, P₁, P₂)
    σ = zeros(Float64, length(ν))
    cia!(σ, ν, x, T, P₁, P₂, Pₐ)
    return σ
end

"""
    cia(ν, x::CIATables, T, Pₐ, g₁::Gas, g₂::Gas)

Compute a collision induced absorption cross-section, retrieving the total absorption coefficient from a [`CIATables`](@ref) object and computing partial pressures from gas objects.

# Arguments
* `ν`: vector of wavenumbers [cm``^{-1}``]
* `x`: [`CIATables`](@ref) object
* `T`: temperature [K]
* `Pₐ`: total air pressure [Pa]
* `g₁`: gas object representing the first component of the CIA pair
* `g₂`: gas object representing the second component of the CIA pair
"""
function cia(ν::Real, x::CIATables, T, Pₐ, g₁::Gas, g₂::Gas)
    P₁ = Pₐ*concentration(g₁, T, Pₐ)
    P₂ = Pₐ*concentration(g₂, T, Pₐ)
    cia(ν, x, T, Pₐ, P₁, P₂)
end

#-------------------------------------------------------------------------------
# computing cross-sections efficiently with known gas objects

export CIA

"""
Container for a [`CIATables`](@ref) object and the two gasses representing the CIA components. Specializes with the type of each gas for fast retreival of absorption cross-sections from CIA data and partial pressures.

| Field | Type | Description |
| ----- | :--- | :---------- |
| `name` | `String` | molecular symbol, i.e. `"CO2-H2"` |
| `formulae` | `Tuple{String,String}` | split molecular formulae, i.e `("CO2", "H2")` |
| `x` | `CIATables` | collision induced absorption tables |
| `g₁` | `<:Gas` | first component of CIA pair |
| `g₁` | `<:Gas` | second component of CIA pair |

# Constructors

    CIA(ciatables::CIATables, g₁::Gas, g₂::Gas)

The name and formulae are taken from `ciatables`.

    CIA(ciatables::CIATables, gases::Tuple)

Using the formulae in `ciatables`, the correct pair of gases is automatically selected from a VarArg collection of gases.

# Example

A `CIA` object is [function-like](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects). Use it like a function, passing it wavenumber, temperature, and pressure arguments to compute an absorption cross-section. Underneath, the [`CIATables`](@ref) object is interpolated and partial pressures are computed using the concentrations stored with the gases.

```julia
#load gases
ν = LinRange(1, 2500, 2500);
Ω = AtmosphericDomain((100,350), 8, (0.9,2e5), 16);
co2 = WellMixedGas("data/par/CO2.par", 0.96, ν, Ω);
ch4 = WellMixedGas("data/par/CH4.par", 1e-6, ν, Ω);

#create CIA object
co2ch4 = CIA(CIATables("data/cia/CO2-CH4_2018.cia"), co2, ch4);

#compute a cross-section
ν = 667;
T = 250;
P = 1e5;
σ = co2ch4(ν, T, P)
```
"""
struct CIA{T,U}
    name::String
    formulae::Tuple{String,String}
    x::CIATables
    g₁::T
    g₂::U
end

function CIA(ciatables::CIATables, g₁::Gas, g₂::Gas)
    @assert g₁.formula in ciatables.formulae "gas $(g₁.formula) not found in $(ciatables.name) CIATables"
    @assert g₂.formula in ciatables.formulae "gas $(g₂.formula) not found in $(ciatables.name) CIATables"
    CIA(ciatables.name, ciatables.formulae, ciatables, g₁, g₂)
end

function findgas(f::String, cianame::String, gases::Gas...)
    idx = findall(g -> g.formula == f, gases)
    @assert length(idx) > 0 "pairing failed for $cianame CIA, gas $f is missing"
    @assert length(idx) == 1 "pairing failed for $cianame CIA, duplicate $f gases found"
    return gases[idx[1]]
end

function CIA(ciatables::CIATables, gases::Tuple{})
    error("no Gas objects provided, cannot create CIA object")
end

function CIA(ciatables::CIATables, gases::Tuple)
    #gas formulae
    f₁, f₂ = ciatables.formulae
    #find matching gases
    g₁, g₂ = findgas(f₁, ciatables.name, gases...), findgas(f₂, ciatables.name, gases...)
    #make a CIA object
    CIA(ciatables, g₁, g₂)
end

(χ::CIA)(ν, T, P) = cia(ν, χ.x, T, P, χ.g₁, χ.g₂)
