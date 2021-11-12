export readpar, SpectralLines

#-------------------------------------------------------------------------------
#mappin between isotopologue labels and integer indices

const ISOINDEX = Dict{Char,Int64}(
    '1'=>1,  '2'=>2,  '3'=>3,  '4'=>4,  '5'=>5,  '6'=>6,
    '7'=>7,  '8'=>8,  '9'=>9,  '0'=>10, 'A'=>11, 'B'=>12,
    'C'=>13, 'D'=>14, 'E'=>15, 'F'=>16, 'G'=>17, 'H'=>18,
    'I'=>19, 'J'=>20, 'K'=>21, 'L'=>22, 'M'=>23, 'N'=>24,
    'O'=>25, 'P'=>26, 'Q'=>27, 'R'=>28, 'S'=>29, 'T'=>30,
    'U'=>31, 'V'=>32, 'W'=>33, 'X'=>34, 'Y'=>35, 'Z'=>36
)

#-------------------------------------------------------------------------------
# struct for molecule parameters

struct MolParam
    #molecule number
    M::Int64
    #molecular formula
    formula::String
    #molecule name
    name::String
    #global isotopologue codes
    I::Vector{Int64}
    #isopologue formulae
    isoform::Vector{String}
    #AFGL isotopologue codes
    AFGL::Vector{Int64}
    #abundance fractions
    A::Vector{Float64}
    #molecular masses [kg/mole]
    μ::Vector{Float64}
    #Qref
    Qref::Vector{Float64}
    #flag, has interpolating chebyshev polynomial
    hascheb::Vector{Bool}
    #length of cheby polys
    ncheb::Vector{Int64}
    #maximum rel err of cheb polys
    maxrelerr::Vector{Float64}
    #chebyshev expansion coefficients
    cheb::Vector{Vector{Float64}}
end

#for empty structs
MolParam() = MolParam(-1, "", "", [], [], [], [], [], [], [], [], [], [])

#-------------------------------------------------------------------------------

"""
    readpar(filename; νmin=0, νmax=Inf, Scut=0, I=[], maxlines=-1, progress=true)

Read an absoption line file from the HITRAN database, which should have the ".par" extension. These files are available at [`https://hitran.org/lbl`](https://hitran.org/lbl) after registering for a free account.

# Keyword Arguments
* `νmin`: smallest line wavenumber to include
* `νmax`: largest line wavenumber to include
* `Scut`: smallest spectral line intensity
* `I`: array of isotopologue numbers to include (excludes all others)
* `maxlines`: maximum number of lines to include (includes only the most intense `maxlines` lines)
* `progress`: whether to display the progress meter

A dictionary of vectors is returned, reflecting the definitions from
1. [HITRAN website](https://hitran.org/docs/definitions-and-units`)
2. [Rothman, Laurence S., et al. "The HITRAN 2004 molecular spectroscopic database." Journal of quantitative spectroscopy and radiative transfer 96.2 (2005): 139-204.](https://www.sciencedirect.com/science/article/abs/pii/S0022407313002859)

| Key | Vector Type | Description |
| --- | :---------- | :---------- |
| `M`    | `Int16`   | [HITRAN molecular identification number](https://hitran.org/docs/molec-meta) |
| `I`    | `Char`    | HITRAN isotopologue identification symbol |
| `ν`    | `Float64` | spectral line wavenumber [cm``^{-1}``] in a vacuum |
| `S`    | `Float64` | spectral line intensity [cm``^{-1}``/(molecule``\\cdot``cm``^{-2}``)] at 296 K |
| `A`    | `Float64` | Einstein-A coefficient (s``^{-1}``) of a transition |
| `γa`   | `Float64` | air-broadened half width at half maximum (HWHM) [cm``^{-1}``/atm] at 296 K and 1 atm |
| `γs`   | `Float64` | self-broadened half width at half maximum (HWHM) [cm``^{-1}``/atm] at 296 K and 1 atm |
| `Epp`  | `Float64` | lower-state energy of the transition [cm``^{-1}``] |
| `na`   | `Float64` | coefficient of temperature dependence of air-broadened half width |
| `δa`   | `Float64` | pressure shift [cm``^{-1}``/atm] at 296 K and 1 atm of the line position with respect to vacuum transition wavenumber |
| `Vp`   | `String`  | upper-state "global" quanta |
| `Vpp`  | `String`  | lower-state "global" quanta |
| `Qp`   | `String`  | upper-state "local" quanta |
| `Qpp`  | `String`  | lower-state "local" quanta |
| `Ierr` | `String`  | uncertainty indices |
| `Iref` | `String`  | reference indices |
| `*`    | `Char`    | flag (?) |
| `gp`   | `String`  | statistical weight of upper state |
| `gpp`  | `String`  | statistical weight of lower state |
"""
function readpar(filename::String;
                 νmin::Real=0,
                 νmax::Real=Inf,
                 Scut::Real=0,
                 I::Vector=[],
                 maxlines::Int=-1,
                 progress::Bool=true)
    @assert filename[end-3:end] == ".par" "expected file with .par extension, downloaded from https://hitran.org/lbl/"
    lines = readlines(filename)
    N = length(lines)
    #hard code the format to make sure nothing slows it down ¯\_(ツ)_/¯
    par = Dict(
        "M"   =>zeros(Int16, N),
        "I"   =>Vector{Char}(undef, N),
        "ν"   =>zeros(Float64, N),
        "S"   =>zeros(Float64, N),
        "A"   =>zeros(Float64, N),
        "γa"  =>zeros(Float64, N),
        "γs"  =>zeros(Float64, N),
        "Epp" =>zeros(Float64, N),
        "na"  =>zeros(Float64, N),
        "δa"  =>zeros(Float64, N),
        "Vp"  =>Vector{String}(undef, N),
        "Vpp" =>Vector{String}(undef, N),
        "Qp"  =>Vector{String}(undef, N),
        "Qpp" =>Vector{String}(undef, N),
        "Ierr"=>Vector{String}(undef, N),
        "Iref"=>Vector{String}(undef, N),
        "*"   =>Vector{Char}(undef, N),
        "gp"  =>Vector{String}(undef, N),
        "gpp" =>Vector{String}(undef, N)
    )
    #parse the values 
    if progress
        prg = Progress(length(lines), 0.1, "Parsing \"$(splitdir(filename)[end])\" ")
    end
    @inbounds for i ∈ eachindex(lines)
        #get the line in question
        line = lines[i]
        #parse all the fields
        par["M"][i]    = parse(Int16, line[1:2])
        par["I"][i]    = line[3]
        par["ν"][i]    = parse(Float64, line[4:15])
        par["S"][i]    = parse(Float64, line[16:25])
        par["A"][i]    = parse(Float64, line[26:35])
        par["γa"][i]   = parse(Float64, line[36:40])
        par["γs"][i]   = parse(Float64, line[41:45])
        par["Epp"][i]  = parse(Float64, line[46:55])
        par["na"][i]   = parse(Float64, line[56:59])
        par["δa"][i]   = parse(Float64, line[60:67])
        par["Vp"][i]   = line[68:82]
        par["Vpp"][i]  = line[83:97]
        par["Qp"][i]   = line[98:112]
        par["Qpp"][i]  = line[113:127]
        par["Ierr"][i] = line[128:133]
        par["Iref"][i] = line[134:145]
        par["*"][i]    = line[146]
        par["gp"][i]   = line[147:153]
        par["gpp"][i]  = line[154:160]
        #update progress meter
        progress && next!(prg)
    end
    #filtering
    mask = ones(Bool, N)
    mask .&= par["ν"] .>= νmin
    mask .&= par["ν"] .<= νmax
    mask .&= par["S"] .>= Scut
    if length(I) > 0
        for j = 1:N
            #isotopologue character
            c = par["I"][j]
            #isotopologue integer
            i = ISOINDEX[c]
            #check if I and its integer counterpart are not in I
            if !(c in I) & !(i in I)
                mask[j] = false
            end
        end
    end
    #check that there will be information left
    @assert any(mask) "par information has been filtered to nothing!"
    #slice all the par arrays
    for (k,v) ∈ par
        par[k] = v[mask]
    end
    #take strongest lines if needed
    if maxlines > 0
        @assert maxlines > 0 "maxlines must be a positive integer"
        #check if there are fewer lines than desired already
        if N > maxlines
            idx = reverse(sortperm(par["S"]))[1:maxlines]
            for (k,v) in par
                par[k] = v[idx]
            end
        end
    end
    #ensure sorted by wavenumber
    idx = sortperm(par["ν"])
    for (k,v) ∈ par
        par[k] = v[idx]
    end
    return par
end

"""
Organizing type for spectral line data of a single gas

| Field | Type | Description |
| ----- | :--- | :---------- |
| `name`    | `String` | gas name |
| `formula` | `String` | gas formula |
| `N`       | `Int64` | number of lines |
| `M`       | `Int16` | see [`readpar`](@ref) |
| `I`       | `Vector{Int16}` | see [`readpar`](@ref) |
| `μ`       | `Vector{Float64}` | molar mass of isotopologues [kg/mole] |
| `A`       | `Vector{Float64}` | isotopologue abundance (Earth) |
| `ν`       | `Vector{Float64}` | see [`readpar`](@ref) |
| `S`       | `Vector{Float64}` | see [`readpar`](@ref) |
| `γa`      | `Vector{Float64}` | see [`readpar`](@ref) |
| `γs`      | `Vector{Float64}` | see [`readpar`](@ref) |
| `Epp`     | `Vector{Float64}` | see [`readpar`](@ref) |
| `na`      | `Vector{Float64}` | see [`readpar`](@ref) |

# Constructors

    SpectralLines(par::Dict)

Construct a `SpectralLines` object from a dictionary of line data. That dictionary can be created with [`readpar`](@ref).

    SpectralLines(filename, νmin=0, νmax=Inf, Scut=0, I=[], maxlines=-1)

Read a `.par` file directly into a `SpectralLines` object. Keyword arguments are passed through to [`readpar`](@ref).
"""
struct SpectralLines
    #gas name
    name::String
    #gas formula
    formula::String
    #number of lines
    N::Int64
    #molecule number (https://hitran.org/docs/molec-meta/)
    M::Int16
    #ordered isotopologue integer (https://hitran.org/docs/iso-meta/)
    I::Vector{Int16}
    #molar mass of isotopologues [kg/mole]
    μ::Vector{Float64}
    #isotopologue abundance in HITRAN
    A::Vector{Float64}
    #wavenumbers of absorption lines [cm^-1]
    ν::Vector{Float64}
    #spectral line intensity [cm^-1/(molecule*cm^-2)]
    S::Vector{Float64}
    #air-broadened half-width [cm^-1/atm]
    γa::Vector{Float64}
    #self-broadened half-width [cm^-1/atm]
    γs::Vector{Float64}
    #lower state energy, E'' [cm^-1]
    Epp::Vector{Float64}
    #temperature-dependence exponent for γair [unitless]
    na::Vector{Float64}
end

function SpectralLines(par::Dict)
    #length
    N = length(par["ν"])
    #make sure there is only one molecule present
    @assert length(unique(par["M"])) == 1 "SpectralLines objects must contain only one molecule's lines"
    #get the molecule name and formula
    M = par["M"][1]
    name = MOLPARAM[M].name
    form = MOLPARAM[M].formula
    #create arrays for isotopologue index, abundance, and molar mass
    I = map(i->ISOINDEX[par["I"][i]], 1:N)
    A = map(i->MOLPARAM[M].A[I[i]], 1:N)
    μ = map(i->MOLPARAM[M].μ[I[i]], 1:N)
    #get sorting indices to make sure ν vector is in ascending order
    idx = sortperm(par["ν"])
    #create a SpectralLines structure
    SpectralLines(
        name,
        form,
        N,
        M,
        I[idx],
        μ[idx],
        A[idx],
        par["ν"][idx],
        par["S"][idx],
        par["γa"][idx],
        par["γs"][idx],
        par["Epp"][idx],
        par["na"][idx]
    )
end

SpectralLines(fn::String; kwargs...) = SpectralLines(readpar(fn; kwargs...))
