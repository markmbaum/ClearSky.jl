pushfirst!(LOAD_PATH, "C:/Users/markm/Documents/code")
println("WOO")
using ClearSky
using Test
using SpecialFunctions

using ClearSky: faddeyeva, MOLPARAM

@testset "Faddeyeva" begin include("test_faddeyeva.jl") end
@testset "Molparam" begin include("test_molparam.jl") end
@testset "Atmospherics" begin include("test_atmospherics.jl") end
@testset "Modeling" begin include("test_modeling.jl") end
