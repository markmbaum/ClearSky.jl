using ClearSky
using Test
using SpecialFunctions

using ClearSky: faddeyeva, MOLPARAM

@testset "Faddeyeva" begin include("test_faddeyeva.jl") end
@testset "Molparam" begin include("test_molparam.jl") end
