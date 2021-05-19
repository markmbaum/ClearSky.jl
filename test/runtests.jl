using ClearSky, Test
using SpecialFunctions

@testset "Faddeeva" begin include("test_faddeeva.jl") end
@testset "Molparam" begin include("test_molparam.jl") end
