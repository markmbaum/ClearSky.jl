using ClearSky
using Test

using ClearSky: MOLPARAM

@testset "Molparam" begin include("test_molparam.jl") end
@testset "Gray" begin include("test_gray.jl") end
