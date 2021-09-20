using ClearSky
using Test

using ClearSky: MOLPARAM

@testset "Molparam" begin include("test_molparam.jl") end
@testset "Atmospherics" begin include("test_atmospherics.jl") end
@testset "Gray" begin include("test_gray.jl") end
#@testset "Modeling" begin include("test_modeling.jl") end
