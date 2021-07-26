using ClearSky
using Test
using SpecialFunctions

using ClearSky: faddeyeva, MOLPARAM

@testset "Faddeyeva" begin include("test_faddeyeva.jl") end
@testset "Molparam" begin include("test_molparam.jl") end

using Downloads

fn = Downloads.download("https://www.dropbox.com/s/qdtyrdlkzcuekoc/CH4.par?dl=0")
ν = LinRange(1, 2500, 2500)
Ω = AtmosphericDomain((100,350), 12, (0.9,2e5), 24)
co2 = WellMixedGas(fn, 270e-6, ν, Ω)
println(co2.name)