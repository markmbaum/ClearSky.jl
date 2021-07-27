pushfirst!(LOAD_PATH, "C:/Users/markm/Documents/code")

using ClearSky
using Test
using SpecialFunctions
using Downloads
using Base.Filesystem: rm

using ClearSky: faddeyeva, MOLPARAM

@testset "Faddeyeva" begin include("test_faddeyeva.jl") end
@testset "Molparam" begin include("test_molparam.jl") end

#------------------------------------------------------------------------------

urlco2 = "https://www.dropbox.com/s/7o61uofl39dl7pg/CO2.par?dl=0"
urlch4 = "https://www.dropbox.com/s/qdtyrdlkzcuekoc/CH4.par?dl=0"
urlh2o = "https://www.dropbox.com/s/xd3xcgqneim5w4c/H2O.par?dl=0"
urlco2co2 = "https://www.dropbox.com/s/83u9j776v13lukr/CO2-CO2_2018.cia?dl=0"
urlco2ch4 = "https://www.dropbox.com/s/755bffs13yap3ry/CO2-CH4_2018.cia?dl=0"

fnco2 = Downloads.download(urlco2)
fnch4 = Downloads.download(urlch4)
fnh2o = Downloads.download(urlh2o)
fnco2co2 = Downloads.download(urlco2co2)
fnco2ch4 = Downloads.download(urlco2ch4)

ν = [1:2500; 2510:10:75000]
Ω = AtmosphericDomain((100,350), 12, (0.9,2e5), 24)
Γ = MoistAdiabat(280, 1e5, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Tstrat=200);

co2 = WellMixedGas(fnco2, 400e-6, ν, Ω)
ch4 = WellMixedGas(fnch4, 1e-5, ν, Ω)
h2o = VariableGas(fnh2o, (T,P)->0.01, ν, Ω)
h2o = condensibleprofile(h2o, Γ, T->0.8*psatH2O(T))

co2co2 = CIATables(fnco2co2)
co2ch4 = CIATables(fnco2ch4)

P = reverse(logrange(1e5, 1, 250))
T = Γ.(P)
A = AcceleratedAbsorber(P, T, co2, ch4, h2o, co2co2, co2ch4, (ν,T,P)->1e-30)

Pₛ = 1e5
Pₜ = 1
g = 9.8
fμ(T,P) = 0.029
θ = acos(2/3)
fstellar(ν) = 1336*normplanck(ν, 5500)
falbedo(ν) = 0.2

τ = opticaldepth(Pₛ, Pₜ, g, Γ, fμ, θ, A)
t = transmittance(Pₛ, Pₜ, g, Γ, fμ, θ, A)
olr = outgoing(Pₛ, g, Γ, fμ, A)
Fₜ⁻, Fₜ⁺ = topfluxes(Pₛ, g, Γ, fμ, fstellar, falbedo, A)
Fₛ⁻, Fₛ⁺ = bottomfluxes(Pₛ, g, Γ, fμ, fstellar, falbedo, A)
