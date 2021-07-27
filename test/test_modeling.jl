Pₛ = 1e5
Pₜ = 1
g = 9.8
fμ(T,P) = 0.029
θ = acos(2/3)
fstellar(ν) = 1336*normplanck(ν, 5500)
falbedo(ν) = 0.2

ν = [1:2500; 2510:10:75000]
Ω = AtmosphericDomain((100,350), 8, (0.9,2e5), 16)

co2 = WellMixedGas("HITRAN/CO2.par", 400e-6, ν, Ω, PHCO2!)
ch4 = WellMixedGas("HITRAN/CH4.par", 1e-5, ν, Ω, lorentz!)
h2o = VariableGas("HITRAN/H2O.par", (T,P)->0.01, ν, Ω, voigt!)

co2co2 = CIATables("HITRAN/CO2-CO2_2018.cia", verbose=false)
co2ch4 = CIATables("HITRAN/CO2-CH4_2018.cia", verbose=false)

Γ = MoistAdiabat(280, Pₛ, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Tstrat=190)
h2o = condensibleprofile(h2o, Γ, T->0.8*psatH2O(T))
P = reverse(logrange(1e5, 1, 250))
T = Γ.(P)
A = AcceleratedAbsorber(P, T, co2, ch4, h2o, co2co2, co2ch4, (ν,T,P)->1e-30)

τ = opticaldepth(Pₛ, Pₜ, g, Γ, fμ, θ, A)
@test trapz(ν, τ) ≈ 65001.19194742341

t = transmittance(Pₛ, Pₜ, g, Γ, fμ, θ, A)
@test trapz(ν, t) ≈ 74670.50132922476

olr = outgoing(Pₛ, g, Γ, fμ, A)
@test trapz(ν, olr) ≈ 336.96022895994196

Fₜ⁻, Fₜ⁺ = topfluxes(Pₛ, g, Γ, fμ, fstellar, falbedo, A)
@test trapz(ν, Fₜ⁻) ≈ 283.5281957875912
@test trapz(ν, Fₜ⁺) ≈ 393.3404384262924

Fₛ⁻, Fₛ⁺ = bottomfluxes(Pₛ, g, Γ, fμ, fstellar, falbedo, A)
@test trapz(ν, Fₛ⁻) ≈ 315.6348654709627
@test trapz(ν, Fₛ⁺) ≈ 405.0246045370683